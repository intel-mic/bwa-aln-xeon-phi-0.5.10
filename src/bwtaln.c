#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"

#include <pthread.h>
#include <semaphore.h>

bwt_t *bwt[2];			// For task parallelism, shared between tasks

// Optimizations
#ifndef HAVE_PTHREAD
#define OMP_OPT			// OpenMP multithread instead of pthread
#define MULTIBUFS_OPT	// Use multiple buffers to pipeline read and calculation
#endif

// Signals of finishing reading and calculation for each task
sem_t finish_read[MAX_TASKS], finish_cal[MAX_TASKS];	

#ifdef OMP_OPT
#include "omp.h"		// Include openmp header file
#define OMP_GRAIN 100	// Dynamic scheduling grain size of OpenMP multithread
#endif

double t_cal[MAX_TASKS] = {0};	// Calculation fime for each task

#ifdef MULTIBUFS_OPT
#define NUM_BUFFERS 10	// Number of buffers
#define START_INDEX 1	// Start calculation index after reading
static int tot_seqs[MAX_TASKS] = {0};
sem_t finish_read_array[MAX_TASKS][NUM_BUFFERS], finish_cal_array[MAX_TASKS][NUM_BUFFERS];
#endif	// MULTIBUFS_OPT

#ifdef USE_MPI
#include "mpi.h"
extern int mpi_id, mpi_numprocs, mpi_name_len;
extern char mpi_name[MPI_MAX_PROCESSOR_NAME];
#endif

// Timing function
double mysecond()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

gap_opt_t *gap_init_opt()
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
	o->max_diff = -1; o->max_gapo = 1; o->max_gape = 6;
	o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32; o->max_seed_diff = 2;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

// width must be filled as zero
static int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = rbwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
			k = rbwt->L2[c] + ok + 1;
			l = rbwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = rbwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int num_tasks, int task_id)
{
	int ncpus = kmp_get_affinity_max_proc();

#ifdef OMP_OPT
#pragma omp parallel num_threads(opt->n_threads)
#endif
{
#ifdef USE_MPI

	if(omp_get_thread_num() == 0)
		fprintf(stderr, "I am on %s, task_id = %d, num_tasks = %d, num_threads = %d\n", mpi_name, task_id, num_tasks, omp_get_num_threads());
		
#endif

	// Binding threads
	int ithr = task_id * (ncpus / num_tasks) + omp_get_thread_num();	
	int icpu = (ithr + 1) % ncpus;								// +1: binding logical cores, [1 ~ 243, 0]
	
	kmp_affinity_mask_t m;										// set kmp affinity
	kmp_create_affinity_mask(&m);
	kmp_set_affinity_mask_proc(icpu, &m);
	kmp_set_affinity(&m);
	kmp_destroy_affinity_mask(&m);

	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;

	// Initialization priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	
#ifdef OMP_OPT
#pragma omp for schedule (dynamic, OMP_GRAIN)
#endif
	for (i = 0;  i < n_seqs; ++i) {
	//for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
#ifdef HAVE_PTHREAD
		if (i % opt->n_threads != tid) continue;
#endif
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		seq[0] = p->seq; seq[1] = p->rseq;
		if (max_l < p->len) {
			max_l = p->len;
			w[0] = (bwt_width_t*)realloc(w[0],(max_l + 1)* sizeof(bwt_width_t));
			w[1] = (bwt_width_t*)realloc(w[1],(max_l + 1)* sizeof(bwt_width_t));
			memset(w[0], 0, (max_l + 1) * sizeof(bwt_width_t));
			memset(w[1], 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], p->len, seq[0], w[0]);
		bwt_cal_width(bwt[1], p->len, seq[1], w[1]);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
		if (p->len > opt->seed_len) {
			bwt_cal_width(bwt[0], opt->seed_len, seq[0] + (p->len - opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len, seq[1] + (p->len - opt->seed_len), seed_w[1]);
		}
		// core function
		p->aln = bwt_match_gap(bwt, p->len, seq, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
		// store the alignment
		free(p->name); free(p->seq); free(p->rseq); free(p->qual);
		p->name = 0; p->seq = p->rseq = p->qual = 0;
	}
	free(seed_w[0]); free(seed_w[1]);
	free(w[0]); free(w[1]);
	gap_destroy_stack(stack);
}
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif

#ifdef MULTIBUFS_OPT	// MULTIBUFS_OPT

typedef struct {
	int tid;
	int join;
	bwt_t *bwt[2];
	
	int n_seqs;
	int *n_seqs_array;
	
	bwa_seq_t *seqs;
	bwa_seq_t **seqs_array;
	
	const gap_opt_t *opt;
	
	// For task parallelism
	const char* res_aln;	// output file name
	int num_tasks;			// number of tasks
	int task_id;			// task id
	FILE* output_file;		// output file name
	
} thread_aux_t;

// Calculation thread
static void *thread_calculate(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	
	int first_flag = 0;
	int while_flag = 0;
	
	// For task parallelism
	const char* res_aln = d->res_aln;
	int num_tasks = d->num_tasks;	
	int task_id = d->task_id;
	FILE* output_file = d->output_file;
	
	// Total sequences
	tot_seqs[task_id] = 0;

	double t, t1;
	while(while_flag == 0){

		if(first_flag == 0){	// firt time wait the signal finish_read to start
			sem_wait(&finish_read[task_id]);
			first_flag = 1;
		}

		bwa_seq_t **seqs_array = d->seqs_array;
		int *n_seqs_array = d->n_seqs_array;

		int i, j;
		// For all data buffers
		for(i = 0; i < NUM_BUFFERS; ++i){
		
			sem_wait(&finish_read_array[task_id][i]);	// wait signal: finish reading data for each buffer

			bwa_seq_t *p_seqs = seqs_array[i];
			int p_n_seqs = n_seqs_array[i];

			// If the data buffer has contant, calculation the data buffer
			if(p_seqs != 0){
				t = mysecond();
				bwa_cal_sa_reg_gap(0, d->bwt, p_n_seqs, p_seqs, d->opt, num_tasks, task_id);
				t1 =  mysecond() - t; 
				t_cal[task_id] += t1;
				
#ifdef USE_MPI

				fprintf(stderr, "I am on %s, task_id = %d, [bwa_aln_core] calculate SA coordinate %.2f sec\n", mpi_name, task_id, t1);

				// Write data after calculation
				t = mysecond();
				for (j = 0; j < p_n_seqs; ++j) {
					bwa_seq_t *p = p_seqs + j;
					err_fwrite(&p->n_aln, 4, 1, output_file);
					if (p->n_aln) err_fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, output_file);					
				}

				// Should free calculation buffer each time
				bwa_free_read_seq(p_n_seqs, p_seqs);
				
				// Summary number of sequences that have been processed
				tot_seqs[task_id] += p_n_seqs;
				
				fprintf(stderr, "I am on %s, task_id = %d, [bwa_aln_core] %d sequences have been processed.\n\n", mpi_name, task_id, tot_seqs[task_id]);

#endif	// USE_MPI
	
			}else{	// p_seqs == 0, finish reading all data, exit while loop
				while_flag = 1;				
				break;	// break NUM_BUFFERS loop
			}

			sem_post(&finish_cal_array[task_id][i]);	// send signal: finish calculation for each buffer
		}

	}// end while
	
	// Finish all calculation, send signal finish_cal
	sem_post(&finish_cal[task_id]);
	
	return 0;
}

#endif	// MULTIBUFS_OPT

bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa, int num_tasks, int task_id)
{
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE) which |= 4;
		if (mode & BWA_MODE_BAM_READ1) which |= 1;
		if (mode & BWA_MODE_BAM_READ2) which |= 2;
		if (which == 0) which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else ks = bwa_seq_open(fn_fa, num_tasks, task_id);
	return ks;
}

void bwa_aln_core(const char *prefix, const char *fn_fa, const char *res_aln, const gap_opt_t *opt, int num_tasks, int task_id)
{
#ifndef MULTIBUFS_OPT
	int i, n_seqs, tot_seqs = 0;
#else
	int i;
#endif
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	double t, t1;
	double start;
	double end;
	//bwt_t *bwt[2];

	// Only task 0 load BWT, shared with other tasks
	if(task_id == 0){
		fprintf(stderr, "Prefix = %s, fn_fa = %s, res_aln = %s, num_tasks = %d\n", prefix, fn_fa, res_aln, num_tasks);
		start = mysecond();
		{ // load BWT
			char *str = (char*)calloc(strlen(prefix) + 10, 1);
			strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
			strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
			free(str);
		}
		end = mysecond() - start;
		fprintf(stderr, "I am on %s, task_id = %d, Loading index time %.2f sec\n", mpi_name, task_id, end);
	}
	
	// Initialization input file
	ks = bwa_open_reads(opt->mode, fn_fa, num_tasks, task_id);	

	// Initialization Output file
	FILE* output_file;
	char res_aln_file[256];
	sprintf(res_aln_file, "%s_%d.sai", res_aln, task_id);
	fprintf(stderr, "I am on %s, task_id = %d, Opening output fie %s\n", mpi_name, task_id, res_aln_file);
	output_file = fopen(res_aln_file, "wb");
	if(output_file == NULL) fprintf(stderr, "I am on %s, task_id = %d, Failed to open output fie %s\n", mpi_name, task_id, res_aln_file);	
	
	// Only task 0 output gap_opt_t
	if(task_id == 0)
		err_fwrite(opt, sizeof(gap_opt_t), 1, /*stdout*/output_file);

	// Core loop
	start = mysecond();
	
#ifdef MULTIBUFS_OPT	// MULTIBUFS_OPT
	
	int first_flag = 0;
	int first_flag_array = 0;
	int while_flag = 0;
	
	// Initialization all signals
	sem_init(&finish_read[task_id], 0, 0);
	sem_init(&finish_cal[task_id], 0, 0);
	for(i = 0; i < NUM_BUFFERS; ++i){
		sem_init(&finish_read_array[task_id][i], 0, 0);
		sem_init(&finish_cal_array[task_id][i], 0, 0);
	}

	// Creat a thread for calculation
	pthread_t tid;
	pthread_attr_t attr;
	thread_aux_t data;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	pthread_create(&tid, &attr, thread_calculate, &data);
	
	// Constant
	data.bwt[0] = bwt[0]; 
	data.bwt[1] = bwt[1];
	data.opt = opt;
	data.join = 0;
	
	// For task parallelism
	data.res_aln = res_aln;
	data.num_tasks = num_tasks;
	data.task_id = task_id;
	data.output_file = output_file;

	// Allocate NUM_BUFFERS buffer pointers
	// 0x40000 = 262144
	bwa_seq_t *seqs_array[NUM_BUFFERS];
	int n_seqs_array[NUM_BUFFERS] = {0};

	// Transfer buffer pointers to calculation thread
	data.n_seqs_array = n_seqs_array; 
	data.seqs_array = seqs_array;
	
	while (while_flag == 0) {

		// For number of buffers
		for(i = 0; i < NUM_BUFFERS; ++i){

			if(first_flag_array != 0){						
				sem_wait(&finish_cal_array[task_id][i]);	// wait sigal: finish calculation for each buffer
				first_flag_array = 1;						// first time needn't to wait
			}

			// Read data to buffers
			seqs_array[i] = bwa_read_seq(ks, 0x40000, &n_seqs_array[i], opt->mode, opt->trim_qual, num_tasks, task_id);

			// Finish reading all data
			if(seqs_array[i] == 0) {
				while_flag = 1;								// exit while loop
				sem_post(&finish_read_array[task_id][i]);	// send last singal to release the last wait in calculation thread
				break;										// break the NUM_BUFFERS
			}

			// When i reach START_INDEX, send signal finish_read to start calclation, only do once
			if(i >= START_INDEX && first_flag == 0){
				sem_post(&finish_read[task_id]);
				first_flag = 1;
			}
			
			sem_post(&finish_read_array[task_id][i]);		// send signal: finish reading data for each buffer
		}
		
	}// end while

	// Wait for all calculation to finish
	sem_wait(&finish_cal[task_id]);

	// Join calculation thread
	pthread_join(tid, 0);

	end = mysecond();
	double timeused = end - start;
	fprintf(stderr, "Total Calculation time %.2f sec, I am on %s, task_id = %d\n", t_cal[task_id], mpi_name, task_id);
	fprintf(stderr, "Total Time %.2f sec, I am on %s, task_id = %d\n\n", timeused, mpi_name, task_id); 	

#endif	// MULTIBUFS_OPT

#ifdef USE_MPI
	fclose(output_file);
#endif

	// All tasks should barrier here
#pragma omp barrier

	// Only task 0 need to destroy BWT
	if(task_id == 0){
		bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}
	bwa_seq_close(ks);
}

int bwa_aln(int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;

	opt = gap_init_opt();
	while ((c = getopt(argc, argv, "n:o:e:i:d:l:k:cLR:m:t:NM:O:E:q:f:b012IYB:")) >= 0) {
		switch (c) {
		case 'n':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'o': opt->max_gapo = atoi(optarg); break;
		case 'e': opte = atoi(optarg); break;
		case 'M': opt->s_mm = atoi(optarg); break;
		case 'O': opt->s_gapo = atoi(optarg); break;
		case 'E': opt->s_gape = atoi(optarg); break;
		case 'd': opt->max_del_occ = atoi(optarg); break;
		case 'i': opt->indel_end_skip = atoi(optarg); break;
		case 'l': opt->seed_len = atoi(optarg); break;
		case 'k': opt->max_seed_diff = atoi(optarg); break;
		case 'm': opt->max_entries = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': opt->max_top2 = atoi(optarg); break;
		case 'q': opt->trim_qual = atoi(optarg); break;
		case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
		case 'f': xreopen(optarg, "wb", stdout); break;
		case 'b': opt->mode |= BWA_MODE_BAM; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': opt->mode |= atoi(optarg) << 24; break;
		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa aln [options] <prefix> <in.fq>\n\n");
		fprintf(stderr, "Options: -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n",
				BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, "         -o INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, "         -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, "         -i INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, "         -d INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);
		fprintf(stderr, "         -m INT    maximum entries in the queue [%d]\n", opt->max_entries);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -M INT    mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, "         -O INT    gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, "         -E INT    gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, "         -R INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
        fprintf(stderr, "         -f FILE   file to write output to instead of stdout\n");
		fprintf(stderr, "         -B INT    length of barcode\n");
		fprintf(stderr, "         -c        input sequences are in the color space\n");
		fprintf(stderr, "         -L        log-scaled gap penalty for long deletions\n");
		fprintf(stderr, "         -N        non-iterative mode: search for all n-difference hits (slooow)\n");
		fprintf(stderr, "         -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
		fprintf(stderr, "         -b        the input read file is in the BAM format\n");
		fprintf(stderr, "         -0        use single-end reads only (effective with -b)\n");
		fprintf(stderr, "         -1        use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr, "         -2        use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr, "         -Y        filter Casava-filtered sequences\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[bwa_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}

	// Task parallelism
	int num_tasks = atoi(argv[optind+3]);	
	omp_set_nested(1);
	
	// Create multiple tasks
#pragma omp parallel num_threads(num_tasks)
{
	int task_id = omp_get_thread_num();
	bwa_aln_core(argv[optind], argv[optind+1], argv[optind+2], opt, num_tasks, task_id);
}

	free(opt);
	return 0;
}

/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
__cigar_op and __cigar_len while keeping stdaln stand alone */
bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar)
{
	uint32_t *cigar32;
	bwa_cigar_t *cigar;
	int i;
	cigar32 = aln_path2cigar32((path_t*) path, path_len, n_cigar);
	cigar = (bwa_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
                cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}

