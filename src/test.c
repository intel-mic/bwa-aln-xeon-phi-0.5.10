

static void *thread_calculate(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	
	while(1){

		sem_wait(&finish_read);
	
		if(d->join) break;
	
		bwa_cal_sa_reg_gap(0, d->bwt, d->n_seqs, d->seqs, d->opt);

		sem_post(&finish_cal);
	}	
	
	return 0;
}



	while (1) {
		seqs0 = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual);	
		if(seqs0 == 0)
			break;
		tot_seqs += n_seqs;
		if(first_flag != 0){
			sem_wait(&finish_cal);	// wait for signal finish_cal
		}
		if(first_flag != 0){
			for (i = 0; i < n_seqs; ++i) {
				bwa_seq_t *p = seqs + i;
				err_fwrite(&p->n_aln, 4, 1, stdout);
				if (p->n_aln) err_fwrite(p->aln, sizeof(bwt_aln1_t), p->n_aln, stdout);
			} 
			bwa_free_read_seq(n_seqs, seqs);
		}
		seqs = seqs0;
		data.n_seqs = n_seqs; 
		data.seqs = seqs; 		
		sem_post(&finish_read);		// send signal finish_read
		first_flag = 1;
	}
	
	
	