bwa-aln-xeon-phi-0.5.10
============================

I. Overview
============================
This article provides a recipe for compiling, optimizing and running bwa-aln-xeon-phi-0.5.10 for Intel® Xeon® processors and Intel® Xeon Phi™ coprocessors.

II. Introduction
============================
BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack(ALN), BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. 

Original bwa-0.5.10 package can be downloaded from:
http://sourceforge.net/projects/bio-bwa/files/

This project bwa-aln-xeon-phi-0.5.10 optimizes bwa-0.5.10 ALN module performance on both Xeon and Xeon Phi platform, and supports symmetric running model on hybrid Xeon and Xeon Phi nodes.
Optimization methods include:

(1) Replace Pthreads with OpenMP for better load balance. 

(2) Overlap File IO and computation to improve thread utilization. 

(3) Use Intel TBB memory allocator for efficient memory management.

(4) Vectorization of performance critical loop. 

(5) Data prefetch Intrinsics to reduce memory latency. 

(6) Add task parallelism besides OpenMP to further improve load balance.

III. Preliminaries
============================
1. Install Intel® MPI Library 4.1 and Intel® C++ Compiler 13.3 or higher on your host system.
2. Download bwa-aln-xeon-phi-0.5.10 from https://github.com/intel-mic/bwa-aln-xeon-phi-0.5.10
3. Install NFS and start NFS service, export /home directory and mount it to Intel Xeon Phi:

	> service nfs start

	> vi /etc/exports

Add:

	/home          172.31.0.0/16(rw,insecure,no_root_squash,async)

Run:
	
	> exportfs -au

	> exportfs -ar

	> showmount -e

	Export list for Host:
	/home        172.31.0.0/16

Mout /home directory to Intel Xeon Phi:

	> service mpss stop

	> micctrl --addnfs=/home --dir=/home

	> service mpss start

4．Download zlib-1.2.8 from http://www.zlib.net/ and compile for Xeon Phi:

	> export CC=icc

	> export CFLAGS=-mmic

	> ./configure --prefix=/home/zlib

	> make

	icc -mmic  -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN -I. -c -o example.o test/example.c

	……

	> make install

zlib that support Xeon Phi will be installed at /home/zlib

Then upload dynamic libraries to Xeon Phi:

	> scp /home/zlib/lib/*.so* mic0:/lib64

5. Download public workload from here:

	Reference file:
	
	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/

	100bp read pair:
	
	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12044/sequence_read/SRR766060_1.filt.fastq.gz
	
	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12044/sequence_read/SRR766060_2.filt.fastq.gz

	then put them to the directory of bwa-aln-xeon-phi-0.5.10/data.
	
IV. Compiling bwa-aln-xeon-phi
============================
1.	Set up the Intel MPI and Intel Compiler environments:

	> source /opt/intel/impi/<version>/bin64/mpivars.sh

	> source /opt/intel/compserxe/bin/compilervars.sh intel64

	> vi ~/.bashrc

Add:

	export I_MPI_MIC=enable

	export I_MPI_MIC_POSTFIX=_mic

	export I_MPI_FABRICS=shm:tcp

	export I_MPI_PIN=enable

Run ~/.bashrc

Upload Intel MPI and Intel Compiler Libraries to Intel Xeon Phi:

	> scp /opt/intel/impi/<version>/mic/bin/mpiexec      mic0:/bin/

	> scp /opt/intel/impi/<version>/mic/bin/pmi_proxy    mic0:/bin/

	> scp /opt/intel/impi/<version>/mic/lib/lib*.so*     mic0:/lib64/

	> scp /opt/intel/composer_xe_<version>/compiler/lib/mic/*.so* mic0:/lib64

	> scp /opt/intel/composer_xe_<version>/tbb/lib/mic/*.so* mic0:/lib64

	> scp /opt/intel/composer_xe_<version>/mkl/lib/mic/lib*.so  mic0:/lib64

2.	Unpack the source code to any directory of /home and build the executables for Intel Xeon and Intel Xeon processor

	> tar –xzvf bwa-aln-xeon-phi-0.5.10.tar.gz

	> cd bwa-aln-xeon-phi-0.5.10/scripts

	> ./build-src

It will compile both executable binaries for Intel Xeon and Intel Xeon Phi.

V. Build bwa reference index file for bwa aln
============================

	> ./build-index

Which runs bwa index -a bwtsw $ref_file to construct the FM-index for the reference genome.

VI. Configure bwa aln for running
============================

	> vi aln-head
	
Modify configuration:

	host1=crt03					# name of host1
	host1_exe=bwa				# executable file name on host1
	host1_num_tasks=2			# number of tasks on host1
	host1_num_threads=16		# number of threads for each task
	host1_ratio=10				# input data file ratio for host1

	host2=${host1}-mic0			# host2: mic0 of host1
	host2_exe=bwa_mic			# executable file name on host2
	host2_num_tasks=3			# number of tasks on host2
	host2_num_threads=80		# number of threads for each task
	host2_ratio=5				# input data file ratio for mic0

	host3=${host1}-mic1			# host3: mic1 of host1
	host3_exe=bwa_mic			# executable file name on host3
	host3_num_tasks=3			# number of tasks on host3
	host3_num_threads=80		# number of threads for each task
	host3_ratio=5				# input data file ratio for mic1
	
	src_directory=../src							# source code

	data_directory=../data/public_workload			# input data directory
	ref_file=${data_directory}/hs37d5.fa			# reference data file
	file1=${data_directory}/SRR766060_1.filt.fastq	# input data file1
	file2=${data_directory}/SRR766060_2.filt.fastq	# input data file2

	res_tmp_file1=tmp_aln_file1						# temporary output file name of each task
	res_tmp_file2=tmp_aln_file2

	res_file1=aln_1.sai								# aln output file name
	res_file2=aln_2.sai

	# for bwa sampe output and picard output
	sampe_output=output.sam
	picard_output=sorted.sam
	
Note that: 

(1) Intel Xeon Phi is considered as a independent host

(2) num_tasks * num_threads must <= number of logical cores on this host

(3) The ratio of host:host-mic0=host1_ratio:host2_ratio according to the speed of their independent running

(4) Set aln output file name by res_file1 and res_file2, needn’t set by -f output file name

VII. Run bwa aln on Intel Xeon
============================

1. Split input data file for multiple tasks

	> ./aln-cpu-split

2. Run 

	> ./aln-cpu-run

Which runs bwa aln on Intel Xeon.

VIII. Run bwa aln on Intel Xeon Phi
============================

1. Split input data file for multiple tasks

	> ./aln-1mic-split

2. Run 

	> ./aln-1mic-run

Which runs bwa aln on Intel Xeon.

IX. Run bwa aln on Intel Xeon and Intel Xeon Phi symmetric
============================

1. Modify host1_ratio and host2_ratio in aln_head according to the performance on Intel Xeon and Intel Xeon Phi

2. Split input data file for multiple tasks

	> ./aln-cpu-and-1mic-split

3. Run

	> ./aln-cpu-and-1mic-run

Which runs bwa aln on Intel Xeon and Intel Xeon Phi.	

Contributors:
============================

You, Liang <liang.you@intel.com>

Congdon, Charles <charles.w.congdon@intel.com>

Ramanujam, Ram <ram.ramanujam@intel.com>






