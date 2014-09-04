Building and Optimizing BWA* ALN 0.5.10 for Intel® Xeon Phi™ Coprocessors
============================

Authors
============================
You, Liang (Intel); Congdon, Charles (Intel)

I. Overview
============================

This article provides a recipe for how to obtain, compile, and run an optimized version of BWA ALN 0.5.10 on Intel® Xeon® processors and Intel® Xeon Phi™ coprocessors.

The source for this version of BWA ALN 0.5.10 can be downloaded from: 

https://github.com/intel-mic/bwa-aln-xeon-phi-0.5.10

II. Introduction
============================

“BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack(ALN), BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp.” (http://bio-bwa.sourceforge.net/bwa.shtml)

This project, bwa-aln-xeon-phi-0.5.10, optimizes the performance of the BWA 0.5.10 ALN module on both Intel® Xeon® processors and Intel® Xeon Phi™ coprocessors, and supports a symmetric execution model that uses both architectures in cooperation for genome mapping.

Optimizations in this package include:

(1) Replace pthreads with OpenMP* for better load balancing. 

(2) Overlap file I/O and computation to improve thread utilization. 

(3) Use the Intel® Threading Building Blocks (Intel® TBB) memory allocator for efficient cross-thread memory management.

(4) Vectorization of performance-critical loops. 

(5) Data prefetch intrinsics to reduce memory latency. 

(6) Add task parallelism in addition to OpenMP to further improve load balancing.

The original bwa-0.5.10 package can be downloaded from:

http://sourceforge.net/projects/bio-bwa/files/

III. Preliminaries
============================

1. To build this package, install the Intel® MPI Library 4.1 and Intel® C++ Composer XE 13.3 or higher products on your host system.  Your host system must also have the Intel® MPSS for Linux* installed (installed separately from Intel® C++ Compiler XE).

2. Download bwa-aln-xeon-phi-0.5.10 from https://github.com/intel-mic/bwa-aln-xeon-phi-0.5.10

3. Install NFS and start the NFS service, export the /home directory, and mount it to the Intel® Xeon Phi™ coprocessor:

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

	Mount the /home directory to the Intel® Xeon Phi™ coprocessor:

	> service mpss stop
	
	> micctrl --addnfs=/home --dir=/home
	
	> service mpss start

4. Set up the Intel® MPI Library and Intel® C++ Compiler environments:	
	
	> source /opt/intel/impi/<version>/bin64/mpivars.sh
	
	> source /opt/intel/ composer_xe_<version>/bin/compilervars.sh intel64
	
	> vi ~/.bashrc

	Add:

	export I_MPI_MIC=enable
	
	export I_MPI_MIC_POSTFIX=_mic
	
	export I_MPI_FABRICS=shm:tcp
	
	export I_MPI_PIN=enable
	
	Run ~/.bashrc
	
5. Download zlib-1.2.8 from http://www.zlib.net/ and compile it for the Intel® Xeon Phi™ coprocessor:

	> export CC=icc
	
	> export CFLAGS=-mmic
	
	> ./configure --prefix=/home/zlib
	
	> make

	icc -mmic  -D_LARGEFILE64_SOURCE=1 -DHAVE_HIDDEN -I. -c -o example.o test/example.c

	……

	> make install

A version of zlib* that supports the Intel® Xeon Phi™ coprocessor will be installed at /home/zlib (note:  this library cannot be used on Intel® Xeon® processors)

Upload the zlib dynamic libraries to the Intel® Xeon Phi™ coprocessor:

	> scp /home/zlib/lib/*.so* mic0:/lib64	

IV. Compile bwa-aln-xeon-phi
============================

1. Upload the Intel® MPI Library and Intel® C Compiler components to the Intel® Xeon Phi™ coprocessor:

	> scp /opt/intel/impi/<version>/mic/bin/mpiexec      mic0:/bin/

	> scp /opt/intel/impi/<version>/mic/bin/pmi_proxy    mic0:/bin/

	> scp /opt/intel/impi/<version>/mic/lib/lib*.so*     mic0:/lib64/

	> scp /opt/intel/composer_xe_<version>/compiler/lib/mic/*.so* mic0:/lib64

	> scp /opt/intel/composer_xe_<version>/tbb/lib/mic/*.so* mic0:/lib64

	> scp /opt/intel/composer_xe_<version>/mkl/lib/mic/lib*.so  mic0:/lib64

2.	Unpack the source code to any directory of /home and build the executables for the Intel® Xeon® processor and the Intel® Xeon Phi™ coprocessor

	> tar –xzvf bwa-aln-xeon-phi-0.5.10.tar.gz

	> cd bwa-aln-xeon-phi-0.5.10/scripts
	
	> chmod +x *

	> ./build-src

This will build executables for both the Intel® Xeon® processor and Intel® Xeon Phi™ coprocessor.

Note:  If you unpacked and built zlib to some location other than /home/zlib, you will need to edit bwa-aln-xeon-phi-0.5.10/src/Makefile.mic to point to your preferred location.

V. Configure bwa aln for running
============================

To test this code, we recommend starting with publically available data:

•	Get the reference genome from ftp://ftp.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/
•	Download and unpack the following read pair:
	o	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12044/sequence_read/SRR766060_1.filt.fastq.gz
	o	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12044/sequence_read/SRR766060_2.filt.fastq.gz

Now prepare for testing by modifying the included scripts to match your test environment:
	
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
	host2_ratio=7				# input data file ratio for mic0

	host3=${host1}-mic1			# host3: mic1 of host1
	host3_exe=bwa_mic			# executable file name on host3
	host3_num_tasks=3			# number of tasks on host3
	host3_num_threads=80		# number of threads for each task
	host3_ratio=7				# input data file ratio for mic1
	
	src_directory=../src							# source code

	data_directory=../data/public_workload			# input data directory
	ref_file=${data_directory}/hs37d5.fa			# reference data file
	file1=${data_directory}/SRR766060_1.filt.fastq	# input data file1
	file2=${data_directory}/SRR766060_2.filt.fastq	# input data file2

	res_tmp_file1=tmp_aln_file1						# temporary output file name of each task
	res_tmp_file2=tmp_aln_file2

	res_file1=aln_1.sai								# ALN output file name
	res_file2=aln_2.sai

	# for bwa sampe output and picard output
	sampe_output=output.sam
	picard_output=sorted.sam
	
Note that: 

(1) Intel Xeon Phi is considered as a independent host

(2) num_tasks * num_threads must <= number of logical cores on this host

(3) Input data file is splited ino seral files according to the ratio for each host. The ratio of host to host-mic0 (equals host1_ratio to host2_ratio) is set according to the speed of their using the speed that each host independently processes the input data files.

(4) Set aln output file name by res_file1 and res_file2, needn’t set by -f output-file-name

VI. Build bwa reference index file for bwa aln
============================

	> ./build-index

Which runs bwa index -a bwtsw $ref_file to construct the FM-index for the reference genome.

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

Which runs bwa aln native on Intel Xeon.

IX. Run bwa aln on Intel Xeon and Intel Xeon Phi in symmetric model
============================

1. Modify host1_ratio and host2_ratio in aln_head according to the performance on Intel Xeon and Intel Xeon Phi

2. Split input data file for multiple tasks

	> ./aln-cpu-and-1mic-split
	
	Split input file for Intel Xeon and 1 Intel Xeon Phi card
	
	Or:

	> ./aln-cpu-and-1mic-split
	
	Split input file for Intel Xeon and 2 Intel Xeon Phi cards	

3. Run

	> ./aln-cpu-and-1mic-run
	
	Which runs bwa aln on both Intel Xeon and 1 Intel Xeon Phi card in symmetric model.

	Or:

	> ./aln-cpu-and-2mic-run

	Which runs bwa aln on both Intel Xeon and 2 Intel Xeon Phi cards in symmetric model.

Contributors:
============================

You, Liang <liang.you@intel.com>

Congdon, Charles <charles.w.congdon@intel.com>

Ramanujam, Ram <ram.ramanujam@intel.com>






