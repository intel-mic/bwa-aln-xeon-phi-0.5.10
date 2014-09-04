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

1. To test this code, we recommend starting with publically available data:

	Get the reference genome from ftp://ftp.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/
	
	Download and unpack the following read pair:

	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12044/sequence_read/SRR766060_1.filt.fastq.gz
	
	ftp://ftp.ncbi.nih.gov/1000genomes/ftp/data/NA12044/sequence_read/SRR766060_2.filt.fastq.gz

2. Make a directory public_workload in bwa-aln-xeon-phi-0.5.10/data and put all the download files to bwa-aln-xeon-phi-0.5.10/data/public_workload.

3. Prepare for testing by modifying the included scripts to match your test environment:
	
	> vi aln-head
	
Modify configuration:

	host1=crt03					# hostname of host1
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

(1) The Intel® Xeon Phi™ coprocessor is considered a host where MPI is concerned by this implementation.

(2) num_tasks * num_threads must be less than or equal to the number of logical cores on each host

(3) Using the provided scripts, the input data files are split into several intermediate files according to the ratio for each host. The ratio of host to host-mic0 (equals host1_ratio to host2_ratio) should be set using the speed that each host independently processes the input data files (which you may need to determine empirically).   For example, if the Intel® Xeon Phi™ coprocessor’s performance is 0.7~0.8x that of the same problem running on the Intel® Xeon® processor (with optimized code), then the ratio is 10 for the Intel® Xeon® processor, and 7 for the Intel® Xeon Phi™ coprocessor.

(4) Set the ALN output file names using the res_file1 and res_file2 variables.

VI. Build the bwa reference index file for bwa aln
============================

	> ./build-index

Which runs “bwa index -a bwtsw $ref_file” to construct the FM-index for the reference genome.

NOTE:  This script runs a version of BWA called “bwa-ori”.  This is the binary you get when you download BWA 0.5.10 from http://sourceforge.net/projects/bio-bwa/files/ and build it with no modifications (beyond renaming the resulting executable to “bwa-ori”).

VII. Run bwa aln on the Intel® Xeon® processor
============================

1. Split the input data file for processing by multiple tasks

	> ./aln-cpu-split

2. Run 

	> ./aln-cpu-run

Which runs bwa aln on the Intel® Xeon® processor.

VIII. Run bwa aln on Intel® Xeon Phi™ coprocessor
============================

1.	Split the input data file for processing by multiple tasks

	> ./aln-1mic-split

2.	Run 

	> ./aln-1mic-run

Which runs bwa aln natively on the Intel® Xeon Phi™ coprocessor.

IX. Run bwa aln on both the Intel® Xeon® processor and Intel® Xeon Phi™ coprocessor in symmetric mode
============================

1. Modify host1_ratio and host2_ratio in aln_head according to the performance seen during the standalone Intel® Xeon® processor and Intel® Xeon Phi™ coprocessor runs above

2. Split the input data file for processing by multiple tasks depending on how many Intel® Xeon Phi™ coprocessors you will use.   

	For an Intel® Xeon® processor and 1 Intel® Xeon Phi™ coprocessor, do this:

	> ./aln-cpu-and-1mic-split
	
	For an Intel® Xeon® processor and 2 Intel® Xeon Phi™ coprocessors, do this:

	> ./aln-cpu-and-2mic-split

3.	Run on the desired number of Intel® Xeon Phi™ coprocessors as well as the Intel® Xeon® host.   

	To run bwa aln on both the Intel® Xeon® processor and 1 Intel® Xeon Phi™ coprocessor at once, do this:

	> ./aln-cpu-and-1mic-run
	
	To run bwa aln on both the Intel® Xeon® processor and 2 Intel® Xeon Phi™ coprocessors at once, do this:
	
	> ./aln-cpu-and-2mic-run

X. Performance gain
============================

For the public workload we described above, the following graph shows the speedup achieved from bwa-aln-xeon-phi-0.5.10. As you can see, we get:

•	Up to a 1.2x speedup can be achieved when processing both sets of reads using this code vs. the original code (bwa-0.5.10 aln, 48 threads) running on a 2-Socket Intel® Xeon® Processor E5-2697 v2. 

•	Up to a 1.76x speedup can be achieved when processing both sets of reads using this code on one Intel® Xeon Phi™ Coprocessor + Intel® Xeon® Processor E5-2697 v2 vs. the original code (bwa-0.5.10 aln, 48 threads) running on a 2-Socket Intel® Xeon® Processor E5-2697 v2.

•	Up to 2.76x speedup can be achieved when processing both sets of reads using this code on two Intel® Xeon Phi™ coprocessors + Intel® Xeon® Processor E5-2697 v2 vs. original code(bwa-0.5.10 aln, 48 threads) running on a 2-Socket Intel® Xeon® Processor E5-2697 v2.

Testing platform configuration:

•	Server with Intel® Xeon® Processors E5-2697 v2: Two sockets, 12-core, 2.7 GHz, 64GB DDR3-1600, 8.0 GT/s, OS version: Red Hat Enterprise Linux Server release 6.5, Intel® Turbo Boost Technology enabled, Intel® Hyper-Threading Technology (Intel® HT Technology) enabled

•	Intel® Xeon Phi™ Coprocessor: 7120A, 61-core 1.238 GHz, 16GB GDDR5-5500,  5.5 GT/s, MPSS 2.1.6720-21, Flash version 2.1.03.0386, uOS version : 2.6.38.8-gefd324e, ECC enabled, Intel® Turbo Boost Technology disabled

Legal Disclaimers
============================

For more information go to http://www.intel.com/performance

Software and workloads used in performance tests may have been optimized for performance only on Intel microprocessors. Performance tests, are measured using specific computer systems, components, software, operations and functions. Any change to any of those factors may cause the results to vary. You should consult other information and performance tests to assist you in fully evaluating your contemplated purchases, including the performance of that product when combined with other products.

Optimization Notice

Intel's compilers may or may not optimize to the same degree for non-Intel microprocessors for optimizations that are not unique to Intel microprocessors. These optimizations include SSE2, SSE3, and SSE3 instruction sets and other optimizations. Intel does not guarantee the availability, functionality, or effectiveness of any optimization on microprocessors not manufactured by Intel. Microprocessor-dependent optimizations in this product are intended for use with Intel microprocessors. Certain optimizations not specific to Intel microarchitecture are reserved for Intel microprocessors. Please refer to the applicable product User and Reference Guides for more information regarding the specific instruction sets covered by this notice.

Notice revision #20110804

By using this document, in addition to any agreements you have with Intel, you accept the terms set forth below.
You may not use or facilitate the use of this document in connection with any infringement or other legal analysis concerning Intel products described herein.  You agree to grant Intel a non-exclusive, royalty-free license to any patent claim thereafter drafted which includes subject matter disclosed herein.

INFORMATION IN THIS DOCUMENT IS PROVIDED IN CONNECTION WITH INTEL PRODUCTS. NO LICENSE, EXPRESS OR IMPLIED, BY ESTOPPEL OR OTHERWISE, TO ANY INTELLECTUAL PROPERTY RIGHTS IS GRANTED BY THIS DOCUMENT. EXCEPT AS PROVIDED IN INTEL'S TERMS AND CONDITIONS OF SALE FOR SUCH PRODUCTS, INTEL ASSUMES NO LIABILITY WHATSOEVER AND INTEL DISCLAIMS ANY EXPRESS OR IMPLIED WARRANTY, RELATING TO SALE AND/OR USE OF INTEL PRODUCTS INCLUDING LIABILITY OR WARRANTIES RELATING TO FITNESS FOR A PARTICULAR PURPOSE, MERCHANTABILITY, OR INFRINGEMENT OF ANY PATENT, COPYRIGHT OR OTHER INTELLECTUAL PROPERTY RIGHT. 

A "Mission Critical Application" is any application in which failure of the Intel Product could result, directly or indirectly, in personal injury or death. SHOULD YOU PURCHASE OR USE INTEL'S PRODUCTS FOR ANY SUCH MISSION CRITICAL APPLICATION, YOU SHALL INDEMNIFY AND HOLD INTEL AND ITS SUBSIDIARIES, SUBCONTRACTORS AND AFFILIATES, AND THE DIRECTORS, OFFICERS, AND EMPLOYEES OF EACH, HARMLESS AGAINST ALL CLAIMS COSTS, DAMAGES, AND EXPENSES AND REASONABLE ATTORNEYS' FEES ARISING OUT OF, DIRECTLY OR INDIRECTLY, ANY CLAIM OF PRODUCT LIABILITY, PERSONAL INJURY, OR DEATH ARISING IN ANY WAY OUT OF SUCH MISSION CRITICAL APPLICATION, WHETHER OR NOT INTEL OR ITS SUBCONTRACTOR WAS NEGLIGENT IN THE DESIGN, MANUFACTURE, OR WARNING OF THE INTEL PRODUCT OR ANY OF ITS PARTS. 

Intel may make changes to specifications and product descriptions at any time, without notice. Designers must not rely on the absence or characteristics of any features or instructions marked "reserved" or "undefined". Intel reserves these for future definition and shall have no responsibility whatsoever for conflicts or incompatibilities arising from future changes to them. The information here is subject to change without notice. Do not finalize a design with this information. 

The products described in this document may contain design defects or errors known as errata which may cause the product to deviate from published specifications. Current characterized errata are available on request. 

Contact your local Intel sales office or your distributor to obtain the latest specifications and before placing your product order. 

Copies of documents which have an order number and are referenced in this document, or other Intel literature, may be obtained by calling 1-800-548-4725, or go to: http://www.intel.com/design/literature.htm

Other names and brands may be claimed as the property of others.
Copyright © 2014 Intel Corporation. All rights reserved









