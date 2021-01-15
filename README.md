QCluster: Extending Alignment-free Measures with Quality Values for Reads Clustering

Introduction:

The data volume generated by next-generation sequencing (NGS) technologies is growing at a pace that is now challenging the storage and data processing capacities of modern computer systems. In this context an important aspect is the reduction of data complexity by collapsing redundant reads in a single cluster to improve the run time, memory requirements, and quality of post-processing steps like assembly. Several alignment-free measures, based on k-mers counts, have been used to cluster reads. Quality scores produced by NGS platforms are fundamental for various analysis of NGS data: reads mapping; error detection; etc. Moreover future-generation sequencing reads promise to contain a large number of erroneous bases (up to 10$\%$). Thus it will be fundamental to exploit quality value information within the alignment-free framework. We have introduced a family of alignment-free measures, called $D^q$-type, that incorporate quality value information and k-mers counts for the comparison of reads data. A set of experiments on simulated and real reads data confirm that the new measures are superior to other classical alignment-free statistics, especially when erroneous reads are considered.

Software:

Here you can find the C++ application QCluster with an example.

The program must be compiled with the command 'make'. QCluster has been developed under the Linux environment, so in other systems you may need to install other software (e.g. CYGWIN).

You can run the program using the following command:

QCluster -d dist_type -c num_cluster -k kmer_length file_fastq

Where dist_type is the distance used for the clustering, num_cluster is the number of partitions in which you want to divide the set of sequences and kmer_length is the length of the words used by the algorithm. The input file (file_fastq) must be in fastq format, as it is defined in Illumina 1.8+ (Phred+33). The program returns on standard output the list of sequences with the number of the corresponding cluster. With the '-w' option the program also outputs the clusters in different files. Use the option '-h' to get the complete list of options.

To run the included examples, type:

QCluster -d a -c 3 -k 3 -t 5 -S 0 -w Example/sequences.fastq

The output will be the same of the content of the directory /Example

Licence

The software is freely available for academic use.

For questions about the tool, please contact Matteo Comin.

Reference

Please cite the following paper:

M.Comin, A. Leoni, M. Schimd
"QCluster: Extending Alignment-free Measures with Quality Values for Reads Clustering",
Proceedings of the 14th Workshop on Algorithms in Bioinformatics (WABI) 2014.
Lecture Notes in BIoinformatics (LNBI) 2014, 8701, pp. 1-13. 
