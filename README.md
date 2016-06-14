# ncsiFCM
Notice:
This is MetaBin2.0 version.
The max read length is 1e4bp. The largest number of reads are 2^31. The input file must be FASTA format. Please read Requirement and How to use this software before running.

 Thanks for reporting bugs and unexpected output.

Requriement:

This software is suitable for all unix-like 64-bit system with gcc installed.

How to use this software:
1. make (then a executive file titled ncsiFCM will be generated)
2. ./ncsiFCM inputfile number_of_bins
	for example: ./ncsiFCM example.fna 5

Output:
result.txt file with the cluster label of each DNA fragment, from Cluster 0 to Cluster k.

Contact us:
laoniu313@qq.com
