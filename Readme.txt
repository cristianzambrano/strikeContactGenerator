
STRIKE v. 1.2


1 Introduction

This program allows the computation of the contact score given an alignment and a PDB (http://www.rcsb.org) structure.
Please write an email to carsten.kemena@crg.es in case you have any questions or problems concerning the program.

2. Installation

The strike program can be compiled with:

@: make


If you have OpenMP on you system you can enable it by compiling with:
<pre>
@: make parallel=y
</pre>


This enables a new option: '--nc' which allows you to specify the number of cores you want to use. This is only used for the computation of contacts which of a long chain chan take some time because each ATOM has to be compared with each other.


3. Using Strike
3.1. Scoring an alignment

To score an alignment two files have to be given to RankAli: The alignment file and a connection file. The connection file consists of three columns. The first column contains the sequence name in the alignment, the second name the complete path to the PDB file which is associated with the sequence in the first column and the third column the chain to use.


The alignment has to be in FASTA, ClustalW, or MSF format. They are recognized automatically.

To Score an alignment the files has to be given like this:


@: strike -a \<alignment file> -c \<connection file>



The output will be a value for each given PDB/Sequence connection and the average value. This output will be printed to standard out. If you want to put it into a file you can use the following parameter set:


@: strike -a \<alignment file> -c \<connection file> -o \<output file>



The connection file is similar to the T-Coffee template file. It consists of the sequence name, the connection to a PDB file and a third column specifying a chain. In case the chain is not specified the first chain will be used. Example line:

my_seq1 ./PDB.file A


3.2 Normalized score

To calculate the normalized score add the parameter '-n' to the strike call.


@: strike -a \<alignment file> -c \<connection file> -o \<output file> -n



3.3 Input formats

Aligments have to be in either MSF, FASTA or CLUSTALW format.

The connectionfiles conists of at least one line in the following format:
\<sequence name> _P_ \<pdb file> \<chain>

The chain is optional. If no chain is given the first chain found in the PDB file will be used.

3.4 Example run


To run an example change into the example directory and run the following command:

<pre>
../bin/strike -c AAA.con -a AAA.msf
</pre>
