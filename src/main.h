


/*! \file main.h
\brief The main contains the main program and the parameter parsing.
*/


/*! \mainpage STRIKE v. 1.2



\section intro_sec Introduction

This program allows the computation of the contact score given an alignment and a <a href="http://www.rcsb.org">PDB</a>   structure.

\section install_sec Installation

The strike program can be compiled with:

<pre>
@: make mode=release
</pre>


If you have OpenMP on you system you can enable it by compiling with:
<pre>
@: make mode=release parallel=y
</pre>


This enables a new option: '--nc' which allows you to specify the number of cores you want to use.


\section manual Using RankAli




\subsection step1 Scoring an alignment

To score an alignment two files have to be given to RankAli: The alignment file and a connection file. The connection file consists of three columns. The first column contains the sequence name in the alignment, the second name the complete path to the PDB file which is associated with the sequence in the first column and the third column the chain to use.


The alignment has to be in FASTA, ClustalW, or MSF format. They are recognized by the ending of the files: \n
FASTA:    .fa/.fasta \n
ClustalW: .aln \n
MSF:      .msf \n

To Score an alignment the files has to be given like this:

<pre>
@: strike -a \<alignment file> -c \<connection file>
</pre>


The output will be a value for each given PDB/Sequence connection and the average value. This output will be printed to standard out. If you want to put it into a file you can use the following parameter set:

<pre>
@: strike -a \<alignment file> -c \<connection file> -o \<output file>
</pre>


The connection file is similar to the T-Coffee template file. It consists of the sequence name, the connection to a PDB file and a third column specifying a chain. In case the chain is not specified the first chain will be used. Example line:

my_seq1 ./PDB.file A


\subsection step2 Normalized score

To calculate the normalized score add the parameter '-n' to the strike call.

<pre>
@: strike -a \<alignment file> -c \<connection file> -o \<output file> -n
</pre>


\subsection formats Input formats

Aligments have to be in either MSF, FASTA or CLUSTALW format.

The connectionfiles conists of at least one line in the following format:\n
\<sequence name> _P_ \<pdb file> \<chain>

The chain is optional. If no chain is given the first chain found in the PDB file will be used.

\subsection example Example run


To run an example change into the example directory and run the following command:

<pre>
../bin/strike -c AAA.con -a AAA.msf
</pre>

*/



