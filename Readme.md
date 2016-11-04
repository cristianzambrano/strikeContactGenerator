
#STRIKE Contacts_File Generator

##Introduction

This utilitary is part of the program STRIKE(http://www.tcoffee.org/Projects/strike/) (Single sTRucture Induced Evaluation) a program to evaluate protein multiple sequence alignments using a single protein structure , and has been developed to support the funcionality of the  jMetalMSA(https://github.com/jMetal/jMetalMSA) a Parallel software tool for Multiple Sequence Alignment with multi-objective metaheuristics .

STRIKE Contacts_File Generator allows the computation of the STRIKE Contact matrix of a sequence using the PDB structural information as a source for amino acid frequencies and contacts. A text file is created with the matrix, called as the sequence with the prefix .contacts
This contact file is requeried to the estimate the STRIKE score given an alignment that includes the sequence.

## Installation

The strike program can be compiled with:
```
@: make
```

## Using Contacts_File Generator

To generate the Contacts file of a given sequence, four parameters have to be given to the program: the sequence name, the path to the PDB file which is associated with the sequence, can be downloaded form Protein Data Bank(http://www.rcsb.org), the chain to use and the output directory where contacts file is saved. 

To run an example change into the program directory and run the following command:
```
@: ./bin/strikeContactGenerator SequenceName PDBPathFile Chain OutputDirectory
```



