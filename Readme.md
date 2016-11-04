
#STRIKE Contact Matrix Generator

##Introduction

This utilitary is part of the program [STRIKE](http://www.tcoffee.org/Projects/strike/) (Single sTRucture Induced Evaluation) a program to evaluate protein multiple sequence alignments using a single protein structure , and has been developed to support the funcionality of the  [jMetalMSA](https://github.com/jMetal/jMetalMSA) a Parallel software tool for Multiple Sequence Alignment with multi-objective metaheuristics .

STRIKE Contact Matrix Generator allows the computation of the STRIKE Contact matrix of a given sequence using its PDB structural information as a source for amino acid frequencies and contacts. A text file is created, called as the sequence with the prefix `.contacts`. This contact file is requeried to estimate the STRIKE score given an alignment that includes the sequence.

## Installation

The STRIKE Contact Matrix Generator program can be compiled with:

```
@: make
```

## Using STRIKE Contact Matrix Generator

To generate the Contacts file of a given sequence, four parameters have to be given to the program: 
* The sequence name
* The path to the PDB file which is associated with the sequence, can be downloaded from the [Protein Data Bank](http://www.rcsb.org), 
* The chain to use
* The output directory where contacts file is saved. 

To run an example change into the program directory and run the following command:
```
@: ./bin/strike_pdbcontactsgenerator SequenceName PDBPathFile Chain OutputDirectory
```



