
# STRIKE Contact Matrix Generator

## Introduction

This utility is part of the program [STRIKE](http://www.tcoffee.org/Projects/strike/) (Single sTRucture Induced Evaluation), a program to evaluate protein multiple sequence alignments using a single protein structure. 

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

## Example 1: Generate the contact file of a protein sequence with PDB ID: 1aab 

To generate the contact file of the protein sequence 1aab, we have to download the PDB file from the [Protein Data Bank](http://www.rcsb.org), with the PDB ID 1aab. 
The PDB file of the sequence (1aab_.pdb) is saved in the directory example and the output directory is the current path of this application.
The command is following:

```
@: ./bin/strike_pdbcontactsgenerator 1aab example/1aab_.pdb A ./
```
The 1aab.contacts is created in the current directory.

## Example 2: Generate the contact file of the protein sequence of Grenache (Vitis vinifera) Polyphenol Oxidase

The PDB File is saved in the example directory (2p3x.pdb). The PDB ID is 2P3X. 
The command is:

```
@: ./bin/strike_pdbcontactsgenerator 2p3x example/2p3x.pdb A ./
```

The 2p3x.contacts is created in the current directory.
