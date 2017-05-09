Hetero2 : A program to simulate heterogeneous multiple alignment sequences

==========================
To compile the programs
==========================

$ tar -zxvf Hetero-2.2.tar.gz
$ cd Hetero-2.2
$ make

An executable files called "Hetero2" will appear

=====================
To run the programs 
=====================

Syntax:
  ./Hetero2 <tree file> <site info file> <param file list> <other options>
  ./Hetero2 -h

  <tree file>           : "Tree file" lists the tree of each site category, the
                          edge lengths and the labels of terminal/internal nodes

  <site info file>      : "Site info file" lists the detailed information of
                          each site category, including the site proportion and
                          the nucleotide distribution at the root

  <param file list>     : "Param file list" shows the name of the parameter file
                          of each variant site category

other options:

  -l <sequence length>  : The length of sequences to be simulated
                          (default: 10,000)

  -f <output format>    : The format of simulated multiple sequence alignment
                          1 - FASTA format
                          2 - Sequential PHYLIP format (default)

  -o <output prefix>    : Prefix for output files
                          (default: <tree file> w/o .ext)

  -h                    : The help page

Output files:
  <output prefix>.out   : The simulated multiple sequence alignment file


=====================
Example files
=====================

The following example files are available for reference:

1. trees.txt            : An example of "tree file" which lists the tree of
                          each site category, the edge lengths and the labels
                          of terminal/internal nodes

2. site_info_file.txt   : An example of "Site info file" lists the detailed
                          information of each site category, including the
                          site proportion and the nucleotide distribution at
                          the root

3. param_file_list.txt  : An example of "parameter list file" showing the name
                          of the parameter file of each variant site category

4. parameter_1/2.txt    : An example of parameter file showing the detailed
                          parameters of the rate matrix of each edge leading
                          to the corresponding node

To try the example, you may use the following command:
$./Hetero2 trees.txt site_info_file.txt param_file_list.txt

The simulated multiple sequence alignment would be in the file: "trees.out".
