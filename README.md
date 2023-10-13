# CGRWDL
CGRWDL is a new alignment-free method to construct the phylogenetic tree. We combine the dynamical language (DL) model and CGR to obtain new sequence information by considering both the frequency and position of k-mers in the sequence.The program is written in python(in Linux).
# Usage
python cgrwdl_main.py protein --files $file_name.fasta --k 3
# Parameters
Sequence type:dna or protein
--files sequences file of species (input file,fasta file format only)
--k  length of k-mer
--savefolder position of saviong the results,default='./output/ '
# Output
The output is the mega format, and can be imported directly in to MEGA X to obtain a phylogenetic tree.
# Example
python cgrwdl_main.py protein --files example.fasta --k 3
# Citation
Ting Wang, Zuguo Yu, Jinyan Li, CGRWDL : Alignment-free phylogeny reconstruction method for viruses based on chaos game representation weighted by dynamical language model
