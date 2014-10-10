
find_pssm_dna is a simple tool for comparing PSSMs against DNA sequences.

1. Compiling

On linux/unix system with gcc installed running "make" should be sufficient. Note find_pssm_dna uses GNU getopt_long to parse command line arguments. Executable named "find_pssm_dna" will be created in the same directory.

2. Running

The command to run the program is of form

"find_pssm_dna [options] p-value sequence matrix1 matrix2 ..."

where p-value is given as decimal, sequence is path to the sequence file and matrix1, matrix2 and so on are paths to weight matrix files describing matrices to be matched against sequence. There can be an arbitrary amount of matrix files. If matrix file path contains wildcards, then all matching files are used.

Symbols other than [acgtnACGTN] are ignored in the sequence, including SNPs. FASTA header lines are ignored entirely. Symbols 'n' and 'N' are randomly interpreted to mean one of 'A', 'C', 'G' or 'T'.

Options are: 


 -h, --help
   Prints help.
 -f, --flat_bg
   Assumes that the sequence has a flat background distribution, instead of calculating it from the sequence.
 -a[number], --algorithm=[number]
   Selects the algorithm used for matching. Algorithms are:
      0: Naive algorithm
      1: Naive superalphabet algorithm
      2: Permutated lookahead search
      3: Lookahead filtration
      4: Multiple matrix lookahead filtration
 -p[number], --parameter=[number]
   Tuning parameter for scanning algorithm. Defines window width for ACFA and LFA, and the size of q-gram for NS.

3. References

Cinzia Pizzi, Pasi Rastas, Esko Ukkonen. Fast search algorithms for position specific scoring matrices. First International Conference on Bioinformatics Research and Development - BIRD 2007. 
   
