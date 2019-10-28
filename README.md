Generate a k-mer-based presence/absence matrix
==============================================

This is not, at present, a general-purpose tool, and will not work with data other than that it was built on.
If you don't believe me, just try it...

To build on an Ubuntu-like system:

    sudo apt install zlib1g-dev
    git clone http://github.com/jwanglab/pam
    cd pam
    git clone http://github.com/attractivechaos/klib
    make
    ./pam

Usage:

    pam <command> [options]
    
    Commands:
      make:      build a presence/absence matrix from a read set
      collapse:  collapse several presence/absence matrices into a pattern set
      count:     compute basic stats from stored presence/absence matrix
      triangle:  compute a Jaccard distance matrix among a set of PAMs
    Options:
      make -f <reference.fa> -r <reads.fq> -k <21>
        -f, --ref:       A reference FASTA/Q[.gz] file
        -r, --reads:     A FASTA/Q[.gz] file with reads
        -c, --coverage:  k-mer coverage to be considered 'present' (default: 5)
        -k:              k-mer length (default: 21)
      collapse <sample.pam> <sample.pam> [...]
        -m, --matrix:    matrix produced by 'pam make'
        -p, --pthrshld:  pattern threshold as a fraction of |matrix| (default: 0.00001)
        -n, --num_kmers: number of k-mers to sample, randomly (default: all)
        --random: probabilistically assign present/absent based on coverage
      count -c <5> -m <sample.pam> [-m <sample.pam> ...]
        -m, --matrix:    matrix produced by 'pam make'
      triangle -c <5> <sample.pam> <sample.pam> [...]
      -o, --output:      output file for whatever
      -t, --threads:     number of processing threads (default: 1) [MULTIPROCESSING IS NOT CURRENTLY SUPPORTED]
      -v, --verbose:     more output
      -h, --help:        print this
    Example usage:
      pam make -f ref.fasta -k 21 -r read_R1.fastq.gz -r read_R2.fastq.gz -o out.pam

Workflow:

1. Generate a reasonable set of reference k-mers for which to count presence/absence. In many cases, this can be a representative reference genome for the clade or taxa you wish the compare.
2. Call 'make' on each individual/line/taxa that will be in your matrix (eventual tree), with a genome, assembly, or sequencing reads.
3. Use 'collapse' on all of these PAMs to merge them, filter, subsample, and generate a matrix in NEXUS format

Alternatively, one can build a matrix with approximate Jaccard distances (a la MASH) from a set of PAMs using 'triangle'.

By default, 'collapse' generates all presence/absence patterns (analogous to SDPs) represented by >0.001% of all sampled k-mers with a minimum of 5 occurrences of each k-mer. Note that, if you want to generate a pam from an assembled genome, you MUST set -c to 0 (or 1), and for sequencing reads, -c should be proportional to your total coverage in such a way that it largely eliminates sequencing errors but is << median read coverage. Other options, including --num\_kmers and --random allow a random subset of k-mers to be considered and/or presence/absence assigned probabilistically based on the k-mer count relative to the expected coverage.
