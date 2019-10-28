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
      -t, --threads:     number of processing threads (default: 1)
      -v, --verbose:     more output
      -h, --help:        print this
    Example usage:
      pam make -f ref.fasta -k 21 -r read_R1.fastq.gz -r read_R2.fastq.gz -o out.pam
