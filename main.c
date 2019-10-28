#include <zlib.h>
#include <math.h>
#include <getopt.h>
#include "klib/kstring.h"
#include "klib/khash.h"
#include "klib/kvec.h"
#include "klib/ksort.h"

#ifndef _kseq_
#define _kseq_

#include "klib/kseq.h"

// init kseq struct
KSEQ_INIT(gzFile, gzread)

#endif

// creates str(k-mer):u32 hash
KHASH_MAP_INIT_STR(kmerHash, uint32_t);

// creates int(p/a pattern):u32 hash
KHASH_MAP_INIT_INT64(patternHash, uint32_t);

void usage() {
  printf("Usage: pam <command> [options]\n");
  printf("Commands:\n");
  printf("  make:      build a presence/absence matrix from a read set\n");
  printf("  collapse:  collapse several presence/absence matrices into a pattern set\n");
  printf("  count:     compute basic stats from stored presence/absence matrix\n");
  printf("  triangle:  compute a Jaccard distance matrix among a set of PAMs\n");
  printf("Options:\n");
  printf("  make -f <reference.fa> -r <reads.fq> -k <21>\n");
  printf("    -f, --ref:       A reference FASTA/Q[.gz] file\n");
  printf("    -r, --reads:     A FASTA/Q[.gz] file with reads\n");
  printf("    -c, --coverage:  k-mer coverage to be considered 'present' (default: 5)\n");
  printf("    -k:              k-mer length (default: 21)\n");
  printf("  collapse <sample.pam> <sample.pam> [...]\n");
  printf("    -m, --matrix:    matrix produced by 'pam make'\n");
  printf("    -p, --pthrshld:  pattern threshold as a fraction of |matrix| (default: 0.00001)\n");
  printf("    -n, --num_kmers: number of k-mers to sample, randomly (default: all)\n");
  printf("    --random: probabilistically assign present/absent based on coverage\n");
  printf("  count -c <5> -m <sample.pam> [-m <sample.pam> ...]\n");
  printf("    -m, --matrix:    matrix produced by 'pam make'\n");
  printf("  triangle -c <5> <sample.pam> <sample.pam> [...]\n");
  printf("  -o, --output:      output file for whatever\n");
  printf("  -t, --threads:     number of processing threads (default: 1)\n");
  printf("  -v, --verbose:     more output\n");
  printf("  -h, --help:        print this\n");
  printf("Example usage:\n");
  printf("  pam make -f ref.fasta -k 21 -r read_R1.fastq.gz -r read_R2.fastq.gz -o out.pam\n");
}

typedef struct {
  uint32_t ct;
  uint64_t pattern;
} pa_ct;

#define pa_lt(a, b) ((a).ct > (b).ct)
KSORT_INIT(pa_ct_decreasing, pa_ct, pa_lt)

static struct option long_options[] = {
// if these are the same as a single-character option, put that character in the 4th field instead of 0
  { "ref",         required_argument, 0, 'f' },
  { "reads",       required_argument, 0, 'r' },
  { "coverage",    required_argument, 0, 'c' },
  { "matrix",      required_argument, 0, 'm' },
  { "output",      required_argument, 0, 'o' },
  { "threads",     required_argument, 0, 't' },
  { "pthrshld",    required_argument, 0, 'p' },
  { "num_kmers",   required_argument, 0, 'n' },
  { "random",      no_argument,       0, 0 },
  { "help",        no_argument,       0, 'h' },
  { "verbose",     no_argument,       0, 'v' },
  { 0, 0, 0, 0}
};

char* tobitstring(uint64_t i) {
  char* s = malloc(65 * sizeof(char));
  s[64] = '\0';
  int j;
  for(j = 0; j < 64; j++) {
    if((i>>(63-j))&1)
      s[j] = '1';
    else
      s[j] = '0';
  }
  return s;
}

// requires a null-terminated string
// the first 2 bytes are the start position, the second 2 are the length
uint32_t basename(char* f) {
  int i;
  uint32_t pos = 0;
  uint16_t len = 0;
  for(i = 0;; i++) {
    if(f[i] == '\0') {
      return (pos<<16) + len;
    } else if(f[i] == '/') {
      pos = i+1;
      len = 0;
    } else if(f[i] != '.' && pos+len == i) {
      len++;
    }
  }
}

float get_covg(char* name) {
  if(strcmp(name, "albom_Ishigaki") == 0)
    return 25.8;
  if(strcmp(name, "albost_Borneo") == 0)
    return 23;
  if(strcmp(name, "albost_Cambodia") == 0)
    return 25.6;
  if(strcmp(name, "albost_Indonesia") == 0)
    return 16.8;
  if(strcmp(name, "albost_Luzon") == 0)
    return 30.8;
  if(strcmp(name, "albostrigata_SriLanka") == 0)
    return 36.2;
  if(strcmp(name, "bilim_Guam") == 0)
    return 25.6;
  if(strcmp(name, "bilim_Oahu") == 0)
    return 38;
  if(strcmp(name, "hypo_Guam") == 0)
    return 18;
  if(strcmp(name, "hypo_Luzon") == 0)
    return 20.4;
  if(strcmp(name, "immigrans") == 0)
    return 30.4;
  if(strcmp(name, "kep_Brunei") == 0)
    return 26.2;
  if(strcmp(name, "kep_Sarawak") == 0)
    return 27.4;
  if(strcmp(name, "koh_Phillipines") == 0)
    return 23.2;
  if(strcmp(name, "koh_Sarawak") == 0)
    return 28;
  if(strcmp(name, "nas_Mombasa") == 0)
    return 27.6;
  if(strcmp(name, "nas_Mysore") == 0)
    return 20.2;
  if(strcmp(name, "neohypo_NGuinea") == 0)
    return 23.8;
  if(strcmp(name, "neon_Mysore") == 0)
    return 27.4;
  if(strcmp(name, "nivei") == 0)
    return 22.4;
  if(strcmp(name, "pallidi") == 0)
    return 13.2;
  if(strcmp(name, "pula_Sarawak") == 0)
    return 14.2;
  if(strcmp(name, "siam_Cambodia") == 0)
    return 30.6;
  if(strcmp(name, "sulf_NGuinea") == 0)
    return 25.2;
  if(strcmp(name, "sulf_NIreland") == 0)
    return 27;
  if(strcmp(name, "taxonF") == 0)
    return 28.2;
  if(strcmp(name, "taxonG") == 0)
    return 31.8;
  if(strcmp(name, "taxonJ") == 0)
    return 20.8;
  fprintf(stderr, "ERROR: no coverage hardcoded for '%s', all present k-mers will be used\n       -- get Jeremy to implement a real covg estimator!\n", name);
  return 0;
}

// check if the pattern at index <i> in <a> is 4-gamete compatible with all those indices in <used>
// m: number of samples (number of bits to evaluate)
int pattern_ok(pa_ct* a, uint8_t* used, uint64_t n, uint32_t i, uint32_t m) {
  uint64_t mask = (1<<m)-1;
  uint64_t p = a[i].pattern; // our target pattern - mask out upper bits (count)
  uint64_t q; // the other pattern
  uint32_t j; // just a counter
  for(j = 0; j < n; j++) {
    if(used[j]) {
      q = a[j].pattern;
      /*
      fprintf(stderr, "comparing %s to %s\n", tobitstring(p), tobitstring(q));
      fprintf(stderr, "p & q == %s\n", tobitstring(p&q));
      fprintf(stderr, "~p & ~q & mask == %s\n", tobitstring(~p&~q&mask));
      fprintf(stderr, "p ^ q & p == %s\n", tobitstring(((p ^ q) & p)));
      fprintf(stderr, "p ^ q & q == %s\n", tobitstring(((p ^ q) & q)));
      */
      if(((p & q) > 0) + ((~p & ~q & mask) > 0) + (((p ^ q) & p) > 0) + (((p ^ q) & q) > 0) == 4) { // fails the 4-gamete test
        return 0;
      }
    }
  }
  return 1;
}

char* get_newick(pa_ct* a, uint8_t* used, uint64_t n) {
  uint32_t i;
  int j = 0; // number of characters printed
  char* newick = malloc(1000 * sizeof(char));
  for(i = 0; i < n; i++) {
    if(used[i]) {
      j = j + sprintf(newick+j, "%u: %s\n", a[i].ct, tobitstring(a[i].pattern));
    }
  }
  return newick;
}

// branch and bound search of k-mer patterns - is guaranteed to find tree with highest total k-mer support
// total: total k-mer occurrences, to help bounding
int bb(pa_ct *a, uint32_t n, uint32_t total, char** top_n, uint32_t *scores, int keep_top_n, uint32_t n_samples) {
  int j;
  uint32_t i = 0; // index into a
  uint32_t used_occ = 0, unused_occ = 0; //used/unused (from the top, to this point) k-mer occurrences, to help bounding
  uint32_t n_used = 0;
  uint8_t *used = calloc(n, sizeof(uint8_t)); // just a boolean marker as to whether each pattern is currently used or not (0 == False, 1 == True)
  uint8_t backtrack = 0; // flag to force backtracking
  while(1) {
    /*
    fprintf(stderr, "i = %u, n_used = %u, used_occ = %u, unused_occ = %u, backtrack = %u\n", i, n_used, used_occ, unused_occ, backtrack);
    fprintf(stderr, "used: ");
    for(j = 0; j < n; j++) {
      fprintf(stderr, "%u ", used[j]);
    }
    fprintf(stderr, "\n");
    */
    if(backtrack || i >= n) { // if we reached the end
      i--; // step back from the edge
      backtrack = 0; // reset backtrack flag

      // check if the current score is in our top set
      // - used_occ is the current "score"
      if(used_occ > scores[0]) {
        if(top_n[0]) // basically, if this isn't the very first score
          free(top_n[0]); // we only ever want free the first element, then just shift the rest
        char* newick_str = get_newick(a, used, n);
        for(j = 1; j < keep_top_n && used_occ > scores[j]; j++) { // find the slot where this score goes, sliding the rest down
          scores[j-1] = scores[j];
          top_n[j-1] = top_n[j];
        }
        // replace this slot with the current tree
        scores[j-1] = used_occ;
        top_n[j-1] = newick_str;
        fprintf(stderr, "with score %u:\n%s\n", used_occ, newick_str);
      }

      while(i < n && !used[i]) { // backtrack to most recent used pattern (we check i < n because it is unsigned and will UNDERFLOW when we want to stop)
        unused_occ = unused_occ - a[i].ct; // since we're moving before this, it's no longer "unused"
        i--;
      }
      //fprintf(stderr, "backtracked to %u\n", i);
      if(i >= n) break; // if we've already reset the first pattern, we're done, get out
      n_used--;
      used_occ = used_occ - a[i].ct;
      unused_occ = unused_occ + a[i].ct;
      used[i++] = 0;
    } else {
      // each item only has 2 states - USE and DON'T USE, so we'll just toggle them in the list
      if(pattern_ok(a, used, n, i, n_samples)) {
        n_used++;
        used_occ = used_occ + a[i].ct;
        used[i++] = 1;
      } else {
        unused_occ = unused_occ + a[i].ct;
        i++;
        if(used_occ + (total - unused_occ) <= scores[0]) { // prune this branch, it can't possibly beat our minimum best
          backtrack = 1; // force backtrack
        }
      }
    }
  }
  //free(used);
  return 0;
}

// TODO: under development!
// build a sort of maximum parsimony tree (assuming a static clock)
// m: number of samples
// n: matrix length
// n_patterns: number of patterns meeting the given threshold
int make_tree(khash_t(patternHash) *h, float pattern_threshold, uint32_t n, uint32_t n_patterns, uint32_t n_samples) {
  pa_ct *a = malloc(n_patterns * sizeof(pa_ct));
  khint_t bin; // hash bin (result of kh_put)
  uint32_t i = 0, j = 0;
  uint32_t tot = 0; // total occurrences
  uint32_t *singletons = calloc(n_samples, sizeof(uint32_t)); // counts of singleton pattern for each sample - they are always compatible, so we don't have to worry about them
  for(bin = kh_begin(h); bin != kh_end(h); ++bin) {
    if(kh_exist(h, bin) && kh_value(h, bin) >= pattern_threshold * n) {
      if(kh_key(h, bin) == 0 || ~kh_key(h, bin) == 0) {
        fprintf(stderr, "%llu k-mers do not vary\n", kh_value(h, bin));
        continue;
      }
      // check for singleton
      for(j = 0; j < n_samples; j++) {
        if(kh_key(h, bin) - (1<<j) == 0 || (~kh_key(h, bin)) - (1<<j) == 0) {
          singletons[j] = kh_value(h, bin);
          fprintf(stderr, "sample %u leaf dist: %llu\n", j, singletons[j]);
          break;
        }
      }
      if(j != n_samples) continue; // it was a singleton, don't include it
      //fprintf(stderr, "%s: %llu\n", tobitstring(kh_key(h, bin)), kh_value(h, bin));
      a[i].pattern = kh_key(h, bin);
      a[i++].ct = kh_value(h, bin);
      tot = tot + kh_value(h, bin);
    }
  }
  fprintf(stderr, "Making tree with %u total patterns\n", tot);
  // sort a (in decreasing order of occurrences)
  ks_mergesort(pa_ct_decreasing, n_patterns, a, 0);
  int keep_top_n = 10;
  char** top_n = calloc(keep_top_n, sizeof(char*));
  uint32_t* scores = calloc(keep_top_n, sizeof(uint32_t));
  if(bb(a, n_patterns, tot, top_n, scores, keep_top_n, n_samples) != 0) {
    fprintf(stderr, "ERROR: Tree construction failed somewhere in branch-and-bound parsimony search.\n");
  }
  for(i = 0; i < keep_top_n; i++) {
    if(!top_n[i]) break; // in case there are fewer than <keep_top_n> answers
    fprintf(stdout, "Score %u: %u\n", i, scores[i]);
    fprintf(stdout, "%s\n", top_n[i]);
    //free(top_n[i]);
  }
  free(scores);
  free(top_n);
  free(a);
  return 0;
}


int main(int argc, char *argv[]) {
  srand(242424); // so that it will always be the same?
  int verbose = 0;
  int coverage_threshold = 5;
  float pattern_threshold = 0.00001;
  int k = 21;
  int threads = 1;
  char* ref_fasta = NULL;
  char* out_file = NULL;
  int num_kmers = 0;
  int random_pa = 0;

  kvec_t(char*) read_fastas;
  kv_init(read_fastas);

  kvec_t(char*) matrices;
  kv_init(matrices);

  int opt, long_idx;
  opterr = 0;
  while ((opt = getopt_long(argc, argv, "f:r:k:c:m:o:p:n:hv", long_options, &long_idx)) != -1) {
    switch (opt) {
      case 'h':
        usage();
        return 0;
        break;
      case 'v':
        verbose = 1;
        break;
      case 'f':
        ref_fasta = optarg;
        break;
      case 'o':
        out_file = optarg;
        break;
      case 'r':
        kv_push(char*, read_fastas, optarg);
        break;
      case 'm':
        kv_push(char*, matrices, optarg);
        break;
      case 'k':
        k = atoi(optarg);
        break;
      case 'c':
        coverage_threshold = atoi(optarg);
        break;
      case 'p':
        pattern_threshold = atof(optarg);
        break;
      case 't':
        threads = atoi(optarg);
        break;
      case 'n':
        num_kmers = atoi(optarg);
        break;
      case '?':
        if (optopt == 'f' || optopt == 'k' || optopt == 'c' || optopt == 'm' || optopt == 'o' || optopt == 't' || optopt == 'p' || optopt == 'n')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      case 0:
        //fprintf(stderr, "option %s\n", long_options[long_idx].name);
        if (long_idx == 8) random_pa = 1; // --random
        //else if (long_idx == 1) fn = atof(optarg); // --fn
        break;
      default:
        usage();
        return 1;
    }
  }

  int index;
  int skip = 0;
  char* command = NULL;
  for (index = optind; index < argc; index++) {
    if(index == optind) {
      command = argv[index];
    } else {
      if(skip) {
        skip = 0;
        continue;
      }
      if(argv[index][0] == '-') {
        if(argv[index][1] == '-') {
          skip = 1;
        } else {
          if(strlen(argv[index]) == 2) {
            skip = 1;
          }
        }
      } else {
        // a positional argument, by default -- will be used as an input matrix file
        kv_push(char*, matrices, argv[index]);
      }
    }
  }
  if(command == NULL) {
    fprintf(stderr, "ERROR: missing command\n\n");
    usage();
    return 1;
  }

  if(strcmp(command, "make") == 0) {
    gzFile f = gzopen(ref_fasta, "r");
    kseq_t* seq = kseq_init(f);
    fprintf(stderr, "Loading FASTA file: %s\n", ref_fasta);

    khash_t(kmerHash) *h = kh_init(kmerHash); // allocate hash table
    kh_resize(kmerHash, h, 150000000); // make expected size so that we don't have to do a lot of copying later as it expands

    int l; // sequence length
    int i; // just a counter
    uint32_t kmer_i = 0; // k-mer index, monotonically increasing
    char *kmer;
    khint_t bin; // hash bin (result of kh_put)
    int absent;
    while ((l = kseq_read(seq)) >= 0) {
      // name: seq->name.s, seq: seq->seq.s, length: l
      //if(verbose)
      //  fprintf(stderr, "Reading sequence '%s' (%i bp).\n", seq->name.s, l);
      for(i = 0; i <= l-k; i++) {
        kmer = malloc((k+1) * sizeof(char));
        kmer[k] = '\0';
        strncpy(kmer, seq->seq.s+i, k);
        bin = kh_put(kmerHash, h, kmer, &absent);
        if(absent) { // bin is empty (unset)
          kh_value(h, bin) = kmer_i++;
        }
      }
    }
    kseq_destroy(seq);
    gzclose(f);

    fprintf(stderr, "%u unique k-mers\n", kmer_i);

    char tmp;
    uint32_t tot = 0;
    uint32_t unique = 0;

    uint8_t *cts = calloc(kmer_i, sizeof(uint8_t));

    int r;
    for(r = 0; r < kv_size(read_fastas); r++) {
      f = gzopen(kv_A(read_fastas, r), "r");
      seq = kseq_init(f);
      fprintf(stderr, "Reading FASTA file: %s\n", kv_A(read_fastas, r));
      while ((l = kseq_read(seq)) >= 0) {
        // name: seq->name.s, seq: seq->seq.s, length: l
        if(verbose)
          fprintf(stderr, "Reading sequence '%s' (%i bp).\n", seq->name.s, l);
        for(i = 0; i <= l-k; i++) {
          tmp = seq->seq.s[i+k];
          *(seq->seq.s+i+k) = '\0';
          bin = kh_get(kmerHash, h, seq->seq.s+i);
          *(seq->seq.s+i+k) = tmp;
          if(bin == kh_end(h)) {
          } else {
            if(!cts[kh_value(h, bin)])
              unique++;
            if(cts[kh_value(h, bin)] < 255)
              cts[kh_value(h, bin)]++;
            tot++;
          }
        }
      }
      kseq_destroy(seq);
      gzclose(f);
    }
    fprintf(stderr, "%u raw k-mers matched (may be counted multiple times)\n", tot);
    fprintf(stderr, "%u unique\n", unique);

    uint32_t met_threshold = 0;
    for(i = 0; i < kmer_i; i++) {
      if(cts[i] >= coverage_threshold)
        met_threshold++;
    }
    fprintf(stderr, "%u with count >= %d\n", met_threshold, coverage_threshold);

    FILE *fout = fopen(out_file, "wb");
    fwrite(&kmer_i, sizeof(uint32_t), 1, fout);
    fwrite(cts, sizeof(uint8_t), kmer_i, fout);
    fclose(fout);

    // free counts array
    free(cts);
    // free k-mer hash keys
    for(bin = kh_begin(h); bin != kh_end(h); ++bin) {
      if(kh_exist(h, bin))
        free((char*)kh_key(h, bin));
    }
    // free k-mer hash itself
    kh_destroy(kmerHash, h);
    // free read_fastas list
    kv_destroy(read_fastas);
  }

  if(strcmp(command, "collapse") == 0) {
    int m; // matrix index
    int o; // other matrix index
    int i; // index into matrix
    uint32_t n; // matrix size
    uint8_t **mm = malloc(kv_size(matrices) * sizeof(uint8_t*)); // array of matrices
    if(kv_size(matrices) > 64) {
      fprintf(stderr, "ERROR: collapse cannot currently accept more than 64 matrices because we are bit-packing them into an 8-byte int\n");
      return 1;
    }
    for(m = 0; m < kv_size(matrices); m++) {
      fprintf(stderr, "Reading '%s'\n", kv_A(matrices, m));
      FILE *f = fopen(kv_A(matrices, m), "rb");
      fread(&n, sizeof(uint32_t), 1, f); // read one uint32_t - this is the size of the matrix
      fprintf(stderr, "Reading matrix of size %u\n", n);
      mm[m] = malloc(n * sizeof(uint8_t));
      fread(mm[m], sizeof(uint8_t), n, f);
      /*
      for(i = 0; i < n; i++) {
        fprintf(stderr, "%u ", mm[m][i]);
      }
      fprintf(stderr, "\n");
      */
      fclose(f);
      uint32_t p_l = basename(kv_A(matrices, m));
      uint16_t p = p_l >> 16 & ((1<<16)-1);
      uint16_t l = p_l & ((1<<16)-1);
      //fprintf(stderr, "%u -> %u, %u\n", p_l, p, l);
      kv_A(matrices, m)[p+l] = '\0';
      kv_A(matrices, m) = kv_A(matrices, m)+p;
      fprintf(stderr, "  name: '%s'\n", kv_A(matrices, m));
    }

    fprintf(stderr, "Building pattern hash... this will take a while.\n");
    khash_t(patternHash) *h = kh_init(patternHash); // allocate hash table
    int absent;
    int idx;
    float *covg = malloc(kv_size(matrices) * sizeof(float));
    for(m = 0; m < kv_size(matrices); m++) covg[m] = get_covg(kv_A(matrices, m));
    khint_t bin; // hash bin (result of kh_put)
    for(i = 0; i < (num_kmers != 0 ? num_kmers : n); i++) {
      idx = i;
      if(num_kmers != 0) {
        idx = (int)(((double)rand()) / RAND_MAX * n);
      }

      // use this k-mer if:
      // 1) at least one sample has 5 <= x
      // 2) x <= 2*cov for EVERY sample
      int good = 0;
      for(m = 0; m < kv_size(matrices); m++) {
        // covg is a float, so it might not be exactly 0...
        if(covg[m] > 0.00001 && mm[m][idx] > 2 * covg[m]) {
          break;
        }
        if(mm[m][idx] >= coverage_threshold) {
          good = 1;
        }
      }
      if(m < kv_size(matrices) || !good) {
        i--;
        continue;
      }

      uint64_t bits = 0;
      for(m = 0; m < kv_size(matrices); m++) {
        if(random_pa && (covg[m] < 0.00001 || (rand()/(double)RAND_MAX < mm[m][idx] / covg[m]))) {
          bits = bits + (1<<m);
        } else if(!random_pa && mm[m][idx] >= coverage_threshold) {
          bits = bits + (1<<m);
        }
      }
      //if(bits & 1 == 1) bits = ~bits; // normalize bits so that the first bit is always 0
      bin = kh_put(patternHash, h, bits, &absent);
      if(absent) { // bin is empty (unset)
        kh_value(h, bin) = 1;
      } else {
        kh_value(h, bin)++;
      }
    }
    free(covg);

    fprintf(stderr, "Producing NEXUS output...\n");
    FILE *fout;
    if(out_file != NULL) {
      uint32_t n_good_patterns = 0;
      for(bin = kh_begin(h); bin != kh_end(h); ++bin) {
        if(kh_exist(h, bin) && kh_value(h, bin) >= pattern_threshold * n) {
          n_good_patterns++;
        }
      }

      fout = fopen(out_file, "w");
      fprintf(fout, "#NEXUS\n");
      fprintf(fout, "Begin data;\n");
      fprintf(fout, "Dimensions ntax=%u nchar=%u;\n", kv_size(matrices), n_good_patterns);
      fprintf(fout, "Format datatype=dna missing=? gap=- interleave=yes;\n");
      fprintf(fout, "Matrix\n");
      char **sample_seqs = malloc(kv_size(matrices) * sizeof(char*));
      int row_len = 80;
      for(i = 0; i < kv_size(matrices); i++) {
        sample_seqs[i] = malloc((row_len+1) * sizeof(char));
        sample_seqs[i][row_len] = '\0';
      }
      // print PA matrix encoded as nucleotides (A: absent, C: present)
      int r = 0; // position in row
      for(bin = kh_begin(h); bin != kh_end(h); ++bin) {
        if(kh_exist(h, bin) && kh_value(h, bin) >= pattern_threshold * n) {
          for(i = 0; i < kv_size(matrices); i++) {
            if((kh_key(h, bin)>>i)&1) {
              sample_seqs[i][r] = 'C';
            } else {
              sample_seqs[i][r] = 'A';
            }
          }
          if(++r >= row_len) {
            fprintf(stderr, "printing batch with length %d\n", r);
            for(i = 0; i < kv_size(matrices); i++) {
              fprintf(fout, "%s  %s\n", kv_A(matrices, i), sample_seqs[i]);
            }
            fprintf(fout, "\n");
            r = 0;
          }
        }
      }
      // last batch (<row_len)
      fprintf(stderr, "printing last batch with length %d\n", r);
      if(r > 0) {
        for(i = 0; i < kv_size(matrices); i++) {
          sample_seqs[i][r] = '\0';
          fprintf(fout, "%s  %s\n", kv_A(matrices, i), sample_seqs[i]);
        }
        fprintf(fout, "\n");
      }
      fprintf(fout, ";\n");
      fprintf(fout, "End;\n");
      fprintf(fout, "BEGIN ASSUMPTIONS; WTSET kmer_counts (VECTOR) =");
      // print weights
      int patterns = 0;
      for(bin = kh_begin(h); bin != kh_end(h); ++bin) {
        if(kh_exist(h, bin) && kh_value(h, bin) >= pattern_threshold * n) {
          fprintf(fout, " %llu", kh_value(h, bin));
          patterns++;
        }
      }
      fprintf(stderr, "total patterns with weights: %d\n", patterns);
      fprintf(fout, "; END;\n");

      // make and print the maximum parsimony tree
      //make_tree(h, pattern_threshold, n, patterns, kv_size(matrices));

      /*
      fprintf(fout, "begin taxa;\n");
      fprintf(fout, "  dimensions ntax=%u;\n", kv_size(matrices));
      fprintf(fout, "  taxlabels\n");
      fprintf(fout, "    ;\n");
      fprintf(fout, "end;\n");
      fprintf(fout, "begin trees;\n");
      fprintf(fout, "   tree quick_pattern_tree = %s\n", "newick");
      fprintf(fout, "end;\n");
      */
      fclose(fout);
    } else {
      // stdout output
      for(bin = kh_begin(h); bin != kh_end(h); ++bin) {
        if(kh_exist(h, bin) && kh_value(h, bin) >= pattern_threshold * n) {
          char *s = tobitstring(kh_key(h, bin));
          fprintf(stdout, "%s\t%llu\n", s, kh_value(h, bin));
          free(s);
        }
      }
    }

    // cleanup
    /*
    // don't need to clean these up - they were originally argvs...
    for(m = 0; m < kv_size(matrices); m++)
      free(mm[m]);
    */
    free(mm);
    // free pattern hash
    kh_destroy(patternHash, h);
  }

  if(strcmp(command, "count") == 0) {
    int m;
    int i;
    for(m = 0; m < kv_size(matrices); m++) {
      fprintf(stderr, "Reading '%s'\n", kv_A(matrices, m));
      FILE *f = fopen(kv_A(matrices, m), "rb");
      uint32_t n;
      fread(&n, sizeof(uint32_t), 1, f); // read one uint32_t - this is the size of the matrix
      fprintf(stderr, "Reading matrix of size %u\n", n);
      uint8_t *cts = malloc(n * sizeof(uint8_t)); // get this from the file itself later
      fread(cts, sizeof(uint8_t), n, f);
      fclose(f);

      uint32_t tot = 0;
      uint32_t unique = 0;
      uint32_t met_threshold = 0;
      for(i = 0; i < n; i++) {
        if(cts[i] > 0)
          unique++;
        if(cts[i] >= coverage_threshold)
          met_threshold++;
      }
      free(cts);

      fprintf(stderr, "%u unique found\n", unique);
      fprintf(stderr, "%u with count >= %d\n", met_threshold, coverage_threshold);
    }
  }

  if(strcmp(command, "triangle") == 0) {
    int m; // matrix index
    int o; // other matrix index
    int i; // index into matrix
    uint32_t n; // matrix size
    uint8_t **mm = malloc(kv_size(matrices) * sizeof(uint8_t*)); // array of matrices
    uint32_t intersect; // intersection of hits
    uint32_t union_; // union of hits
    fprintf(stdout, "\t%u\n", kv_size(matrices));
    for(m = 0; m < kv_size(matrices); m++) {
      //fprintf(stderr, "Reading '%s'\n", kv_A(matrices, m));
      fprintf(stdout, "%s", kv_A(matrices, m));
      FILE *f = fopen(kv_A(matrices, m), "rb");
      fread(&n, sizeof(uint32_t), 1, f); // read one uint32_t - this is the size of the matrix
      //fprintf(stderr, "Reading matrix of size %u\n", n);
      mm[m] = malloc(n * sizeof(uint8_t)); // get this from the file itself later
      fread(mm[m], sizeof(uint8_t), n, f);
      fclose(f);

      for(o = 0; o < m; o++) {
        intersect = 0;
        union_ = 0;
        for(i = 0; i < n; i++) {
          if(mm[m][i] >= coverage_threshold || mm[o][i] >= coverage_threshold) {
            union_++;
            if(mm[m][i] >= coverage_threshold && mm[o][i] >= coverage_threshold)
              intersect++;
          }
        }
        // do MASH2-type distance adjustment (based on k) WARNING: this will use the value of k given for this run, not necessarily that used to compute the matrices
        float j = (float)intersect / union_; // jaccard estimate
        float dist = (j == 0) ? 1 : (-1.0/k * log(2*j / (1+j)));
        fprintf(stdout, "\t%f", dist);
      }
      fprintf(stdout, "\n");
    }

    // cleanup
    for(m = 0; m < kv_size(matrices); m++)
      free(mm[m]);
    free(mm);
  }
}
