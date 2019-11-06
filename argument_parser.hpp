#ifndef _ARGUMENT_PARSER_HPP_
#define _ARGUMENT_PARSER_HPP_

#include <iostream>
#include <sstream>
#include <getopt.h>

static const char *USAGE_MESSAGE =
  "Usage: shark [-v] -r <references> -1 <sample1> [-2 <sample2>] [-k <kmer size>] [-c <confidence>] [-b <filter size>] [-q <min base quality>] [-s]\n"
  "Top notch description of this tool\n"
  "\n"
  "      -h, --help                        display this help and exit\n"
  "      -r, --reference                   reference sequences in FASTA format (can be gzipped)\n"
  "      -1, --sample1                     sample in FASTA/Q (can be gzipped)\n"
  "      -2, --sample2                     second sample in FASTA/Q (optional, can be gzipped)\n"
  "      -k, --kmer-size                   size of the kmers to index (default:17)\n"
  "      -c, --confidence                  confidence for associating a read to a gene (default:0.6)\n"
  "      -b, --bf-size                     bloom filter size in GB (default:1)\n"
  "      -q, --min-base-quality            minimum base quality (assume FASTQ Illumina 1.8+ Phred scale, default:0, i.e., no filtering)\n"
  "      -s, --single                      report an association only if a single gene is found\n"
  "      -t, --threads                     number of threads (default:1)\n"
  "      -v, --verbose                     verbose mode\n"
  "\n";

namespace opt {
  static std::string fasta_path = "";
  static std::string sample1_path = "";
  static std::string sample2_path = "";
  static std::string out1_path = "";
  static std::string out2_path = "";
  static bool paired_flag = false;
  static uint k = 17;
  static double c = 0.6;
  static uint64_t bf_size = ((uint64_t)0b1 << 33);
  static char min_quality = 0;
  static bool single = false;
  static bool verbose = false;
  static int nThreads = 1;
}

static const char *shortopts = "t:r:1:2:o:p:k:c:b:q:svh";

static const struct option longopts[] = {
  {"reference", required_argument, NULL, 'r'},
  {"threads", required_argument, NULL, 't'},
  {"sample1", required_argument, NULL, '1'},
  {"sample2", required_argument, NULL, '2'},
  {"out1", required_argument, NULL, 'o'},
  {"out2", required_argument, NULL, 'p'},
  {"kmer-size", required_argument, NULL, 'k'},
  {"confidence", required_argument, NULL, 'c'},
  {"bf-size", required_argument, NULL, 'b'},
  {"min-base-quality", required_argument, NULL, 'q'},
  {"single", no_argument, NULL, 's'},
  {"verbose", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, 'h'},
  {NULL, 0, NULL, 0}
};

void parse_arguments(int argc, char **argv) {
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1; ) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'r':
      arg >> opt::fasta_path;
      break;
    case 't':
      arg >> opt::nThreads;
      break;
    case '1':
      arg >> opt::sample1_path;
      break;
    case '2':
      arg >> opt::sample2_path;
      opt::paired_flag = true;
      break;
    case 'o':
      arg >> opt::out1_path;
      break;
    case 'p':
      arg >> opt::out2_path;
      break;
    case 'k':
      arg >> opt::k;
      break;
    case 'c':
      arg >> opt::c;
      break;
    case 'b':
      // Let's consider this as GB
      arg >> opt::bf_size;
      opt::bf_size = opt::bf_size * ((uint64_t)0b1 << 33);
      break;
    case 'q':
      int mq;
      arg >> mq;
      opt::min_quality = static_cast<char>(mq);
      break;
    case 's':
      opt::single = true;
      break;
    case 'v':
      opt::verbose = true;
      break;
    case 'h':
      std::cerr << USAGE_MESSAGE;
      exit(EXIT_SUCCESS);
    default:
      std::cerr << "shark : unknown argument" << std::endl;
      std::cerr << "\n" << USAGE_MESSAGE;
      exit(EXIT_FAILURE);
    }
  }

  if (opt::fasta_path == "" || opt::sample1_path == "") {
    std::cerr << "shark : missing required files" << std::endl;
    std::cerr << "\n" << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }
}

#endif
