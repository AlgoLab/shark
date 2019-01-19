#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

#include <zlib.h>

#include "kseq.h"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/util.hpp"

#include "argument_parser.hpp"
#include "bloomfilter.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

auto start_t = chrono::high_resolution_clock::now();

void pelapsed(const string &s = "") {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[read-filter/" << s << "] Time elapsed "
       << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000
       << endl;
}

/*** METHODS TO READ AND PARSE A GTF FILE *************************************/
/**
 * !!! I'm assuming a GTF from ensembl !!!
 * Given a gtf feature, returns the third (type) and the last
 * (attributes) columns as a pair of strings.
 **/
pair<string, string> parse_line(string line) {
  int i = 0;
  string delimiter = "\t";
  string token;
  size_t pos;
  string type = "";
  string attributes = "";
  while ((pos = line.find(delimiter)) != string::npos) {
    ++i;
    token = line.substr(0, pos);
    line.erase(0, pos + delimiter.length());
    if(i == 3)
      type = token;
  }
  attributes = line;
  return make_pair(type, attributes);
}

/**
 * Given the attributes of a feature, returns the gene_id and the
 * transcript_id as a pair of strings.
 **/
pair<string, string> parse_attributes(string line) {
  string gene_idx = "";
  string transcript_idx = "";

  string gene_key = "gene_id";
  string transcript_key = "transcript_id";

  string delimiter = "; ";
  string token;
  size_t pos;

  while ((pos = line.find(delimiter)) != string::npos) {
    token = line.substr(0, pos);
    line.erase(0, pos + delimiter.length());
    if(token.compare(0, gene_key.size(), gene_key) == 0)
      gene_idx = token.substr(gene_key.size()+2, token.size()-gene_key.size()-3);
    else if(token.compare(0, transcript_key.size(), transcript_key) == 0)
      transcript_idx = token.substr(transcript_key.size()+2, token.size()-transcript_key.size()-3);
  }
  if(line.compare(0, gene_key.size(), gene_key) == 0)
    gene_idx = line.substr(gene_key.size()+2, line.size()-gene_key.size()-3);
  else if(line.compare(0, transcript_key.size(), transcript_key) == 0)
    transcript_idx = line.substr(transcript_key.size()+2, line.size()-transcript_key.size()-3);

  return make_pair(gene_idx, transcript_idx);
}

/**
 * Given a gtf file, extract a map <transcript_idx : gene_idx>
 **/
unordered_map<string, string> read_gtf(string gtf_path) {
  unordered_map<string, string> transcriptidx2geneidx;

  std::ifstream gtf;
  std::string line;
  gtf.open(gtf_path);
  if(gtf.is_open()) {
    while(getline(gtf,line)) {
      pair<string, string> type_attributes = parse_line(line);
      if(type_attributes.first.compare("transcript") == 0) {
        pair<string, string> gene_transcript = parse_attributes(type_attributes.second);
        transcriptidx2geneidx[gene_transcript.second] = gene_transcript.first;
      }
    }
  }
  gtf.close();
  return transcriptidx2geneidx;
}

/******************************************************************************/

void analyze_read(const kseq_t *seq, BF &bloom, vector<string> legend_ID, const uint &k, const uint &c) {
  map<int, int> classification_id;
  string read_seq = seq->seq.s;

  if (read_seq.size() >= k) {
    string kmer (read_seq, 0, k);
    transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
    // cout << kmer << endl;
    IDView id_kmer = bloom.get_index(kmer);
    while (id_kmer.has_next()) {
      int x = id_kmer.get_next();
      // cout << "- " << x << endl;
      ++classification_id[x];
    }

    for (uint p = k; p < read_seq.size(); ++p) {
      char c = toupper(read_seq[p]);
      kmer.erase(0, 1);
      kmer += c;
      // cout << kmer << endl;
      id_kmer = bloom.get_index(kmer);
      while (id_kmer.has_next()) {
        int x = id_kmer.get_next();
        // cout << "- " << x << endl;
        ++classification_id[x];
      }
    }
  }

  int max = std::max(max_element(begin(classification_id),
                                 end(classification_id),
                                 [](const pair<int, int> &a,
                                    const pair<int, int> &b) {
                                   return a.second < b.second;
                                 })
                     ->second,
                     (int)3);
  // save in file all headers of transcripts probably associated to the read
  // search and store in a file all elements of classification with
  // max(classification[i])
  for (auto it_class = classification_id.cbegin(); it_class != classification_id.cend(); ++it_class) {
    if (it_class->second == max && it_class->second >= opt::c) {
      // legendid[it_class->first] is the name of the transcript, mapped
      /// with index it_class->first
      cout << seq->name.s << "\t" << legend_ID[it_class->first] << "\t" << read_seq << endl;
    }
  }
}

/*****************************************
 * Main
 *****************************************/
int main(int argc, char *argv[]) {
  parse_arguments(argc, argv);

  /*** 0. Check input files and initialize variables **************************/
  // Transcripts
  gzFile transcript_file = gzopen(opt::fasta_path.c_str(), "r");
  kseq_t *seq = kseq_init(transcript_file);

  // Annotation
  std::ifstream gtf;
  gtf.open(opt::gtf_path);
  gtf.close();

  // Sample 1
  gzFile read1_file = gzopen(opt::sample1_path.c_str(), "r");
  seq = kseq_init(read1_file);

  // Sample 2
  gzFile read2_file;
  if(opt::paired_flag) {
    read2_file = gzopen(opt::sample2_path.c_str(), "r");
    seq = kseq_init(read2_file);
  }

  BF bloom(opt::bf_size);
  vector<string> legend_ID;
  int file_line;

  if(opt::verbose) {
    cerr << "Transcripts: " << opt::fasta_path << endl;
    cerr << "Annotation: " << opt::gtf_path << endl;
    cerr << "Sample 1: " << opt::sample1_path << endl;
    if(opt::paired_flag)
      cerr << "Sample 2: " << opt::sample2_path << endl;
    cerr << "K-mer length: " << opt::k << endl;
    cerr << "Threshold value: " << opt::c << endl << endl;
  }

  // Map <transcript_idx : gene_idx>
  unordered_map<string, string> transcriptidx2geneidx = read_gtf(opt::gtf_path);

  /****************************************************************************/

  /*** 1. First iteration over transcripts ************************************/
  seq = kseq_init(transcript_file);
  while ((file_line = kseq_read(seq)) >= 0) {
    string input_seq = seq->seq.s;

    if (input_seq.size() >= opt::k) {
      // Build kmers and add them to bf
      string kmer (input_seq, 0, opt::k);
      transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      bloom.add_kmer(kmer);
      for (uint p = opt::k; p < input_seq.size(); ++p) {
        char c = toupper(input_seq[p]);
        kmer.erase(0, 1);
        kmer += c;
        bloom.add_kmer(kmer);
      }
    }
  }

  kseq_destroy(seq);
  gzclose(transcript_file);

  pelapsed("Transcript file processed");

  bloom.switch_mode(1);

  /****************************************************************************/

  /*** 2. Second iteration over transcripts ***********************************/
  transcript_file = gzopen(opt::fasta_path.c_str(), "r");
  seq = kseq_init(transcript_file);
  string last_gene_idx = "";
  int nidx = -1;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((file_line = kseq_read(seq)) >= 0) {
    string input_name = seq->name.s;
    /**
     * ! this is needed to remove inconsistency between transcript ids
     * in the ensembl cdna (ending with .number) and ensembl annotation !
     **/
    size_t pos;
    while ((pos = input_name.find(".")) != string::npos) {
      input_name = input_name.substr(0, pos);
      break;
    }
    // input_name = input_name.substr(0, input_name.size()-2); // this doesn't work if the number is >9

    string input_seq = seq->seq.s;

    if (input_seq.size() >= opt::k) {
      // Storing a numeric ID for each gene
      if(transcriptidx2geneidx[input_name] != last_gene_idx) {
        ++nidx;
        legend_ID.push_back(transcriptidx2geneidx[input_name]);
        last_gene_idx = transcriptidx2geneidx[input_name];
      }

      // Build kmers and store their nidx in the bf
      string kmer (input_seq, 0, opt::k);
      transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      // cout << "# " << kmer << " " << nidx << endl;
      bloom.add_to_kmer(kmer, nidx);
      for (uint p = opt::k; p < input_seq.size(); ++p) {
        char c = toupper(input_seq[p]);
        kmer.erase(0, 1);
        kmer += c;
        // cout << "# " << kmer << " " << nidx << endl;
        bloom.add_to_kmer(kmer, nidx);
      }
    }
  }

  kseq_destroy(seq);
  gzclose(transcript_file);

  pelapsed("BF created from transcripts (" + to_string(nidx+1) + " genes)");

  bloom.switch_mode(2);

  /****************************************************************************/

  /*** 3. Iteration over first sample *****************************************/
  seq = kseq_init(read1_file);
  while ((file_line = kseq_read(seq)) >= 0) {
    analyze_read(seq, bloom, legend_ID, opt::k, opt::c);
  }

  kseq_destroy(seq);
  gzclose(read1_file);

  pelapsed("First sample completed");

  /****************************************************************************/

  /*** 4. Iteration over second sample ****************************************/
  if(opt::paired_flag){
    read2_file = gzopen(opt::sample2_path.c_str(), "r");

    seq = kseq_init(read2_file);
    while ((file_line = kseq_read(seq)) >= 0) {
      analyze_read(seq, bloom, legend_ID, opt::k, opt::c);
    }

    kseq_destroy(seq);
    gzclose(read2_file);

    pelapsed("Second sample completed");
  }

  /****************************************************************************/

  pelapsed("Association done");
  return 0;
}
