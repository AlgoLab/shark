#include "bloomfilter.h" 
#include "kseq.h" 
#include "sdsl/int_vector.hpp" 
#include "sdsl/int_vector.hpp" // for the bit_vector class 
#include "sdsl/util.hpp" 
#include <algorithm> 
#include <fstream> 
#include <iostream> 
#include <stdio.h> 
#include <string> 
#include <vector> // std::vector 
#include <zlib.h>

static size_t sizebloom = (10e8);



// function search, returns indexes in a vector
IDView search(const string &kmer, BF &bloomfilter) {
  return bloomfilter.get_index(kmer);
}
// using namespace seqan;

void print_menu() {
  cout << endl;
  cout << "1. Choose k value - default: k length = 31" << endl;
  cout << "2. Change threshold value:" << endl;
  cout << "\t Default value: reads with 0 error (read's length - k + 1) "
       << endl;
  cout << "3. Launch tool" << endl;
  cout << endl;
  cout << "Choice: ";
}

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)


/*****************************************
 * Main
 *****************************************/

int main(int argc, char *argv[]) {
  gzFile transcript_file;
  kseq_t *seq;
  int file_line;
  map<int, string> legend_ID;
  int mapped_ID = 0;
  size_t kmer_length = 31;
  vector<string> transcript_kmers_vec;
  BF bloom(sizebloom);
  int threshold;
  bool double_argument = false;
  string transcript_name;
  string read_name;
  string read_name2;


  time_t start = time(0);


    if (argc > 1) {
    	for(int i =0; i < argc; i++){
    	string option = argv[i];
    		if (option == "-t")
    			transcript_name = argv[i+1];
    		if (option == "-r1")
    			read_name = argv[i+1];
    		if (option == "-r2"){
    			double_argument = true;
    			read_name2 = argv[i+1];
    		}
    	}
    } else {
      std::cerr << "Error in input" << std::endl;
      return 1;
    }

	cout << transcript_name << endl;
	cout << read_name << endl;
	cout << read_name2 << endl;

  threshold = 33;
  cout << "K-mer length: " << kmer_length << endl;
  cout << "Threshold value: " << threshold << endl;
  
  time_t read_time = time(0);


  transcript_file =
      gzopen(transcript_name.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(transcript_file);         // STEP 3: initialize seq
  string input_seq;
  // open and read the .fa
  while ((file_line = kseq_read(seq)) >=
         0) { // STEP 4: read sequence of transcript

    // seq name is the key map for the transcript, and has an assigned int
    input_seq = seq->seq.s;
    //  std::cerr << name_transcript << std::endl;

    if (input_seq.size() >= kmer_length) {
      legend_ID[mapped_ID] = seq->name.s;
      // split each sequence of a transcritp in k-mer with k=n

      transcript_kmers_vec.resize(input_seq.size() - (kmer_length - 1));
      transform(input_seq.cbegin(), input_seq.cend() - (kmer_length - 1),
                transcript_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // add all k-mers to BF

      for (auto &kmer : transcript_kmers_vec) {
        bloom.add_kmer(kmer);
      }
      mapped_ID++;

      transcript_kmers_vec.clear();
    }
    input_seq.clear();
  }

  std::cerr << "return value: " << file_line << std::endl;
  kseq_destroy(seq);        // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler
  

  std::cerr << "Transcript file processed" << std::endl;
  cout << mapped_ID << endl;

  int seconds = 0;
  int minutes = 0;
  int hour = 0;
  time_t end_read_time = time(0);
  int time_reading = difftime(end_read_time, read_time) * 1000.0;

  seconds = (time_reading / 1000) % 60;
  minutes = (int)((time_reading / (1000 * 60)) % 60);

  std::cerr << "Time to process transcripts:  " << minutes << ":" << seconds
            << std::endl;

  bloom.switch_mode(1);

  time_t add_time = time(0);

  transcript_file = gzopen(transcript_name.c_str(), "r");
  seq = kseq_init(transcript_file);
  int idx = 0;
  // open and read the .fa, every time a kmer is found the relative index is
  // added to BF
  while ((file_line = kseq_read(seq)) >= 0) {
    input_seq = seq->seq.s;

    if (input_seq.size() >= kmer_length) {

      // split each sequence of a transcritp in k-mer with k=n

      transcript_kmers_vec.resize(input_seq.size() - (kmer_length - 1));
      transform(input_seq.cbegin(), input_seq.cend() - (kmer_length - 1),
                transcript_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // add for each k-mer its id to BF
      for (auto &kmer : transcript_kmers_vec)
        bloom.add_to_kmer(kmer, idx);

      idx++;
    }
  }

  std::cerr << "return value: " << file_line << std::endl;
  kseq_destroy(seq);        // STEP 5: destroy seq
  gzclose(transcript_file); // STEP 6: close the file handler

  
  std::cerr << "Transcript indexes added to Bloom filter" << std::endl;
  time_t end_add_time = time(0);
  int time_adding = difftime(end_add_time, add_time) * 1000.0;

  seconds = (int)((time_adding / 1000) % 60);
  minutes = (int)((time_adding / (1000 * 60)) % 60);

  std::cerr << "Time to add indexes:  " << minutes << ":" << seconds
            << std::endl;

  bloom.switch_mode(2);

  time_t association_time = time(0);
  
  FILE * pFile;
  string output_file_name = "id_results_test2.tsv";
  pFile = fopen (output_file_name.c_str(), "w");
  //ofstream file;
  //file.open(output_file_name);
  gzFile read_file;
  string read_seq;
  vector<string> read_kmers_vec;
  map<int, int> classification_id;
  IDView id_kmer;

// open .fq file that contains the reads
  read_file = gzopen(read_name.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(read_file);                 // STEP 3: initialize seq

  while ((file_line = kseq_read(seq)) >= 0) { // STEP 4: read sequence of reads

    read_seq = seq->seq.s;

    if (read_seq.size() >= kmer_length) {

      read_kmers_vec.resize(read_seq.size() - (kmer_length - 1));
      transform(read_seq.cbegin(), read_seq.cend() - (kmer_length - 1),
                read_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // search all kmers of each read in the BF and store relative indexes

      for (auto &kmer : read_kmers_vec) {
        id_kmer = bloom.get_index(kmer);
        //		int* buffer = &id_kmer[0];
        //		fwrite (buffer, sizeof(int), sizeof(buffer), pFile);
        while (id_kmer.has_next())
          classification_id[id_kmer.get_next()]++;
      }
    }

    int max = std::max(std::max_element(begin(classification_id),
                                        end(classification_id),
                                        [](const std::pair<int, int> &a,
                                           const std::pair<int, int> &b) {
                                          return a.second < b.second;
                                        })
                           ->second,
                       (int)3);
    // save in file all headers of transcripts probably associated to the read
    // search and store in a file all elements of classification with
    // max(classification[i])
    for (auto it_class = classification_id.cbegin();
         it_class != classification_id.cend(); it_class++) {
      if (it_class->second == max && it_class->second >= threshold) {
        // legendid[it_class->first] is the name of the transcript, mapped
        /// with index it_class->first
        //file << seq->name.s << "\t" << (legend_ID.at(it_class->first)).c_str()
        //     << "\n";
               fwrite(seq->name.s, 1, seq->name.l, pFile);
               fwrite("\t", 1, sizeof("\t"), pFile);
               fwrite((legend_ID.at(it_class->first)).c_str(), 1,
               strlen((legend_ID.at(it_class->first)).c_str()), pFile);
               fwrite("\n", 1, sizeof("\n"), pFile);
               fflush(pFile);
      }
    }
    classification_id.clear();
    id_kmer.clear();
    read_kmers_vec.clear();
  }

  //  fclose(pFile);

  std::cerr << "return value: " << file_line << std::endl;
  kseq_destroy(seq);  // STEP 5: destroy seq
  gzclose(read_file); // STEP 6: close the file handler

  if(double_argument){
  string read_seq;
  vector<string> read_kmers_vec;
  map<int, int> classification_id;
  IDView id_kmer;

// open .fq file that contains the reads
  read_file = gzopen(read_name2.c_str(), "r"); // STEP 2: open the file handler
  seq = kseq_init(read_file);                 // STEP 3: initialize seq

  while ((file_line = kseq_read(seq)) >= 0) { // STEP 4: read sequence of reads

    read_seq = seq->seq.s;

    if (read_seq.size() >= kmer_length) {

      read_kmers_vec.resize(read_seq.size() - (kmer_length - 1));
      transform(read_seq.cbegin(), read_seq.cend() - (kmer_length - 1),
                read_kmers_vec.begin(), [kmer_length](const auto &i) {
                  return string(&i, kmer_length);
                });

      // search all kmers of each read in the BF and store relative indexes

      for (auto &kmer : read_kmers_vec) {
        id_kmer = bloom.get_index(kmer);
        //		int* buffer = &id_kmer[0];
        //		fwrite (buffer, sizeof(int), sizeof(buffer), pFile);
        while (id_kmer.has_next())
          classification_id[id_kmer.get_next()]++;
      }
    }

    int max = std::max(std::max_element(begin(classification_id),
                                        end(classification_id),
                                        [](const std::pair<int, int> &a,
                                           const std::pair<int, int> &b) {
                                          return a.second < b.second;
                                        })
                           ->second,
                       (int)3);
    // save in file all headers of transcripts probably associated to the read
    // search and store in a file all elements of classification with
    // max(classification[i])
    for (auto it_class = classification_id.cbegin();
         it_class != classification_id.cend(); it_class++) {
      if (it_class->second == max && it_class->second >= threshold) {
        // legendid[it_class->first] is the name of the transcript, mapped
        /// with index it_class->first
        //file << seq->name.s << "\t" << (legend_ID.at(it_class->first)).c_str()
        //     << "\n";
               fwrite(seq->name.s, 1, seq->name.l, pFile);
               fwrite("\t", 1, sizeof("\t"), pFile);
               fwrite((legend_ID.at(it_class->first)).c_str(), 1,
               strlen((legend_ID.at(it_class->first)).c_str()), pFile);
               fwrite("\n", 1, sizeof("\n"), pFile);
               fflush(pFile);
      }
    }
    classification_id.clear();
    id_kmer.clear();
    read_kmers_vec.clear();
  }

  fclose(pFile);

  std::cerr << "return value: " << file_line << std::endl;
  kseq_destroy(seq);  // STEP 5: destroy seq
  gzclose(read_file); // STEP 6: close the file handler

    }
  //file.close();
  
  std::cerr << "Output file created:" << output_file_name << endl; 
  std::cerr << "Association done." << std::endl;
  time_t end_assoc_time = time(0);
  int time_assoc = difftime(end_assoc_time, association_time) * 1000.0;

  seconds = (int)((time_assoc / 1000) % 60);
  minutes = (int)((time_assoc / (1000 * 60)) % 60);
  hour = (int) (minutes / 60);
  std::cerr << "Time to associate reads to transcript:  " << minutes << ":"
            << seconds << std::endl;

  time_t end = time(0);
  int time_proc = difftime(end, start) * 1000.0;

  seconds = (int)((time_proc / 1000) % 60);
  minutes = (int)((time_proc / (1000 * 60)) % 60);

  std::cerr << "Time of entire process:  " << hour << ":" << minutes << ":" << seconds
            << std::endl;

  std::cerr << "Reading transcripts: " << (int)(time_reading / time_proc) * 100
            << std::endl;
  std::cerr << "Adding indexes to BF: " << (int)(time_adding / time_proc) * 100
            << std::endl;
  std::cerr << "Associating reads to transcript: "
            << (int)(time_assoc / time_proc) * 100 << std::endl;

  return 0;
}
