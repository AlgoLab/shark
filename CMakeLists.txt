cmake_minimum_required(VERSION 3.9)
project(read-filter)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)
include_directories(htslib)
include_directories(htslib/cram)
include_directories(htslib/htslib)
include_directories(htslib/os)
include_directories(KMC/kmc_api)
include_directories(KMC/kmc_dump)
include_directories(KMC/kmc_dump_sample)
include_directories(KMC/kmc_tools)
include_directories(KMC/kmc_tools/libs)
include_directories(KMC/kmer_counter)
include_directories(KMC/kmer_counter/libs)
include_directories(KMC/kmer_counter/libs/vectorclass)

include_directories("/usr/include")
link_directories("/usr/lib")

include_directories("/home/tamara/include")
link_directories("/home/tamara/lib")

add_executable(read-filter main.cpp bloomfilter.h htslib KMC MurmurHash3.hpp MurmurHash3.cpp)

target_link_libraries(read-filter sdsl divsufsort divsufsort64)