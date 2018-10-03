#!/bin/bash

cd sdsl-lite/build
./build.sh
cd ../../KMC
make
cd ../htslib
make
cd ..
