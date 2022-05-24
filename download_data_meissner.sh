#!/bin/bash

# Get Mouse atlas data from Marioni

mkdir Analysis/Meissner/data
mkdir Analysis/Meissner/data/raw

curl https://cloud.mpi-cbg.de/index.php/s/6eTwcA6Bhrw9l0f/download > ./Analysis/Meissner/data/raw/atlas_data.tar.gz

tar -xf ./Analysis/Meissner/data/raw/atlas_data.tar.gz atlas/genes.tsv atlas/barcodes.tsv atlas/meta.csv atlas/matrix.mtx -C ./Analysis/Meissner/data/raw
mv -f ./atlas/* ./Analysis/Meissner/data/raw
rm -r ./atlas