#!/bin/bash

# Get Mouse atlas data from Marioni

mkdir Analysis/Pijuan/data
mkdir Analysis/Pijuan/data/raw

curl https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz > ./Analysis/Pijuan/data/raw/atlas_data.tar.gz

tar -xf ./Analysis/Pijuan/data/raw/atlas_data.tar.gz atlas/genes.tsv atlas/barcodes.tsv atlas/meta.tab atlas/raw_counts.mtx -C ./Analysis/Pijuan/data/raw
mv -f ./atlas/* ./Analysis/Pijuan/data/raw
rm -r ./atlas