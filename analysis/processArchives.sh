#!/bin/bash
cd ../data/data_raw
for archive in ./*
do
    unzip "${archive}"
done
mkdir eQTL_finemapping_results
mv *.gz eQTL_finemapping_results
rm -f *.tbi

cd eQTL_finemapping_results
for folder in ./*.gz
do
    if test -e "${folder}"; then
        gunzip "${folder}"
    fi
done
cd ..
cd ..
mkdir data_clean/eQTL_finemapping
mv data_raw/README data_clean/eQTL_finemapping
mv data_raw/eQTL_finemapping_results/*.txt data_clean/eQTL_finemapping
rm -r -f data_raw/eQTL_finemapping_results

cd data_raw/Sprime\ results\ for\ 1000\ Genomes\ non-African\ populations\ and\ SGDP\ Papuans
for folder in ./*.gz
do
    if test -e "${folder}"; then
        gunzip "${folder}"
    fi
done

for tar in ./*.tar
do
    tar -xvf ${tar}
done
cd ..
cd ..
mkdir data_clean/Sprime_results
mv data_raw/Sprime\ results\ for\ 1000\ Genomes\ non-African\ populations\ and\ SGDP\ Papuans/README data_clean/Sprime_results
mv data_raw/Sprime\ results\ for\ 1000\ Genomes\ non-African\ populations\ and\ SGDP\ Papuans/mendeley_data/* data_clean/Sprime_results
rm -r -f data_raw/Sprime\ results\ for\ 1000\ Genomes\ non-African\ populations\ and\ SGDP\ Papuans

