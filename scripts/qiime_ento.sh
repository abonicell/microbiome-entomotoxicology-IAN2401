#!/bin/bash

# Import paired-end reads using a manifest file
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Summarize demultiplexed sequences
qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization paired-end-demux.qzv

# Denoise with DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-chimera-method pooled \
  --p-trim-left-f 20 \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 160

# Visualize denoising stats
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# Summarize feature table
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

# Tabulate representative sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# Generate phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree-rep-seqs.qza \
  --o-rooted-tree rooted-tree-rep-seqs.qza

# Taxonomic classification with SILVA 138.2 (V4 region)
qiime feature-classifier classify-sklearn \
  --i-classifier /Users/andreabonicelli/Documents/forensOMICS/silva_2024.10/silva-138.2-ssu-nr99-seqs-515f-806r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy-silva-v4.qza \
  --p-n-jobs 0 \
  --p-reads-per-batch auto \
  --verbose

# Visualize taxonomy classification
qiime metadata tabulate \
  --m-input-file taxonomy-silva-v4.qza \
  --o-visualization taxonomy-silva-v4.qzv

# Create taxonomy bar plots
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy-silva-v4.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization taxa-bar-plots-silva-v4.qzv
