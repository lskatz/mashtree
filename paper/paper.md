---
title: 'Mashtree: a rapid comparison of whole genome sequence files'
authors:
- affiliation: "1, 2"
  name: Lee S. Katz
  orcid: 0000-0002-2533-9161
- affiliation: 1
  name: Taylor Griswold
- affiliation: 3
  name: Shatavia Morrison
  orcid: 0000-0002-4658-5951
- affiliation: 3
  name: Jason Caravas
  orcid: 0000-0001-9111-406X
- affiliation: 2
  name: Shaokang Zhang
  orcid: 0000-0003-0874-2212
- affiliation: 2
  name: Henk C. den Bakker
  orcid: 0000-0002-4086-1580
- affiliation: 2
  name: Xiangyu Deng
- affiliation: 2
  name: Heather A. Carleton
date: "11 September 2019"
output:
  html_document:
    df_print: paged
  pdf_document: default
  word_document: default
bibliography: paper.bib
tags:
- dendrogram
- mash
- sketch
- tree
- rapid
affiliations:
- index: 1
  name: Enteric Diseases Laboratory Branch, Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 2
  name: Center for Food Safety, University of Georgia, Griffin, GA, USA
- index: 3
  name: Respiratory Diseases Laboratory Branch, Centers for Disease Control and Prevention,
    Atlanta, GA, USA
---

# Summary

In the past decade, the number of publicly available bacterial genomes has increased dramatically.
These genomes have been generated for impactful initiatives, especially in the field of genomic epidemiology [@Timme:2017; @Brown:2019].
Genomes are sequenced, shared publicly, and subsequently analyzed for phylogenetic relatedness.
If two genomes of epidemiological interest are found to be related, further investigation might be prompted.
However, comparing the multitudes of genomes for phylogenetic relatedness is computationally expensive and, with large numbers, laborious.

Consequently, there are many strategies to reduce the complexity of the data for downstream analysis.
One strategy is to create a fast and accurate *de novo* assembly with software such as SKESA [@Souvorov:2018].
The sizes of the resulting assemblies are usually in the megabases, whereas the raw data are usually in the hundreds of megabases. Therefore, downstream analyses on the compact assemblies can be much faster.

One strategy is to create multilocus sequence type (MLST) profiles from each genome to be analyzed.  In MLST analysis, only the allelic variants are analyzed, either via reference-based mapping or via *de novo* genome assemblies. This reduction of data from raw sequence reads or assemblies, down to a subset of seven to thousands of genes, reduces the complexity of the analysis from millions of base pairs to a subset of genes and their allele calls [@Moura:2017;@Maiden:2013; @Alikhan:2018]. The allele calls can then be used to compute a distance matrix rapidly.

Another fast downstream analysis that is routinely performed is on the NCBI Pathogens site, where the SKESA assemblies are used to identify clusters of genomes in a whole-genome MLST (wgMLST) analysis.
For each cluster of genomes as determined by wgMLST, single nucleotide polymorphisms (SNPs) are identified. These SNPs are used to inter a phylogenetic signal, and those results are displayed online [@NCBI_Pathogens].

Yet another classic strategy is reducing each genome down to split kmers. With split kmer analysis, kmers on both sides of a variable site are recorded, and the variable nucleotide is identified. When comparing two or more genomes, the variable sites are compared. Split kmers have been implemented in many software packages including KSNP and SKA [@Gardner:2015; @Harris:2018].
This strategy essentially eliminates the assembly steps required in a genome-wide approach and allows for comparative whole genome analysis based on short fragments, usually from 13 to 31 base pairs at a time.

Several bioinformatics applications have been released which utilize the kmer-based strategy to convert next generation sequencing data into manageable datasets, usually called sketches [@Ondov:2016; @Baker:2019; @Zhao:2018].
Most notably, an algorithm called min-hash was implemented in the Mash package [@Ondov:2016].  In the min-hash algorithm, all nucleotides of length _k_ (kmers) are recorded and transformed into integers using an algorithm called hashing.  These hashed kmers are sorted and only the first several kmers are retained.  The kmers that appear at the top of the sorted list are collectively called the sketch.  Sketches are a representation of an organism's genome, and comparing multiple sketches can be used to estimate relatedness.
In other words, two genomes that have been sketched with the min-hash algorithm can be compared by counting how many hashed kmers they have in common. 

A dendrogram can be calculated from pairwise genome distances by several algorithms [@Kuhner:1994]. One widely used algorithm is neighbor-joining (NJ) [@Saitou:1987].  The NJ algorithm reads a matrix of pairwise genomic distances. Each genome starts as a leaf node in a star-like tree.  In each iteration of NJ, the most related genomes or nodes (neighbors) are identified by seeking the lowest distance value in the matrix. Neighbors are given a common node on the growing tree.  Then, new distances are calculated between all nodes and the new node.  The product of the NJ algorithm is a binary tree, where each node except the root has one parent node and zero or two children.

Because min-hash creates distances between any two genomes and NJ uses distances to create trees, min-hash values can be used to rapidly cluster genomes into trees.  We implemented this idea in software called Mashtree, which quickly and efficiently generates large trees that would be computationally intensive using other methods.

# Implementation

## Base workflow

Mashtree builds on two major algorithms that are already implemented in other software packages. The first is the min-hash algorithm, which is implemented in the software Mash [@Ondov:2016].  Mashtree uses Mash to create sketches of the genomes with the function `mash sketch`.  We elected to keep most default Mash parameters but increased the sketch size (number of hashed kmers) from 1,000 to 10,000 to increase discriminatory power.  Then, Mash is used to calculate the distances between genomes with `mash dist`.  Mashtree records these distances into a pairwise distance matrix. Next, Mashtree calls the neighbor-joining (NJ) algorithm which is implemented in the software QuickTree [@Howe:2002].  The Mash distance matrix is used with QuickTree with default options to generate a dendrogram.  The workflow is depicted in [Figure 1].

## Fast or accurate mode

While Mashtree has been developed for speed, some might want higher accuracy.  In response, we have implemented an optional minimum abundance filter for raw read input files.  This filter creates a histogram of kmers using Mash and detects the valleys between peaks.  The first (left-most) peak of a kmer histogram represents kmers that occur only once or rarely.  These "rare" kmers usually represent sequencing errors.  Therefore, the minimum abundance filter can optionally be used to see where the first valley is and automatically set that value for the minimum kmer coverage parameter.  Using this filter helps avoid noise in the input data.  Running in accurate mode is shorthand for running Mashtree with the minimum abundance filter turned on.  Otherwise by default, a user is running in fast mode.

## Confidence values

Most phylogenetic trees have confidence scores associated with each hypothetical ancestor node.
These support scores are typically generated by resampling the input data many times, analyzing each replicate in context of the original phylogeny.
Strong  support for a node in the input data will result in consistent recovery of that node in the majority of resampled replicate trees.  Conversely, nodes which have low support in the input data will be sensitive to stochastic changes in the resampled sketches which can result in alternate topologies of involved isolates.  Therefore, the support value correlates with the robustness and reproducibility of a recovered relationship.  As the NJ tree recovered from a single tree reconstruction will always be a fully resolved binary tree, the addition of these resampled support values provides a valuable tool for assessing the confidence and reliability of a node.  This can be especially helpful when trying to assess the exact topology of clades containing  three or more very closely related genomes.  

Although Mashtree does not infer phylogeny, we have borrowed the ideas behind phylogenetic confidence values to yield confidence values for each parent node in the tree. There are two resampling methods implemented in Mashtree to assign support values to infernal nodes, bootstrapping and jackknifing. Initially, both methods create a tree as depicted in [Figure 1]. 

The bootstrapping method infers confidence on the tree by creating one new tree per replication ([Figure 2]). However, instead of accepting the default seed value for Mash, a new one is generated.
Changing this seed value will affect the value of hashed kmers and how they sort, therefore providing a different sampling of kmers in the sketch.
Due to changes in sketch sampling, distances between genomes can differ and the resulting tree in the replication can differ. Then, nodes in the original tree are assigned support values based on how frequently the same bipartition was recovered in the replicate trees.

In the jackknife method, for each replication, each genome is picked one at a time as the query genome ([Figure 3]).  Its hashes are sampled without replacement for half of the original hashes.  Then, a Mash distance is recomputed between the query genome and all other genomes in the analysis.
When the distances are recorded in the replicate's pairwise distance matrix, the average distance between each comparison is used of each query-reference genome pair, and vice versa. This new distance matrix is used for computing a NJ tree. Nodes in the original tree are assigned support values based on how frequently the same bipartition was recovered in replicate trees.

## Other features

Mashtree has several other useful features.  First, Mashtree can read any common sequence file type and can read gzip-compressed files (e.g., fastq, fastq.gz, fasta).  This is a major advantage in being compatible with a wide variety of databases and with space-saving file compression.  Second, Mashtree takes advantage of multithreading.  The number of requested threads is used to determine how many genomes are sketched at the same time and how many sketches can be compared at the same time.  When the number of threads requested outnumbers the number of operations that it can parallelize, Mashtree uses the multithreading already encoded in Mash sketches and distances.  Third, Mashtree uses an SQLite database which can be used to cache results between runs.

# Installation

The Mashtree package is programmed in Perl, and is available in the CPAN repository.  Documentation can be found at https://github.com/lskatz/mashtree.  

# Figures

![Figure 1. The Mashtree workflow.  Step 1) Sketch genomes with Mash. In this schematic, there is a green circle representing each genome in the analysis.  Filled-in brown circles indicate the presence of a kmer.  Missing circles represent true absence.  After hashing with a sketch size of six (after the arrow), some kmers are not represented in the Mash sketch either because they are not present in the original genome or because only a finite number of kmers are sketched (e.g., six in this example).  Henceforth, truly missing hashes or hashes not included in the Mash sketch are represented by empty circles.  Step 2) Calculate distances with `Mash dist`.  Distances in the figure are represented by Jaccard distances, which are calculated as the intersection divided by the union.  In this example, the genomes are separated by Jaccard distances of 5/9, 4/9, and 3/9.  These jaccard distances are internally transformed into Mash distances [@Ondov:2016].  Step 3) Create dendrogram with Quicktree using the Mash distance matrix.  ](Mashtree_workflow.png)  


![Figure 2. The Mashtree bootstrap workflow.  Step 1) Generate a tree with the normal workflow [Figure 1](Figure 1). This is the main tree.  Step 2) Run the normal workflow once per replicate but with a different random seed. In this example, the top right replicate differs from the main tree.  All ten of these trees are the bootstrap tree replicates.  Step 3) For each parent node in the main tree, quantify how many bootstrap tree replicates have the same node with the same children. Record that percentage next to each parent node. This percentage quantifies how confident the Mashtree cluster is, controlling for the random seed in the Mash program.  ](Bootstrap_workflow.png)  

![Figure 3.  The Mashtree jackknife workflow.  Step 1) Generate a tree with the normal workflow [Figure 1](Figure 1). This is the main tree.  Step 2) For each replicate, sample the half hashes without replacement for each query genome. Recalculate the Mash distance between the query genome and all other genomes, reducing the denominator to one half, rounding up, to reflect the smaller pool of hashes.  After all genomes have been selected for query genomes, average the distances to create a new distance matrix.  Create the dendrogram from the new distance matrix.  For brevity, only one detailed replicate is shown.  Step 3) For each replication, calculate the new tree from the new distance matrix. In this example, the top right replication differs from the main tree. All ten of these trees are the jack knife tree replicates.  Step 4) For each parent node in the main tree, quantify how many jack knife tree replicates have the same node with the same children. Record that percentage next to each parent node.  This percentage quantifies how confident Mashtree is at clustering, controlling for stochasticity in hashes.  ](Jackknife_workflow.png)

# Acknowledgements

This work was made possible through support from the Advanced Molecular Detection (AMD) Initiative at the Centers for Disease Control and Prevention.  Thank you Sam Minot, Andrew Page, Brian Raphael, and Torsten Seemann for helpful discussions.
The findings and conclusions in this report are those of the authors and do not necessarily represent the official position of the Centers for Disease Control and Prevention.

# References
