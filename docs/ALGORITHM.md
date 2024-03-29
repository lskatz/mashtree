# Algorithms

## Fast or accurate mode

While Mashtree has been developed for speed, some might want higher accuracy.
In response, we have implemented an optional minimum abundance filter for raw read input files.
This filter creates a histogram of kmers using Mash and detects the valleys between peaks.
The first (left-most) peak of a kmer histogram represents kmers that occur only once or rarely.
These "rare" kmers usually represent sequencing errors.
Therefore, the minimum abundance filter can optionally be used to see where the first valley is and automatically set that value for the minimum kmer coverage parameter.
Using this filter helps avoid noise in the input data.
Running in accurate mode is shorthand for running Mashtree with the minimum abundance filter turned on.
Otherwise by default, a user is running in fast mode.

To invoke accurate mode, run with `--mindepth 0`, e.g.,

```bash
mashtree --mindepth 0 --numcpus 12 *.fastq.gz [*.fasta] > mashtree.dnd
```

## Confidence values

Most phylogenetic trees have confidence scores associated with each hypothetical ancestor node.
These support scores are typically generated by resampling the input data many times, analyzing each replicate in context of the original phylogeny.
Strong  support for a node in the input data will result in consistent recovery of that node in the majority of resampled replicate trees.  Conversely, nodes which have low support in the input data will be sensitive to stochastic changes in the resampled sketches which can result in alternate topologies of involved isolates.  Therefore, the support value correlates with the robustness and reproducibility of a recovered relationship.  As the NJ tree recovered from a single tree reconstruction will always be a fully resolved binary tree, the addition of these resampled support values provides a valuable tool for assessing the confidence and reliability of a node.  This can be especially helpful when trying to assess the exact topology of clades containing  three or more very closely related genomes.

Although Mashtree does not infer phylogeny, we have borrowed the ideas behind phylogenetic confidence values to yield confidence values for each parent node in the tree. There are two resampling methods implemented in Mashtree to assign support values to infernal nodes, bootstrapping and jackknifing. Initially, both methods create a tree as depicted in [Figure 1](../paper/Mashtree_workflow.png).

The bootstrapping method infers confidence on the tree by creating one new tree per replication [Figure 2](../paper/Bootstrap_workflow.png). However, instead of accepting the default seed value for Mash, a new one is generated.
Changing this seed value will affect the value of hashed kmers and how they sort, therefore providing a different sampling of kmers in the sketch.
Due to changes in sketch sampling, distances between genomes can differ and the resulting tree in the replication can differ. Then, nodes in the original tree are assigned support values based on how frequently the same bipartition was recovered in the replicate trees.

In the jackknife method, for each replication, each genome is picked one at a time as the query genome [Figure 3](../paper/Jackknife_workflow.png).  Its hashes are sampled without replacement for half of the original hashes.  Then, a Mash distance is recomputed between the query genome and all other genomes in the analysis.
When the distances are recorded in the replicate's pairwise distance matrix, the average distance between each comparison is used of each query-reference genome pair, and vice versa. This new distance matrix is used for computing a NJ tree. Nodes in the original tree are assigned support values based on how frequently the same bipartition was recovered in replicate trees.

# Figures

![Figure 1](../paper/Mashtree_workflow.png)  

**Figure 1.** The Mashtree workflow.  Step 1) Sketch genomes with Mash. In this schematic, there is a green circle representing each genome in the analysis.  Filled-in brown circles indicate the presence of a kmer.  Missing circles represent true absence.  After hashing with a sketch size of six (after the arrow), some kmers are not represented in the Mash sketch either because they are not present in the original genome or because only a finite number of kmers are sketched (e.g., six in this example).  Henceforth, truly missing hashes or hashes not included in the Mash sketch are represented by empty circles.  Step 2) Calculate distances with `Mash dist`.  Distances in the figure are represented by Jaccard distances, which are calculated as the intersection divided by the union.  In this example, the genomes are separated by Jaccard distances of 5/9, 4/9, and 3/9.  These jaccard distances are internally transformed into Mash distances [@Ondov:2016].  Step 3) Create dendrogram with Quicktree using the Mash distance matrix.

![Figure 2](../paper/Bootstrap_workflow.png)  

**Figure 2.** The Mashtree bootstrap workflow.  Step 1) Generate a tree with the normal workflow as in Figure 1. This is the main tree.  Step 2) Run the normal workflow once per replicate but with a different random seed. In this example, the top right replicate differs from the main tree.  All ten of these trees are the bootstrap tree replicates.  Step 3) For each parent node in the main tree, quantify how many bootstrap tree replicates have the same node with the same children. Record that percentage next to each parent node. This percentage quantifies how confident the Mashtree cluster is, controlling for the random seed in the Mash program.

![Figure 3](../paper/Jackknife_workflow.png)  

**FIgure 3.** The Mashtree jackknife workflow.  Step 1) Generate a tree with the normal workflow as in Figure 1. This is the main tree.  Step 2) For each replicate, sample the half hashes without replacement for each query genome. Recalculate the Mash distance between the query genome and all other genomes, reducing the denominator to one half, rounding up, to reflect the smaller pool of hashes.  After all genomes have been selected for query genomes, average the distances to create a new distance matrix.  Create the dendrogram from the new distance matrix.  For brevity, only one detailed replicate is shown.  Step 3) For each replication, calculate the new tree from the new distance matrix. In this example, the top right replication differs from the main tree. All ten of these trees are the jack knife tree replicates.  Step 4) For each parent node in the main tree, quantify how many jack knife tree replicates have the same node with the same children. Record that percentage next to each parent node.  This percentage quantifies how confident Mashtree is at clustering, controlling for stochasticity in hashes.


