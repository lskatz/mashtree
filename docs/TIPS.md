# Tips and tricks

## Leveraging a high performance computer for bootstrapping

Bootstrapping is simply running mashtree several times and then quantifying how often a parent node
occurs for each group of genomes.  This frequency is usually displayed as a percentage
at each node.

### Bootstrap trees

Here is a bash loop to make 100 bootstrap trees

    mkdir bootstrapTrees
    export FASTA=$(ls *.fasta)
    # Run the first tree without using --seed, to form your actual tree
    # without bootstraps, using many threads so that it can be done before
    # any of the bootstrap trees.
    echo 'mashtree --numcpus $NSLOTS --mindepth 0 $FASTA' | \
      qsub -V -cwd -pe smp 1-36 -N mashtree -o mashtree.dnd -e mashtree.log
    
    # Make 100 bootstrap trees
    for i in `seq 1 100`; do 
      # The echo statement makes a short bash script, which
      # is piped to qsub
      echo '
        mashtree --numcpus $NSLOTS --mindepth 0 --seed $RANDOM $FASTA' 
      |\
        qsub -V -cwd -pe smp 1-12 -N mashtree$i -o bootstrapTrees/$i.dnd -e bootstrapTrees/$i.log; 
    done;
    
### Combine bootstrap trees

#### Bootstrapping with RAxML

RAxML is not a prerequisite for Mashtree, and so it might not be installed on your computer even if Mashtree is working.

    cat bootstrapTrees/*.dnd > bootstrapTrees.dnd
    raxmlHPC -f b -t mashtree.dnd -z bootstrapTrees.dnd -m GTRCAT -n TEST
    # output file: RAxML_bipartitions.TEST
    ln -sv RAxML_bipartitions.TEST bsTree.dnd # in case you want to have a sane file extension

#### Bootstrapping with BioPerl

This script requires bioperl, which is a prerequisite for Mashtree.

    perl -MBio::TreeIO -MBio::Tree::Statistics -e '
      # Gather all the trees into an array
      for my $file(glob("bootstrapTrees/*.dnd")){
        next if(! -s $file);
        my $in=Bio::TreeIO->new(-file=>$file);
        while($tree = $in->next_tree){
          push(@tree, $tree);
        }
      }
      
      # Combine trees
      my $baseTree = Bio::TreeIO->new(-file=>"mashtree.dnd")->next_tree;
      $stats=Bio::Tree::Statistics->new; 
      my $bsTree = $stats->assess_bootstrap(\@tree,$baseTree);
      # BioPerl is very respectful about which values are IDs vs which are
      # bootstraps. However, most tree drawing programs look at IDs, so
      # we have to alter the ID
      for ($bsTree->get_nodes){
        next if($_->is_Leaf); # Do not alter leaves
        $_->id($_->bootstrap);
      }
      
      print $bsTree->as_text("newick") ."\n";
    ' > bsTree.dnd
    
