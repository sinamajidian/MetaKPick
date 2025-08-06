

import argparse
import sys
import os
import logging
import pandas as pd
from _utils import parse_tree_file, find_tax2root



def get_tax_index(tree_file, tax_genome_file):

    tree_file="/vast/blangme2/smajidi5/metagenomics/changek/kraken1/kraken_db/standard/k31/taxonomy/nodes.dmp"
    tax_genome_file="/vast/blangme2/smajidi5/metagenomics/changek/kraken1/kraken_db/standard/seqid2taxid.map_tax_uniq"
    # tax_genome_file="/vast/blangme2/smajidi5/metagenomics/changek/kraken1/kraken_db/48genomes/added_header_tax_uniq"
    # tree_file= "/vast/blangme2/smajidi5/metagenomics/changek/kraken1/kraken_db/48genomes/k31/taxonomy/nodes.dmp"

    print("Tree file: ", tree_file)
    tree_df = pd.read_csv(tree_file, sep="\t", header=None)
    print("Read the tree file") 
    children = tree_df[0]
    parents = tree_df[2]
    info, Tree = parse_tree_file(tree_df,children,parents)
    print("Parsed the tree file. Size of info: ", len(info),"Size of Tree: ", len(Tree))

    tax_genome_f = open(tax_genome_file,'r') # genomes in the kraken index  
    tax_genome = set()
    for line in tax_genome_f:
        tax_genome.add(int(line.strip()))
    print('Number of genomes/strain ',len(tax_genome),"in kraken index")

    tax_index_=[]
    tax_notfound=[]
    for tax_ in tax_genome:
        list_tax2root = find_tax2root(info,parents, tax_)
        if list_tax2root!=-1:
            tax_index_ +=list_tax2root
        else:
            tax_notfound.append(tax_)
    tax_index =set(tax_index_) # all the tax in genomes of kraken, plus all up to root 
    print('Number of tax ids',len(tax_index),len(tax_index_)) # 
    print("Genome tax not found in kraken tree",tax_notfound)
    
    return tax_index



def main():
    parser = argparse.ArgumentParser(description='Metakpick: for metagenomic classification of reads')
    
    # Main operation mode arguments
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--train', action='store_true', help='Run the training of the model')
    mode_group.add_argument('--classify', action='store_true', help='Run the classification of the model')


    parser.add_argument('--reads', type=str, help='Input reads file')
    parser.add_argument('--output', type=str, help='Output file')
    #parser.add_argument('--model', type=str, help='Model file')
   # parser.add_argument('--threads', type=int, help='Number of threads')



if __name__ == "__main__":
    main()
