import pandas as pd

def parse_tree_file(tree_df,children,parents):
    # create a dictionary from taxonomic id to index of corresponding row
    info = {}
    for i in range(len(children)):
        info[int(children[i])] = int(i)

    # Convert to tree structure (dictionary of sets)
    Tree = {}
    for i in range(len(tree_df[0])):
        parent = int(parents[i])
        child = int(children[i])
        if parent in Tree:
            Tree[parent].add(child)
        else:
            Tree[parent] = {child}
    return info, Tree 


def find_tax2root(info,parents, tid):
    if tid not in info:
        #print(tid)
        return -1
    list_tax2root=[]
    while tid != 1:    
        list_tax2root.append(tid)
        tid = int(parents[info[tid]])
    list_tax2root.append(tid)
    return list_tax2root



def get_tax_info(tree_file,tax_genome_file):
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

    return info, Tree, tax_index, tax_genome, parents


def get_tax2path(tax_genome, info, parents):

    tax2root_all=[]
    tax2root_all_dic={} # paths endint at all levels 
    # paths endint at leaf level (no Vague positive )
    tax_index_list= list(tax_genome)  # tax_genome_specieslevel 
    for tax in tax_index_list:   
        tax2root= find_tax2root(info, parents, tax)
        if tax2root==-1:
            print("tax",tax)
        else:    #elif 2 in tax2root: # selecting only bacteria
            tax2root_all.append(tax2root) # # [655353, 655352, 655351, 356, 28211, 1224, 3379134, 2, 131567, 1]
            tax2root_all_dic[tax]=tax2root
    print("number of paths",len(tax2root_all))
    tax2path = {} # a tax is present in which paths
    tax2depth = {} 
    for tax2root_path_idx, tax2root_path in enumerate(tax2root_all):
        for tax_idx,tax in enumerate(tax2root_path[::-1]):
            tax2depth[tax]=tax_idx
            if tax in tax2path:
                tax2path[tax].append(tax2root_path_idx) #tax2root_path[0]
            else:
                tax2path[tax]=[tax2root_path_idx]
    print(len(tax2path),len(tax2depth))

    return tax2path, tax2depth