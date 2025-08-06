
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
    
def find_tax_level(info,tree_df,parents, tid, tax_level="species"):
    if tid not in info:
        #print(tid)
        return -1
    while tid != 1:    
        tid_row = info[tid]
        if tree_df[4][tid_row] == tax_level:
            break
        else:
            tid = parents[info[tid]] #find_parent_node(tid)
            # if retrun 1, it means that the tid is more towaards root that the requested level
    return tid


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



def lca_finder_list(info,parents, tax_list):
    lca=tax_list[0]
    
    for tax2 in tax_list[1:]:
        lca = lca_finder(info,parents, lca,tax2)
        
    return lca

def lca_finder(info,parents, tax1,tax2):
    list_tax2root_1 = find_tax2root(info,parents, tax1)
    list_tax2root_2 = find_tax2root(info,parents, tax2)
    
    list_tax2root1 =list_tax2root_1 [::-1] # reverse from root to leaf
    list_tax2root2 =list_tax2root_2 [::-1]
    
    depth_i = 0
    while(depth_i < len(list_tax2root1) and depth_i < len(list_tax2root2)):
        #print(depth_i)
        if list_tax2root1[depth_i] != list_tax2root2[depth_i]:
            break
        depth_i += 1
    assert list_tax2root1[depth_i-1] == list_tax2root2[depth_i-1], "issue in LCA finding"
    lca=list_tax2root1[depth_i-1]
    #print(lca)
    return lca

    
def first_nonzero(lst):
    for num in lst:
        if num != 0:
            return num
    return None  # In case all elements are zero


def first_nonzero_idx(lst):
    for idx, num in enumerate(lst):
        if num != 0:
            return idx
    return None  # In case all elements are zero    

def get_tax_index(tree_file, tax_genome_file):
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