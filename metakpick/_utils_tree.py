
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



def get_tax2root_all(info, parents, tax_genome):
    tax2root_all=[]
    tax2root_all_dic={}
    # paths endint at all levels 
    #tax_index_list= list(tax_index) 
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
    return tax2root_all, tax2root_all_dic


def get_tax2path_all(tax2root_all):
        # the following is for all tax at all levels
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




def get_read_tax_depth(merged, info, parents, cases):
        
    read_tax_depth={}
    for case in cases: #'k19'
        read_tax_depth[case]={}
        tocheck_=""
        for index, row in merged.iterrows():
            taxid=row['taxon_tool_'+case]
            tax2root_list = find_tax2root(info,parents, taxid) # [2711156, 1649486, 41297, 204457, 28211, 1224, 3379134, 2, 131567, 1]
            if tax2root_list!=-1:
                depth=len(tax2root_list)
            else:
                depth=0
            if depth<3:
                tocheck_=taxid
            read_tax_depth[case][row['read_name']]=depth
            
    print(tocheck_)
    print(len(read_tax_depth),len(read_tax_depth[case]))
    #Counter(read_tax_depth[case].values())
    return read_tax_depth
    