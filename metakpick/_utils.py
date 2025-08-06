
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
