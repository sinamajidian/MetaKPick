


import logging
import numpy as np
import _utils_tree
import _utils_kraken


def get_features_read(tax_kmer_num_dic, num_nodes_tree,kmer_reported_tax,rlen):
    
    # tax_kmer_num_dic[tax]+=num
    num_kmer_nonzero=np.array(list(tax_kmer_num_dic.values()),dtype=np.int32)
    num_kmer_norm=num_kmer_nonzero/rlen # read_length_dic[read_id]

    mean_exc_tax = (np.sum(num_kmer_nonzero)-kmer_reported_tax)/num_nodes_tree
    
    #feature_names=['mean_all','mean_nonzero','max','sum', 'mean_exc_tax']
    features_read=[np.mean(num_kmer_nonzero), np.mean(num_kmer_nonzero)*len(num_kmer_nonzero)/num_nodes_tree, np.max(num_kmer_nonzero),np.sum(num_kmer_nonzero),mean_exc_tax]
    #feature_names+=['mean_all/rlen','mean_nonzero/rlen','max/rlen','sum/rlen', 'mean_exc_tax/rlen'] # /readlength
    features_read+=[np.mean(num_kmer_norm)/rlen, np.mean(num_kmer_norm)*len(num_kmer_norm)/num_nodes_tree/rlen, 
                np.max(num_kmer_norm)/rlen,np.sum(num_kmer_norm)/rlen, mean_exc_tax/rlen]
    return features_read



def get_features_depth(case, read_id, read_tax_depth, tax2depth,reported_tax,tax_kmer_dic):
    if case in read_tax_depth and read_id in read_tax_depth[case]:
        txdepth_norm= read_tax_depth[case][read_id]
    else:
        txdepth_norm=0
    
    depth_reported= tax2depth[reported_tax]
    
    numkmer_consecutive=[]
    weighted_depth_numkmer =[]
    for tax, pos_num in tax_kmer_dic.items():
        numkmer_consecutive += [i[1] for i in pos_num] # how many kmer are one after each other (~ length of math)
        if tax in tax2depth:
            tax_depth=tax2depth[tax] 
            weighted_depth_numkmer += [i[1]*tax_depth for i in pos_num]
    if weighted_depth_numkmer:
        numkmer_consecutive_mean=np.mean(numkmer_consecutive)
        weighted_depth_numkmer_mean=np.mean(weighted_depth_numkmer)
    else:
        numkmer_consecutive_mean= 0 
        weighted_depth_numkmer_mean=0
        
    features_depth=[depth_reported,numkmer_consecutive_mean,weighted_depth_numkmer_mean ]
    #feature_names+=['depth_reported_tax','Avg_kmer_consecutive', 'WAvg_kmer_consecutive']
    return features_depth


def get_features_tax(reported_tax, tax2root_dic, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen):
    kmer_reported_uppertax=0
    kmer_reported_belowtax=0
    kmer_other_reported_tax=0
    kmer_other_reported_uppertax=0 # towards root (smaller depth value)
    kmer_other_reported_belowtax=0 # towards leaves
    
    tax_krak_2root=tax2root_dic.get(reported_tax, [reported_tax,1]) #_utils_tree.find_tax2root(info, parents, reported_tax)
    if tax_krak_2root==-1:
        kmer_reported_uppertax=0
    else:   
        for tax_upper in tax_krak_2root[1:]: #excluding tax_krak towards the root  tax_krak_2root = [201174, 1783272, 2, 131567, 1] where tax_krak=201174 
            if tax_upper in tax_kmer_num_dic: 
                kmer_reported_uppertax +=tax_kmer_num_dic[tax_upper]
        
    below_tax_all_perread=[]
    for taxi in tax_kmer_num_dic:
        # to get the below for the reported 
        tax2root_t= tax2root_dic.get(taxi, [taxi,1]) #_utils_tree.find_tax2root(info, parents, taxi)
        if tax2root_t!=-1 and reported_tax in tax2root_t: # tax_krak is lca of taxi
            kmer_reported_belowtax += tax_kmer_num_dic[taxi]   
            below_tax_all_perread.append(taxi)
    
        # to get the below/above/at the reported tax across all (not just max path)
        taxi_depth = tax2depth.get(taxi)
        if taxi_depth is None:
            continue
        if depth_reported == taxi_depth:
            kmer_other_reported_tax += tax_kmer_num_dic[taxi]
        elif depth_reported < taxi_depth:  # towards leaves (deeper)
            kmer_other_reported_belowtax += tax_kmer_num_dic[taxi]
        else:  # depth_reported > taxi_depth: towards root (shallower)
            kmer_other_reported_uppertax += tax_kmer_num_dic[taxi]
    
    features_tax =[kmer_reported_tax, kmer_reported_uppertax,kmer_reported_belowtax,
        kmer_reported_tax/rlen, kmer_reported_uppertax/rlen,kmer_reported_belowtax/rlen,
        kmer_other_reported_tax, kmer_other_reported_uppertax,kmer_other_reported_belowtax,
        kmer_other_reported_tax/rlen, kmer_other_reported_uppertax/rlen,kmer_other_reported_belowtax/rlen]
    
    #feature_names+=['kmer_reported_tax','kmer_tax_above','kmer_tax_below','kmer_tax/rlen','kmer_tax_above/rlen','kmer_tax_below/rlen','kmer_othertax','kmer_othertax_above','kmer_othertax_below','kmer_othertax/rlen','kmer_othertax_above/rlen','kmer_othertax_below/rlen',]
    return features_tax, below_tax_all_perread

def get_features_half_read(reported_tax, tax_kmer_dic, rlen, below_tax_all_perread):
    tax_krak_touse=reported_tax
    if reported_tax not in tax_kmer_dic: # when lca is reprorted, use onf its lower ones. 
        tax_krak_touse= below_tax_all_perread[0]
    segment_num=2 # 3
    segment_len=int(rlen/segment_num)
    pos_num_reportedtax = tax_kmer_dic[tax_krak_touse] # [(77, 1), (136, 4), (139, 2), (144, 4), (148, 2), (162, 2), (233, 1)]
    cnt_perbin=np.zeros((segment_num,1))
    for pos,numkmer in pos_num_reportedtax:
        bin_idx= int(pos/segment_len)
        cnt_perbin[bin_idx]+=numkmer
    diff_fromnext_seg = [ np.abs(cnt_perbin[segment_i+1]-cnt_perbin[segment_i])[0] for segment_i in range(segment_num-1)]
    
    feature_half_read = [np.mean(diff_fromnext_seg)]
    return feature_half_read




# def traverse_from_root(current, tree, num_kmers_path, tax_kmer_num_dic, child2parent): # current is a node of tree
#     #print(current)  
#     num_kmers_node = tax_kmer_num_dic.get(current,0)
#     if current in child2parent:
#         parent = child2parent[current]
#         num_kmers_path_parent= num_kmers_path[parent]
#         if num_kmers_node:
#             last_nonzero_parent= parent
#             if num_kmers_path_parent[0]!=0: # if parent kmer is nonzero 
#                 last_nonzero_parent= num_kmers_path_parent[1]
#             num_kmers_path[current] = (num_kmers_path[last_nonzero_parent][0] + num_kmers_node, current)

#             # if num_kmers_path[parent][0]: # if parent kmer is nonzero 
#             #     num_kmers_path_upto_parent = num_kmers_path[parent][0]
#             # else:
#             #     last_nonzero_parent= num_kmers_path[parent][1]
#             #     num_kmers_path_upto_parent= num_kmers_path[last_nonzero_parent][0] #print('here',current, parent, last_nonzero_parent)
#             # num_kmers_path[current] = (num_kmers_path_upto_parent + num_kmers_node, current)
#             #print('current:', current,' parent:',parent,' num_parrent:', num_kmers_path[parent][0],' num_current_path:', num_kmers_path[current])
#         else:
#             last_nonzero_parent=num_kmers_path_parent[1]
#             num_kmers_path[current] = (0, last_nonzero_parent)
        
#     else:  # current node is root
#         num_kmers_path[current] = (num_kmers_node, current)

#     if current in tree: 
#         children = tree[current]
#         for child in children:
#             traverse_from_root(child, tree, num_kmers_path,tax_kmer_num_dic, child2parent)
#     # else: current is a leaf, no child 
#     return 1



def traverse_from_root(current, tree, num_kmers_path, tax_kmer_num_dic, child2parent): # current is a node of tree
    num_kmers_node = tax_kmer_num_dic.get(current,0)
    parent = child2parent[current]
    num_kmers_path_parent0,num_kmers_path_parent1= num_kmers_path[parent]
    if num_kmers_node:  # current node has kmers assinged
        if num_kmers_path_parent0: # if kmers of parent and current are nonzero 
            num_kmers_path[current] = (num_kmers_path_parent0 + num_kmers_node, current)
        else: #  kmers of current is nonzero and parent is zero 
            num_kmers_path[current] = (num_kmers_path[num_kmers_path_parent1][0]+ num_kmers_node, current)
    else: #  kmers of current and parent are zero 
        num_kmers_path[current] = (0, num_kmers_path_parent1) # the first item is arbitary, as a way to be checked for zero nodes
    if current in tree: 
        children = tree[current]
        for child in children:
            traverse_from_root(child, tree, num_kmers_path,tax_kmer_num_dic, child2parent)
    # else: current is a leaf, no child 
    #return 1


def calculate_num_kmers_path(Tree_updated, tax_kmer_num_dic,child2parent):
    
    # logging.debug('number of leaves',len(child2parent))
    root=1
    num_kmers_path={root:(tax_kmer_num_dic.get(root,0), root)}
    num_kmers_path[0]=(0,0)
    child2parent[1]=0
    traverse_from_root(root, Tree_updated, num_kmers_path,tax_kmer_num_dic, child2parent)
    #print('number of nodes',len(num_kmers_path))

    num_kmers_path_updated ={node:accum for node, (accum,nn) in num_kmers_path.items() if accum!=0}

    #len(Tree), len(tax_kmer_num_dic), len(num_kmers_path), len(num_kmers_path_updated)
    # Tree2={1:{2,3,8},2:{4,5},3:{6,7},8:{11,12}}
    # child2parent={}
    # for parent, children in Tree2.items():
    #     for child in children:
    #         child2parent[child]=parent
    # print('number of leaves',len(child2parent))

    # num_kmers_path={}
    # tax_kmer_num_dic={i:i*10 for i in range(1,13)}
    # tax_kmer_num_dic[11]=0
    # tax_kmer_num_dic[12]=0

    return num_kmers_path_updated




def get_features_path(path_feature, Tree_updated, reported_tax, tax_kmer_num_dic, rlen,num_nodes_tree,info, parents,child2parent):
    if not path_feature:
        return [0]*9
    num_kmers_path_updated= calculate_num_kmers_path(Tree_updated, tax_kmer_num_dic,child2parent)
    num_kmer_path_all=np.array(list(num_kmers_path_updated.values()),dtype=np.int32)
    num_kmer_norm=num_kmer_path_all/rlen # read_length_dic[read_id]
    if reported_tax in num_kmers_path_updated:
        kmer_reported_tax =  num_kmers_path_updated[reported_tax]
    else:
        kmer_reported_tax= 0 # todo uncomment the below
        #path_reported_tax= _utils_tree.find_tax2root(info, parents, reported_tax)
        #kmer_reported_tax = np.sum([ tax_kmer_num_dic.get(tax, 0) for tax in path_reported_tax])
        
    mean_exc_tax = (np.sum(num_kmer_path_all)-kmer_reported_tax)/num_nodes_tree # todo is num_nodes_tree different than numebr of paths
    #feature_names=['mean_nonzero','mean_all','max','sum', 'mean_exc_tax']
    features_read_path=[np.mean(num_kmer_path_all), np.mean(num_kmer_path_all)*len(num_kmer_path_all)/num_nodes_tree,
                       np.max(num_kmer_path_all), np.sum(num_kmer_path_all), mean_exc_tax]
    #feature_names+=['mean_nonzero/rlen','max/rlen','sum/rlen', 'mean_exc_tax/rlen'] # /readlength
    features_read_path+=[np.mean(num_kmer_norm), np.mean(num_kmer_norm)*len(num_kmer_norm)/num_nodes_tree, 
                np.max(num_kmer_norm), np.sum(num_kmer_norm)]

    return features_read_path




def parse_fastq_reads(fastq_file):
    reads_quality = {}
    handle=open(fastq_file, 'r')
    line_count = 0
    for line in handle:
        line = line.strip()
        line_count += 1
        if line_count % 4 == 1:
            read_id=line[1:].strip()
            if not line.startswith('@'):
                raise ValueError(f"Invalid FASTQ format: Expected header line starting with '@', got: {line}")
            #current_read = {'header': line[1:], 'sequence': '', 'quality_header': '', 'quality': ''}                    
        elif line_count % 4 == 0:
            quality_scores = [ord(q) - 33 for q in line]  # Convert to Phred scores
            reads_quality[read_id]= quality_scores
       
    return reads_quality

def get_features_readq(reads_quality):
    # reads_quality is a dictionary of read_name and the read_quality list (integer values )
    reads_quality_features={}
    if reads_quality:
        for read_name, read_quality in reads_quality.items():
            reads_quality_features[read_name]= [np.mean(read_quality), np.median(read_quality), np.percentile(read_quality, 10), np.percentile(read_quality, 90)]
    
    return reads_quality_features




def get_features_all(read_names_list, tax2path, kraken_kmers_cases, read_tax_depth, tax2depth, info, parents, Tree, tax_index, reads_quality_features,path_feature):
    
    Tree_updated ={node:node_set.intersection(tax_index) for node,node_set in  Tree.items() if  node in tax_index}
    if 1 in Tree_updated and 1 in Tree_updated[1]:
        Tree_updated[1].remove(1)
    child2parent = { }
    for parent, children in Tree_updated.items():
        for child in children:
            child2parent[child]=parent

    tax2root_dic={}
    for taxii in tax_index:
        tax2root_dic[taxii] = _utils_tree.find_tax2root(info, parents, taxii)



    feature_names = ['mean_nonzero', 'mean_all', 'max', 'sum', 'mean_exc_tax', 'mean_all/rlen', 'mean_nonzero/rlen', 'max/rlen', 'sum/rlen', 'mean_exc_tax/rlen', 
                'depth_reported_tax', 'Avg_kmer_consecutive', 'WAvg_kmer_consecutive', 'kmer_reported_tax', 'kmer_tax_above', 'kmer_tax_below', 
                'kmer_tax/rlen', 'kmer_tax_above/rlen', 'kmer_tax_below/rlen', 'kmer_othertax', 'kmer_othertax_above', 'kmer_othertax_below', 
                'kmer_othertax/rlen', 'kmer_othertax_above/rlen', 'kmer_othertax_below/rlen', 'diff#kmers_halfRead','rlen',
                'path_mean_nonzero','path_mean_all','path_max','path_sum', 'path_mean_exc_tax',
                'path_mean_nonzero/rlen','path_mean_all/rlen','path_max/rlen','path_sum/rlen',
                'mean_readq','median_readq','readq_10p','readq_90p']
    logging.debug("feature_names are : "+str(feature_names))
    
    num_nodes_tree = len(tax2path)    # 52229
    cases=list(kraken_kmers_cases.keys())
    features_cases={}
    not_found_in_tax2depth=[]
    for case in cases: 
        logging.debug("Extracting features for case: "+case)
        kmer_size= int(case[1:])
        X3=[]
        for read_idx,read_name in enumerate(read_names_list):
            if read_idx%10000==0:
                logging.debug("Extracting features for read idx: "+str(read_idx)+"  out of "+str(len(read_names_list)))

            reported_tax, rlen, tax_kmer_dic, tax_kmer_num_dic = kraken_kmers_cases[case].get(read_name, (0,0,{},{}))  # read length
            #if read_name not in kraken_kmers_cases[case]:
            if reported_tax==0:
                features=np.zeros(len(feature_names))
                X3.append(features)
                continue 
            
            kmer_reported_tax = tax_kmer_num_dic.get(reported_tax, 0)    # kraken tax may be LCA of others

            features_read = get_features_read(tax_kmer_num_dic, num_nodes_tree, kmer_reported_tax, rlen)
            
            if reported_tax in tax2depth:
                features_depth= get_features_depth(case, read_name, read_tax_depth, tax2depth,reported_tax,tax_kmer_dic)
            else:
                not_found_in_tax2depth.append((read_name,case,reported_tax))
                features_depth = [0, 0, 0]
            depth_reported= features_depth[0]
            features_tax, below_tax_all_perread = get_features_tax(reported_tax, tax2root_dic, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen)

            feature_half_read = get_features_half_read(reported_tax, tax_kmer_dic, rlen, below_tax_all_perread)
       
            features_path=get_features_path(path_feature, Tree_updated, reported_tax, tax_kmer_num_dic, rlen,num_nodes_tree,info, parents,child2parent)
            if reads_quality_features and read_name in reads_quality_features:
                features_readq = reads_quality_features[read_name]
            else:
                features_readq = [0]*4
            


            features=features_read + features_depth+ features_tax+ feature_half_read + [rlen] + features_path + features_readq
            #print(sum(features))




            X3.append(features)
            #print(sum(sum(X3)))
        X3=np.array(X3)
        features_cases[case]=X3    

    logging.debug("number of kmer sizes "+str(len(features_cases))+" number of reads "+str(len(features_cases[case]))+" number of features "+str(len(features_cases[case][0])))
    logging.debug("number of reads not found in tax2depth "+str(len(not_found_in_tax2depth))+" a few examples: "+str(not_found_in_tax2depth[:10]))
    return features_cases, feature_names


def read_read_names(read_name_file):
    read_names_list = []
    with open(read_name_file, 'r') as f:
        for line in f:
            read_names_list.append(line.strip())
    return read_names_list