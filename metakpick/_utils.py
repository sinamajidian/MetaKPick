


import logging
import numpy as np
import _utils_tree



def get_features_read(tax_kmer_num_dic, num_nodes_tree,kmer_reported_tax,rlen):
    
    
    num_kmer_all=np.array(list(tax_kmer_num_dic.keys()),dtype=np.int32)
    num_kmer_norm=num_kmer_all/rlen # read_length_dic[read_id]
    


    mean_exc_tax = (np.sum(num_kmer_all)-kmer_reported_tax)/num_nodes_tree
    
    #feature_names=['mean_all','mean_nonzero','max','sum', 'mean_exc_tax']
    features_read=[np.mean(num_kmer_all), np.mean(num_kmer_all)*len(num_kmer_all)/num_nodes_tree, np.max(num_kmer_all),np.sum(num_kmer_all),mean_exc_tax]
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


def get_features_tax(reported_tax, info, parents, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen):
    kmer_reported_uppertax=0
    kmer_reported_belowtax=0
    kmer_other_reported_tax=0
    kmer_other_reported_uppertax=0 # towards root (smaller depth value)
    kmer_other_reported_belowtax=0 # towards leaves
    
    tax_krak_2root= _utils_tree.find_tax2root(info, parents, reported_tax)
    if tax_krak_2root==-1:
        kmer_reported_uppertax=0
    else:   
        for tax_upper in tax_krak_2root[1:]: #excluding tax_krak towards the root  tax_krak_2root = [201174, 1783272, 2, 131567, 1] where tax_krak=201174 
            if tax_upper in tax_kmer_num_dic: 
                kmer_reported_uppertax +=tax_kmer_num_dic[tax_upper]
        
    below_tax_all_perread=[]
    for taxi in tax_kmer_num_dic:
        # to get the below for the reported 
        tax2root_t= _utils_tree.find_tax2root(info, parents, taxi)
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






def get_features_all(read_names_list, tax2path, kraken_kmers_cases, read_tax_depth, tax2depth, info, parents):
        
    feature_names = ['mean_all', 'mean_nonzero', 'max', 'sum', 'mean_exc_tax', 'mean_all/rlen', 'mean_nonzero/rlen', 'max/rlen', 'sum/rlen', 'mean_exc_tax/rlen', 
                'depth_reported_tax', 'Avg_kmer_consecutive', 'WAvg_kmer_consecutive', 'kmer_reported_tax', 'kmer_tax_above', 'kmer_tax_below', 
                'kmer_tax/rlen', 'kmer_tax_above/rlen', 'kmer_tax_below/rlen', 'kmer_othertax', 'kmer_othertax_above', 'kmer_othertax_below', 
                'kmer_othertax/rlen', 'kmer_othertax_above/rlen', 'kmer_othertax_below/rlen', 'diff#kmers_halfRead']
    logging.debug("feature_names are : "+str(feature_names))
    
    num_nodes_tree = len(tax2path)  # 52229
    cases=list(kraken_kmers_cases.keys())
    features_cases={}
    not_found_in_tax2depth=[]
    for case in cases: 
        logging.debug("Extracting features for case: "+case)
        kmer_size= int(case[1:])
        X3=[]
        for read_idx,read_name in enumerate(read_names_list):
            if read_idx%20000==0:
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
            features_tax, below_tax_all_perread = get_features_tax(reported_tax, info, parents, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen)

            feature_half_read = get_features_half_read(reported_tax, tax_kmer_dic, rlen, below_tax_all_perread)
            features=features_read + features_depth+ features_tax+ feature_half_read 
            #print(sum(features))

            X3.append(features)
            #print(sum(sum(X3)))
        X3=np.array(X3)
        features_cases[case]=X3    

    logging.debug("number of kmer sizes "+str(len(features_cases))+" number of reads "+str(len(features_cases[case]))+" number of features "+str(len(features_cases[case][0])))
    logging.debug("number of reads not found in tax2depth "+str(len(not_found_in_tax2depth))+" a few examples: "+str(not_found_in_tax2depth[:10]))
    return features_cases, feature_names

