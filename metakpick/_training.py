import _utils_kraken
import numpy as np

def read_kraken_all(cases, folder, readids_max):
    print("read kraken's k-mer  count per tax")
    #folder="/vast/blangme2/smajidi5/metagenomics/changek/simulatation/classification/" # 
    #cases=['k19','k25','k31'] # ,'k21'
    dic_cases={}
    for case in reversed(cases): # 
        print(case)
        try:
            #dic_cases[case]=read_kraken_limited(folder+case+"_out",10000)
            #dic_cases[case]=read_kraken(folder+case+"_out")
            dic_cases[case]=  _utils_kraken.read_kraken_set(folder+case+"_out",readids_max)
        except:
            print("n",case)
        print(case,len(dic_cases[case]))
    
    cases_readids=set()
    for case_k, case in  enumerate(cases): 
        read_ids_k=set(dic_cases[case].keys())
        cases_readids |= read_ids_k
    len(cases_readids)
    return dic_cases, cases_readids






def get_features(dic_cases, case, read_id, num_nodes_tree):
    tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic = dic_cases[case][read_id] # read length
        
    num_kmer_all=np.array(list(tax_kmer_num_dic.keys()),dtype=np.int32)
    num_kmer_norm=num_kmer_all/rlen # read_length_dic[read_id]
    

    if tax_krak in tax_kmer_num_dic:
        kmer_reported_tax = tax_kmer_num_dic[tax_krak] 
    else:
        kmer_reported_tax=0 # probably this kraken tax is the lca of a few tax closer to leaves
    
    mean_exc_tax = (np.sum(num_kmer_all)-kmer_reported_tax)/num_nodes_tree
    
    #feature_names=['mean_all','mean_nonzero','max','sum', 'mean_exc_tax']
    features=[np.mean(num_kmer_all), np.mean(num_kmer_all)*len(num_kmer_all)/num_nodes_tree, np.max(num_kmer_all),np.sum(num_kmer_all),mean_exc_tax]
    #feature_names+=['mean_all/rlen','mean_nonzero/rlen','max/rlen','sum/rlen', 'mean_exc_tax/rlen'] # /readlength
    features+=[np.mean(num_kmer_norm)/rlen, np.mean(num_kmer_norm)*len(num_kmer_norm)/num_nodes_tree/rlen, 
                np.max(num_kmer_norm)/rlen,np.sum(num_kmer_norm)/rlen, mean_exc_tax/rlen]
    return features



def get_features_depth(dic_cases, case, read_id, read_tax_depth, tax2depth):
    if read_id in read_tax_depth:
        txdepth_norm= read_tax_depth[case][read_id]
    else:
        txdepth_norm=0
    
    depth_reported= tax2depth[tax_krak]
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
        
    features+=[depth_reported,numkmer_consecutive_mean,weighted_depth_numkmer_mean ]
    #feature_names+=['depth_reported_tax','Avg_kmer_consecutive', 'WAvg_kmer_consecutive']
    return features


def get_features_tax(dic_cases, case, read_id, read_tax_depth, tax2depth, info, parents, tax_kmer_num_dic):
    kmer_reported_uppertax=0
    kmer_reported_belowtax=0
    kmer_other_reported_tax=0
    kmer_other_reported_uppertax=0 # towards root (smaller depth value)
    kmer_other_reported_belowtax=0 # towards leaves
    
    tax_krak_2root= _utils_tree.find_tax2root(info, parents, tax_krak)
    for tax_upper in tax_krak_2root[1:]: #excluding tax_krak towards the root  tax_krak_2root = [201174, 1783272, 2, 131567, 1] where tax_krak=201174 
        if tax_upper in tax_kmer_num_dic: 
            kmer_reported_uppertax +=tax_kmer_num_dic[tax_upper]
    
    below_tax_all_perread=[]
    for taxi in tax_kmer_num_dic:
        # to get the below for the reported 
        tax2root_t= find_tax2root(info, parents, taxi)
        if tax2root_t!=-1 and tax_krak in tax2root_t: # tax_krak is lca of taxi
            kmer_reported_belowtax += tax_kmer_num_dic[taxi]   
            below_tax_all_perread.append(taxi)
    
        # to get the below/above/at  the reported tax, accross all but not the max path 
        if depth_reported==tax_depth:
            kmer_other_reported_tax+=tax_kmer_num_dic[taxi]   
        elif depth_reported<tax_depth: # towards root (smaller depth value)
            kmer_other_reported_uppertax+=tax_kmer_num_dic[taxi]   
        elif depth_reported>tax_depth:
            kmer_other_reported_belowtax+=tax_kmer_num_dic[taxi]   
    
    features+=[kmer_reported_tax, kmer_reported_uppertax,kmer_reported_belowtax,
        kmer_reported_tax/rlen, kmer_reported_uppertax/rlen,kmer_reported_belowtax/rlen,
        kmer_other_reported_tax, kmer_other_reported_uppertax,kmer_other_reported_belowtax,
        kmer_other_reported_tax/rlen, kmer_other_reported_uppertax/rlen,kmer_other_reported_belowtax/rlen]
    
    #feature_names+=['kmer_reported_tax','kmer_tax_above','kmer_tax_below','kmer_tax/rlen','kmer_tax_above/rlen','kmer_tax_below/rlen','kmer_othertax','kmer_othertax_above','kmer_othertax_below','kmer_othertax/rlen','kmer_othertax_above/rlen','kmer_othertax_below/rlen',]
    return features