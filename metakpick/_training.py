

import numpy as np
from collections import Counter
from sklearn.ensemble import RandomForestRegressor
import logging
from sklearn.tree import _tree
from sklearn import tree
import matplotlib.pyplot as plt


import _utils_kraken
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

    
    num_nodes_tree = len(tax2path)  # 52229
    cases=list(kraken_kmers_cases.keys())
    features_cases={}
    for case in cases: 
        logging.debug("working on case: "+case)
        kmer_size= int(case[1:])
        X3=[]
        for read_idx,read_name in enumerate(read_names_list):
            if read_idx%10000==0:
                logging.debug("Working on read idx: "+str(read_idx)+"  out of "+str(len(read_names_list)))

            reported_tax, rlen, tax_kmer_dic, tax_kmer_num_dic = kraken_kmers_cases[case][read_name]  # read length
            #if read_name not in kraken_kmers_cases[case]:
            if reported_tax==0:
                features=np.zeros(len(feature_names))
                X3.append(features)
                continue 
            
            kmer_reported_tax = tax_kmer_num_dic.get(reported_tax, 0)    # kraken tax may be LCA of others

            features_read = get_features_read(tax_kmer_num_dic, num_nodes_tree, kmer_reported_tax, rlen)
            features_depth= get_features_depth(case, read_name, read_tax_depth, tax2depth,reported_tax,tax_kmer_dic)
            depth_reported = tax2depth[reported_tax]
            features_tax, below_tax_all_perread = get_features_tax(reported_tax, info, parents, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen)

            feature_half_read = get_features_half_read(reported_tax, tax_kmer_dic, rlen, below_tax_all_perread)
            features=features_read + features_depth+ features_tax+ feature_half_read 
            #print(sum(features))

            X3.append(features)
            #print(sum(sum(X3)))
        X3=np.array(X3)
        features_cases[case]=X3    

    logging.debug("number of kmer sizes"+str(len(features_cases))+" number of reads"+str(len(features_cases[case]))+" number of features"+str(len(features_cases[case][0])))
    return features_cases, feature_names


def train_RF_model(X_input,Y_input, n_estimators=1000, max_features=float(0.8), max_leaf_nodes=50, random_state=14, n_jobs=1):
    regr = RandomForestRegressor(n_estimators=n_estimators,  
                                    max_features=max_features, #1.0 all 
                                    max_leaf_nodes=max_leaf_nodes,
                                    random_state=random_state, verbose=True, n_jobs=n_jobs)  # “log2” sqrt # max_depth=10,  'log2' min_samples_leaf=10)  # 
    regr.fit(X_input, Y_input)
    return regr



def train_RF_model_all(features_cases, tp_binary_reads_cases,read_names_list,n_estimators, max_features, max_leaf_nodes, random_state, n_jobs=1):
    
    num_reads=len(read_names_list)
    regr_dic = {}
    for case, features_case in features_cases.items():
   
        X_input = features_cases[case]
        Y_input = tp_binary_reads_cases[case]

        regr_dic[case] = train_RF_model(X_input, Y_input, n_estimators, max_features, max_leaf_nodes, random_state, n_jobs)

        logging.debug("X_input.shape "+str(X_input.shape)+" len(Y_input) "+str(len(Y_input)))


    return regr_dic



def get_tp_binary_reads_cases(read_names_list, reads_tp_cases):

    tp_binary_reads_cases={}
    cases=list(reads_tp_cases.keys())
    for case in cases:
        tp_binary_reads=[]    
        for read_name in read_names_list:
            if read_name in reads_tp_cases[case]:
                tp_binary_reads.append(1) 
            else:
                tp_binary_reads.append(0)
        tp_binary_reads_cases[case] = tp_binary_reads
        
    logging.debug("len(tp_binary_reads_cases) "+str(len(tp_binary_reads_cases))+" len(tp_binary_reads_cases[case]) "+str(len(tp_binary_reads_cases[case]))+" len(read_names_list) "+str(len(read_names_list))+" sum(tp_binary_reads_cases[case]) "+str(sum(tp_binary_reads_cases[case])))
    for case in cases:
        logging.debug("case"+str(case)+" Counter(tp_binary_reads_cases[case])"+str(Counter(tp_binary_reads_cases[case])))
    return tp_binary_reads_cases




def write_estimated_tax(estimated_tax_dict,output_file_name="estimated_tax.csv"):
    output_file =open(output_file_name,"w")

    for read_name, estimated_tax in estimated_tax_dict.items():
        output_file.write(read_name+","+str(estimated_tax)+"\n")
    output_file.close()

    return output_file_name

def read_estimated_tax(output_file_name="estimated_tax.csv"):
    estimated_tax_dict = {}
    with open(output_file_name, "r") as file:
        for line in file:
            read_name, estimated_tax = line.strip().split(",")
            estimated_tax_dict[read_name] = int(estimated_tax)
    return estimated_tax_dict


def decision_tree_to_code(tree, feature_names):

    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]
    code_txt ="def tree({}):".format(", ".join(feature_names))+"\n"
    print(code_txt)
    def recurse_func(node, depth):
        #code_txt=""
        indent = "  " * depth
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]
            code_txt = "{}if {} <= {}:".format(indent, name, threshold)+"\n"
            print(code_txt)
            recurse_func(tree_.children_left[node], depth + 1)
            code_txt = "{}else:  # if {} > {}".format(indent, name, threshold)+"\n"
            print(code_txt)
            recurse_func(tree_.children_right[node], depth + 1)
        else:
            code_txt = "{}return {}".format(indent, tree_.value[node])+"\n"
            print(code_txt)
        #print(code_txt)
        return 1
    
    recurse_func(0, 1)

    # with open(workingdir+"../results/tree_code.txt", "w") as file:
    #     for line in code_txt.split("\n"):
    #         file.write(line+"\n")
    # logging.debug("Tree code saved in "+workingdir+"../results/tree_code.txt")
    return 1




def plot_tree(regr_dic, feature_names, workingdir,num_trees=1):
        
    for case,regr in regr_dic.items():
        for tree_idx in range(num_trees):
            individual_tree = regr.estimators_[tree_idx]  # Get the first tree (you can choose any index)
            decision_tree_to_code(individual_tree,feature_names)
            plt.figure(figsize=(30, 30))
            tree.plot_tree(individual_tree, filled=True ,rounded=True,fontsize=14, feature_names=feature_names) # , class_names=class_names,
            plt.savefig(workingdir+"../results/tree_"+case+"_"+str(tree_idx)+".png",dpi=100)

    return 1



