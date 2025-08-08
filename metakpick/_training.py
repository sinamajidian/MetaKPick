import _utils_kraken
import _utils_tree

import numpy as np
from collections import Counter
from sklearn.ensemble import RandomForestRegressor


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
            #dic_cases[case]=  _utils_kraken.read_kraken_set(folder+case+"_out",readids_max)
            # todo: read the full file 
            dic_cases[case]=  _utils_kraken.read_kraken_num(folder+case+"_out",1000)
        except:
            print("n",case)
        print(case,len(dic_cases[case]))
    
    cases_readids=set()
    for case_k, case in  enumerate(cases): 
        read_ids_k=set(dic_cases[case].keys())
        cases_readids |= read_ids_k
    len(cases_readids)
    return dic_cases, cases_readids






def get_features_read(tax_kmer_num_dic, num_nodes_tree,kmer_reported_tax,rlen):
    
    
    num_kmer_all=np.array(list(tax_kmer_num_dic.keys()),dtype=np.int32)
    num_kmer_norm=num_kmer_all/rlen # read_length_dic[read_id]
    


    mean_exc_tax = (np.sum(num_kmer_all)-kmer_reported_tax)/num_nodes_tree
    
    #feature_names=['mean_all','mean_nonzero','max','sum', 'mean_exc_tax']
    features=[np.mean(num_kmer_all), np.mean(num_kmer_all)*len(num_kmer_all)/num_nodes_tree, np.max(num_kmer_all),np.sum(num_kmer_all),mean_exc_tax]
    #feature_names+=['mean_all/rlen','mean_nonzero/rlen','max/rlen','sum/rlen', 'mean_exc_tax/rlen'] # /readlength
    features+=[np.mean(num_kmer_norm)/rlen, np.mean(num_kmer_norm)*len(num_kmer_norm)/num_nodes_tree/rlen, 
                np.max(num_kmer_norm)/rlen,np.sum(num_kmer_norm)/rlen, mean_exc_tax/rlen]
    return features



def get_features_depth(case, read_id, read_tax_depth, tax2depth,tax_krak,tax_kmer_dic):
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
        
    features=[depth_reported,numkmer_consecutive_mean,weighted_depth_numkmer_mean ]
    #feature_names+=['depth_reported_tax','Avg_kmer_consecutive', 'WAvg_kmer_consecutive']
    return features


def get_features_tax(tax_krak, info, parents, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen):
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
        tax2root_t= _utils_tree.find_tax2root(info, parents, taxi)
        if tax2root_t!=-1 and tax_krak in tax2root_t: # tax_krak is lca of taxi
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
    
    features+=[kmer_reported_tax, kmer_reported_uppertax,kmer_reported_belowtax,
        kmer_reported_tax/rlen, kmer_reported_uppertax/rlen,kmer_reported_belowtax/rlen,
        kmer_other_reported_tax, kmer_other_reported_uppertax,kmer_other_reported_belowtax,
        kmer_other_reported_tax/rlen, kmer_other_reported_uppertax/rlen,kmer_other_reported_belowtax/rlen]
    
    #feature_names+=['kmer_reported_tax','kmer_tax_above','kmer_tax_below','kmer_tax/rlen','kmer_tax_above/rlen','kmer_tax_below/rlen','kmer_othertax','kmer_othertax_above','kmer_othertax_below','kmer_othertax/rlen','kmer_othertax_above/rlen','kmer_othertax_below/rlen',]
    return features






def get_features_all(read_names_list, tax2path, dic_cases, cases, read_tax_depth, tax2depth, info, parents, classification_folder):
        
    feature_names = ['mean_all', 'mean_nonzero', 'max', 'sum', 'mean_exc_tax', 'mean_all/rlen', 'mean_nonzero/rlen', 'max/rlen', 'sum/rlen', 'mean_exc_tax/rlen', 
                'depth_reported_tax', 'Avg_kmer_consecutive', 'WAvg_kmer_consecutive', 'kmer_reported_tax', 'kmer_tax_above', 'kmer_tax_below', 
                'kmer_tax/rlen', 'kmer_tax_above/rlen', 'kmer_tax_below/rlen', 'kmer_othertax', 'kmer_othertax_above', 'kmer_othertax_below', 
                'kmer_othertax/rlen', 'kmer_othertax_above/rlen', 'kmer_othertax_below/rlen', 'diff#kmers_halfRead']

    
    num_nodes_tree = len(tax2path)  # 52229
    read_names_ = list(read_names_list)
    dic_matrix2={}
    for case in cases[::-1]: 
        print(case)
        kmer_size= int(case[1:])
        X3=[]
        for read_idx,read_id in enumerate(read_names_[:100]):
            if read_idx%10000==0:
                print(read_idx,len(read_names_))
            if read_id not in dic_cases[case]:
                features=np.zeros(len(feature_names))
                X3.append(features)
                continue 

            tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic = dic_cases[case][read_id]  # read length
            kmer_reported_tax = tax_kmer_num_dic.get(tax_krak, 0)  # kraken tax may be LCA of others

            features = get_features_read(tax_kmer_num_dic, num_nodes_tree, kmer_reported_tax, rlen)
            features_depth= get_features_depth(case, read_id, read_tax_depth, tax2depth,tax_krak,tax_kmer_dic)
            depth_reported = tax2depth[tax_krak]
            features_tax = get_features_tax(tax_krak, info, parents, tax_kmer_num_dic, depth_reported, tax2depth, kmer_reported_tax, rlen)
            features=np.concatenate([features, features_depth, features_tax])

            X3.append(features)
        X3=np.array(X3)
        dic_matrix2[case]=X3    

    print("number of kmer sizes",len(dic_matrix2),"number of reads",len(dic_matrix2[case]),len(dic_matrix2[case][0]))
    return dic_matrix2


def train_RF_model(X_input,Y_input):
    regr = RandomForestRegressor(n_estimators=1000,  
                                    max_features=float(0.8), #1.0 all 
                                    max_leaf_nodes=50,
                                    random_state=14, verbose=True, n_jobs=20)  # “log2” sqrt # max_depth=10,  'log2' min_samples_leaf=10)  # 
    regr.fit(X_input, Y_input)
    return regr



def train_RF_model_all(cases, features, true_k,read_names_list):

    num_reads=len(read_names_list)
    read_k_prob={}
    for read_name in read_names_list:
        read_k_prob[read_name]=np.zeros(len(cases))    
    regr_dic = {}
    for case_idx, case in enumerate(cases):
        X_input = features[case_idx]
        Y_input = true_k[case_idx]

        regr_dic[case] = train_RF_model(X_input, Y_input)

        print(X_input.shape, Y_input.shape)

        y_pred = regr_dic[case].predict(X_input)
        y_pred_binary = np.round(y_pred)  # round(0.55)=1
        print(regr_dic[case])
        print(sum([1 for i in range(len(y_pred)) if y_pred_binary[i] == Y_input[i]]) / len(Y_input))
        # sklearn.metrics.confusion_matrix(y_train, y_traing_pred)

    for read_name_idx,read_name in enumerate(read_names_list):        
        read_k_prob[read_name][case_idx]= y_pred_binary[read_name_idx]


    return regr_dic



def get_best_tax(read_k_prob,cases,read_names_list,merged,thr_minprob=0.5):
    best_k_dic={}
    for read_name in read_names_list:
        best_k= cases[np.argmax(read_k_prob[read_name])]
        max_val = np.max(read_k_prob[read_name])

        if max_val >thr_minprob:
            best_k = cases[np.argmax(read_k_prob[read_name])]
        else:
            best_k=-1

        best_k_dic[read_name]=best_k                
    print(Counter(best_k_dic.values()))
    estimated_tax_dict={}
    for index, row in merged.iterrows():
        read=row['read_name']   
        if  read in  best_k_dic and best_k_dic[read]!=-1:
            estimated_tax= row['taxon_tool_'+best_k_dic[read]]
        else:
            estimated_tax=0
        estimated_tax_dict[read]=estimated_tax

    return best_k_dic,estimated_tax_dict
