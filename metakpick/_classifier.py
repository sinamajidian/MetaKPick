import logging
import numpy as np
from collections import Counter

import _utils_tree

def apply_RF_model(cases_classify_intersect,features_cases,read_names_list,loaded_regression_dic):
    #cases=list(loaded_regression_dic.keys())
    read_k_prob={}
    for read_name in read_names_list:
        read_k_prob[read_name]=np.zeros(len(cases_classify_intersect))    
    

    for case_idx, case in enumerate(cases_classify_intersect):
        X_input = features_cases[case]
        
        y_pred = loaded_regression_dic[case].predict(X_input)

        logging.debug("X_input.shape "+str(X_input.shape)+" y_pred.shape "+str(y_pred.shape))
        #y_pred_binary = np.round(y_pred)  # round(0.55)=1
        #logging.debug("loaded_regression_dic[case]"+str(loaded_regression_dic[case]))

        for read_name_idx,read_name in enumerate(read_names_list):        
            read_k_prob[read_name][case_idx]= y_pred[read_name_idx]

    return read_k_prob


def calculate_accuracy(y_pred_prob, y_true):
    y_pred_binary=np.round(y_pred_prob)
    #print("y_pred_prob.shape", y_pred_prob.shape)
    if len(y_pred_prob.shape)==1:
        accuracy=sum([1 for i in range(len(y_true)) if y_pred_binary[i] == y_true[i]]) / len(y_true)
    elif len(y_pred_prob.shape)==2 and y_pred_prob.shape[1]>1:
        accuracy=np.mean(np.abs(y_pred_prob-y_true))/len(y_true)
    logging.debug(" Accuracy of the model is "+str(accuracy))
    return accuracy



def get_best_tax_old(read_k_prob,read_names_list,kraken_kmers_cases,thr_minprob=0.5):

    cases=list(kraken_kmers_cases.keys())
    best_k_dic={}
    estimated_tax_dict={}
    for read_name in read_names_list:
        max_val = np.max(read_k_prob[read_name])
        if max_val >thr_minprob:
            best_k = cases[np.argmax(read_k_prob[read_name])]
        else:
            best_k=-1
        best_k_dic[read_name]=best_k                
        if  best_k!=-1:
            estimated_tax= kraken_kmers_cases[best_k][read_name][0]
        else:
            estimated_tax=0
        estimated_tax_dict[read_name ]=estimated_tax

    logging.debug("** best k: "+str(Counter(best_k_dic.values())))
    

    return best_k_dic,estimated_tax_dict

def get_genus_tax(estimated_tax_raw,info,parents):
    list_tax2root = _utils_tree.find_tax2root(info,parents,estimated_tax_raw)
    idx= list_tax2root.index(estimated_tax_raw)
    if idx-1 > 0:
        estimated_tax=list_tax2root[idx-1]
    else:
        estimated_tax=estimated_tax_raw

    return estimated_tax

def get_best_tax(read_k_prob,read_names_list,kraken_kmers_cases,thr_minprob,thr_highprob_lca,thr_minprob_genus,info,parents,version_decision):
    
    cases=list(kraken_kmers_cases.keys())
    best_k_dic={}
    estimated_tax_dict={}
    for read_name in read_names_list:
        #print(read_name)
        estimated_tax=0

        max_val = np.max(read_k_prob[read_name])
        if version_decision=="1":
            #print("** Version decision is 1")
            if max_val >thr_minprob:
                best_k = cases[np.argmax(read_k_prob[read_name])]
            else:
                best_k=-1
            best_k_dic[read_name]=best_k                
            if  best_k!=-1:
                estimated_tax= kraken_kmers_cases[best_k][read_name][0]
        elif version_decision=="2":
            #print("** Version decision is 2")
            if max_val >thr_minprob:
                k_high_prob=[]
                for case_idx, case in enumerate(cases):
                    prob_k= read_k_prob[read_name][case_idx]
                    if prob_k >thr_highprob_lca:
                        k_high_prob.append(case)
                
                if k_high_prob:
                    best_k_dic[read_name]= str(k_high_prob)
                    tax_list=[kraken_kmers_cases[k][read_name][0] for k in k_high_prob]
                    estimated_tax=  _utils_tree.lca_finder_list(info,parents, tax_list)
                else: # no tax with prob > high_prob but prob > thr_minprob
                    best_k = cases[np.argmax(read_k_prob[read_name])]
                    best_k_dic[read_name]=best_k  
                    if best_k!=-1:    
                        estimated_tax= kraken_kmers_cases[best_k][read_name][0]
        elif version_decision=="3": 
            #print("** Version decision is 1")
            if max_val >thr_minprob:
                best_k = cases[np.argmax(read_k_prob[read_name])]
                if  best_k!=-1:
                    estimated_tax= kraken_kmers_cases[best_k][read_name][0]
            elif max_val >thr_minprob_genus:
                best_k = cases[np.argmax(read_k_prob[read_name])]
                if  best_k != -1:
                    estimated_tax_raw= kraken_kmers_cases[best_k][read_name][0]
                    estimated_tax= get_genus_tax(estimated_tax_raw,info,parents)
            else:
                best_k=-1
            best_k_dic[read_name]=best_k 
            
                
            
        estimated_tax_dict[read_name ]=estimated_tax

    logging.debug("** best k: "+str(Counter(best_k_dic.values())) + " version decision: "+version_decision + " high_prob: "+str(thr_highprob_lca) + " min_prob: "+str(thr_minprob_genus) + " min_prob: "+str(thr_minprob))
    return best_k_dic,estimated_tax_dict