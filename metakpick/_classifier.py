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


def calculate_accuracy(y_pred_prob,y_true):
    y_pred_binary=np.round(y_pred_prob)
    accuracy=sum([1 for i in range(len(y_true)) if y_pred_binary[i] == y_true[i]]) / len(y_true)
    logging.debug(" Accuracy of regression model is "+str(accuracy))
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


def get_best_tax(read_k_prob,read_names_list,kraken_kmers_cases,thr_minprob,info,parents,version_decision='1'):

    high_prob= 0.7
    version_decision=1
    cases=list(kraken_kmers_cases.keys())
    best_k_dic={}
    estimated_tax_dict={}
    for read_name in read_names_list:
        max_val = np.max(read_k_prob[read_name])
        if version_decision=="1":
            if max_val >thr_minprob:
                best_k = cases[np.argmax(read_k_prob[read_name])]
            else:
                best_k=-1
            best_k_dic[read_name]=best_k                
            if  best_k!=-1:
                estimated_tax= kraken_kmers_cases[best_k][read_name][0]
            else:
                estimated_tax=0
        if version_decision=="2":
            estimated_tax=0
            if max_val >thr_minprob:
                k_high_prob=[]
                for case_idx, case in enumerate(cases):
                    prob_k= read_k_prob[read_name][case_idx]
                    if prob_k >high_prob:
                        k_high_prob.append(case)
                if k_high_prob:
                    tax_list=[kraken_kmers_cases[k][read_name][0] for k in k_high_prob]
                    estimated_tax=  _utils_tree.lca_finder_list(info,parents, tax_list)
                else: # no tax with prob > high_prob
                    best_k = cases[np.argmax(read_k_prob[read_name])]
                    if best_k!=-1:
                        estimated_tax= kraken_kmers_cases[best_k][read_name][0]
                    else:
                        estimated_tax=0

        estimated_tax_dict[read_name ]=estimated_tax

    logging.debug("** best k: "+str(Counter(best_k_dic.values())))
    return best_k_dic,estimated_tax_dict