import logging
import numpy as np
from collections import Counter
import random
random.seed(42)
import numpy as np
np.random.seed(42)
import _utils_tree
import _utils_kraken

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





def print_statiscs(estimated_tax_dict,cases_classify_intersect,kraken_kmers_cases,read_names_list,dic_tax_truth,info,tree_df,parents,tax_index):
    read_name_set=set(read_names_list)
    reads_tp_cases_all={}
    for tax_level in ['species','genus','family','order','class']:
        reads_tp_cases_all[tax_level] = _utils_kraken.calculate_true_k(kraken_kmers_cases,dic_tax_truth,info,tree_df,parents,tax_level,tax_index,read_names_list)


    logging.info("cleaning reported tax for all kmer sizes")
    kraken_reportedtax_cases=dict()
    kraken_reportedtax_cases['RF']=estimated_tax_dict
    for case in cases_classify_intersect:
        kraken_reportedtax_cases[case]={}
        for read_name, kraken_info in kraken_kmers_cases[case].items():
            reported_tax = kraken_info[0]
            kraken_reportedtax_cases[case][read_name]=reported_tax

    kraken_reportedtax_cases["Random"]={} 
    for read_name in read_names_list:
        case_random = cases_classify_intersect[random.randint(0, len(cases_classify_intersect)-1)]
        kraken_reportedtax_cases["Random"][read_name]  = kraken_reportedtax_cases[case_random][read_name]      	     

    kraken_reportedtax_cases[ "Oracle"]={}
    for tax_level in ['species','genus','family','order','class']:
        kraken_reportedtax_cases[ "Oracle"][tax_level]={}
        for read_name in read_names_list:
            true_k_set= reads_tp_cases_all[tax_level][read_name]
            found_true_tax= 0
            if true_k_set:
                case_true_k = 'k'+str(list(true_k_set)[0])
                found_true_tax = kraken_reportedtax_cases[case_true_k][read_name]
            kraken_reportedtax_cases[ "Oracle"][tax_level][read_name]  = found_true_tax           	     

    logging.info("Calculating the tp fp")

    print("case\tF1\tprecision\trecall\tTP\tFP\tVP")
    for level in ['species','genus','family','order','class']:
        print("level: "+level)
        for case in list(cases_classify_intersect)+["RF","Random","Oracle"]:
            logging.info("Calculating the tp fp for case: "+case+" level: "+level)
            if case=='Oracle':
                read_tpfp_dic = _utils_kraken.calculate_tp_fp('predicted_tax',kraken_reportedtax_cases[case][level],dic_tax_truth,info,tree_df,parents,level,tax_index,read_name_set)
            else:
                read_tpfp_dic = _utils_kraken.calculate_tp_fp('predicted_tax',kraken_reportedtax_cases[case],dic_tax_truth,info,tree_df,parents,level,tax_index,read_name_set)
        
            logging.info("Number of reads in the TP: "+str(len(read_tpfp_dic['TP'])) + " case: "+case+ " level: "+level)

            FP=len(read_tpfp_dic['FP-level-index'])+len(read_tpfp_dic['FP-higher-index'])+len(read_tpfp_dic['FP-level-notindex'])+len(read_tpfp_dic['FP-higher-notindex'])
            recall=len(read_tpfp_dic['TP'])/(len(read_tpfp_dic['TP']) + len(read_tpfp_dic['VP']) + len(read_tpfp_dic['FN']) +FP )
            if len(read_tpfp_dic['TP'])!=0:
                precision=len(read_tpfp_dic['TP'])/(len(read_tpfp_dic['TP']) + FP)
                F1= 2* precision* recall/(precision+recall)
            else:
                precision=0
                F1=0
            
            print('=='+level+'_'+case+','+str(round(F1,4))+","+str(round(precision,4))+","+str(round(recall,4))+","+str(len(read_tpfp_dic['TP']))+","+str(FP)+","+str(len(read_tpfp_dic['VP'])))

    return F1 