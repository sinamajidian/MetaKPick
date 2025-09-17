

import numpy as np
from collections import Counter
from sklearn.ensemble import RandomForestRegressor    
from sklearn.ensemble import RandomForestClassifier # RandomForestRegressor

import logging
from sklearn.tree import _tree
from sklearn import tree
import matplotlib.pyplot as plt

#import _utils_kraken
#import _utils_tree
import _classifier

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
   
        X_input = features_case #features_cases[case]
        Y_input = tp_binary_reads_cases[case]

        regr_dic[case] = train_RF_model(X_input, Y_input, n_estimators, max_features, max_leaf_nodes, random_state, n_jobs)

        logging.debug("case: "+str(case)+" X_input.shape "+str(X_input.shape)+" len(Y_input) "+str(len(Y_input)))
        
        y_pred_prob=regr_dic[case].predict(X_input)
        accuracy= _classifier.calculate_accuracy(y_pred_prob,Y_input)
        logging.debug("Accuracy for case "+str(case)+" is "+str(accuracy))


    return regr_dic



def get_tp_binary_reads_cases(cases, read_names_list, reads_tp_cases):

    tp_binary_reads_cases={}
    for case in cases:
        logging.debug("TP for case "+str(case))
        tp_binary_reads=[]    
        for read_name in read_names_list:
            true_k_set= reads_tp_cases[read_name]
            kmer= int(case[1:])
            if kmer in true_k_set:
                tp_binary_reads.append(1) 
            else:
                tp_binary_reads.append(0)
        tp_binary_reads_cases[case] = tp_binary_reads
        
    logging.debug("len(tp_binary_reads_cases): "+str(len(tp_binary_reads_cases))+", len(tp_binary_reads_cases[case]): "+str(len(tp_binary_reads_cases[case]))+", len(read_names_list): "+str(len(read_names_list))+", sum(tp_binary_reads_cases[case]): "+str(sum(tp_binary_reads_cases[case]))+" for case: "+case)
    for case in cases:
        logging.debug("case: "+str(case)+", Counter(tp_binary_reads_cases[case]): "+str(Counter(tp_binary_reads_cases[case])))
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




def plot_tree(regr_dic, feature_names, model_folder,num_trees=1):
        
    for case,regr in regr_dic.items():
        for tree_idx in range(num_trees):
            individual_tree = regr.estimators_[tree_idx]  # Get the first tree (you can choose any index)
            decision_tree_to_code(individual_tree,feature_names)
            plt.figure(figsize=(30, 30))
            tree.plot_tree(individual_tree, filled=True ,rounded=True,fontsize=14, feature_names=feature_names) # , class_names=class_names,
            plt.savefig(model_folder+"/tree_"+case+"_"+str(tree_idx)+".png",dpi=200)
            
    return 1



def topRF_model(read_k_prob,read_names_list,tp_binary_reads_cases, kraken_kmers_cases, cases):
    #print(sum(sum(X3)))
    X4=[]
    for read in read_names_list:
        features=read_k_prob[read]
        X4.append(features)
        #X_input = features_cases[case]
    X4=np.array(X4)

    #np.column_stack((a,b))
    Y4=[tp_binary_reads_cases[case] for case in cases]
    Y4=np.transpose(np.array(Y4))

    n_estimators=100
    max_leaf_nodes=10
    random_state=14
    n_jobs=1

    regr = RandomForestClassifier(n_estimators=n_estimators,  
                                    max_features=1, #1.0 all 
                                    max_leaf_nodes=max_leaf_nodes,
                                    random_state=random_state, verbose=True, n_jobs=n_jobs)  # “log2” sqrt # max_depth=10,  'log2' min_samples_leaf=10)  # 
    regr.fit(X4, Y4)


    y_pred2=regr.predict(X4)
    y_pred2.shape[1]
    accuracy= _classifier.calculate_accuracy(y_pred2,Y4)
    logging.info("accuracy of top model: "+str(accuracy))


    estimated_tax_dict={}
    best_k_dic={}
    ks=[int(case[1:]) for case in cases]
    true_ks={}
    for read_idx, read_name in enumerate(read_names_list):
        estimated_tax=0
        y_pred2_read=y_pred2[read_idx,:]
        true_ks[read_name] =[cases[predict_idx] for predict_idx, precict in enumerate(y_pred2_read)  if precict==  1 ]
        if true_ks[read_name]:
            best_k=int(np.median([ int(k[1:]) for k in true_ks[read_name] ])) 
            closest_k=min(ks, key=lambda x:abs(x-best_k)) # find the closest kmer size to median of true kmer sizes
            best_k= 'k'+str(closest_k)
        else:
            best_k=-1
        if best_k!=-1:
            estimated_tax= kraken_kmers_cases[best_k][read_name][0]
        estimated_tax_dict[read_name ]=estimated_tax
        best_k_dic[read_name]=best_k
    print(len(best_k_dic))
    return best_k_dic, estimated_tax_dict, regr
