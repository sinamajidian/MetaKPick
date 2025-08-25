

import numpy as np
from collections import Counter
from sklearn.ensemble import RandomForestRegressor
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
   
        X_input = features_cases[case]
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



