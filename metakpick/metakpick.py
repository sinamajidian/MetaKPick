

import argparse
import sys
import os
import logging
import pandas as pd
import numpy as np
import random
import pickle
from datetime import datetime
import time 
np.random.seed(42)

import __init__
import _utils
import _training
import _classifier
import _utils_tree
import _utils_kraken

# import matplotlib.pyplot as plt
# import seaborn as sns
# from collections import Counter

def main():
    start_time = time.time()    
    version_info = str(__init__.__packagename__) + str(__init__.__version__)
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(description='Metakpick: for metagenomic classification of reads')
    parser.add_argument('--mode', type=str, default='classify', help='train or classify')
    
    # Input arguments
    # parser.add_argument('--workingdir', type=str, help='Working directory')
    parser.add_argument('--kraken_out', type=str, help='Kraken output folder')
    parser.add_argument('--truth_file', type=str, help='Truth file')
    parser.add_argument('--model_file', type=str, help='Model file')
    parser.add_argument('--tree_file', type=str, help='Tree file')
    parser.add_argument('--tax_genome_file', type=str, help='Tax genome file')
    parser.add_argument('--fastq', type=str, help='Use fastq file, read quality, to generate features')  
    parser.add_argument('--read_name_file', type=str, help='to use in calssify mode,  a file with each line a read name, to limit kraken output to the read names in the file')

    parser.add_argument('--output_file_name', type=str, help='Output file name')

    parser.add_argument('--plot_tree', action='store_true', help='Plot the tree')
    parser.add_argument("--version", action="version", version=version_info)

    # Algorithm arguments
    parser.add_argument('--kmer_list', type=str, help='K-mer list file comma separated to use during training/testing') # 19,21,23 
    parser.add_argument('--path_feature', action='store_true', help=' generate and use path feature, slower')
    parser.add_argument('--version_decision', type=str, default='1', help='Version decision')
    parser.add_argument('--thr_minprob', type=str, default='0.5', help='with version_decision 1/2/3, threshold for minimum probability to decide unclassified reads')
    parser.add_argument('--thr_highprob_lca', type=str, default='0.5', help='with version_decision 2, threshold for high probability to decide lca')
    parser.add_argument('--thr_minprob_genus', type=str, default='0.5', help='with version_decision 3,threshold for minimum probability to decide genus')
    parser.add_argument('--topRF', action='store_true', help='Use a RF on top of all the RFs')
    parser.add_argument('--training_level', type=str, default='species', help='taxanomic level/rank for training')
    #parser.add_argument('--help', action='store_true', help='show this help message and exit')
    #parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)
    
    # Print version information
    logging.info(f"MetaKPick version: {version_info}")
    logging.info("Args: "+str(args))
    mode=args.mode # train or classify
    kraken_output_folder=args.kraken_out
    truth_file=args.truth_file
    output_file_name=args.output_file_name
    tree_file=args.tree_file
    tax_genome_file=args.tax_genome_file
    plot_tree=args.plot_tree
    model_file=args.model_file
    path_feature=args.path_feature
    version_decision=args.version_decision
    topRF=args.topRF
    tax_level_training=args.training_level
    thr_minprob= float(args.thr_minprob)
    thr_highprob_lca= float(args.thr_highprob_lca)
    thr_minprob_genus= float(args.thr_minprob_genus)
    fastq_file=args.fastq
    read_name_file=args.read_name_file
    kmer_list=[] if args.kmer_list is None else [int(kmer) for kmer in args.kmer_list.split(',')]
    logging.info("Input kmer list: "+str(kmer_list))


    #  For either training or classification
    logging.info("Starting the program and loading the tax info")
    info, Tree, tax_index, tax_genome, parents, tree_df = _utils_tree.get_tax_info(tree_file,tax_genome_file)    
    tax_level_path='species' # the lowest resolution taxanomic level/rank to use for the path/ features or taxnomic assignments   
    tax2path, tax2depth, tax2root_all_dic, tax_genome_specieslevel = _utils_tree.get_tax2path(tax_genome, info, parents,tree_df,tax_level_path)

    if mode=="train": 
        logging.info("Reading the truth file: "+truth_file)
        dic_tax_truth = _utils_kraken.read_truth_file(truth_file)
        logging.info("Number of reads in the truth file: "+str(len(dic_tax_truth)))
        cases=[i.split("_")[0] for i in os.listdir(kraken_output_folder) if i.endswith('_out')]
        logging.info("List of kraken indexes  aka cases: "+str(cases))
        if kmer_list:
            cases=[case for case in cases if int(case[1:]) in kmer_list] #'k19'
        logging.info("List of kraken indexes  intersected with the input kmer list: "+str(cases))

        logging.info("Reading the kraken kmers from: "+kraken_output_folder)
        read_names_list, kraken_kmers_cases = _utils_kraken.read_kraken_all(cases, kraken_output_folder)
        logging.info("Getting the tax depth and true k-mer size at the training level: "+tax_level_training)
        read_tax_depth = _utils_kraken.get_tax_depth(kraken_kmers_cases, info, parents)
        reads_tp_cases = _utils_kraken.calculate_true_k(kraken_kmers_cases,dic_tax_truth,info,tree_df,parents,tax_level_training,tax_index,read_names_list)
        logging.info("Getting the True positive read  names for each k-mer size case")
        tp_binary_reads_cases = _training.get_tp_binary_reads_cases(cases, read_names_list, reads_tp_cases)
        reads_quality_features={}
        if fastq_file:
            logging.info("Reading the reads quality from fastq file: "+fastq_file)
            reads_quality = _utils.parse_fastq_reads(fastq_file)
            logging.info("Getting the reads quality features")
            reads_quality_features = _utils.get_features_readq(reads_quality)
        features_cases, feature_names = _utils.get_features_all(read_names_list, tax2path, kraken_kmers_cases, read_tax_depth, tax2depth, info, parents, Tree, tax_index, reads_quality_features, path_feature)
        logging.info("Cases in features: "+str(features_cases.keys()))
        
        logging.info("Training the RF model for all k-mer sizes")
        regr_dic = _training.train_RF_model_all(features_cases, tp_binary_reads_cases, read_names_list, n_estimators=1000, max_features=float(0.8), max_leaf_nodes=50, random_state=14, n_jobs=20) 
        logging.info("Saving the model")
        logging.info(regr_dic)


        if topRF:
            logging.info("Working on mode TopRF, using a RF on top of all the RFs")

            logging.info("Applying the model")
            read_k_prob= _classifier.apply_RF_model(cases, features_cases,read_names_list,regr_dic)
            logging.info("Number of reads in the read_k_prob: "+str(len(read_k_prob)))
            #print("\n\n\n*&&__&&  read_k_prob",list(read_k_prob.keys())[0], read_k_prob[list(read_k_prob.keys())[0]])
            best_k_dic, estimated_tax_dict, regr_topRF = _training.topRF_model(read_k_prob,read_names_list,tp_binary_reads_cases, kraken_kmers_cases, cases)
            regr_dic['topRF']=regr_topRF

        logging.info("Saving all models in: "+model_file) 
        pickle.dump(regr_dic, open(model_file, "wb"))
        if plot_tree:
            folder_to_save=  "/".join(model_file.split("/")[:-1])
            _training.plot_tree(regr_dic, feature_names, folder_to_save, num_trees=1)
            logging.info("A few example decision trees plotted and saved in: "+folder_to_save)
    

    elif mode=="classify":

        logging.info("Loading the model: "+model_file)
        loaded_regression_dic= pickle.load(open(model_file, "rb"))
        logging.info("Model loaded: "+str(loaded_regression_dic))
        cases_classify =[i.split("_")[0] for i in os.listdir(kraken_output_folder) if i.endswith('_out')]
        cases_model= list(loaded_regression_dic.keys())
        cases_classify_intersect=sorted(list(set(cases_classify).intersection(set(cases_model))))
        if kmer_list:   
            kmer_list_set=set(['k'+str(kmer) for kmer in kmer_list])
            cases_classify_intersect=sorted(list(set(cases_classify_intersect).intersection(set(kmer_list_set))))
        logging.info("Cases in the input for classification intersected with the kmer list: "+str(cases_classify_intersect))    
            
        logging.info("Cases in the input for classification: "+str(cases_classify))
        logging.info("Cases in the model: "+str(cases_model))
        logging.warning("Working on the intersection of cases in the input for classification and in the model: "+str(cases_classify_intersect))        
        logging.info("Reading the kraken kmers for these cases")
        read_names_list_all, kraken_kmers_cases = _utils_kraken.read_kraken_all(cases_classify_intersect, kraken_output_folder)
        logging.info("Number of reads in the kraken kmers: "+str(len(kraken_kmers_cases)))

        logging.info("Getting the tax depth")
        read_tax_depth = _utils_kraken.get_tax_depth(kraken_kmers_cases, info,parents)
        logging.info("Getting the features")

        reads_quality_features={}
        if fastq_file:
            logging.info("Reading the reads quality from: "+fastq_file)
            reads_quality = _utils.parse_fastq_reads(fastq_file)
            logging.info("Getting the reads quality features")
            reads_quality_features = _utils.get_features_readq(reads_quality)

        read_names_list_=read_names_list_all
        if read_name_file:
            logging.info("Reading the read names from: "+read_name_file)
            read_names_list_input = _utils.read_read_names(read_name_file)
            logging.info("Number of reads in the read names file: "+str(len(read_names_list_input)))
            read_names_list_ = [ i for i in read_names_list_all if i in   read_names_list_input ]
            logging.info("Number of reads in the read names list intersected with the read names file: "+str(len(read_names_list_)))

        features_cases, feature_names = _utils.get_features_all(read_names_list_, tax2path, kraken_kmers_cases, read_tax_depth, tax2depth, info, parents, Tree, tax_index, reads_quality_features, path_feature)
        logging.info("Cases in features: "+str(features_cases.keys()))
        logging.info("Applying the model")
        read_k_prob= _classifier.apply_RF_model(cases_classify_intersect, features_cases,read_names_list_,loaded_regression_dic)
        logging.info("Number of reads in the read_k_prob: "+str(len(read_k_prob)))

        if topRF:
            cases=cases_classify_intersect
            logging.info("workingmode TopRF")

            X4=_training.get_topRF_features(read_k_prob,read_names_list_)
            y_pred2=loaded_regression_dic['topRF'].predict(X4)
           
            logging.info("Getting the best tax for TopRF")
            estimated_tax_dict, best_k_dic = _training.get_best_tax_topRF(read_k_prob,read_names_list_,y_pred2,cases,kraken_kmers_cases)
            
        else:
            logging.info("Getting the best tax")
            best_k_dic, estimated_tax_dict = _classifier.get_best_tax(read_k_prob,read_names_list_,kraken_kmers_cases,thr_minprob,thr_highprob_lca,thr_minprob_genus,info,parents,version_decision) 
        logging.info("Writing the estimated tax")
        output_file_name= _training.write_estimated_tax(estimated_tax_dict,output_file_name)
        logging.info("Number of reads in the dic_tax_estimated: "+str(len(estimated_tax_dict)))

        
        if truth_file:

            logging.info("Reading the truth file: "+truth_file)
            dic_tax_truth = _utils_kraken.read_truth_file(truth_file)
            logging.info("Number of reads in the truth file: "+str(len(dic_tax_truth)))
            reads_tp_cases_all={}
            for tax_level in ['species','genus','family','order','class']:
                reads_tp_cases_all[tax_level] = _utils_kraken.calculate_true_k(kraken_kmers_cases,dic_tax_truth,info,tree_df,parents,tax_level,tax_index,read_names_list_)
        logging.info("Calculating and printing the statistics")
        _classifier.print_statiscs(estimated_tax_dict,cases_classify_intersect,kraken_kmers_cases,read_names_list_,dic_tax_truth,info,tree_df,parents,tax_index,read_names_list_)

    
    # Calculate and print total runtime
    end_time = time.time()
    total_runtime = end_time - start_time
    hours = int(total_runtime // 3600)
    minutes = int((total_runtime % 3600) // 60)
    seconds = total_runtime % 60
    runtime_str = f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"
    logging.info(f"Total runtime: {runtime_str} (HH:MM:SS.sss)")

if __name__ == "__main__":
    main()
    print("Done!")
    logging.info("Done!")
