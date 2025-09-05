
import os
import pandas as pd
import _utils_tree
from collections import Counter
import logging
import copy

def read_truth_file(truth_file):
    merged = pd.read_csv(truth_file, names=["read_name", "taxon"])
    dic_tax_truth={}
    for index, row in merged.iterrows():
        dic_tax_truth[row['read_name']]=row['taxon']
    return dic_tax_truth


def get_tax_depth(kraken_kmers_cases, info, parents):
    logging.debug("Getting the tax depth")
    read_tax_depth={}
    for case, kraken_kmers in kraken_kmers_cases.items(): 
        read_tax_depth[case]={}
        for read_name, kraken_read_info in kraken_kmers.items():
            reported_tax=kraken_read_info[0]
            tax2root_list =  _utils_tree.find_tax2root(info,parents, reported_tax) # [2711156, 1649486, 41297, 204457, 28211, 1224, 3379134, 2, 131567, 1]
            if tax2root_list!=-1:
                depth=len(tax2root_list)
            else:
                depth=0
            read_tax_depth[case][read_name]=depth
            
    logging.debug("Number of reads in read_tax_depth: "+str(len(read_tax_depth))+" "+str(len(read_tax_depth[case])))
    logging.debug("Counter of read_tax_depth: "+str(Counter(read_tax_depth[case].values())))

    return read_tax_depth  


def calculate_true_k(kraken_kmers_cases,dic_tax_truth,info,tree_df,parents,tax_level,tax_index,read_names_list):
    tp_cases_dic={}
    for case, kraken_kmers in kraken_kmers_cases.items():
        logging.debug('finding true k for case '+case)
        #merged  = _utils_kraken.read_kraken_classification(cases, truth_file, classification_folder )    
        
        
        read_tpfp_dic_case= calculate_tp_fp("raw_kraken",kraken_kmers,dic_tax_truth,info,tree_df,parents,tax_level,tax_index)
        logging.debug("Number of reads in the case TP : "+str(len(read_tpfp_dic_case['TP'])))
        if len(read_tpfp_dic_case['notruth']):
            logging.debug("Number of reads in the case not truth level : "+str(len(read_tpfp_dic_case['notruth']))+" a few examples: "+str(list(read_tpfp_dic_case['notruth'])[:10]))
        tp_cases_dic[case]=read_tpfp_dic_case['TP']

    cases=list(kraken_kmers_cases.keys())
    reads_tp_cases={}
    for read_name in read_names_list:
        reads_tp_cases[read_name]=set()
        for case_k, case in  enumerate(cases): 
            k = int(case[1:3])
            if read_name in tp_cases_dic[case]:
                reads_tp_cases[read_name].add(k)
    logging.debug("Number of reads in the case TP : "+str(len(reads_tp_cases)))
    logging.debug("Number of kmer sizes in the TP for read"+read_name+" : "+str(len(reads_tp_cases[read_name])))
    return reads_tp_cases





# def limit_read_per_genome(merged,num_max):

#     counter_tax={}
#     readids_max=set()
#     for index, row in merged.iterrows():
#         tax_true= row['taxon']
#         read_name=row['read_name']

#         to_add=True
#         if not tax_true in counter_tax:
#             counter_tax[tax_true]=1
#         else:
#             if counter_tax[tax_true]>num_max:
#                 to_add=False 
#             counter_tax[tax_true]+=1
                
#         if to_add:
#             readids_max.add(read_name)
#     logging.debug("limiting reads per genome: "+str(len(merged))+" "+str(len(readids_max))+" "+str(len(counter_tax))+" "+str(counter_tax[tax_true]))
#     merged_limited=merged.copy()
#     merged_limited= merged[merged['read_name'].isin(readids_max)]
#     merged_limited.reset_index(drop=True, inplace=True)
#     logging.debug("Number of reads in merged_limited: "+str(len(merged_limited))+" "+str(len(merged)))

#     return merged_limited, readids_max


# def write_max_pergenome_calssificaiton(dic_tax_truth,num_max):

#     num_max=10
#     counter_tax={}
#     readids_max=set()
#     for read_name, tax_true in dic_tax_truth.items():

#         to_add=True
#         if not tax_true in counter_tax:
#             counter_tax[tax_true]=1
#         else:
#             if counter_tax[tax_true]>num_max:
#                 to_add=False 
#             counter_tax[tax_true]+=1
                
#         if to_add:
#             readids_max.add(read_name)
#     print(len(dic_tax_truth),len(readids_max))


#     folder_clas="simulatation/classification/"
#     for case in ['k17', 'k19', 'k21', 'k23', 'k25', 'k27', 'k29', 'k31']:
#         print(case)
#         kraken_file=folder_clas+case+"_out"
#         kraken_file_out=folder_clas+"/max15/"+case+"_out"
#         file_in=open(kraken_file,'r')
#         file_out=open(kraken_file_out,'w')
#         for line in file_in: 
#             line_split=line.strip().split("\t")
#             #classified ,read_id, tax_krak, read_len, dic = line_split
#             read_id=line_split[1]
#             if read_id in readids_max:
#                 file_out.write(line.strip()+"\n")
#         file_in.close()
#         file_out.close()

# file_in=open(folder_clas+"true_tax.csv",'r')
# file_out=open(folder_clas+"true_tax_max15.csv",'w')
# for line in file_in: 
#     line_split=line.strip().split(",")
#     #classified ,read_id, tax_krak, read_len, dic = line_split
#     read_id=line_split[0]
#     if read_id in readids_max:
#         file_out.write(line.strip()+"\n")
# file_in.close()
# file_out.close()

#     return readids_max




def calculate_tp_fp(mode,input_dic,dic_tax_truth,info,tree_df,parents,tax_level,tax_index):

    """
    input_dic:
        raw_kraken: kraken_kmers_case     mode:"raw_kraken"
        predicted_tax: predicted_tax_dic  mode:"predicted_tax"
    """
    read_tpfp_dic={"TP":set(),"FP-level-notindex":set(), 'FP-higher-notindex':set(),"FP-level-index":set(), 'FP-higher-index':set(),
               'inconsistent-tool':set(), 'TN':set(),
               'FN':set(),'VP':set(), "truth-level":set(), "no-truth-level":set(),"total_unclassified":set(),'notruth':[]} #  at the level ,'taxNotfoundintool':[]
    
    tid_true_taxlevel_dic={}
    tax_true_all = set(dic_tax_truth.values())
    for tid_true in tax_true_all:
        tid_true_taxlevel_dic[tid_true] = int(_utils_tree.find_tax_level(info,tree_df,parents, tid_true, tax_level))
    tax_predicted_level_dic={}
    tax_predicted_raw=set()    
    for read_name, kraken_read_info in input_dic.items():
        if mode=="raw_kraken":
            tax_predicted_raw.add(kraken_read_info[0])
        elif mode=="predicted_tax":
            tax_predicted_raw.add(kraken_read_info)
    tax_predicted_raw.add(2) # todo check the other place  where we added 2  
    for tax_predicted in tax_predicted_raw:
        tax_predicted_level_dic[tax_predicted] = int(_utils_tree.find_tax_level(info,tree_df,parents, tax_predicted, tax_level))
    tax_predicted_list_tocheck=[]
    for read_name, kraken_read_info in input_dic.items():
        if read_name not in dic_tax_truth:
            read_tpfp_dic['notruth'].append(read_name)
            continue
        tax_true = dic_tax_truth[read_name]
        if mode=="raw_kraken":
            tax_predicted = kraken_read_info[0]
        elif mode=="predicted_tax":
            tax_predicted = kraken_read_info
        else:
            raise ValueError("mode not supported")

        if tax_predicted==1:
            tax_predicted=2 # todo probably 1 has another meaning in find tax level, that's fine the level should be much deeper 
    
        tax_true_level = tid_true_taxlevel_dic[tax_true]
        tax_predicted_level = tax_predicted_level_dic[tax_predicted]
    
        if tax_predicted_level==-1 and tax_predicted!=0:
            read_tpfp_dic['inconsistent-tool'].add(read_name)
            # tid_tool_taxlevel==-1 we are ignoring these reads if the tax is not in the node.dump 
            # we think these cases are limited and probably happened due to different node.dump in build vs evaluation 
    
        if tax_predicted==0 :
            read_tpfp_dic["total_unclassified"].add(read_name)
        
        if tax_predicted==-1 or tax_true_level < 1:
            tax_predicted_list_tocheck.append(tax_predicted)
        

    
        if tax_true_level ==1 or  tax_true_level==-1 or (tax_predicted_level==-1 and tax_predicted!=0): # true is not at the level, is higher than the level     or true doest not exist in info
            #including 'inconsistent-tool'
            read_tpfp_dic["no-truth-level"].add(read_name)
            
        else:
            read_tpfp_dic["truth-level"].add(read_name)
    
        if tax_true_level >1 and tax_predicted==0:
            if tax_true in tax_index:
                read_tpfp_dic['FN'].add(read_name) # unclassified, among reads with true label at the level 
            else:
                read_tpfp_dic['TN'].add(read_name) # unclassified, among reads with true label at the level 
            
        if tax_true_level>1 and tax_predicted_level==tax_true_level :
            read_tpfp_dic['TP'].add(read_name) # correct
    
        if tax_true_level>1 and tax_predicted>1 and tax_predicted_level==1 :  # tool labled at higher than taxlevel 
            tool_level= tree_df[4][info[tax_predicted]]
            tax_true_toollevel = _utils_tree.find_tax_level(info,tree_df,parents, tax_true, tool_level)
            if tax_true_toollevel == tax_predicted: # at the tool level is correct 
                read_tpfp_dic['VP'].add(read_name) # vague positive, one deeper lvel (ancestror)  match    
            if tax_true_toollevel != tax_predicted:
                if tax_true in tax_index:
                    read_tpfp_dic['FP-higher-index'].add(read_name) # tool labeld at higher level and worng (not a VP) 
                else:
                    read_tpfp_dic['FP-higher-notindex'].add(read_name) # tool labeld at higher level and worng (not a VP) 
    
        if tax_true_level > 1 and tax_predicted_level>1 and  tax_predicted_level != tax_true_level :
            if tax_true in tax_index:
                read_tpfp_dic['FP-level-index'].add(read_name) # tool labeled at the tax level (or lower) and wrong 
            else:
                read_tpfp_dic['FP-level-notindex'].add(read_name) # tool labeled at the tax level (or lower) and wrong             
    if tax_predicted_list_tocheck:
        logging.debug("to check tax_predicted==-1 or tax_true_level < 1: "+str(len(tax_predicted_list_tocheck))+" a few of them are: "+str(tax_predicted_list_tocheck[:20]))        
    return read_tpfp_dic

# def get_tax_depth(merged, cases, info,parents):
#     logging.debug("Getting the tax depth")
#     read_tax_depth={}
#     for case in cases: 
#         read_tax_depth[case]={}
#         tocheck_=""
#         for index, row in merged.iterrows():
#             taxid=row['taxon_tool_'+case]
#             tax2root_list = _utils_tree.find_tax2root(info,parents, taxid) # [2711156, 1649486, 41297, 204457, 28211, 1224, 3379134, 2, 131567, 1]
#             if tax2root_list!=-1:
#                 depth=len(tax2root_list)
#             else:
#                 depth=0
#             if depth<3:
#                 tocheck_=taxid
#             read_tax_depth[case][row['read_name']]=depth
            
    
#     logging.debug("Number of reads in read_tax_depth: "+str(len(read_tax_depth))+" "+str(len(read_tax_depth[case])))
#     logging.debug("Counter of read_tax_depth: "+str(Counter(read_tax_depth[case].values())))

#     return read_tax_depth   


# def read_kraken_classification(cases, truth_file, classification_folder):
    
#     #cases=cases[:-1]; for k in range(20,26): cases.append("k"+str(k)+"l15")
#     logging.debug("Number of cases of kmer sizes: "+str(len(cases))+" "+str(cases))
#     tool_res_files=[]
#     logging.debug(classification_folder)
#     for case in cases: 
#         f=classification_folder+case+"_out.csv"
#         tool_res_files.append(f)
#     len(tool_res_files)
#     merged = pd.read_csv(truth_file, names=["read_name", "taxon"])
#     for f in tool_res_files:
#         c=f.split("/")[-1][1:3]
#         print(c)
#         merged_ = pd.read_csv(f, names=["read_name", "taxon_tool_k"+c],sep=",")
#         print(f, len(merged_))
#         merged = pd.merge(merged, merged_, on="read_name")
#     merged.taxon = merged.taxon.astype(int)
#     # merging only gets the intersection of all (which is fine, files should have all)

#     column_names = list(merged.columns)
#     logging.debug("Number of reads in merged: "+str(len(merged))+" "+str(column_names))
#     return merged 




# def calculate_tp_fp(dfmerged_taxa,info,tree_df,parents,tax_level,col_name,tax_index):

#     #print(tax_level)
#     cases_dic={"TP":set(),"FP-level-notindex":set(), 'FP-higher-notindex':set(),"FP-level-index":set(), 'FP-higher-index':set(),
#                'inconsistent-tool':set(), 'TN':set(),
#                'FN':set(),'VP':set(), "truth-level":set(), "no-truth-level":set(),"total_unclassified":set()} #  at the level ,'taxNotfoundintool':[]
#     num_reads=len(dfmerged_taxa["taxon"])
#     for j in range(num_reads):
#         #if j%200000==0: print("working at row",j," out of",num_reads, "total reads")
        
#         tid_true = dfmerged_taxa["taxon"][j]
#         tid_tool = dfmerged_taxa[col_name][j]
#         if tid_tool==1:
#             tid_tool=2
        
#         tid_true_taxlevel = _utils_tree.find_tax_level(info,tree_df,parents, tid_true, tax_level)
#         tid_tool_taxlevel = _utils_tree.find_tax_level(info,tree_df,parents, tid_tool, tax_level)
    
#         if tid_tool_taxlevel==-1 and tid_tool!=0:
#             cases_dic['inconsistent-tool'].add(dfmerged_taxa['read_name'][j])
#             # tid_tool_taxlevel==-1 we are ignoring these reads if the tax is not in the node.dump 
#             # we think these cases are limited and probably happened due to different node.dump in build vs evaluation 
    
        
#         if tid_tool==0 :
#             cases_dic["total_unclassified"].add(dfmerged_taxa['read_name'][j])
#         if tid_tool==-1 or tid_true < 1:
#             print("to check ",j,tid_tool, tid_true, tid_true_taxlevel)

#             #to check  9400 632 0 -1
    
#         if tid_true_taxlevel ==1 or  tid_true_taxlevel==-1 or (tid_tool_taxlevel==-1 and tid_tool!=0): # true is not at the level, is higher than the level     or true doest not exist in info
#             #including 'inconsistent-tool'
#             cases_dic["no-truth-level"].add(dfmerged_taxa['read_name'][j])
            
#         else:
#             cases_dic["truth-level"].add(dfmerged_taxa['read_name'][j])
    
#         if tid_true_taxlevel >1 and tid_tool==0:
#             if tid_true in tax_index:
#                 cases_dic['FN'].add(dfmerged_taxa['read_name'][j]) # unclassified, among reads with true label at the level 
#             else:
#                 cases_dic['TN'].add(dfmerged_taxa['read_name'][j]) # unclassified, among reads with true label at the level 
            
#         if tid_true_taxlevel>1 and tid_tool_taxlevel==tid_true_taxlevel :
#             cases_dic['TP'].add(dfmerged_taxa['read_name'][j]) # correct
    
#         if tid_true_taxlevel>1 and tid_tool>1 and tid_tool_taxlevel==1 :  # tool labled at higher than taxlevel 
#             tool_level= tree_df[4][info[tid_tool]]
#             tid_true_toollevel = _utils_tree.find_tax_level(info,tree_df,parents, tid_true, tool_level)
#             if tid_true_toollevel == tid_tool: # at the tool level is correct 
#                 cases_dic['VP'].add(dfmerged_taxa['read_name'][j]) # vague positive, one deeper lvel (ancestror)  match    
#             if tid_true_toollevel != tid_tool:
#                 if tid_true in tax_index:
#                     cases_dic['FP-higher-index'].add(dfmerged_taxa['read_name'][j]) # tool labeld at higher level and worng (not a VP) 
#                 else:
#                     cases_dic['FP-higher-notindex'].add(dfmerged_taxa['read_name'][j]) # tool labeld at higher level and worng (not a VP) 
    
#         if tid_true_taxlevel > 1 and tid_tool_taxlevel>1 and  tid_tool_taxlevel != tid_true_taxlevel :
#             if tid_true in tax_index:
#                 cases_dic['FP-level-index'].add(dfmerged_taxa['read_name'][j]) # tool labeled at the tax level (or lower) and wrong 
#             else:
#                 cases_dic['FP-level-notindex'].add(dfmerged_taxa['read_name'][j]) # tool labeled at the tax level (or lower) and wrong             
               
#     return cases_dic



# def find_true_ksize(cases,merged,info,tree_df,parents,tax_index, readids_max, level='species'):
#     tp_cases_dic={}
#     for case in cases:
#         case1= 'taxon_tool_'+case
#         tp_cases_dic[case]= calculate_tp_fp(merged,info,tree_df,parents,level,case1,tax_index) #+case
#     print(len(tp_cases_dic),len(tp_cases_dic[case]))

#     true_k={}
#     for read_idx, read_id in enumerate(readids_max):
#         true_k[read_id]=set()
#         for case_k, case in  enumerate(cases): # reversed(
#             k = int(case[1:3])
#             if read_id in tp_cases_dic[case]['TP']:
#                 true_k[read_id].add(k)
#     print(len(true_k),len(true_k[read_id]))

#     return tp_cases_dic, true_k


def read_kraken_file(kraken_file):

    # "C"/"U": a one letter code indicating that the sequence was either classified or unclassified.
    # The sequence ID, obtained from the FASTA/FASTQ header.
    # The taxonomy ID Kraken 2 used to label the sequence; this is 0 if the sequence is unclassified.
    # The length of the sequence in bp. In the case of paired read data, this will be a string containing the lengths of the two sequences in bp, separated by a pipe character, e.g. "98|94".
    # A space-delimited list indicating the LCA mapping of each  k-mer in the sequence(s).
    # For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
    # the first 13 k-mers mapped to taxonomy ID #562 ... 31  ambiguous nucleotide (A)
    file_o=open(kraken_file,'r')
    dic_kraken={}
    for line in file_o:
        line_split=line.strip().split("\t")
        classified ,read_id, tax_krak, read_len, tax_kmers = line_split
        if classified=="C":
            
            tax_kmer_dic = {}
            tax_kmer_num_dic={}
            pos=-1
            for tax_kmer in tax_kmers.split(" "):
                tax_, num_raw =  tax_kmer.split(":")
                num=int(num_raw)
                pos+=num     # how many consecutive kmer has this taxid # 0:6 543:1 91347:1 six unclassified, 
                if tax_ =='0' or tax_ =='A':
                    continue
                tax=int(tax_)
                if tax in tax_kmer_dic:
                    tax_kmer_dic[tax].append( (pos,num))
                    tax_kmer_num_dic[tax]+=num
                else:    
                    tax_kmer_dic[tax]=[(pos,num)]
                    tax_kmer_num_dic[tax]=num
            kraken_read_info = int(tax_krak), int(read_len), tax_kmer_dic, tax_kmer_num_dic
        elif classified=="U":
            
            kraken_read_info = 0, int(read_len), {}, {}
        else:
            print("Error: ",line_split)
            kraken_read_info = 0, 0, {}, {}
        
        dic_kraken[read_id] = kraken_read_info 
    file_o.close()
    return dic_kraken


# def read_kraken_set(file,read_names):
#     file_o=open(file,'r')
#     dic_krak={}
#     for line in file_o:
#         line_split=line.strip().split("\t")
#         if line_split[0]=="C":
#             read_id, tax_krak, read_len, tax_kmers = line_split[1:]
#             if read_id in read_names:
#                 tax_kmer_dic={}
#                 tax_kmer_num_dic={}
#                 pos=-1
#                 for tax_kmer in tax_kmers.split(" "):
                    
#                     tax_, num_raw =  tax_kmer.split(":")
                    
#                     num=int(num_raw)
#                     pos+=num     # how many consecutive kmer has this taxid # 0:6 543:1 91347:1 six unclassified, 
#                     if tax_=='0' or tax_ =='A':
#                         continue

#                     tax=int(tax_)
#                     if tax in tax_kmer_dic:
#                         tax_kmer_dic[tax].append( (pos,num))
#                         tax_kmer_num_dic[tax]+=num
#                     else:    
#                         tax_kmer_dic[tax]=[(pos,num)]
#                         tax_kmer_num_dic[tax]=num
#                 dic_krak[read_id] = int(tax_krak), int(read_len), tax_kmer_dic, tax_kmer_num_dic
#     return dic_krak



# def read_kraken_num(file,num_reads=1000):
#     file_o=open(file,'r')
#     dic_krak={}
#     read_id_cnt=0
#     for line in file_o:
#         line_split=line.strip().split("\t")
#         if line_split[0]=="C":
#             read_id, tax_krak, read_len, tax_kmers = line_split[1:]
#             read_id_cnt+=1
#             if read_id_cnt <num_reads:
#                 tax_kmer_dic={}
#                 tax_kmer_num_dic={}
#                 pos=-1
#                 for tax_kmer in tax_kmers.split(" "):
                    
#                     tax_, num_raw =  tax_kmer.split(":")
                    
#                     num=int(num_raw)
#                     pos+=num     # how many consecutive kmer has this taxid # 0:6 543:1 91347:1 six unclassified, 
#                     if tax_=='0' or tax_ =='A':
#                         continue

#                     tax=int(tax_)
#                     if tax in tax_kmer_dic:
#                         tax_kmer_dic[tax].append( (pos,num))
#                         tax_kmer_num_dic[tax]+=num
#                     else:    
#                         tax_kmer_dic[tax]=[(pos,num)]
#                         tax_kmer_num_dic[tax]=num
#                 dic_krak[read_id] = int(tax_krak), int(read_len), tax_kmer_dic, tax_kmer_num_dic
#     return dic_krak



def read_kraken_all(cases, classification_folder): # , readids_max, num_reads=10000
    logging.info("Reading kraken's k-mer count per tax from: "+classification_folder)
    #folder="/vast/blangme2/smajidi5/metagenomics/changek/simulatation/classification/" # 
    #cases=['k19','k25','k31'] # ,'k21'
    kraken_kmers_cases={}
    for case in cases: 
        logging.info('Reading kraken file for case: '+case)
        kraken_kmers_cases[case]=  read_kraken_file(classification_folder+case+"_out")
        logging.info("Number of reads for case: "+case+" is: "+str(len(kraken_kmers_cases[case])))
    
    for case_k, case in  enumerate(cases): 
        read_names_case=set(kraken_kmers_cases[case].keys())
        if case_k==0:
            read_names_case0= read_names_case
        else:
            if read_names_case0 != read_names_case:
                logging.warning("we expect the same reads in all cases. The "+case+" is not the same as the case "+case+str(len(read_names_case0))+" "+str(len(read_names_case))) 

    return list(read_names_case0), kraken_kmers_cases

def traverse_from_root(current, tree, num_kmers_path, tax_kmer_num_dic, child2parent):
    #print(current)  
    num_kmers_node = tax_kmer_num_dic.get(current,0)
    if current in child2parent:
        parent = child2parent[current]
        if num_kmers_node:
            if num_kmers_path[parent][0]: # if parent kmer is nonzero 
                num_kmers_path_upto_parent = num_kmers_path[parent][0]
            else:
                last_nonzero_parent= num_kmers_path[parent][1]
                num_kmers_path_upto_parent= num_kmers_path[last_nonzero_parent][0] #print('here',current, parent, last_nonzero_parent)
            num_kmers_path[current] = (num_kmers_path_upto_parent + num_kmers_node, current)
            #print('current:', current,' parent:',parent,' num_parrent:', num_kmers_path[parent][0],' num_current_path:', num_kmers_path[current])
        else:
            last_nonzero_parent=num_kmers_path[parent][1]
            num_kmers_path[current] = (0, last_nonzero_parent)
        
    else:  # this is root
        num_kmers_path[current] = (num_kmers_node, current)
    if current in tree: 
        children = tree[current]
        for child in children:
            traverse_from_root(child, tree, num_kmers_path,tax_kmer_num_dic, child2parent)
    # else: current is a leaf, no child 
    return 1



def calculate_num_kmers_path(Tree, tax_kmer_num_dic, tax_index):
    Tree_updated ={node:node_set.intersection(tax_index) for node,node_set in  Tree.items() if  node in tax_index}

    if 1 in Tree_updated and 1 in Tree_updated[1]:
        Tree_updated[1].remove(1)

    root=1
    child2parent = { }
    for parent, children in Tree_updated.items():
        for child in children:
            child2parent[child]=parent
    # logging.debug('number of leaves',len(child2parent))

    num_kmers_path={}
    traverse_from_root(root, Tree_updated, num_kmers_path,tax_kmer_num_dic, child2parent)
    #print('number of nodes',len(num_kmers_path))

    traverse_from_root(1, Tree_updated,num_kmers_path,tax_kmer_num_dic, child2parent) #, num_kmers_path, tax_kmer_num_dic)
    num_kmers_path_updated ={node:accum for node, (accum,nn) in num_kmers_path.items() if accum!=0}

    #len(Tree), len(tax_kmer_num_dic), len(num_kmers_path), len(num_kmers_path_updated)
    # Tree2={1:{2,3,8},2:{4,5},3:{6,7},8:{11,12}}
    # child2parent={}
    # for parent, children in Tree2.items():
    #     for child in children:
    #         child2parent[child]=parent
    # print('number of leaves',len(child2parent))

    # num_kmers_path={}
    # tax_kmer_num_dic={i:i*10 for i in range(1,13)}
    # tax_kmer_num_dic[11]=0
    # tax_kmer_num_dic[12]=0

    return num_kmers_path_updated


