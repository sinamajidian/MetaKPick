import os
import sys
import pandas as pd
import sys
import numpy as np

"""
python eval_kraken  kraken_out.1st4th.csv  true.csv nodes.dmp taxa.ids > out .stat 
For comparing kraken output with truth and calculate F1 when true reads defined at species level
"""



#folder="/vast/blangme2/smajidi5/metagenomics/movi_classif/eval/"
#tax_index_file=folder+"taxa.ids"
#tree_file =folder+"nodes.dmp" #sys.argv[3] #"/vast/blangme2/mzakeri1/Move/steven/classification_data_pseudomonadota/kraken2_library/taxonomy/nodes.dmp"
#tax_genome_index_file= folder+"bacteria_kraken2_library.matching_ids_taxuniq"

#truth_file=folder+"reads_mapping_0_long.csv"
#input1= "CAMI_pseudo_out_bigmem_cutoff1_1threads_slurm_bigmem_rep0.csv"
#tool_res_file=folder+input1

tool_res_file=sys.argv[1]
truth_file=sys.argv[2]
tree_file =sys.argv[3] #"/vast/blangme2/mzakeri1/Move/steven/classification_data_pseudomonadota/kraken2_library/taxonomy/nodes.dmp"



tax_index=set()
if len(sys.argv)>4:
    tax_index_file=sys.argv[4] #"/vast/blangme2/mzakeri1/Move/steven/classification_data_pseudomonadota/taxa.ids"
    tax_index_f=open(tax_index_file,'r')
    tax_index=set()
    for line in tax_index_f:
        tax_index.add(int(line.strip()))
    print(len(tax_index))



import pandas as pd
import sys
import numpy as np




def parse_tree_file(tree_df,children,parents):
    # create a dictionary from taxonomic id to index of corresponding row
    info = {}
    for i in range(len(children)):
        info[children[i]] = i

    # Convert to tree structure (dictionary of sets)
    Tree = {}
    for i in range(len(tree_df[0])):
        parent = parents[i]
        child = children[i]
        if parent in Tree:
            Tree[parent].add(child)
        else:
            Tree[parent] = {child}
    return info, Tree  
def find_tax_level(info,tree_df,parents, tid, tax_level="species"):
    if tid not in info:
        #print(tid)
        return -1
    while tid != 1:    
        tid_row = info[tid]
        if tree_df[4][tid_row] == tax_level:
            break
        else:
            tid = parents[info[tid]] #find_parent_node(tid)
    return tid



print("Tree file: ", tree_file)
tree_df = pd.read_csv(tree_file, sep="\t", header=None)
print("Read the tree file")
children = tree_df[0]
parents = tree_df[2]
info, Tree = parse_tree_file(tree_df,children,parents)
print("Parsed the tree file")
print("Size of info: ", len(info))
print("Size of Tree: ", len(Tree))

truth = pd.read_csv(truth_file, names=["read_name", "taxon"])




#tax_genome_index_f_=open(tax_genome_index_file,'r')
#tax_geneme=set()
#for line in tax_genome_index_f_:
#    tax_geneme.add(int(line.strip()))
#print(len(tax_geneme))
#tax_geneme_species = set()
#for tax_geneme_i in  tax_geneme:
#    tax_geneme_species.add(int(find_tax_level(info,tree_df,parents, tax_geneme_i, 'species')))
#len(tax_geneme),len(tax_geneme_species)



tool_res = pd.read_csv(tool_res_file, names=["read_name", "taxon_tool"])

merged = pd.merge(truth, tool_res, on="read_name")
total_count = len(merged["taxon_tool"])
merged

merged.taxon_tool = merged.taxon_tool.astype(int)
merged.taxon = merged.taxon.astype(int)


tid_true_taxlevel_dic={}
for tax_level in  [ "species","genus", "family", "order", "class"]:
    tid_true_taxlevel_dic[tax_level]={}
    for tid_true in  set(truth['taxon']):
        tid_true_taxlevel_dic[tax_level][tid_true] = int(find_tax_level(info,tree_df,parents, tid_true, tax_level))





tid_tool_taxlevel_dic={}
for tax_level in  [ "species","genus", "family", "order", "class"]:
    tid_tool_taxlevel_dic[tax_level]={}
    for tid_tax in set(list(tool_res['taxon_tool'])+[2]):
        if tid_tax=='lca': continue # first row
        tid_tax=int(tid_tax)
        tid_tool_taxlevel_dic[tax_level][tid_tax] = int(find_tax_level(info,tree_df,parents, tid_tax, tax_level))

dic_if_species={}
for tax_true in set(truth['taxon']):
    tid_true_specieslevel = find_tax_level(info,tree_df,parents, tax_true, 'species')
    if tid_true_specieslevel==1 or tid_true_specieslevel==-1:
        dic_if_species[tax_true]='to_ignore'
    else:
        dic_if_species[tax_true]='species_or_lower'
len(dic_if_species)        

#dic_if_species_notinindex={}
#for tid_true in set(truth['taxon']):
    #tid_true_specieslevel = tid_true_taxlevel_dic['species'][tid_true]
    #dic_if_species_notinindex[tid_true]="A"
    #if  tid_true_specieslevel==1 or tid_true_specieslevel==-1:        
    #if tid_true_specieslevel not in  tax_geneme_species:
    #    dic_if_species_notinindex[tid_true]='to_ignore'

#len(dic_if_species_inindex), sum([1 for i in list(dic_if_species_notinindex.values()) if i=='A']),sum([1 for i in list(dic_if_species_notinindex.values()) if i=='to_ignore'])




for tax_level in  [ "species","genus", "family", "order", "class"]:
    #tax_level="genus" # species

    print("h",tax_level)
    cases_dic={"truth-species-not-in-index":[], "truth-species":[],"no-truth-species":[] ,"TP":[],"FP-level-notindex":[], 'FP-higher-notindex':[],"FP-level-index":[], 'FP-higher-index':[],
               'inconsistent-tool':[], 'TN':[],
               'FN':[],'VP':[], "truth-level":[] , "no-truth-level":[],"total_unclassified":[]} #  at the level ,'taxNotfoundintool':[]
    num_reads=len(merged["taxon_tool"])
    for j in range(num_reads):
        #if j%200000==0: print("working at row",j," out of",num_reads, "total reads")

        tid_true = merged["taxon"][j]
        tid_tool = merged["taxon_tool"][j]
        if tid_tool==1:
            tid_tool=2

        tid_true_taxlevel =  tid_true_taxlevel_dic[tax_level][int(tid_true)] #find_tax_level(info,tree_df,parents, tid_true, tax_level)
        tid_tool_taxlevel =  tid_tool_taxlevel_dic[tax_level][int(tid_tool)] #find_tax_level(info,tree_df,parents, tid_tool, tax_level)

        if tid_tool_taxlevel==-1 and tid_tool!=0:
            cases_dic['inconsistent-tool'].append(j)
            # tid_tool_taxlevel==-1 we are ignoring these reads if the tax is not in the node.dump
            # we think these cases are limited and probably happened due to different node.dump in build vs evaluation


        if tid_tool==0 :
            cases_dic["total_unclassified"].append(j)
        if tid_tool==-1 or tid_true < 1:
            print("to check ",j,tid_tool, tid_true, tid_true_taxlevel)

        if dic_if_species[tid_true]=='to_ignore': #tid_true_taxlevel ==1 or  tid_true_taxlevel==-1 or (tid_tool_taxlevel==-1 and tid_tool!=0): # true is not at the level, is higher than the level     or true doest not exist in info
            #including 'inconsistent-tool'
            cases_dic["no-truth-species"].append(j)
            continue             
        else:
            # if dic_if_species_notinindex[tid_true]=='to_ignore':
            #     cases_dic["truth-species-not-in-index"].append(j)
            #     continue

            cases_dic["truth-species"].append(j)

        if tid_true_taxlevel >1 and tid_tool==0:
            if tid_true in tax_index:
                cases_dic['FN'].append(j) # unclassified, among reads with true label at the level
            else:
                cases_dic['TN'].append(j) # unclassified, among reads with true label at the level

        if tid_true_taxlevel>1 and tid_tool_taxlevel==tid_true_taxlevel :
            cases_dic['TP'].append(j) # correct

        if tid_true_taxlevel>1 and tid_tool>1 and tid_tool_taxlevel==1 :  # tool labled at higher than taxlevel
            tool_level= tree_df[4][info[tid_tool]]
            tid_true_toollevel = find_tax_level(info,tree_df,parents, tid_true, tool_level)
            if tid_true_toollevel == tid_tool: # at the tool level is correct
                cases_dic['VP'].append(j) # vague positive, one deeper lvel (ancestror)  match
            if tid_true_toollevel != tid_tool:
                if tid_true in tax_index:
                    cases_dic['FP-higher-index'].append(j) # tool labeld at higher level and worng (not a VP)
                else:
                    cases_dic['FP-higher-notindex'].append(j) # tool labeld at higher level and worng (not a VP)

        if tid_true_taxlevel > 1 and tid_tool_taxlevel>1 and  tid_tool_taxlevel != tid_true_taxlevel :
            if tid_true in tax_index:
                cases_dic['FP-level-index'].append(j) # tool labeled at the tax level (or lower) and wrong
            else:
                cases_dic['FP-level-notindex'].append(j) # tool labeled at the tax level (or lower) and wrong

    print("total reads:,"+str(num_reads))
    # 'no-truth-species',
    for case  in ['truth-species', 'TN','FN','TP', 'VP','FP-level-index','FP-higher-index','FP-level-notindex','FP-higher-notindex']:
        list1= cases_dic[case]
        print(case +","+str(len(list1)))

    FP=len(cases_dic['FP-level-index'])+len(cases_dic['FP-higher-index']+cases_dic['FP-level-notindex'])+len(cases_dic['FP-higher-notindex'])

    
    recall=len(cases_dic['TP'])/(len(cases_dic['TP']) + len(cases_dic['VP']) + len(cases_dic['FN']) +FP )
    print("recall"+","+str(np.round(recall,4))) # sensitivity
    precision=len(cases_dic['TP'])/(len(cases_dic['TP']) + FP)
    print("precision"+","+str(np.round(precision,4))) # sensitivity
    if precision+recall:
        F1= 2* precision* recall/(precision+recall)
        print("F1"+","+str(np.round(F1,4))) # sensitivity

    print("TP+FP+VP+TN+FN"+","+str(len(cases_dic['TN']) +len(cases_dic['TP'])+FP+len(cases_dic['FN']) +len(cases_dic['VP'])))
    for case  in ['total_unclassified', 'inconsistent-tool']:
        list1= cases_dic[case]
        print(case +","+str(len(list1)))


print("Note that when there last argument is empty, FN/TN,  not index should be considered the opposite ")
