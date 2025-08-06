



def def_cal(dfmerged_taxa,info,tree_df,parents,tax_level,col_name,tax_index):

    #print(tax_level)
    cases_dic={"TP":set(),"FP-level-notindex":set(), 'FP-higher-notindex':set(),"FP-level-index":set(), 'FP-higher-index':set(),
               'inconsistent-tool':set(), 'TN':set(),
               'FN':set(),'VP':set(), "truth-level":set(), "no-truth-level":set(),"total_unclassified":set()} #  at the level ,'taxNotfoundintool':[]
    num_reads=len(dfmerged_taxa["taxon"])
    for j in range(num_reads):
        #if j%200000==0: print("working at row",j," out of",num_reads, "total reads")
        
        tid_true = dfmerged_taxa["taxon"][j]
        tid_tool = dfmerged_taxa[col_name][j]
        if tid_tool==1:
            tid_tool=2
        
        tid_true_taxlevel = find_tax_level(info,tree_df,parents, tid_true, tax_level)
        tid_tool_taxlevel = find_tax_level(info,tree_df,parents, tid_tool, tax_level)
    
        if tid_tool_taxlevel==-1 and tid_tool!=0:
            cases_dic['inconsistent-tool'].add(dfmerged_taxa['read_name'][j])
            # tid_tool_taxlevel==-1 we are ignoring these reads if the tax is not in the node.dump 
            # we think these cases are limited and probably happened due to different node.dump in build vs evaluation 
    
        
        if tid_tool==0 :
            cases_dic["total_unclassified"].add(dfmerged_taxa['read_name'][j])
        if tid_tool==-1 or tid_true < 1:
            print("to check ",j,tid_tool, tid_true, tid_true_taxlevel)

            #to check  9400 632 0 -1
    
        if tid_true_taxlevel ==1 or  tid_true_taxlevel==-1 or (tid_tool_taxlevel==-1 and tid_tool!=0): # true is not at the level, is higher than the level     or true doest not exist in info
            #including 'inconsistent-tool'
            cases_dic["no-truth-level"].add(dfmerged_taxa['read_name'][j])
            
        else:
            cases_dic["truth-level"].add(dfmerged_taxa['read_name'][j])
    
        if tid_true_taxlevel >1 and tid_tool==0:
            if tid_true in tax_index:
                cases_dic['FN'].add(dfmerged_taxa['read_name'][j]) # unclassified, among reads with true label at the level 
            else:
                cases_dic['TN'].add(dfmerged_taxa['read_name'][j]) # unclassified, among reads with true label at the level 
            
        if tid_true_taxlevel>1 and tid_tool_taxlevel==tid_true_taxlevel :
            cases_dic['TP'].add(dfmerged_taxa['read_name'][j]) # correct
    
        if tid_true_taxlevel>1 and tid_tool>1 and tid_tool_taxlevel==1 :  # tool labled at higher than taxlevel 
            tool_level= tree_df[4][info[tid_tool]]
            tid_true_toollevel = find_tax_level(info,tree_df,parents, tid_true, tool_level)
            if tid_true_toollevel == tid_tool: # at the tool level is correct 
                cases_dic['VP'].add(dfmerged_taxa['read_name'][j]) # vague positive, one deeper lvel (ancestror)  match    
            if tid_true_toollevel != tid_tool:
                if tid_true in tax_index:
                    cases_dic['FP-higher-index'].add(dfmerged_taxa['read_name'][j]) # tool labeld at higher level and worng (not a VP) 
                else:
                    cases_dic['FP-higher-notindex'].add(dfmerged_taxa['read_name'][j]) # tool labeld at higher level and worng (not a VP) 
    
        if tid_true_taxlevel > 1 and tid_tool_taxlevel>1 and  tid_tool_taxlevel != tid_true_taxlevel :
            if tid_true in tax_index:
                cases_dic['FP-level-index'].add(dfmerged_taxa['read_name'][j]) # tool labeled at the tax level (or lower) and wrong 
            else:
                cases_dic['FP-level-notindex'].add(dfmerged_taxa['read_name'][j]) # tool labeled at the tax level (or lower) and wrong             
               
    return cases_dic




from sklearn.tree import _tree

def tree_to_code(tree, feature_names):
    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]
    print("def tree({}):".format(", ".join(feature_names)))

    def recurse(node, depth):
        indent = "  " * depth
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]
            print("{}if {} <= {}:".format(indent, name, threshold))
            recurse(tree_.children_left[node], depth + 1)
            print("{}else:  # if {} > {}".format(indent, name, threshold))
            recurse(tree_.children_right[node], depth + 1)
        else:
            print("{}return {}".format(indent, tree_.value[node]))

    recurse(0, 1)



def merge_kraken_results(folder, truth_file):

    cases=[i.split("_")[0] for i in os.listdir(folder) if i.endswith('_out')] #["k"+str(k) for k in range(15,32)]
    #cases=cases[:-1]; for k in range(20,26): cases.append("k"+str(k)+"l15")
    print(len(cases),cases)
    tool_res_files=[]
    print(folder)
    for case in cases: 
        f=folder+case+"_out.csv"
        tool_res_files.append(f)
    len(tool_res_files)
    merged = pd.read_csv(truth_file, names=["read_name", "taxon"])
    for f in tool_res_files:
        c=f.split("/")[-1][1:3]
        print(c)
        merged_ = pd.read_csv(f, names=["read_name", "taxon_tool_k"+c],sep=",")
        print(f, len(merged_))
        merged = pd.merge(merged, merged_, on="read_name")
    merged.taxon = merged.taxon.astype(int)
    # merging only gets the intersection of all (which is fine, files should have all)
    return merged




def get_readids_max(merged, num_max):
    num_max=15
    counter_tax={}
    readids_max=set()
    for index, row in merged.iterrows():
        tax_true= row['taxon']
        read_name=row['read_name']

        to_add=True
        if not tax_true in counter_tax:
            counter_tax[tax_true]=1
        else:
            if counter_tax[tax_true]>num_max:
                to_add=False 
            counter_tax[tax_true]+=1
            
                
        if to_add:
            readids_max.add(read_name)
    print(len(merged), len(readids_max),len(counter_tax), counter_tax[tax_true])
    return readids_max



def get_y_label(merged, cases):
    #read_names=read_name_list#[:1000
    #read_names=merged2a['read_name']
    #true_k={}
    dic_y2={}
    #case='k25'
    for case in cases:
        y2=[]    
        for read_idx, read_id in enumerate(read_names):
            true_k_list_read=[]
            #for case_k, case in  enumerate(cases): # reversed(
                #k = int(case[1:3])
            if read_id in case_dic_all[case]['TP'] :
                y2.append(1) 
            #elif  read_id in case_dic_all[case]['VP']:
            #    y2.append(2)
            else:
                y2.append(0)
        dic_y2[case] = y2
        
    print(len(dic_y2),len(y2),len(read_names),sum(y2))
    # for t, v in case_dic_all[case].items():
    #     print(t,len(v))
    for case in cases:
        print(case ,Counter(dic_y2[case]))
    return dic_y2

def get_features_kmer(tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic, num_nodes_tree, tax2depth):

            num_kmer_all=np.array(list(tax_kmer_num_dic.keys()),dtype=np.int32)
            num_kmer_norm=num_kmer_all/rlen # read_length_dic[read_id]
            
            if tax_krak in tax_kmer_num_dic:
                kmer_reported_tax = tax_kmer_num_dic[tax_krak] 
            else:
                kmer_reported_tax=0 # probably this kraken tax is the lca of a few tax closer to leaves
            
            mean_exc_tax = (np.sum(num_kmer_all)-kmer_reported_tax)/num_nodes_tree
            
            #feature_names=['mean_all','mean_nonzero','max','sum', 'mean_exc_tax']
            features=[np.mean(num_kmer_all), np.mean(num_kmer_all)*len(num_kmer_all)/num_nodes_tree, np.max(num_kmer_all),np.sum(num_kmer_all),mean_exc_tax]
            #feature_names+=['mean_all/rlen','mean_nonzero/rlen','max/rlen','sum/rlen', 'mean_exc_tax/rlen'] # /readlength
            features+=[np.mean(num_kmer_norm)/rlen, np.mean(num_kmer_norm)*len(num_kmer_norm)/num_nodes_tree/rlen, 
                    np.max(num_kmer_norm)/rlen,np.sum(num_kmer_norm)/rlen, mean_exc_tax/rlen]
            return features
        


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
            
        features_tax = [depth_reported,numkmer_consecutive_mean,weighted_depth_numkmer_mean ]
        #feature_names+=['depth_reported_tax','Avg_kmer_consecutive', 'WAvg_kmer_consecutive']
        return features_tax


def get_features_depth(tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic, num_nodes_tree, tax2depth):

    kmer_reported_uppertax=0
    kmer_reported_belowtax=0
    kmer_other_reported_tax=0
    kmer_other_reported_uppertax=0 # towards root (smaller depth value)
    kmer_other_reported_belowtax=0 # towards leaves
    
    tax_krak_2root= find_tax2root(info, parents, tax_krak)
    for tax_upper in tax_krak_2root[1:]: #excluding tax_krak towards the root  tax_krak_2root = [201174, 1783272, 2, 131567, 1] where tax_krak=201174 
        if tax_upper in tax_kmer_num_dic: 
            kmer_reported_uppertax +=tax_kmer_num_dic[tax_upper]
    
    below_tax_all_perread=[]
    for taxi in tax_kmer_num_dic:
        # to get the below for the reported 
        tax2root_t= find_tax2root(info, parents, taxi)
        if tax2root_t!=-1 and tax_krak in tax2root_t: # tax_krak is lca of taxi
            kmer_reported_belowtax += tax_kmer_num_dic[taxi]   
            below_tax_all_perread.append(taxi)
    
        # to get the below/above/at  the reported tax, accross all but not the max path 
        if depth_reported==tax_depth:
            kmer_other_reported_tax+=tax_kmer_num_dic[taxi]   
        elif depth_reported<tax_depth: # towards root (smaller depth value)
            kmer_other_reported_uppertax+=tax_kmer_num_dic[taxi]   
        elif depth_reported>tax_depth:
            kmer_other_reported_belowtax+=tax_kmer_num_dic[taxi]   
    
    features_tax2 =[kmer_reported_tax, kmer_reported_uppertax,kmer_reported_belowtax,
        kmer_reported_tax/rlen, kmer_reported_uppertax/rlen,kmer_reported_belowtax/rlen,
        kmer_other_reported_tax, kmer_other_reported_uppertax,kmer_other_reported_belowtax,
        kmer_other_reported_tax/rlen, kmer_other_reported_uppertax/rlen,kmer_other_reported_belowtax/rlen]
    return features_tax2
    

def get_features(read_id, dic_cases, kmer_size, tax2path, tax2depth, read_tax_depth, tax_index):
       tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic = dic_cases[case][read_id] # read length
        
       features = get_features_kmer(tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic, num_nodes_tree)

       features_tax = get_features_depth(tax_krak, rlen, tax_kmer_dic, tax_kmer_num_dic, num_nodes_tree, tax2depth):

       
       
        #feature_names+=['kmer_reported_tax','kmer_tax_above','kmer_tax_below','kmer_tax/rlen','kmer_tax_above/rlen','kmer_tax_below/rlen','kmer_othertax','kmer_othertax_above','kmer_othertax_below','kmer_othertax/rlen','kmer_othertax_above/rlen','kmer_othertax_below/rlen',]
        
        
        tax_krak_touse=tax_krak
        if tax_krak not in tax_kmer_dic: # when lca is reprorted, use onf its lower ones. 
            tax_krak_touse= below_tax_all_perread[0]
        
        segment_num=2 # 3
        segment_len=int(rlen/segment_num)
        pos_num_reportedtax = tax_kmer_dic[tax_krak_touse] # [(77, 1), (136, 4), (139, 2), (144, 4), (148, 2), (162, 2), (233, 1)]
        cnt_perbin=np.zeros((segment_num,1))
        for pos,numkmer in pos_num_reportedtax:
            bin_idx= int(pos/segment_len)
            cnt_perbin[bin_idx]+=numkmer
        diff_fromnext_seg = [ np.abs(cnt_perbin[segment_i+1]-cnt_perbin[segment_i])[0] for segment_i in range(segment_num-1)]
        
        features += [np.mean(diff_fromnext_seg)]
        #feature_names += ['diff#kmers_halfRead']

    return features
