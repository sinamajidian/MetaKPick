
import os
import pandas as pd

def read_kraken_classification(cases , truth_file, classification_folder):
    
    #cases=cases[:-1]; for k in range(20,26): cases.append("k"+str(k)+"l15")
    print(len(cases),cases)
    tool_res_files=[]
    print(classification_folder)
    for case in cases: 
        f=classification_folder+case+"_out.csv"
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

    column_names = list(merged.columns)
    print(len(merged),column_names )
    return merged 


def get_tax_depth(merged, cases, info,parents):
    read_tax_depth={}
    for case in cases: 
        read_tax_depth[case]={}
        tocheck_=""
        for index, row in merged.iterrows():
            taxid=row['taxon_tool_'+case]
            tax2root_list = find_tax2root(info,parents, taxid) # [2711156, 1649486, 41297, 204457, 28211, 1224, 3379134, 2, 131567, 1]
            if tax2root_list!=-1:
                depth=len(tax2root_list)
            else:
                depth=0
            if depth<3:
                tocheck_=taxid
            read_tax_depth[case][row['read_name']]=depth
            
    print(tocheck_)
    print(len(read_tax_depth),len(read_tax_depth[case]))
    print(Counter(read_tax_depth[case].values()))

    return read_tax_depth   
