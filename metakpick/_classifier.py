

def apply_RF_model(cases, features_cases,read_names_list,regr_dic,n_estimators, max_features, max_leaf_nodes, random_state, n_jobs=1):
    read_k_prob={}
    for read_name in read_names_list:
        read_k_prob[read_name]=np.zeros(len(cases))    
    regr_dic = {}
    y_pred_binary_={}
    for case_idx, case in enumerate(cases):
        X_input = features_cases[case]
        Y_input = tp_binary_reads_cases[case]

        regr_dic[case] = train_RF_model(X_input, Y_input, n_estimators, max_features, max_leaf_nodes, random_state, n_jobs)

        logging.debug("X_input.shape "+str(X_input.shape)+" len(Y_input) "+str(len(Y_input)))

        y_pred = regr_dic[case].predict(X_input)
        #y_pred_binary = np.round(y_pred)  # round(0.55)=1
        logging.debug("regr_dic[case]"+str(regr_dic[case]))

    for read_name_idx,read_name in enumerate(read_names_list):        
        read_k_prob[read_name][case_idx]= y_pred[read_name_idx]

    return read_k_prob


def calculate_accuracy(y_pred,tp_binary_reads_cases):

    accuracy=sum([1 for i in range(len(y_pred)) if y_pred_binary[i] == Y_input[i]]) / len(Y_input)
    logging.debug(" Accuracy of regression model for case "+case+" is "+str(accuracy))




def get_best_tax(read_k_prob,cases,read_names_list,merged,thr_minprob=0.5):
    best_k_dic={}
    for read_name_idx, read_name in enumerate(read_names_list):
        
        max_val = np.max(read_k_prob[read_name])

        if max_val >thr_minprob:
            best_k = cases[np.argmax(read_k_prob[read_name])]
        else:
            best_k=-1

        best_k_dic[read_name]=best_k                
    logging.debug("** best k"+str(Counter(best_k_dic.values())))
    estimated_tax_dict={}
    for index, row in merged.iterrows():
        read_name   =row['read_name']   
        if  read_name in  best_k_dic and best_k_dic[read_name]!=-1:
            estimated_tax= row['taxon_tool_'+best_k_dic[read_name]]
        else:
            estimated_tax=0
        estimated_tax_dict[read_name ]=estimated_tax

    return best_k_dic,estimated_tax_dict


