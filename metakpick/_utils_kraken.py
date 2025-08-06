

def read_kraken_set(file,read_names):

    # "C"/"U": a one letter code indicating that the sequence was either classified or unclassified.
    # The sequence ID, obtained from the FASTA/FASTQ header.
    # The taxonomy ID Kraken 2 used to label the sequence; this is 0 if the sequence is unclassified.
    # The length of the sequence in bp. In the case of paired read data, this will be a string containing the lengths of the two sequences in bp, separated by a pipe character, e.g. "98|94".
    # A space-delimited list indicating the LCA mapping of each  k-mer in the sequence(s).
    # For example, "562:13 561:4 A:31 0:1 562:3" would indicate that:
    # the first 13 k-mers mapped to taxonomy ID #562 ... 31  ambiguous nucleotide (A)
    file_o=open(file,'r')
    dic_krak={}
    for line in file_o:
        line_split=line.strip().split("\t")
        if line_split[0]=="C":
            read_id, tax_krak, read_len, tax_kmers = line_split[1:]
            if read_id in read_names:
                tax_kmer_dic={}
                tax_kmer_num_dic={}
                pos=-1
                for tax_kmer in tax_kmers.split(" "):
                    
                    tax_, num_raw =  tax_kmer.split(":")
                    
                    num=int(num_raw)
                    pos+=num     # how many consecutive kmer has this taxid # 0:6 543:1 91347:1 six unclassified, 
                    if tax_=='0' or tax_ =='A':
                        continue

                    tax=int(tax_)
                    if tax in tax_kmer_dic:
                        tax_kmer_dic[tax].append( (pos,num))
                        tax_kmer_num_dic[tax]+=num
                    else:    
                        tax_kmer_dic[tax]=[(pos,num)]
                        tax_kmer_num_dic[tax]=num
                dic_krak[read_id] = int(tax_krak), int(read_len), tax_kmer_dic, tax_kmer_num_dic
    return dic_krak



def read_kraken_set(folder, cases, readids_max):
    print("read kraken's k-mer  count per tax")
    #folder="/vast/blangme2/smajidi5/metagenomics/changek/simulatation/classification/" # 
    #cases=['k19','k25','k31'] # ,'k21'
    dic_cases={}
    for case in reversed(cases): # 
        print(case)
        try:
            #dic_cases[case]=read_kraken_limited(folder+case+"_out",10000)
            #dic_cases[case]=read_kraken(folder+case+"_out")
            dic_cases[case]=read_kraken_set(folder+case+"_out",readids_max)
        except:
            print("n",case)
        print(case,len(dic_cases[case]))

    cases_readids=set()
    for case_k, case in  enumerate(cases): 
        read_ids_k=set(dic_cases[case].keys())
        cases_readids |= read_ids_k
    len(cases_readids)
    return dic_cases, cases_readids
    