
def read_len_dic(fastq_file):
    def process(lines=None):
        ks = ['name', 'sequence', 'optional', 'quality']
        return {k: v for k, v in zip(ks, lines)}
    n = 4
    #records_seq_len=[]
    read_length_dic={}
    with open(fastq_file, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == n:
                record = process(lines)
                #sys.stderr.write("Record: %s\n" % (str(record)))
                #records_seq_len.append(len(record['sequence']))
                read_length_dic[record['name'][1:]]=len(record['sequence'])
                lines = []
    #print(fn,sum(records_seq_len)/len(records_seq_len))
    print(len(read_length_dic),' reads')
    return read_length_dic


