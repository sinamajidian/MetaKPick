



def write_estimated_tax(estimated_tax_dict,output_file_name="estimated_tax.csv"):
    output_file =open(output_file_name,"w")
    for read_id, estimated_tax in estimated_tax_dict.items():
        output_file.write(read_id+","+str(estimated_tax)+"\n")
    output_file.close()

    return output_file_name


