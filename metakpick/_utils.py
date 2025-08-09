



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


