# MetaKPick


# Work In Progress


A tool for metagenomic classifciation


## installtion



```
conda create -n metak python=3.10
conda activate metak

conda install bioconda::kraken
conda install anaconda::scikit-learn
conda install numpy 


```

## usage

First run Kraken1 with a few k-mer sizes:

```
$ python MetaKPick/metakpick/metakpick.py --help
usage: metakpick.py [-h] [--mode MODE] [--kraken_output_folder KRAKEN_OUTPUT_FOLDER] [--truth_file TRUTH_FILE]
                    [--model_folder MODEL_FOLDER] [--output_file_name OUTPUT_FILE_NAME] [--tree_file TREE_FILE]
                    [--tax_genome_file TAX_GENOME_FILE]

Metakpick: for metagenomic classification of reads

options:
  -h, --help            show this help message and exit
  --mode MODE           train or classify
  --kraken_output_folder KRAKEN_OUTPUT_FOLDER
                        Kraken output folder
  --truth_file TRUTH_FILE
                        Truth file
  --model_folder MODEL_FOLDER
                        Model folder
  --output_file_name OUTPUT_FILE_NAME
                        Output file name
  --tree_file TREE_FILE
                        Tree file
  --tax_genome_file TAX_GENOME_FILE
                        Tax genome file
```



### Classification

Metagenomics classification using pre-trained Random Forest model (from simulated genomes)
```
python metakpick.py --mode classify  --kraken_output_folder $krk  --truth_file $truth  --model_folder $files   --output_file_name $files/estimated_tax_cami.csv  --tree_file  $files/nodes.dmp   --tax_genome_file $files/seqid2taxid.map_tax_uniq
```


### Train

```
python metakpick.py --mode train --kraken_output_folder $krk --truth_file $truth  --model_folder $files --tree_file  $files/nodes.dmp --tax_genome_file $files/seqid2taxid.map_tax_uniq

```