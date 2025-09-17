# MetaKPick


# Work In Progress


A tool for metagenomic classifciation


## installtion

You can install the dependencies using conda: 
```
conda create -n metak python=3.10
conda activate metak
conda install bioconda::kraken
conda install anaconda::scikit-learn
conda install numpy 
```

## Step 1. Kraken build

First download the genomes and the taxanomy from Kraken index zone:
```
wget  https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250714.tar.gz
```
Then you can build krakn indexes using the following
```
for g in $(ls -1 genome); do 
kraken-build --add-to-library genome/${g} --db db
done 

threads=10
for k in 19 23 27 31; do
cp -r db db_${k} 
kraken-build --build --db db_${k} --threads ${threads} --kmer-len ${k}
done
```

## Step 2. Kraken classification

reads=reads.fq
threads=10
mkdir kraken_out
for k in 19 23 27 31; do 
    kraken_db=db_${k}
    kraken --db $kraken_db $read --threads ${threads} --output kraken_out/k${k}_out 
echo "kraken classification is done"



## Step 3. Training the model
Using our previously simulated data, we know the truth taxanomic classification: 

```
kraken_out_dir="kraken_out"
truth="truth.csv"
model="metakpick_k19-k31.pkl"
node_file=nodes.dmp 
tax_genome=seqid2taxid.map
python metakpick.pk --mode train --kraken_out $kraken_out_dir  --truth_file $truth --model_file ${model}.pkl --tree_file ${node_file} --tax_genome_file ${tax_genome} --kmer_list 19,23,27,31
```


# Classifying sequencing reads 
We will soon also provide the kraken indexes and the trained model files to be used for classification.

First, you need to classify reads using with four kraken indexes as described in Step 2. Then, use the provided model to infer the final classification:

```
kraken_out_dir="kraken_out"
model="metakpick_k19-k31.pkl"
node_file=nodes.dmp 
tax_genome=seqid2taxid.map
python metakpick.pk --mode classify --kraken_out $kraken_out_dir --model_file ${model}.pkl --tree_file ${node_file} --tax_genome_file ${tax_genome} --kmer_list 19,23,27,31 --output_file_name classification.csv
```









### All arguments

```
usage: metakpick.py [-h] [--mode MODE] [--kraken_out KRAKEN_OUT] [--truth_file TRUTH_FILE] [--model_file MODEL_FILE]
                    [--output_file_name OUTPUT_FILE_NAME] [--tree_file TREE_FILE] [--tax_genome_file TAX_GENOME_FILE]
                    [--kmer_list KMER_LIST] [--plot_tree] [--path_feature] [--version] [--version_decision VERSION_DECISION]
                    [--thr_minprob THR_MINPROB] [--thr_highprob_lca THR_HIGHPROB_LCA] [--thr_minprob_genus THR_MINPROB_GENUS] [--topRF]

Metakpick: for metagenomic classification of reads

options:
  -h, --help            show this help message and exit
  --mode MODE           train or classify
  --kraken_out KRAKEN_OUT
                        Kraken output folder
  --truth_file TRUTH_FILE
                        Truth file
  --model_file MODEL_FILE
                        Model file
  --output_file_name OUTPUT_FILE_NAME
                        Output file name
  --tree_file TREE_FILE
                        Tree file
  --tax_genome_file TAX_GENOME_FILE
                        Tax genome file
  --kmer_list KMER_LIST
                        K-mer list file comma separated
  --plot_tree           Plot the tree
  --path_feature        To use path feature, slower
  --version             show program's version number and exit
  --version_decision VERSION_DECISION
                        Version decision
  --thr_minprob THR_MINPROB
                        with version_decision 1/2/3, threshold for minimum probability to decide unclassified reads
  --thr_highprob_lca THR_HIGHPROB_LCA
                        with version_decision 2, threshold for high probability to decide lca
  --thr_minprob_genus THR_MINPROB_GENUS
                        with version_decision 3,threshold for minimum probability to decide genus
  --topRF               TopRF mode
```