echo "step 1 downlaoding the kraken indexes"
# step1 download the  kraken indexes
# wget 

# alterantively build the kraken indexes from the genome files
#for g in $(ls -1 genome); do 
#kraken-build --add-to-library genome/${g} --db db  >> kraken_build.log 2>&1
#done 

# kraken-build --build --db $db --threads ${thr} --kmer-len ${k} --minimizer-len ${l}  > kraken_build_${db}.log 2>&1 

kraken_folder="kraken_db"



echo "step 2 running kraken classification"

read=$2 
threads= 5 #$2
for k in 21 23 25 27 29 31; do 
    kraken_db=${kraken_folder}/k${k}
    kraken --db $kraken_db $read --threads ${threads} --output kraken.${k}.classification.out > kraken.${k}.classification.log 2>&1 
echo "kraken classification is done"



# step 3 
echo "step 3 running metakpick"
python metakpick.py --mode classify --read $read --threads $threads






