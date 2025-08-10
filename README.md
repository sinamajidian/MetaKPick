# MetaKPick



A tool for metagenomic classifciation


## installtion

WIP

```
conda create -n metak python=3.10
conda activate metak

conda install bioconda::kraken
conda install anaconda::scikit-learn
conda install numpy 

python install setup.py
 
```

## usage




### Classification

using trained RF model


```
bash metakpick.sh classify 
```



### train the model

```
metakpick  --mode train 
```


