pofotools=/home/jhu/jkim789/scratch_aphill27/software/Pooch/pofo_tools/src/pofo_tools
filter_model=$(realpath ../*pkl)
rob=//path/to/pooch/test/rob/methylation_analysis
xist=$(realpath $rob/chrX_promoter_methylation_rob.tsv)

mkdir -p prediction_xist

head $xist
head $rob/sex.txt
head $rob/rob.fofn

python $pofotools/createMethylDataset_xist.py \
    --chr chrX \
    --input-fofn $rob/rob.fofn \
    --output-file prediction_xist/chrX.createMethylDataset.out \
    --log-file prediction_xist/chrX.createMethylDataset.log \
    --sex-mapping-file $rob/sex.txt \
    --filter-model $filter_model

python $pofotools/createDiffDatasetFromMethylDataset.py \
    --input-methyl-dataset prediction_xist/chrX.createMethylDataset.out \
    -L Hap1 -R Hap2 \
    -c chrX \
    -l prediction_xist/chrX.createDiffDatasetFromMethylDataset.log \
    --output-file prediction_xist/chrX.createDiffDatasetFromMethylDataset.out
    

python $pofotools/predict_parent_origin_xist.py \
    -i prediction_xist/chrX.createDiffDatasetFromMethylDataset.out \
    -M $filter_model \
    --xist $xist\
    -o chrX.predictions.tsv 


#swap: default=['HG02451', 'HG02717', 'HG03579'],