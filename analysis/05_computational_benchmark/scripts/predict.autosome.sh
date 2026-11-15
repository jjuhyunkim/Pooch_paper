pofotools=/home/jhu/jkim789/scratch_aphill27/software/Pooch/pofo_tools/src/pofo_tools

washu=/path/to/pooch/test/washu-pedigree/methylation_analysis
xist=$(realpath xist_level.tsv)
sex=/path/to/pooch/pooch_data/05.dss_using_123_pairing/sex.trio.tsv
chr=$1

if [[ $chr == "chrX" ]];then
fofn=HG01081.fofn
filter_model=$(realpath ../chrX*pkl)
else
fofn=HG002.fofn
filter_model=$(realpath ../elasticnet_bestmodel_${chr}.*.pkl)
fi

echo -e "fofn : $fofn"
echo -e "model : $filter_model"

mkdir -p prediction_xist

python $pofotools/createMethylDataset_xist.py \
    --chr $chr \
    --input-fofn $fofn \
    --output-file $chr.createMethylDataset.out \
    --log-file $chr.createMethylDataset.log \
    --sex-mapping-file $sex \
    --filter-model $filter_model

python $pofotools/createDiffDatasetFromMethylDataset.py \
    --input-methyl-dataset $chr.createMethylDataset.out \
    -L Maternal -R Paternal \
    -c $chr \
    -l $chr.createDiffDatasetFromMethylDataset.log \
    --output-file $chr.createDiffDatasetFromMethylDataset.out
    
if [[ $chr == "chrX" ]]; then
prediction_script=$pofotools/predict_parent_origin_xist.py 
python $prediction_script \
    -i $chr.createDiffDatasetFromMethylDataset.out \
    -M $filter_model \
    --xist $xist \
    -o $chr.predictions.tsv
else 
prediction_script=$pofotools/predict.py
python $prediction_script \
    -i $chr.createDiffDatasetFromMethylDataset.out \
    -M $filter_model \
    -o $chr.predictions.tsv 
fi