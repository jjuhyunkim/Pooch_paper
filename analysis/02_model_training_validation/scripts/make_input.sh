#!/bin/bash

#SBATCH --partition=largemem
#SBATCH --cpus-per-task=50
#SBATCH --mem=3000g
#SBATCH --array=0,8
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/make_input.%j.log
#SBATCH --error=logs/make_input.%j.log
#SBATCH --job-name=make_input

chrs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
i=${chrs[$SLURM_ARRAY_TASK_ID]}

if [ ! -f ./chr${i}.createMethylDataset.out ]; then
  echo "Creating methyl dataset for chr${i}"
  python /data/Phillippy/projects/pickettbd_data/parent-of-origin/pofo_tools/src/pofo_tools/createMethylDataset.py \
  --chr chr${i} \
  --input-fofn /data/Phillippy/projects/pooch_snakemake/pooch_data/05.dss/dss.fofn \
  --output-file chr${i}.createMethylDataset.out \
  --log chr${i}.createMethylDataset.log \
  --sex-mapping-file /data/Phillippy/projects/pooch_snakemake/model_validation/model_juhyun/test_data/sex.123.tsv
else
  echo "Methyl dataset for chr${i} already exists, skipping."
fi

if [ ! -f ./chr${i}.missingrate_5.createMethylDataset.out ]; then
  echo "Filtering methyl dataset for chr${i} with missing rate 5"
  python filter_missing_rate.py chr${i} 5
else
  echo "Filtered methyl dataset for chr${i} with missing rate 5 already exists, skipping."
fi

if [ ! -f ./chr${i}.missingrate_5.createDiffDatasetFromMethylDataset.out ]; then
  echo "Creating differential dataset for chr${i} with missing rate 5"
  python /data/Phillippy/projects/pickettbd_data/parent-of-origin/pofo_tools/src/pofo_tools/createDiffDatasetFromMethylDataset.py \
  --input-methyl-dataset ./chr${i}.missingrate_5.createMethylDataset.out \
  -L Maternal -R Paternal -c chr${i} \
  -l ./chr${i}.createDiffDatasetFromMethylDataset.log \
  --output-file ./chr${i}.missingrate_5.createDiffDatasetFromMethylDataset.out
else
  echo "Differential dataset for chr${i} with missing rate 5 already exists, skipping."
fi


# if [ ! -f elasticnet_selected_cpgs_chr${i}.csv ]; then
#   echo "Running elastic net cross-validation for chr${i} with missing rate 5"
#   python elasticnet.cv.py chr${i}
# else
#   echo "Elastic net cross-validation for chr${i} with missing rate 5 already exists, skipping."
# fi

#### #SBATCH --array=2,16,22
