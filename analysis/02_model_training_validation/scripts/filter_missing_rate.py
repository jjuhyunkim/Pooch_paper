import pandas as pd
import seaborn as sns
import glob as glob
import sys

# get arguments from outside
chr=sys.argv[1]
missing_rate = int(sys.argv[2])  # percentage of missing values allowed in each column


print(f"chr : {chr}")
print(f"missing_rate : {missing_rate}")

# calculate variance for each column
files = glob.glob(f"{chr}.createMethylDataset.out")
if len(files) == 0:
    print(f"No files found for {chr}.createMethylDataset.out")
    sys.exit(1)


# files = glob.glob(f"/data/Phillippy/projects/pooch_snakemake/pooch_data/06.pooch_input/{chr}.N05_filt.createMethylDataset.out")
# files = glob.glob(f"/path/to/pooch/pooch_data/07.models_using_80_filt_regions/chr3.createMethylDataset.out")
print(files)
for file in files:    
    df = pd.read_csv(file, sep = '\t', header = 0 )
    # chr = file.split("/")[-1].split(".")[0]
    if chr == "chr19":
        df = df.loc[df['Sample'] != 'HG02602']

    # calculate number of NA values in each column
    # filter columns with less than 10 NAs
    #na_counts = df.isna().sum()
    #filtered_df = df.loc[:, na_counts < NA_per * len(df)]
    # print(filtered_df.shape)

    # filtered_df.to_csv(file.replace(".createMethylDataset.out", ".createMethylDataset.filtered_NA0.out"), sep='\t', header = True, index = False)

if chr == "chrX" :
    testset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/createDiffDatasetFromMethylDataset/testX.samples", sep ='\t', header = None)[0].to_list()
    trainingset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/createDiffDatasetFromMethylDataset/trainingX.samples", sep ='\t', header = None)[0].to_list()
else :
    testset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/testSet.samplelist", sep ='\t', header = None)[0].to_list()
    trainingset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/trainSet.samplelist", sep ='\t', header = None)[0].to_list()
# testset[0:6]
# test_data = pd.read_csv(f"/vf/users/Phillippy/projects/pickettbd_data/parent-of-origin/classifier/logistic-regression/diff_NA0_xform_V4/final-117.{chr}.methylDs.filt-N90.diffDs.tsv", sep='\t', header = 0)
# test_data = pd.read_csv(f"/vf/users/Phillippy/projects/pickettbd_data/parent-of-origin/classifier/logistic-regression/diff_NA0_xform_V4/final-117.{chr}.methylDs.filt-N90.diffDs.tsv", sep='\t', header = 0)

# Rebuild training subset each run to avoid cumulative filtering side-effects
df_train_curr = df.loc[df['Sample'].isin(trainingset)].copy()
print(f"before filtering : {df_train_curr.shape}")

# Keep metadata columns unchanged; filter only feature columns by NA rate
meta_cols = df_train_curr.columns[:3]
feature_cols = df_train_curr.columns[3:]

na_counts_column = df_train_curr[feature_cols].isna().sum(axis=0)
max_na = df_train_curr.shape[0] * missing_rate / 100

print(f"number of feature columns above threshold ({max_na:.2f} NAs): {(na_counts_column > max_na).sum()}")
keep_feature_cols = na_counts_column[na_counts_column <= max_na].index
df_train = pd.concat([df[meta_cols], df[keep_feature_cols]], axis=1)

print(f"after filtering : {df_train.shape}")
df_train.to_csv(f"{chr}.missingrate_{missing_rate}.createMethylDataset.out", sep='\t', header = True, index = False)