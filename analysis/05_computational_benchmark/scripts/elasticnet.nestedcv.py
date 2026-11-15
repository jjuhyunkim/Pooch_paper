import sys
import pickle
import joblib
import warnings
import pandas as pd
import numpy as np
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, confusion_matrix
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedKFold, GridSearchCV, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    roc_auc_score,
    classification_report,
    confusion_matrix
)
from sklearn.exceptions import ConvergenceWarning, UndefinedMetricWarning
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=UndefinedMetricWarning)
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, balanced_accuracy_score, confusion_matrix, classification_report, brier_score_loss
from sklearn.linear_model import LogisticRegression


thread = 50
max_iter = 100000 # maximum solver iterations

print("All arguments:", sys.argv)
if len(sys.argv) > 1:
    chr = sys.argv[1]
    print("First argument:", chr)

prop = sys.argv[2] if len(sys.argv) > 2 else None

par_filer=False

if chr == "chrX" :
    testset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/createDiffDatasetFromMethylDataset/testX.samples", sep ='\t', header = None)[0].to_list()
    trainingset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/createDiffDatasetFromMethylDataset/trainingX.samples", sep ='\t', header = None)[0].to_list()
else :
    testset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/testSet.samplelist", sep ='\t', header = None)[0].to_list()
    trainingset = pd.read_csv("/path/to/pooch/model_validation/model_juhyun/test_data/trainSet.samplelist", sep ='\t', header = None)[0].to_list()

diff_methyl = pd.read_csv(f"{chr}.missingrate_5.createDiffDatasetFromMethylDataset.out", sep='\t', header = 0)
print(f"left : {diff_methyl['LeftOpHapLabel'].unique()}")

if chr == "chrX" and par_filer == True :
    print("set to filter in only par regions")
    par_bounds = globals().get('x_par')
    if par_bounds is None:
        raise ValueError("x_par is not defined; set par_filer=False or load x_par before using chrX PAR filtering")
    diff_methyl_T = diff_methyl.T
    diff_methyl_T_names = diff_methyl_T.iloc[0:6, :]
    diff_methyl_T_loc = diff_methyl_T.iloc[6:, :]
    diff_methyl_T_loc['loc'] = diff_methyl_T_loc.index.astype(int)
    par1 = diff_methyl_T_loc.loc[diff_methyl_T_loc['loc'].between(par_bounds.iloc[0,1],par_bounds.iloc[0,2])]
    par2 = diff_methyl_T_loc.loc[diff_methyl_T_loc['loc'].between(par_bounds.iloc[1,1],par_bounds.iloc[1,2])]
    par = pd.concat([par1, par2])
    par = par.drop('loc', axis = 1)
    diff_methyl_T_filt = pd.concat([diff_methyl_T_names, par])
    diff_methyl = diff_methyl_T_filt.T
else :
    print(f"par_filter : {str(par_filer)}")

# SWAPSAMPLES
sample_list=['HG02451', 'HG02717', 'HG03579']
diff_methyl_sub = diff_methyl.loc[~diff_methyl['Sample'].isin(sample_list),:]
diff_methyl_sub_a = diff_methyl.loc[diff_methyl["Sample"].isin(sample_list),diff_methyl.columns[:3]].copy()
diff_methyl_sub_b = diff_methyl.loc[diff_methyl["Sample"].isin(sample_list),diff_methyl.columns[3:]].copy()
diff_methyl_sub_b= diff_methyl_sub_b *-1
diff_methyl_sub_a_b = pd.concat([diff_methyl_sub_a, diff_methyl_sub_b], axis=1)
diff_methyl = pd.concat([diff_methyl_sub_a_b, diff_methyl_sub], axis = 0)

# Training set
training_data = diff_methyl.loc[diff_methyl['Sample'].isin(trainingset)]
training_data = training_data.iloc[:int(prop), :].copy() if prop is not None else training_data

training_data = training_data.reset_index().drop("index", axis = 1)
training_data_a = training_data.iloc[:,:3].copy()
training_data_a['LeftOpHapLabel'] = "Paternal"
training_data_a['RightOpHapLabel'] = "Maternal"
training_data_b = training_data.iloc[:,3:].copy()
training_data_b= training_data_b *-1
training_data_flip = pd.concat([training_data_a, training_data_b], axis=1).copy()
training_data = pd.concat([training_data, training_data_flip], axis = 0).copy()
training_data = training_data.reset_index().drop("index", axis = 1)
group=training_data[["Sample"]].to_numpy().flatten()
train_X = training_data.iloc[:, 3:]
yd_train = pd.get_dummies(training_data, columns=['RightOpHapLabel'], prefix='', prefix_sep='')
yd_train = yd_train['Paternal'].to_numpy()
print(f"group : {group}")
print(f"training_data : {training_data.shape}")
print(f"train : {yd_train}")


# TEST
test_data = diff_methyl.loc[diff_methyl['Sample'].isin(testset)]
test_data = test_data.reset_index().drop("index", axis = 1)
test_data_a = test_data.iloc[:,:3].copy()
test_data_b = test_data.iloc[:,3:].copy()
test_data_b= test_data_b *-1
test_data_a['LeftOpHapLabel'] = "Paternal"
test_data_a['RightOpHapLabel'] = "Maternal"
test_data_flip = pd.concat([test_data_a, test_data_b], axis=1).copy()
test_data = pd.concat([test_data, test_data_flip], axis = 0).copy()
test_data = test_data.reset_index().drop("index", axis = 1)
yd_test = pd.get_dummies(test_data, columns=['RightOpHapLabel'], prefix='', prefix_sep='')
yd_test = yd_test['Paternal'].to_numpy()
test_X = test_data.iloc[:, 3:]
print(f"test_data : {test_data.shape}")
print(f"test : {yd_test}")

if training_data.shape[0] == 0:
    raise ValueError("No training rows matched trainingset after filtering; check sample lists and input file")

if test_data.shape[0] == 0:
    raise ValueError("No test rows matched testset after filtering; check sample lists and input file")

if train_X.shape[0] == 0:
    raise ValueError("train_X is empty after preprocessing; cannot fit LogisticRegression")
if len(training_data['Sample'].to_numpy()) != len(train_X):
    raise ValueError(f"Group length mismatch: train_groups={len(training_data['Sample'].to_numpy())} train_X={len(train_X)}")

from sklearn.feature_selection import SelectFdr, f_classif
## Fitting
pipe = Pipeline([
    ("select", SelectFdr(score_func=f_classif, alpha=0.05)),
    ("model", LogisticRegression(
        penalty="elasticnet",
        solver="saga",
        max_iter=max_iter,
        random_state=42,
        fit_intercept=False
    ))
])

# ----------------------------
# Hyperparameter grid
# ----------------------------

param_grid = {
    "model__C" : [0.01, 0.1, 1, 10, 100, 1000],   # inverse of regularization strength
    "model__l1_ratio": [0, 0.25, 0.5, 0.75, 1]  # 0=ridge, 1=lasso
}

inner_cv = GroupKFold(n_splits=5)
outer_cv = GroupKFold(n_splits=5)

search = GridSearchCV(
    estimator=pipe,
    param_grid=param_grid,
    scoring=["neg_brier_score", "accuracy"],
    cv=inner_cv,
    n_jobs=-1,
    refit="neg_brier_score"
)

# ----------------------------
# Nested CV estimate
# ----------------------------
# This gives a less biased estimate of performance
cv_pred = np.zeros(len(train_X))

for train_idx, valid_idx in outer_cv.split(
    train_X,
    yd_train,
    groups=group
):

    X_train = train_X.iloc[train_idx]
    X_valid = train_X.iloc[valid_idx]

    y_train = yd_train[train_idx]

    group_train = group[train_idx]


    search.fit(
        X_train,
        y_train,
        groups=group_train
    )


    cv_pred[valid_idx] = (
        search.best_estimator_
        .predict_proba(X_valid)[:,1]
    )


cv_auc = roc_auc_score(
    yd_train,
    cv_pred
)

print(
    f"Nested CV AUC: {cv_auc:.4f}"
)

# cv_auc = roc_auc_score(yd_train, cv_pred)
print(f"Nested CV AUC (train only): {cv_auc:.4f}")




# ----------------------------
# Fit final model on full training data
# ----------------------------
search.fit(train_X, yd_train, groups=group)
best_model = search.best_estimator_

mask = best_model.named_steps["select"].get_support()
selected_features = train_X.columns[mask]

select = best_model['select']
model = best_model['model']
best_model['model'].feature_names_in_ = selected_features
k = len(selected_features)
l1_ratio = best_model['model'].l1_ratio
c = best_model['model'].C

# save the best model
import joblib
import pickle

# Save best model
with open(f"elasticnet_bestmodel_{chr}.c{c}_l{l1_ratio}.{k}.{prop}.pkl", "wb") as f:
    pickle.dump(best_model, f)

with open(f"elasticnet_model_{chr}.c{c}_l{l1_ratio}.{k}.{prop}.pkl", "wb") as f:
    pickle.dump(model, f)

# Save full grid search (for cv_results_)
joblib.dump(search, f"elasticnet_model_gridsearch_{chr}.{prop}.pkl")

# ----------------------------
# Evaluate on test set
# ----------------------------
test_prob = best_model.predict_proba(test_X)[:, 1]
test_pred = (test_prob >= 0.5).astype(int)

print("\nTest performance")
print(f"Accuracy: {accuracy_score(yd_test, test_pred):.4f}")
print(f"Balanced accuracy: {balanced_accuracy_score(yd_test, test_pred):.4f}")
# print(f"AUC: {roc_auc_score(yd_test, test_prob):.4f}")
print("\nConfusion matrix:")
print(confusion_matrix(yd_test, test_pred))
print("\nClassification report:")
print(classification_report(yd_test, test_pred))
## add brier score
brier_score = brier_score_loss(yd_test, test_prob)
print(f"Brier score: {brier_score:.4f}")

# ----------------------------
# Extract selected CpGs
# ----------------------------
coefs = best_model.named_steps["model"].coef_.ravel()

# Use fitted feature names when available to avoid shape/name mismatches.
# feature_names = getattr(best_model, "feature_names_in_", train_X.columns)
# feature_names = getattr(best_model.named_steps["model"], "feature_names_in_", train_X.columns)
feature_names = selected_features
if len(feature_names) != len(coefs):
    print(
        f"Warning: feature name count ({len(feature_names)}) != coef count ({len(coefs)}). "
        "Falling back to train_X columns."
    )
    feature_names = train_X.columns

selected = pd.Series(coefs, index=feature_names)
selected = selected[selected != 0].sort_values(key=np.abs, ascending=False)

print(f"\nOriginal CpGs: {len(coefs)}")
print(f"\nSelected CpGs: {len(selected)}")
print(selected.head(20))

# make one table
selected_df = pd.DataFrame([
    {   "original_CpG": train_X.shape[1],
        "train_num" : prop,
        "total_CpG": len(coefs),
        "selected_CpG": len(selected),
        "cv_auc": f"{cv_auc:.4f}",
        "best_c": c,
        "best_l1_ratio": l1_ratio,
        "Test_Accuracy": accuracy_score(yd_test, test_pred),
        "Test_Balanced_Accuracy": balanced_accuracy_score(yd_test, test_pred),
        "Test_Brier_Score": brier_score_loss(yd_test, test_prob),
    }
])
selected_df.to_csv(f"elasticnet_selected_cpgs_{chr}.{prop}.csv", index=False)