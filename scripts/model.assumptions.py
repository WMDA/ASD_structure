import nibabel as nb
import numpy as np
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
import os
from decouple import config

data = config('data')
os.chdir(data) 

img = nb.load('eres.mgh')
df = pd.read_csv('behavioural_results.csv')
data = img.get_fdata()
voxel_data = np.array(data[10000].flatten(),dtype=float) 
voxel_df = pd.DataFrame(voxel_data)

for covariate in ['BMI_baseline', 'Age']:
    model = pd.concat([df[['G-Number','age_adjusted_group', covariate]],voxel_df],axis=1).rename(columns={0:'voxel'})
    groups = pd.get_dummies(model['age_adjusted_group'])
    model = pd.concat([model,groups], axis=1)
    X_val = model[['voxel','AAN','HC','WR', covariate]].dropna()  
    vif_data = pd.DataFrame()
    vif_data["feature"] = X_val.columns
    vif_data["VIF"] = [variance_inflation_factor(X_val.values, i) for i in range(len(X_val.columns))]
    print(vif_data, '\n')
    print(X_val.corr(), '\n')

X_no_cov = model[['voxel','AAN','HC','WR']]  
vif_data_mincov = pd.DataFrame()
vif_data_mincov["feature"] = X_no_cov.columns
vif_data_mincov["VIF"] = [variance_inflation_factor(X_no_cov.values, i) for i in range(len(X_no_cov.columns))]
print(vif_data_mincov,'\n')
print(X_no_cov.corr())