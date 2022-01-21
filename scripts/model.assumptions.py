import nibabel as nb
import numpy as np
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor
import os
from scipy import stats

os.chdir() #change this to your local file path

img = nb.load('eres.mgh')

df = pd.read_csv('behavioural_results.csv')

data = img.get_fdata()

voxel_data = np.array(data[10000].flatten(),dtype=float) 

voxel_df = pd.DataFrame(voxel_data)

model = pd.concat([df[['G-Number','age_adjusted_group','Age']],voxel_df],axis=1).rename(columns={0:'voxel'})

groups = pd.get_dummies(model['age_adjusted_group'])

model = pd.concat([model,groups], axis=1)

X_age = model[['voxel','AAN','HC','WR','Age']]  

vif_data_age = pd.DataFrame()

vif_data_age["feature"] = X_age.columns
  
vif_data_age["VIF"] = [variance_inflation_factor(X_age.values, i) for i in range(len(X_age.columns))]

print(vif_data_age)

print(X_age.corr())

X_minage = model[['voxel','AAN','HC','WR']]  

vif_data_minage = pd.DataFrame()

vif_data_minage["feature"] = X_minage.columns
  
vif_data_minage["VIF"] = [variance_inflation_factor(X_minage.values, i) for i in range(len(X_minage.columns))]

print(vif_data_minage)

print(X_minage.corr())

print(stats.normaltest(voxel_data))