'''
Script to create the cortical_measures.csv needed for other scripts.
Input data needs to be in set format (see README.md)
'''

import pandas as pd
import os
from decouple import config

data = config('data')

os.chdir(data) 

#Read in behavioural df
behaviour_df = pd.read_csv('behavioural_results.csv')

#This reads in the aseg_volume.dat table from freesurfer and calculates total global white matter

volumes = pd.read_csv('aseg_volume.dat',sep='\t')
cerbellum_white_matter = volumes['Right-Cerebellum-White-Matter'] + volumes['Left-Cerebellum-White-Matter']
volumes['Total_white_matter'] = cerbellum_white_matter + volumes['CerebralWhiteMatterVol']

#area
left_area = pd.read_csv('lh_area.dat',sep='\t')
right_area = pd.read_csv('rh_area.dat',sep='\t')
area_com = pd.concat([left_area, right_area],axis=1)
mean_area = area_com[['lh_WhiteSurfArea_area','rh_WhiteSurfArea_area']].mean(axis=1)
area = pd.concat([area_com,mean_area],axis=1).rename(columns={0:'mean_area'})

#curv
left_curv = pd.read_csv('lh_curv.dat',sep='\t')
right_curv = pd.read_csv('rh_curv.dat',sep='\t')
curv_com = pd.concat([left_curv,right_curv],axis=1).drop(columns=[ 'BrainSegVolNotVent', 'eTIV','rh.aparc.meancurv'])
mean_curv = curv_com[curv_com.columns[1:68]].mean(axis=1)
curv = pd.concat([curv_com,mean_curv],axis=1).rename(columns={0:'mean_curv'})


#lgi
left_lgi = pd.read_csv('lh_lgi.dat',sep='\t')
right_lgi = pd.read_csv('rh_lgi.dat',sep='\t')
lgi_com = pd.concat([left_lgi,right_lgi], axis=1).drop(columns=[ 'BrainSegVolNotVent', 'eTIV','rh.aparc.pial_lgi.thickness'])
mean_lgi = lgi_com[lgi_com.columns[1:68]].mean(axis=1)
lgi = pd.concat([lgi_com,mean_lgi],axis=1).rename(columns={0:'mean_lgi'})

#thickness
left_thickness = pd.read_csv('lh_thickness.dat')
right_thickness = pd.read_csv('rh_thickness.dat', sep='\t')
thickness = pd.concat([left_thickness,right_thickness], axis=1)
mean_thickness = thickness[['rh_MeanThickness_thickness','lh_MeanThickness_thickness']].mean(axis=1)
thickness = pd.concat([thickness,mean_thickness], axis=1).rename(columns={0:'mean_thickness'})

#Any cluster
cluster = pd.read_csv('/home/wmda/Documents/Documents_from_old_comp/BEACON/Write_up/Cortical/results/rh_lgi_cluster_parameters.csv', sep='\s+')

#Concats the dataframes into one useful dataframe
cortical_measures = pd.concat([behaviour_df[['G-Number','age_adjusted_group' ,'Age']],curv['mean_curv'], lgi['mean_lgi'], 
                             thickness['mean_thickness'], area['mean_area'], volumes['TotalGrayVol'], volumes['Total_white_matter'],
                             cluster],axis=1)
                             
cortical_measures['Group'] = cortical_measures['age_adjusted_group'].apply(lambda x: 1 if x =='AAN'  else (2 if x=='WR'  else 0))

cortical_measures.to_csv('cortical_measures.csv')