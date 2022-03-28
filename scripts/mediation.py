'''
Old script to test a mediation effect. However not used in the actual paper.
'''

import statsmodels.api as sm
from statsmodels.stats.mediation import Mediation
import pandas as pd
import os
from decouple import config

data = config('data')

os.chdir(data) 

cortical = pd.read_csv('cortical_measures.csv')
b_df = pd.read_csv('behavioural_results.csv')
plotting_df = pd.concat([cortical, b_df[['BMI_at_scan','ADOS_sterotyped_and_repetititve']]], axis=1)

group = plotting_df.groupby('age_adjusted_group')
hc = group.get_group('HC')
aan = group.get_group('AAN')
wr = group.get_group('WR')

aan_model = aan[['mean_lgi','Age' , 'BMI_at_scan']].dropna() 

outcome_model = sm.OLS.from_formula("mean_lgi ~ Age + BMI_at_scan",
                                     aan_model)

mediator_model = sm.OLS.from_formula("BMI_at_scan ~ Age + mean_lgi", 
aan_model)

print(outcome_model.fit().summary())
print(mediator_model.fit().summary())

med = Mediation(outcome_model, mediator_model, exposure="Age", mediator='BMI_at_scan').fit()

print(med.summary())
