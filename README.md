# Examining the Relationship Between Autistic Spectrum Disorder Characteristics and Structural Brain Differences seen in Anorexia Nervosa.

This is a repo for all the code used in this paper.


Usage for all scripts in scripts folder is: 

~~~
python3 behavioural_differences.py
~~~

Except for plotting scripts which need to be run in a ipython notebook.


## Files needed


behavioural_results csv
	
	- A CSV with all the behavioural results.
	- Ideally columns named:  BMI_at_scan, Initial_EDE_Q_Total, Age, IQ, ADOS_com_soc, ADOS_Creativity, ADOS_sterotyped_and_repetititve and age_adjusted_group
	- If columns not named like this then code needs to be changed.

cortical_measures.csv

	- A CSV with mean_curv, mean_lgi, mean_area, mean_thickness, TotalGrayVol, Total_white_matter, parameter estimates.
	- Created from dataframe_creation.py

Output from Freesurfer:

asegstats2table command:

	- aseg_volume.dat

aparcstats2table command:

	- lh_area.dat
	- rh_area.dat
	- lh_curv.dat
	- rh_curv.dat
	- lh_lgi.dat
	- rh_lgi.dat
	- lh_thickness.dat
	- rh_thickness.dat

## Enviornment

The scripts in this repo get the eniormental variables from a .env (not included in the repo). Create a .env file with the following set up:

```
data=<path to data folder>
```

to allow for scripts to be ran correctly.