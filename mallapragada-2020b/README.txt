This folder contains the pyomo implementation of the design and operations optimization model used in the paper and described in S.1.

The folder includes the following files:
Electrolyzer_v11.py - python file with functions to build optimization model object and store relevant outputs
example.ipynb - Jupyter notebook file that reads input data, solves an instance of the model and stores outputs.

INPUT DATA:
5796_23.65_68.75_tmy.csv - hourly solar resource availability at location specified by Latitude 23.65 degrees and Longitude 68.75 degrees. Resource availability based on typical meteorological data. Data obtained from National Solar Radiation Database (nsrdb.nrel.gov)

CostScenarios - folder that store possible set of cost scenarios to evaluate the model over. Each subfolder corresponds to a particular cost scenario and consists of the following files:
ElyData.xlsx -  Electrolyzer cost and performance parameters
H2StData.xlsx - H2 storage, compressor cost and performance parameters
PVData.xlsx - PV system cost and performance parameters
StorageData.xlsx - Battery storage cost and performance parameters

