from pyutilib.misc.timing import tic, toc
import pyomo.environ as en
import pandas as pd
import numpy as np
import collections
pd.options.mode.chained_assignment = None 
import os
from model_refactor import build_model

cost_scenario_folder = os.path.join('CostScenarios', 
                                    '2020_AG')

StorageData = pd.read_excel(os.path.join(cost_scenario_folder, 
                                         'StorageData.xlsx'),
                            index_col=[0])

Discount_rate = 0.054  
Lifetime = 20.0  
CCF_val = 1/float((Discount_rate+1)/float(Discount_rate)*(1-1/(1+Discount_rate)**Lifetime))  

PVData = pd.read_excel(os.path.join(cost_scenario_folder, 
                                    'PVData.xlsx'),
                       'Data',
                       index_col=[0]) 

ElyData = pd.read_excel(os.path.join(cost_scenario_folder, 
                                     'ElyData.xlsx'),
                        'Data',
                        index_col=[0]) 

H2StData = pd.read_excel(os.path.join(cost_scenario_folder, 
                                      'H2StData.xlsx'),
                         'Data',
                         index_col=[0]) 

cf_file = '5796_23.65_68.75_tmy.csv' 

PVAvail_tmy = pd.read_csv(cf_file,
                          index_col=0, 
                          parse_dates=True,
                          header=None, 
                          squeeze=True)

productionCommitmentLB = int(np.floor(len(PVAvail_tmy)*0.95)) 

minimumProductionShutdownLength = 12 

P_Electricity = 120.0 

LMPData = pd.Series(P_Electricity, 
                    index=range(len(PVAvail_tmy))) 


model = build_model(pvavailarray = PVAvail_tmy,
                    pricearray = LMPData,
                    dfPVData = PVData,
                    dfStorData = StorageData,
                    dfElyData = ElyData,
                    dfH2StData = H2StData,
                    CCF = CCF_val,
                    productionCommitmentLB = productionCommitmentLB,
                    minimumProductionShutdownLength = minimumProductionShutdownLength, 
                    genLP=True, 
                    output_file_name='model_file.lp')
