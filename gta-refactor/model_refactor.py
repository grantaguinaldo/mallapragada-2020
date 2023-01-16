from pyutilib.misc.timing import tic, toc
import pyomo.environ as en
import pandas as pd
import numpy as np
import collections
pd.options.mode.chained_assignment = None 
import os

def build_model(pvavailarray=None,
                pricearray=None,
                dfPVData=None,
                dfStorData=None,
                dfElyData=None,
                dfH2StData=None,
                CCF=CCF_val,
                productionCommitmentLB=None,
                minimumProductionShutdownLength=None, 
                genLP=None, 
                output_file_name=None):
    
    m = en.ConcreteModel()
    
    tval = range(1, len(pvavailarray)+1, 1)
    m.t = en.Set(initialize = tval)

    esval = dfStorData.index.values.tolist()
    m.bes = en.Set(initialize = esval)
    
    ##########################################
    # Model Parameters
    ##########################################
    
    m.pCCF = en.Param(within = en.NonNegativeReals, 
                      initialize = CCF, 
                      mutable = False)
    
    m.pH2DesignFlowRate = en.Param(within = en.NonNegativeReals, 
                                   initialize = 4166.667, 
                                   mutable = False)
    
    m.pCFPlantLB = en.Param(within = en.NonNegativeReals,
                            initialize = 0.9,
                            mutable = False)
    
    m.pH2LHV = en.Param(within = en.NonNegativeReals, 
                        initialize = 120.1, 
                        mutable = False)
    
    m.pProductionSlackCost = en.Param(within = en.NonNegativeReals,
                                      initialize = 0, 
                                      mutable=True)
    
    m.pCapCostPV = en.Param(within = en.NonNegativeReals,
                            initialize = dfPVData.CapCost_dkW.values[0], 
                            mutable = False)
    
    m.pFOMCostPV = en.Param(within = en.NonNegativeReals,
                            initialize = dfPVData.FOM_pct_CAPEX.values[0]*dfPVData.CapCost_dkW.values[0]*1000, 
                            mutable = False)
    
    m.pVOMCostPV = en.Param(within = en.NonNegativeReals,
                            initialize = dfPVData.VOM_dMWh.values[0],
                            mutable = False)
    
    m.pPVCapFactor = en.Param(m.t,
                              within = en.NonNegativeReals, 
                              initialize = {tval[j]:pvavailarray[j] for j in range(len(tval))}, 
                              mutable = True)
    
    m.pGridElecPrice = en.Param(m.t, 
                                initialize={tval[j]:pricearray[j] for j in range(len(tval))}, 
                                mutable=True)
    
    m.pInvEff = en.Param(within = en.NonNegativeReals, 
                         initialize = dfPVData.InvEff.values[0], 
                         mutable = False)
    
    m.pCapCostInv = en.Param(within = en.NonNegativeReals, 
                             initialize = dfPVData.InvCapCost_dkW.values[0], 
                             mutable = False)
    
    m.pStEffChg = en.Param(m.bes, 
                           within = en.NonNegativeReals, 
                           initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                                 'St_eff_chg'] for j in range(len(esval))}, 
                           mutable = False)
    
    m.pStEffDischg = en.Param(m.bes, 
                              within = en.NonNegativeReals, 
                              initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                                    'St_eff_dischg'] for j in range(len(esval))}, 
                              mutable = False)
    
    m.pCapCostPowSt = en.Param(m.bes, 
                               within = en.NonNegativeReals,
                               initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                                     'Power_capex_dpkW'] for j in range(len(esval))}, 
                               mutable = False)
    
    m.pCapCostEnergySt = en.Param(m.bes, 
                                  within = en.NonNegativeReals,
                                  initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                                        'Energy_capex_dpkWh'] for j in range(len(esval))}, 
                                  mutable = False)
    
    m.pFOMCostSt = en.Param(m.bes, 
                            within = en.NonNegativeReals, 
                            initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                                  'FOM_dMWyr'] for j in range(len(esval))}, 
                            mutable = False)
    
    m.pVOMCostSt = en.Param(m.bes, 
                            within = en.NonNegativeReals,
                            initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                                  'VOM_dMWh'] for j in range(len(esval))}, 
                            mutable = False)
    
    m.pDur_UB = en.Param(m.bes,
                         within = en.NonNegativeReals,
                         initialize = {esval[j]:dfStorData.loc[esval[j], 
                                                               'MaxDur_hrs'] for j in range(len(esval))}, 
                         mutable = False)
    
    m.pElySpecPower = en.Param(within = en.NonNegativeReals, 
                               initialize = dfElyData.ElySpecPower_kWhkg.values[0], 
                               mutable = False)
    
    m.pCapCostEly = en.Param(within = en.NonNegativeReals, 
                             initialize = dfElyData.CapCost_dkW.values[0], 
                             mutable = False)
    
    m.pFOMCostEly = en.Param(within = en.NonNegativeReals,
                             initialize = dfElyData.FOM_pct_CAPEX.values[0]*dfElyData.CapCost_dkW.values[0]*1000, 
                             mutable = False)
    
    m.pVOMCostEly = en.Param(within = en.NonNegativeReals,
                             initialize = dfElyData.VOMCost_dMWh.values[0], 
                             mutable = False)

    m.pFeedH2OCostEly = en.Param(within = en.NonNegativeReals,
                                 initialize = dfElyData.Water_cost_d_per_kg_H2.values[0], 
                                 mutable = False)

    m.pCapCostH2Comp = en.Param(within = en.NonNegativeReals, 
                                initialize = dfH2StData.CapCostComp_dkW.values[0], 
                                mutable = False)

    m.pFOMCostH2Comp = en.Param(within = en.NonNegativeReals,
                                initialize = dfH2StData.CompFOM_pct_CAPEX.values[0]*dfH2StData.CapCostComp_dkW.values[0]*1000, 
                                mutable = False)

    m.pCompSpecPower = en.Param(within = en.NonNegativeReals, 
                                initialize = dfH2StData.CompSpecPower_kWhpkg.values[0], 
                                mutable = False)

    m.pCapCostH2st = en.Param(within = en.NonNegativeReals,
                              initialize = dfH2StData.CapCostst_dkg.values[0], 
                              mutable = False)    
    
    m.pFOMCostH2st = en.Param(within = en.NonNegativeReals,
                             initialize = dfH2StData.StFOM_pct_CAPEX.values[0]*dfH2StData.CapCostst_dkg.values[0], 
                              mutable = False)
    
    m.pH2kgpertank = en.Param(within = en.NonNegativeReals, 
                              initialize = dfH2StData.mass_stored_kg.values[0], 
                              mutable = False)
    
    ##########################################
    # Decision Variables                     
    ##########################################
    
    m.vPVInstalledMW = en.Var(within = en.NonNegativeReals)
    
    m.vInvInstalledMW = en.Var(within = en.NonNegativeReals)
    
    m.vPVtoDCPow = en.Var(m.t, 
                          within = en.NonNegativeReals)
    
    m.vPVOutput = en.Var(m.t, 
                         within = en.NonNegativeReals)
    
    m.vStDischarge = en.Var(m.t, 
                            m.bes, 
                            within = en.NonNegativeReals)

    m.vStCharge = en.Var(m.t, 
                         m.bes, 
                         within = en.NonNegativeReals)

    m.vStSoC = en.Var(m.t, 
                      m.bes, 
                      within = en.NonNegativeReals)

    m.vStInstalledMW = en.Var(m.bes, 
                              within = en.NonNegativeReals)

    m.vStInstalledMWh = en.Var(m.bes, 
                               within = en.NonNegativeReals)
    
    m.vElyInstalledMW = en.Var(within = en.NonNegativeReals)

    m.vDCPowtoEly = en.Var(m.t,
                           within = en.NonNegativeReals)
    
    m.vCompInstalledMW = en.Var(within = en.NonNegativeReals)
    
    m.vH2StInstalledNumber = en.Var(within = en.NonNegativeReals)

    m.vH2flowProd = en.Var(m.t, 
                           within = en.NonNegativeReals)
    
    m.vH2flowStorage = en.Var(m.t,
                              within = en.NonNegativeReals)
        
    m.vDCPowtoInv = en.Var(m.t,
                           within = en.NonNegativeReals) 
    
    m.vACPowtoComp = en.Var(m.t, 
                            within = en.NonNegativeReals)
    
    m.vACPowtoGrid = en.Var(m.t, 
                            within = en.NonNegativeReals)
    
    m.vH2StStatekg = en.Var(m.t, 
                            within = en.NonNegativeReals)
    
    m.vH2StflowProd = en.Var(m.t, 
                             within = en.NonNegativeReals)

    m.vTotalH2Output = en.Var(m.t, 
                              within = en.NonNegativeReals)

    # Absolute difference outside of production bands -kg/hr
    # Set this to zero in case grid interactions are ignored
    m.vH2PlantOutputSlack = en.Var(m.t, 
                                   within = en.NonNegativeReals,
                                   initialize={tval[j]:0 for j in range(len(tval))})
    
    ##########################################
    # PV Constraints                    
    ##########################################
    
    # Equation S3
    # Units: MW
    # LHS: Total Amount to Power Block (MW) and 
    #      Total Amount Sent to Storage (MW)
    # RHS: Total output from the PV array (MW)
    # (Legacy Comment): PV energy balance - sum of storage charge and DC power sent to Power block less than PV output (DC)
    def PVEnergyBalance(m, t):
        return sum(m.vStCharge[t, st] for st in m.bes)+m.vPVtoDCPow[t]  == m.vPVOutput[t]
    m.cPVEnergyBal = en.Constraint(m.t, 
                                   rule = PVEnergyBalance)

    # TODO: Look into `pPVCapFactor`
    # Equation S4
    # Units: MW
    # LHS: Total output from the PV array (MW)
    # RHS: Nameplate capacity of PV array (MW) * Availability (unitless)
    def PVCapLim(m,t):
        return m.vPVOutput[t] <= (m.vPVInstalledMW*m.pPVCapFactor[t])
    m.cPVCapLim = en.Constraint(m.t, 
                                rule = PVCapLim)
    
    ##########################################
    # Storage Constraints                    
    ##########################################
    
    # TODO: Double check units to see if this is an energy or power balance. 
    # Equation S5 and S6
    # Units: MW-Hr
    # Balance storage capacity at each time step (MWh)
    def StSoCBal(m, t, st):  # wrapping storage capacity to ensure first and last period are matching
        # Equation S6
        if t in [1]:       # First hour of the year
                           #  first hour specific constraints- wrapping constraints across the year

            # LHS1: Energy Charge State of battery at last time step (MW-Hr) + 
            # LHS2: Battery Charge Efficiency (unitless) + Battery Power Charge (MW) (`vStCharge_{t=1}` (MW)) -
            # LHS3:  Battery Power discharge (MW) (`vStDischarge_{t=1}` (MW)) / Battery Discharge Efficiency (unitless)
            # RHS: Energy Charge State of battery at time `t` (MW-Hr)
            return m.vStSoC[len(pvavailarray),st]+(m.pStEffChg[st]*m.vStCharge[t,st])-(m.vStDischarge[t,st]/m.pStEffDischg[st]) == m.vStSoC[t,st]
        # Equation S5
        else:
            # LHS1: Energy Charge State of battery at last time step (MW-Hr) + 
            # LHS2: Battery Charge Efficiency (unitless) + Battery Power Charge (MW) (`vStCharge_{t>1}` (MW)) -
            # LHS3:  Battery Power discharge (MW) (`vStDischarge_{t>1}` (MW)) / Battery Discharge Efficiency (unitless)
            # RHS: Energy charge state of battery at time `t>1` (MW-Hr)
            return m.vStSoC[t-1, st]+(m.pStEffChg[st]*m.vStCharge[t,st])-(m.vStDischarge[t,st]/m.pStEffDischg[st]) == m.vStSoC[t,st]
    m.cStSoCBal = en.Constraint(m.t, 
                                m.bes, 
                                rule=StSoCBal)

    # Equation S7
    # Units: MW
    # LHS: Nameplate power capacity of the Battery System (MW)
    # RHS: Power Charge efficiency (unitless) * Power Charge State at time `t` (MW)
    # (Legacy Comment): Upper limit on power charge rate into the battery
    def StLimChargeUB(m, t, st):
        return (m.pStEffChg[st]*m.vStCharge[t,st]) <= m.vStInstalledMW[st] 
    m.cStChargeUB = en.Constraint(m.t,
                                  m.bes,
                                  rule = StLimChargeUB)

    # Equation S8
    # Units: MW
    # LHS: Power Charge State at time `t` (MW)
    # RHS: Nameplate power capacity of the Battery System (MW) * Power Discharge efficiency (unitless)
    # (Legacy Comment): Upper limit on discharge rate into the battery
    def StLimDischargeUB(m, t, st):
        return m.vStDischarge[t,st] <= m.pStEffDischg[st]*m.vStInstalledMW[st]
    m.cStLimDischargeUB = en.Constraint(m.t, 
                                        m.bes, 
                                        rule = StLimDischargeUB)

    # Equation S9
    # Units: MW-Hr
    # LHS: Energy Charge State at time `t` (MW-Hr)
    # RHS: Nameplate energy capacity of the Battery System (MW-Hr)  
    # (Legacy Comment): Storage capacity cannot exceed purchased energy capacity
    def StCap_rule(m, t, st):
        return m.vStSoC[t,st] <= m.vStInstalledMWh[st]
    m.cSt_Cap = en.Constraint(m.t, 
                              m.bes, 
                              rule = StCap_rule)

    # Equation S10
    # RHS: Nameplate power capacity of the Battery System (MW) * Max Discharge Duration (hrs)
    # LHS: Nameplate energy capacity of the Battery System (MW-Hr)  
    # (Legacy Comment): Storage duration upper bounded by specified parameter value hours
    def StDur_UB(m, st):
        return m.vStInstalledMWh[st] <= m.pDur_UB[st]*m.vStInstalledMW[st]
    m.cStDur_UB = en.Constraint(m.bes, 
                                rule = StDur_UB)
    
    ##########################################
    # Inverter Constraints                    
    ##########################################
    
    # Equation S11
    # Units: MW
    # RHS: Power Block to Inverter (MW) + Power Block to Electrolyzer (MW)
    # LHS: Power from Battery Discharge (MW) + Power from PV Array (MW)
    # (Legacy Comment): DC Power splitter block balance 
    def DCPowerBlockBal(m, t):
        return sum(m.vStDischarge[t,st] for st in m.bes)+m.vPVtoDCPow[t] == m.vDCPowtoEly[t]+m.vDCPowtoInv[t]
    m.cDCPowerBal = en.Constraint(m.t,
                                  rule =DCPowerBlockBal)

    ##########################################
    # Electrolyzer Constraints                    
    ##########################################

    # Equation S12
    # Units: kg/hr
    # RHS: H2 flow to Production (kg/hr) + H2 flow to Compressor (Storage) (kg/hr)
    # LHS: (1,000 KW / MW) * Power to electrolyzer (MW) / Hydrogen Production Factor ( (KW-Hr) / (kg) ) 
    # (Legacy Comment): electrolyzer energy balance 
    # (Legacy Comment): Multiplying DC power by 1000 to convert to kW - 1 hour resolution unis are kg/hr left and right hand side
    def ElyEnergyBal(m, t):
        return (m.vDCPowtoEly[t]*1000)/m.pElySpecPower == m.vH2flowProd[t]+m.vH2flowStorage[t]
    m.cElyEnergyBal = en.Constraint(m.t, 
                                    rule = ElyEnergyBal)

    # Equation S13
    # Units: MW
    # RHS: Nameplate power capacity of electrolyzer (MW)
    # LHS: Total power supplied to electrolyzer from power block (MW)
    # (Legacy Comment): capacity balance on electrolyzer
    def ElyCapLim(m, t):
        return m.vDCPowtoEly[t] <= m.vElyInstalledMW
    m.cElyCapLim = en.Constraint(m.t, 
                                 rule = ElyCapLim)
    
    # Storage, compressor and inverter balance   
    # Equation S14
    # Units: MW (Inbound), AC (Outbound)
    # RHS: Power from Inverter to Compressor (AC) + Power from Inverter to Grid (AC) 
    # LHS: Power from Power Block (DC) to Inverter (MW) * Inverter Efficiency (unitless)
    # (Legacy Comment): Inverter energy balance in MW
    def InvEnergyBal(m, t):
        return (m.vDCPowtoInv[t]*m.pInvEff) == (m.vACPowtoComp[t]+m.vACPowtoGrid[t])
    m.cInvEnergyBal = en.Constraint(m.t, 
                                    rule = InvEnergyBal)

    # Equation S15
    # Units: MW
    # RHS: Installed Inverter Nameplate Capacity (MW)
    # LHS: Power from Power Block (DC) to Inverter (MW)
    # (Legacy Comment): Inverter capacity rating in MW
    def InvCapLim(m, t):
        return m.vDCPowtoInv[t] <= m.vInvInstalledMW
    m.cInvCapLim = en.Constraint(m.t,
                                 rule = InvCapLim)

    # Equation S16 (Constraint on the Power Requirements for the Compressor)
    # Units: MW
    # RHS: Installed Compressor Nameplate Capacity (MW)
    # LHS: Power from Inverter to Compressor (AC) (MW)
    # (Legacy Comment): Compressor capacity rating inMW
    def CompCapLim(m, t):
        return m.vACPowtoComp[t] <= m.vCompInstalledMW
    m.cCompCapLim = en.Constraint(m.t, 
                                  rule=CompCapLim)

    # Equation S17 (Constraint on the Production of Compressed H2 from Compressor)
    # Units: MW
    # RHS: H2 Production to Storage from Comprssor (kg/hr) * Compressor Production Factor (KW-Hr / kg) * (1 MW / 1000 KW)
    # LHS: Power from Inverter to Compressor (AC) (MW)
    # (Legacy Comment): Compressor energy requirement -CompSpecPower unit MJ/kg of H2, units of equation MW
    def CompPowReq(m, t):
        return m.vACPowtoComp[t] == m.vH2flowStorage[t]*(m.pCompSpecPower/1000)
    m.cCompPowReq = en.Constraint(m.t, 
                                  rule = CompPowReq)
    
    #TODO: Double check units for these constraints.
    # Equation S18 and S19
    def H2StSoCBalance(m, t): 
        if t in [1]:
            # Units: kg/hr
            # RHS: Total amount of H2 being stored at time `t_{8760}` (or last time step) (kg/hr)
            # LHS1: Total amount of H2 being stored at time `t=1` (kg/hr) +
            # LHS2: Total amount of H2 sent to storage at time `t=1` (kg/hr) - 
            # LHS3: Total amount of H2 sent to production at time `t=1` (kg/hr)     
            return m.vH2StStatekg[len(pvavailarray)]+m.vH2flowStorage[t]-m.vH2StflowProd[t] == m.vH2StStatekg[t]
        else:
            
            # Units: kg/hr 
            # RHS: Total amount of H2 being stored at time `t`, such that `t > 1` (kg/hr)
            # LHS1: Total amount of H2 being stored at time `t-1`, such that `t > 1` (kg/hr) +
            # LHS2: Total amount of H2 sent to storage at time `t`, such that `t > 1` (kg/hr) - 
            # LHS3: Total amount of H2 sent to production at time `t`, such that `t > 1` (kg/hr) 
            return m.vH2StStatekg[t-1]+m.vH2flowStorage[t]-m.vH2StflowProd[t] == m.vH2StStatekg[t]
    m.cH2StBal = en.Constraint(m.t, 
                               rule = H2StSoCBalance)

    # Equation S20
    # Units: kg
    # RHS: Number of Installed H2 Tanks (num. units) * Unitized Capacity to store H2 (kg/unit)
    # LHS: Total amount of H2 being stored at time `t` (kg)
    # (Legacy Comment): Capacity rating for H2 storage cannot exceed total storage
    def H2StCapLim(m, t):
        return m.vH2StStatekg[t] <= m.vH2StInstalledNumber*m.pH2kgpertank
    m.cH2StCapLim = en.Constraint(m.t, 
                                  rule =H2StCapLim)
    
    # Total plant H2 production in kg/hr
    def TotalPlantH2bal(m, t):
        return m.vH2flowProd[t]+m.vH2StflowProd[t] == m.vTotalH2Output[t]
    m.cTotalPlantH2bal =en.Constraint(m.t, 
                                      rule = TotalPlantH2bal)

    # Binary variable for each time period
    m.vProductionCommitment = en.Var(m.t, 
                                     within = en.Binary)

    # Equation S21
    # Units: hrs
    # RHS: Total Commitment (hr)
    # LHS: Sum of Indicator Variables for Commitment (hr)
    # (Legacy Comment): Lower bound on the # of hours shutdown
    def ProductionCommitmentLB(m):
        return sum(m.vProductionCommitment[t] for t in m.t) >= productionCommitmentLB #m.pProductionCommitmentLB
    m.cProductionCommitmentLB = en.Constraint(rule = ProductionCommitmentLB)

    # TODO: Look at where does `offset` get introduced into the model. 
    # Unsure of Equations S22, S23, and S24.
    def ProductionCommitmentContiguity(m, t, offset):
        def ind(i):
            return ((i-1) % len(pvavailarray)) + 1
        k = t - offset
        return m.vProductionCommitment[ind(t)]-m.vProductionCommitment[ind(t-1)] <= 1-m.vProductionCommitment[ind(k)]

    m.sContiguityOffset = en.RangeSet(1, 
                                      minimumProductionShutdownLength)

    m.cProductionCommitmentContiguity = en.Constraint(m.t, 
                                                      m.sContiguityOffset,
                                                      rule = ProductionCommitmentContiguity)

    # Equation 2 (Lower Bound)
    # Units: kg/hr
    # RHS: Total H2 Production, from Storage and Electrolyzer (kg/hr)
    # LHS: 
    # (Legacy Comment): H2 plant output generation requirements 
    def H2PlantOutputLB_rule(m, t):
        if en.value(m.pProductionSlackCost) == 0: #- in the absence of slacks
            return (m.pH2DesignFlowRate*m.pCFPlantLB*m.vProductionCommitment[t]) <= m.vTotalH2Output[t] 

        else:# H2 plant output generation requirements with slacks
            #return - m.vH2PlantOutputSlack[t] <= m.vTotalH2Output[t] - m.pH2DesignFlowRate * m.pCFPlantLB * m.vProductionCommitment[t]
            return (m.pH2DesignFlowRate*m.pCFPlantLB*m.vProductionCommitment[t])-m.vH2PlantOutputSlack[t] <= m.vTotalH2Output[t] 
    m.cH2PlantOutputLB = en.Constraint(m.t, 
                                       rule = H2PlantOutputLB_rule)

    # Equation 2 (Upper Bound)
    # Units: kg/hr
    # RHS: ( Design Flow Rate (kg/hr) * Binary Variable (1/0) ) 
    # LHS: Total H2 Production, from Storage and Electrolyzer (kg/hr)
    def H2PlantOutputUB_rule(m, t):
        if en.value(m.pProductionSlackCost) == 0: #- in the absence of slacks
            return m.vTotalH2Output[t] <= (m.pH2DesignFlowRate*m.vProductionCommitment[t])

        else: # H2 plant output generation requirements with slacks
            # Units: kg/hr
            # RHS: ( Design Flow Rate (kg/hr) * Binary Variable (1/0) ) + Plant Output Slack (kg/hr)
            # LHS: Total H2 Production, from Storage and Electrolyzer (kg/hr)
            return m.vH2PlantOutputSlack[t] >= m.vTotalH2Output[t]-(m.pH2DesignFlowRate*m.vProductionCommitment[t])
            #return m.vTotalH2Output[t] <= (m.pH2DesignFlowRate * m.vProductionCommitment[t]) + m.vH2PlantOutputSlack[t]
    m.cH2PlantOutputUB = en.Constraint(m.t, 
                                       rule = H2PlantOutputUB_rule)

    ##########################################
    # Objective Function Terms                   
    ##########################################
    
    def PVFixCost_rule(m):
        val = (m.vPVInstalledMW*m.pCCF*m.pCapCostPV*1000)+(m.vPVInstalledMW*m.pFOMCostPV)
        return val
    m.ePVFixCost = en.Expression(rule =PVFixCost_rule)
    
    def StorFixCost_rule(m):
        val = sum(1000*m.pCapCostPowSt[st]*m.vStInstalledMW[st]*m.pCCF+1000*m.pCapCostEnergySt[st]*m.vStInstalledMWh[st]*m.pCCF+m.pFOMCostSt[st]*m.vStInstalledMW[st] for st in m.bes)
        return val
    m.eStFixCost = en.Expression(rule = StorFixCost_rule)
    
    def ElyFixCost_rule(m):
        val = (m.pCapCostEly*m.vElyInstalledMW*m.pCCF*1000)+(m.pFOMCostEly*m.vElyInstalledMW)
        return val
    m.eElyFixCost = en.Expression(rule = ElyFixCost_rule)
    
    def H2StFixCost_rule(m):
        val= (m.vCompInstalledMW*m.pCapCostH2Comp*m.pCCF*1000) + \
             (m.vCompInstalledMW*m.pFOMCostH2Comp) + \
             (m.vInvInstalledMW*m.pCapCostInv*m.pCCF*1000) + \
             (m.vH2StInstalledNumber*m.pH2kgpertank*m.pCapCostH2st*m.pCCF) +\
             (m.vH2StInstalledNumber*m.pH2kgpertank*m.pFOMCostH2st)
        return val
    m.eH2StFixCost = en.Expression(rule = H2StFixCost_rule)
    
    def SysVOMCost_rule(m):
        val = sum(sum(m.vStCharge[t,st] + m.vStDischarge[t,st] for t in m.t)*m.pVOMCostSt[st] for st in m.bes) + sum((m.vH2flowProd[t]*m.pFeedH2OCostEly)+(m.vH2flowStorage[t]*m.pFeedH2OCostEly) for t in m.t)
        return val
    m.eSysVOMCost = en.Expression(rule = SysVOMCost_rule)
    
    def RevenueGridSale_rule(m):
        val = sum(m.vACPowtoGrid[t]*m.pGridElecPrice[t] for t in m.t)
        return val
    m.eGridRevenue = en.Expression(rule = RevenueGridSale_rule)

    def ProductionSlack_rule(m):
        val = m.pProductionSlackCost*sum(m.vH2PlantOutputSlack[t] for t in m.t)
        return val
    m.eProductionSlack = en.Expression(rule = ProductionSlack_rule)
    
    def SysTotalCost(m):
        val = m.ePVFixCost+m.eStFixCost+m.eElyFixCost+m.eH2StFixCost+m.eSysVOMCost-m.eGridRevenue
        return val
    m.eSysTotalCost = en.Expression(rule = SysTotalCost)

    def Obj_fn(m):
        val = m.eSysTotalCost+m.eProductionSlack
        return val
    m.oObjective = en.Objective(rule = Obj_fn, sense = en.minimize)
    
    if genLP:
        m.write(output_file_name)
    else:
        pass
    
    return m

