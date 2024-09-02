import numpy as np
from types import SimpleNamespace

# ##############################################################
# # Temporary calculator for rolling back to original sty-BA experiments
# ##############################################################
# def calculateFlowrates(method, parameters, numExp=1, strategy="standard"):
#     expIndex = np.linspace(0, numExp - 1, numExp)
# ###### For generalisation, take an argument for the type of polymerisation used and run the calculation for each type ######## e.g. if type=conventionalEmulsion: do these calcs...
    
#     vars = SimpleNamespace(**parameters) # Convert all dict items from variables argument to variables (as lists so need indexing when used)
#     # calculate a list of emulsion concentrations
#     w_emulsion = []
#     for r, f in zip(vars.Seed_fraction, vars.w_f):
#         w_emulsion_temp = (f*vars.w_s[0]*(1 - r))/(vars.w_s[0] - f*r)
#         w_emulsion.append(w_emulsion_temp)

#     # subsequently calculate the emulsion densities (weighted average of monomer and aqueous phase densities)
#     densityEmulsion = [] 
#     for e, m in zip(w_emulsion, vars.densityMonomer):
#         densityEmulsion_temp = m*e + (1-e)*vars.densityAq[0]
#         densityEmulsion.append(densityEmulsion_temp)
   
#     v_seed = []
#     v_emulsion = []
#     v_total = []
#     v_monomers = []
#     v_monomerA = []
#     v_monomerB = []
#     v_Aq = []
#     v_Aq1 = []
#     v_Aq2 = []

#     for i in expIndex.astype(int):
#         if strategy=="optimisation": # If optimisation experiment, we want to keep all variables in the provided 'parameters' dictionary and only find the new flowrates for the last set of conditions
#             i=-1
#         emulsionFeedsIndex = np.linspace(1, vars.numFeeds[i], vars.numFeeds[i])
#         summation = 0
#         for f in emulsionFeedsIndex:
#             summation_temp = 1/(1 + ((vars.densitySeed[i]/densityEmulsion[i])*((vars.w_s[i]/(vars.w_f[i]*vars.Seed_fraction[i])) - 1))*(1-((vars.numFeeds[i]-f)/vars.numFeeds[i])))
#             summation = summation + summation_temp
#         A = (vars.densitySeed[i]/densityEmulsion[i])*((vars.w_s[i]/(vars.w_f[i]*vars.Seed_fraction[i])) - 1)
#         B = (vars.numberCSTRs[i] - vars.numFeeds[i])/(1 + A)
#         v_seed_temp = round(((vars.volumeCSTR[i]/vars.tau[i])*(summation + B)), 4)
#         v_seed.append(v_seed_temp)
#         v_emulsion_temp = ((vars.densitySeed[i]*v_seed_temp)/densityEmulsion[i])*((vars.w_s[i]/(vars.w_f[i]*vars.Seed_fraction[i])) - 1)
#         v_emulsion.append(v_emulsion_temp)
#         v_total_temp = v_seed_temp + v_emulsion_temp
#         v_total.append(v_total_temp)
        
#         v_monomers_temp = (v_emulsion[i]*densityEmulsion[i]*(vars.Surfactant_concentration[i] - w_emulsion[i]))/(vars.densityMonomer[i]*(vars.Surfactant_concentration[i] - 1))
#         v_monomers.append(v_monomers_temp)

#         v_monomerA_temp = round((vars.densityMonomer[i]*v_monomers[i]*vars.monomerAFraction[i]/vars.densityMonomerA[i]), 4)
#         v_monomerA.append(v_monomerA_temp)
        
#         if vars.densityMonomerB[i] == 0: # Conditional to avoid divide by zero error if monomer ratio not varied
#             v_monomerB_temp = 0
#             v_monomerB.append(v_monomerB_temp)
#         else:
#             v_monomerB_temp = round((vars.densityMonomer[i]*v_monomers[i]*(1-vars.monomerAFraction[i])/vars.densityMonomerB[i]), 4)
#             v_monomerB.append(v_monomerB_temp)
        
#         v_Aq_temp = round((((densityEmulsion[i]*v_emulsion[i])-(vars.densityMonomer[i]*v_monomers[i]))/(vars.densityAq[i])), 4)
#         v_Aq.append(v_Aq_temp)

#         v_Aq1_temp = round(((vars.densityAq[i]*v_Aq[i]*(vars.Surfactant_concentration[i] - vars.w_Aq2[i]))/(vars.densityAq[i]*(vars.w_Aq1[i] - vars.w_Aq2[i]))), 4)
#         v_Aq1.append(v_Aq1_temp)

#         v_Aq2_temp = round(((vars.densityAq[i]*v_Aq[i]*(vars.Surfactant_concentration[i] - vars.w_Aq1[i]))/(vars.densityAq[i]*(vars.w_Aq2[i] - vars.w_Aq1[i]))), 4)
#         v_Aq2.append(v_Aq2_temp)

#     flowRatesDict = {
#         'v_seed' : v_seed,
#         'v_Aq1' : v_Aq1,
#         'v_Aq2' : v_Aq2,
#         'v_monomerA' : v_monomerA,
#         'v_monomerB' : v_monomerB
#     }
#     return flowRatesDict
# ##############################################################
# # Temporary calculator for rolling back to original sty-BA experiments
# ##############################################################

def calculateFlowrates(method, parameters, numExp=1, strategy=None, seeded=False):

    expIndex = np.linspace(0, numExp - 1, numExp)
    vars = SimpleNamespace(**parameters)

    ### If method passed is for conventional emulsion polymerisation, calculate the following ###
    if method == 'conventional':

        # calculate a list of emulsion concentrations
        w_emulsion = []
        for r, f in zip(vars.Seed_fraction, vars.w_f):
            w_emulsion_temp = (f*vars.w_s[0]*(1 - r))/(vars.w_s[0] - f*r)
            w_emulsion.append(w_emulsion_temp)
        
        # subsequently calculate the emulsion densities (weighted average of monomer and aqueous phase densities)
        densityEmulsion = [] 
        for e, m in zip(w_emulsion, vars.densityMonomer):
            densityEmulsion_temp = m*e + (1-e)*vars.densityAq[0]
            densityEmulsion.append(densityEmulsion_temp)

        v_seed = []
        v_emulsion = []
        v_total = []
        v_initiator = []
        v_surfactant = []
        v_monomer = []
        v_monomerA = []
        v_monomerB = []
        v_water = []

        for i in expIndex.astype(int):
            if strategy=="optimisation": # If optimisation experiment, we want to keep all variables in the provided 'parameters' dictionary and only find the new flowrates for the last set of conditions
                i=-1
            emulsionFeedsIndex = np.linspace(1, vars.numFeeds[i], vars.numFeeds[i])
            summation = 0
            for f in emulsionFeedsIndex:
                summation_temp = 1/(1 + ((vars.densitySeed[i]/densityEmulsion[i])*((vars.w_s[i]/(vars.w_f[i]*vars.Seed_fraction[i])) - 1))*(1-((vars.numFeeds[i]-f)/vars.numFeeds[i])))
                summation += summation_temp
            A = (vars.densitySeed[i]/densityEmulsion[i])*((vars.w_s[i]/(vars.w_f[i]*vars.Seed_fraction[i])) - 1)
            B = (vars.numberCSTRs[i] - vars.numFeeds[i])/(1 + A)
            v_seed_temp = round(((vars.volumeCSTR[i]/vars.tau[i])*(summation + B)), 4)
            v_seed.append(v_seed_temp)
            v_emulsion_temp = ((vars.densitySeed[i]*v_seed_temp)/densityEmulsion[i])*((vars.w_s[i]/(vars.w_f[i]*vars.Seed_fraction[i])) - 1)
            v_emulsion.append(v_emulsion_temp)
            v_total_temp = round((v_seed_temp + v_emulsion_temp), 4)
            v_total.append(v_total_temp)

            Sterm = (vars.Surfactant_ratio[i]*(1-(vars.w_surfactantStock[i]/vars.w_makeupWater[i])))/(100*vars.surfactantStockConc[i])
            Iterm = (vars.Initiator_ratio[i]*(1-(vars.w_initiatorStock[i]/vars.w_makeupWater[i])))/(100*vars.initiatorStockConc[i])
            x = 1 - (1/vars.w_makeupWater[i]) + Sterm + Iterm
            v_monomer_temp = round((v_emulsion[i]*densityEmulsion[i]*(1-(w_emulsion[i]/vars.w_makeupWater[i])))/(vars.densityMonomer[i]*x), 4)
            v_monomer.append(v_monomer_temp)

            try:
                v_monomerA_temp = round(((v_monomer[i]*vars.densityMonomer[i])/(vars.densityMonomerA[i]*(1 + ((1 - vars.monomerAFraction[i])/vars.monomerAFraction[i])))), 4)
                v_monomerA.append(v_monomerA_temp)
                v_monomerB_temp = round(((v_monomer[i]*vars.densityMonomer[i])/(vars.densityMonomerB[i]))*(1 - (1/(1 + ((1 - vars.monomerAFraction[i])/vars.monomerAFraction[i])))), 4)
                v_monomerB.append(v_monomerB_temp)
            except(AttributeError): 
                v_monomerB.append(0)
                pass

            v_initiator_temp = round((vars.Initiator_ratio[i]*v_monomer_temp*vars.densityMonomer[i])/(100*vars.densityInitiator[i]*vars.initiatorStockConc[i]), 4)
            v_initiator.append(v_initiator_temp)

            v_surfactant_temp = round((vars.Surfactant_ratio[i]*v_monomer_temp*vars.densityMonomer[i])/(100*vars.densitySurfactant[i]*vars.surfactantStockConc[i]), 4)
            v_surfactant.append(v_surfactant_temp)

            v_water_temp = round((v_emulsion[i] - v_monomer[i] - v_surfactant[i] - v_initiator[i]), 4)
            v_water.append(v_water_temp)
        
        flowratesDict = {
            'v_seed' : v_seed,
            'v_monomer' : v_monomer,
            'v_monomerA' : v_monomerA,
            'v_monomerB' : v_monomerB,
            'v_initiator' : v_initiator,
            'v_surfactant' : v_surfactant,
            'v_water' : v_water
        }

    ### If method passed is for RAFT emulsion polymerisation, calculate the following ###
    if method == 'RAFT':

        v_total = []
        v_initiator = []
        v_macroCTA = []
        v_monomer = []
        v_water = []

        if seeded:
            for i in expIndex.astype(int):
                if strategy=="optimisation":
                    i=-1
                v_total_temp = ((vars.volumeCSTR[i]*vars.numberCSTRs[i]) + vars.volumeReactor2[i])/vars.tau[i]
                v_total.append(v_total_temp)

                X = (vars.macroCTAPolymerConc[i] + vars.macroCTAMonomerConc[i])*vars.molarMassMacroCTA[i] + ((vars.macroCTAPolymerConc[i]*vars.molarMassInitiator[i])/vars.ctaToInitiator[i]) + ((vars.DP[i]*vars.macroCTAPolymerConc[i]*vars.molarMassMonomer[i]) - (vars.macroCTAMonomerConc[i]*vars.molarMassMacroCTA[i]))

                v_initiator_temp = round(((v_total[i]*vars.densityProduct[i]*vars.w_f[i]*vars.macroCTAPolymerConc[i]*vars.molarMassInitiator[i])/(vars.ctaToInitiator[i]*vars.initiatorStockConc[i]*vars.densityInitiator[i]*X)), 4)
                v_initiator.append(v_initiator_temp)

                v_macroCTA_temp = round((v_total[i]*vars.densityProduct[i]*vars.w_f[i]*vars.molarMassMacroCTA[i])/(vars.densityMacroCTA[i]*X), 4)
                v_macroCTA.append(v_macroCTA_temp)

                v_monomer_temp = round((v_total[i]*vars.densityProduct[i]*vars.w_f[i]*((vars.DP[i]*vars.macroCTAPolymerConc[i]*vars.molarMassMonomer[i]) - (vars.macroCTAMonomerConc[i]*vars.molarMassMacroCTA[i])))/(vars.densityMonomer[i]*X), 4)
                v_monomer.append(v_monomer_temp)

                v_water_temp = round(((v_total[i]*vars.densityProduct[i])/vars.densityAq[i])*(1 - ((vars.w_f[i]/X)*(vars.molarMassMacroCTA[i] + ((vars.macroCTAPolymerConc[i]*vars.molarMassInitiator[i])/(vars.ctaToInitiator[i]*vars.initiatorStockConc[i])) + ((vars.DP[i]*vars.macroCTAPolymerConc[i]*vars.molarMassMonomer[i]) - vars.macroCTAMonomerConc[i]*vars.molarMassMacroCTA[i])))), 4)
                v_water.append(v_water_temp)
        
        else:

            for i in expIndex.astype(int):
                if strategy=="optimisation": # If optimisation experiment, we want to keep all variables in the provided 'parameters' dictionary and only find the new flowrates for the last set of conditions
                    i=-1
                v_total_temp = ((vars.volumeCSTR[i]*vars.numberCSTRs[i]) + vars.volumeReactor2[i])/vars.tau[i]
                v_total.append(v_total_temp)

                v_initiator_temp = round((vars.molarMassInitiator[i]*v_total[i]*vars.densityProduct[i]*vars.w_f[i])/((vars.initiatorStockConc[i]*vars.densityInitiator[i]*vars.ctaToInitiator[i])*(vars.molarMassMacroCTA[i] + (vars.molarMassInitiator[i]/vars.ctaToInitiator[i]) + (vars.DP[i]*vars.molarMassMonomer[i]))), 4)
                v_initiator.append(v_initiator_temp)

                v_macroCTA_temp = round((vars.molarMassMacroCTA[i]*v_total[i]*vars.densityProduct[i]*vars.w_f[i])/((vars.macroCTAStockConc[i]*vars.densityMacroCTA[i])*(vars.molarMassMacroCTA[i] + (vars.molarMassInitiator[i]/vars.ctaToInitiator[i]) + (vars.DP[i]*vars.molarMassMonomer[i]))), 4)
                v_macroCTA.append(v_macroCTA_temp)

                v_monomer_temp = round((vars.molarMassMonomer[i]*v_total[i]*vars.densityProduct[i]*vars.w_f[i]*vars.DP[i])/((vars.densityMonomer[i])*(vars.molarMassMacroCTA[i] + (vars.molarMassInitiator[i]/vars.ctaToInitiator[i]) + (vars.DP[i]*vars.molarMassMonomer[i]))), 4)
                v_monomer.append(v_monomer_temp)

                v_water_temp = round(((v_total[i]*vars.densityProduct[i])/(vars.densityAq[i]))*(1-(vars.w_f[i]/(vars.molarMassMacroCTA[i] + (vars.molarMassInitiator[i]/vars.ctaToInitiator[i]) + (vars.DP[i]*vars.molarMassMonomer[i])))*((vars.molarMassMacroCTA[i]/vars.macroCTAStockConc[i]) + (vars.molarMassInitiator[i]/(vars.initiatorStockConc[i]*vars.ctaToInitiator[i])) + (vars.DP[i]*vars.molarMassMonomer[i]))), 4)
                v_water.append(v_water_temp)

        flowratesDict = {
            'v_initiator' : v_initiator,
            'v_macroCTA' : v_macroCTA,
            'v_monomer' : v_monomer,
            'v_water' : v_water
        }

        Vol_macroCTA = 0
        Vol_monomer = 0
        Vol_initiator = 0
        Vol_water = 0

        for i in expIndex.astype(int):

            time_SS = vars.tau[i]*3
            time_prep = vars.volumeCSTR[i]*3/(10*(v_total[i] - v_initiator[i]))
            time_sample = vars.volumeSample[i]/v_total[i]
            time_total = time_SS + time_sample

            Vol_macroCTA_temp = (v_macroCTA[i]*time_total) + (10*v_macroCTA[i]*time_prep)
            Vol_macroCTA += Vol_macroCTA_temp

            Vol_monomer_temp = (v_monomer[i]*time_total) + (10*v_monomer[i]*time_prep)
            Vol_monomer += Vol_monomer_temp

            Vol_initiator_temp = v_initiator[i]*time_total
            Vol_initiator += Vol_initiator_temp

            Vol_water_temp = (v_water[i]*time_total) + (10*v_water[i]*time_prep)
            Vol_water += Vol_water_temp

        print(f'Volume macroCTA = {Vol_macroCTA}')
        print(f'Volume monomer = {Vol_monomer}')
        print(f'Volume initiator = {Vol_initiator}')
        print(f'Volume water = {Vol_water}')
    return flowratesDict