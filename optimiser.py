from summit.domain import *
from summit.strategies import SNOBFIT, SOBO, TSEMO, MultitoSingleObjective, LogSpaceObjectives, Chimera
from summit.benchmarks import SnarBenchmark
from summit.utils.dataset import DataSet
import pandas as pd
import time
import latinHypercubeSampling
# from emulator import func
import random
import glob
from pymoo.indicators.hv import HV

def func(X, a, b, c, d, e, f, g, h, i, m):
    x, y, z = X
    return a*x* + b*y + c*z +d*x**2 + e*y**2 + f*z**2 + g*x*y + h*x*z + i*y*z + m

def getNextExperiment(inputs, objectives, iteration, current_dataSet=None):
    '''
    Function takes input variables and objectives as dict of {"parameter name" : [minBound, maxBound]}
    
    For objective variables, the dict keys must be tuples indicating whether the objective is to be maximised
    or minimised, i.e. {("parameter name", "MIN") : [minBound, maxBound]}

    current_dataSet is the previous experimental data passed as a dict of {"parameter name" : [value1, value2, ... valuen]}
    
    iteration (usually used to count the number of iterations for a given optimisation), lets the objective be 'FLIPFLOPed'
    alternately to switch an objective, i.e. in the case where a maximisation and minimisation maps out an objective space
    
    '''
    if iteration%2==0:
        maximize_FLIPFLOP = 1
    else:
        maximize_FLIPFLOP = 0

    domain = Domain()
    for variable in inputs:
        bound = inputs[variable]
        domain += ContinuousVariable(name=str(variable), description='', bounds=bound)

    for objective in objectives:
        bound = objectives[objective]
        if objective[1] == "MIN":
            maximize=False
        elif objective[1] == "MAX":
            maximize=True
        elif objective[1] == "FLIPFLOP":
            maximize=maximize_FLIPFLOP
        domain += ContinuousVariable(name=str(objective[0]), description='', bounds=bound, is_objective=True, maximize=maximize)

    print("Getting new reaction conditions for experiment " + str(iteration) + "...")
    df = pd.DataFrame(data=current_dataSet)
    previous_data = DataSet.from_df(df)
    # transform = LogSpaceObjectives(domain)
    strategy = TSEMO(domain)
    next_experiments = strategy.suggest_experiments(1, previous_data)
    next_experiments_df = pd.DataFrame(next_experiments)
    parameter_names = list(inputs.keys())
    params_dict = {}
    for i, name in zip(range(len(next_experiments_df.columns) - 1), parameter_names):
        params_temp = round(next_experiments_df.iloc[:, i], 5).to_list()
        params_dict[name] = params_temp
    return params_dict

# ''' 
# 1. TODO: Generalise algorithm parameter name and bounds inputs

# '''
# a = -27.782
# b = -26.921
# c = -131.300
# d = -25187.368
# e = -7.612
# f = -1229.133
# g = -660.330
# h = 9018.266
# i = 32.775
# m = 131.737

# expString = r'standard-targeting-min-surfactant-and-seed'
# inputs = {'surfactantConc' : [0.01, 0.05],
#             'seedFrac': [0.02, 0.2086]}
# objectives = {('particleSizeFunc', "MIN") : [0.00001, 50],
#               ('surfactantConcMin', "MIN") : [0.01, 0.05],
#               ('seedFracMin', "MIN") : [0.02, 0.2086]}
# objectivesList = []
# for objective in objectives:
#     objectivesList.append(objective[0])
# variablesList = list(inputs.keys())
# variablesList.extend(objectivesList)

# sizeList = [80, 90, 100, 110, 120]
# n = 5 # number of optimisations to run (repeats)
# LHCs = 10

# for sizeTarget in sizeList:
#     path = str(r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Optimisation\Conventional\Simulated_experiments\CONSOLIDATED\STANDARD\SIZE-TARGETING\MINIMISE-SURFACTANT-AND-SEED' + r'\_' + str(sizeTarget) + 'nm')

#     for exp in range(n):  
#         monomerFrac = []  
#         parameterNames = list(inputs.keys())
#         parameterBounds = list(inputs.values())
#         trainingExps = latinHypercubeSampling.getHypercubeSamplingParams(LHCs, parameter_names=parameterNames, parameter_bounds=parameterBounds)
#         surfactantConc = trainingExps.get('surfactantConc')
#         seedFrac = trainingExps.get('seedFrac')
#         # monomerFrac = trainingExps.get('monomerFrac')
#         monomerFracTemp = 0.5
#         for i in surfactantConc:
#             monomerFrac.append(monomerFracTemp)
#         trainingExps["surfactantConcMin"] = surfactantConc
#         trainingExps["seedFracMin"] = seedFrac
 

#         particleSize = []
#         for x, y, z, in zip(surfactantConc,
#                             monomerFrac,
#                             seedFrac):
#             particleSizeTemp = round(func((x, y, z), a, b, c, d, e, f, g, h, i, m), 2)
#             particleSize.append(particleSizeTemp)
#         trainingExps["particleSize"] = particleSize

#         particleSizeFunc = []
#         for size in particleSize:
#             particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#             particleSizeFunc.append(particleSizeFuncTemp)
#         trainingExps["particleSizeFunc"] = particleSizeFunc

#         LHCRandomStateValue = trainingExps.pop("Random state value")
#         current_dataSet = trainingExps
#         current_dataSet = pd.DataFrame.from_dict(current_dataSet)

#         iterations = len(particleSize)
#         while iterations < 30:

#             new_conditions = getNextExperiment(inputs, objectives, current_dataSet, iteration=iterations)
#             iterations += 1

#             x = new_conditions.get('surfactantConc')[0]
#             # y = new_conditions.get('monomerFrac')[0]
#             y = 0.5
#             z = new_conditions.get('seedFrac')[0]
#             X = x, y, z
#             newSize = round(func(X, a, b, c, d, e, f, g, h, i, m), 2)
            
#             ### Update data set ###
#             surfactantConc.append(x)
#             monomerFrac.append(monomerFracTemp)
#             seedFrac.append(z)
#             particleSize.append(newSize)
#             particleSizeFunc = []
#             for size in particleSize:
#                 particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#                 particleSizeFunc.append(particleSizeFuncTemp)

#             values = [surfactantConc, seedFrac, particleSizeFunc, surfactantConc, seedFrac]
#             current_dataSet = {}
#             for variable, value in zip(variablesList, values):
#                 current_dataSet[variable] = value
#             current_dataSet = pd.DataFrame.from_dict(current_dataSet)
#             current_dataSet = DataSet(current_dataSet)

#             print(particleSize)
#         print(str(iterations) + ' reactions completed')
            
#         LHC_randomStateExp = []
#         for i in particleSize:
#             LHC_randomStateExp.append(LHCRandomStateValue)

#         optimisationSeries = {
#             "Surfactant concentration" : surfactantConc,
#             "Monomer fraction" : monomerFrac,
#             "Seed fraction" : seedFrac,
#             "Particle sizes" : particleSize,
#             "Particle size function" : particleSizeFunc,
#             "LHC random state" : LHC_randomStateExp
#         }
#         optimisationdf = pd.DataFrame(optimisationSeries)
#         optimisationdf.to_excel((path + r'\TSEMO-STANDARD' + str(exp) + '.xlsx'))

#     minParticleSizeFuncDict = {}
#     surfactantDict = {}
#     monomerDict = {}
#     seedDict = {}
#     particleSizeDict = {}
#     particleSizeFuncDict = {}
#     fileCounter = 0

#     for file in glob.glob(path + r"\*.xlsx"):
#         fileCounter += 1

#         data = pd.read_excel(file)
#         numExps = data.shape[0]

#         minParticleSizeFunc = []
#         minParticleSizeFuncTemp = 99999
#         for exp in range(numExps):
#             if data.iloc[exp, 5] < minParticleSizeFuncTemp:
#                 minParticleSizeFuncTemp = data.iloc[exp, 5]
#             minParticleSizeFunc.append(minParticleSizeFuncTemp)
        
#         minParticleSizeFuncDict[str(fileCounter)] = minParticleSizeFunc
        
#         surfactantConc = []
#         monomerFrac = []
#         seedFrac = []
#         particleSize = []
#         particleSizeFunc = []

#         for exp in range(numExps): # get experimental conditions from each 'reaction' of the experiments
#             surfactantConc.append(data.iloc[exp, 1])
#             monomerFrac.append(data.iloc[exp, 2])
#             seedFrac.append(data.iloc[exp, 3])
#             particleSize.append(data.iloc[exp, 4])
#             particleSizeFunc.append(data.iloc[exp, 5])


#         surfactantDict[str(fileCounter)] = surfactantConc
#         monomerDict[str(fileCounter)] = monomerFrac
#         seedDict[str(fileCounter)] = seedFrac
#         particleSizeDict[str(fileCounter)] = particleSize
#         particleSizeFuncDict[str(fileCounter)] = particleSizeFunc

#     surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
#     monomerSummary = pd.DataFrame.from_dict(monomerDict)
#     seedSummary = pd.DataFrame.from_dict(seedDict)
#     particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
#     particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
#     minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

#     surfactantTot = 0
#     surfactantAveList = []
#     monomerTot = 0
#     monomerAveList = []
#     seedTot = 0
#     seedAveList = []
#     particleSizeTot = 0
#     particleSizeAveList = []
#     particleSizeFuncTot = 0
#     particleSizeFuncAveList = []
#     minParticleSizeFuncTot = 0
#     minParticleSizeFuncAveList = []

#     for exp in range(numExps):
#         surfactantTot = sum(surfactantSummary.iloc[exp, :])
#         surfactantAveList.append(surfactantTot/len(surfactantSummary.columns))
#         monomerTot = sum(monomerSummary.iloc[exp, :])
#         monomerAveList.append(monomerTot/len(monomerSummary.columns))
#         seedTot = sum(seedSummary.iloc[exp, :])
#         seedAveList.append(seedTot/len(seedSummary.columns))
#         particleSizeTot = sum(particleSizeSummary.iloc[exp, :])
#         particleSizeAveList.append(particleSizeTot/len(particleSizeSummary.columns))
#         particleSizeFuncTot = sum(particleSizeFuncSummary.iloc[exp, :])
#         particleSizeFuncAveList.append(particleSizeFuncTot/len(particleSizeFuncSummary.columns))
#         minParticleSizeFuncTot = sum(minParticleSizeFuncSummary.iloc[exp, :])
#         minParticleSizeFuncAveList.append(minParticleSizeFuncTot/len(minParticleSizeFuncSummary.columns))
#     surfactantDict["Average"] = surfactantAveList
#     monomerDict["Average"] = monomerAveList
#     seedDict["Average"] = seedAveList
#     particleSizeDict["Average"] = particleSizeAveList
#     particleSizeFuncDict["Average"] = particleSizeFuncAveList
#     minParticleSizeFuncDict["Average"] = minParticleSizeFuncAveList

#     averagesDict = {
#         "Surfactant concentration" : surfactantAveList,
#         "Monomer fraction" : monomerAveList,
#         "Seed fraction" : seedAveList,
#         "Particle size" : particleSizeAveList,
#         "Particle size function" : particleSizeFuncAveList,
#         "Minimised particle size function" : minParticleSizeFuncAveList
#     }

#     surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
#     monomerSummary = pd.DataFrame.from_dict(monomerDict)
#     seedSummary = pd.DataFrame.from_dict(seedDict)
#     particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
#     particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
#     averagesSummary = pd.DataFrame.from_dict(averagesDict)
#     minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

#     surfactantSummary.to_excel(path + r'\surfactantSummary-' + expString + '.xlsx')
#     monomerSummary.to_excel(path + r'\monomerSummary-' + expString + '.xlsx')
#     seedSummary.to_excel(path + r'\seedSummary-' + expString + '.xlsx')
#     particleSizeSummary.to_excel(path + r'\particleSizeSummary-' + expString + '.xlsx')
#     particleSizeFuncSummary.to_excel(path + r'\particleSizeFuncSummary-' + expString + '.xlsx')
#     averagesSummary.to_excel(path + r'\averagesSummary-' + expString + '.xlsx')
#     minParticleSizeFuncSummary.to_excel(path + r'\minimisationSummary-' + expString + '.xlsx')

# #####################################################################################################################
# '''

# NEXT RUN

# '''
# #####################################################################################################################


# expString = r'standard-targeting-min-surfactant'
# inputs = {'surfactantConc' : [0.01, 0.05],
#             'seedFrac': [0.02, 0.2086]}
# objectives = {('particleSizeFunc', "MIN") : [0.00001, 50],
#               ('surfactantConcMin', "MIN") : [0.01, 0.05]
#             }
# objectivesList = []
# for objective in objectives:
#     objectivesList.append(objective[0])
# variablesList = list(inputs.keys())
# variablesList.extend(objectivesList)

# sizeList = [80, 90, 100, 110, 120]
# n = 5 # number of optimisations to run (repeats)
# LHCs = 10

# for sizeTarget in sizeList:
#     path = str(r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Optimisation\Conventional\Simulated_experiments\CONSOLIDATED\STANDARD\SIZE-TARGETING\MINIMISE-SURFACTANT' + r'\_' + str(sizeTarget) + 'nm')

#     for exp in range(n):  
#         monomerFrac = []  
#         parameterNames = list(inputs.keys())
#         parameterBounds = list(inputs.values())
#         trainingExps = latinHypercubeSampling.getHypercubeSamplingParams(LHCs, parameter_names=parameterNames, parameter_bounds=parameterBounds)
#         surfactantConc = trainingExps.get('surfactantConc')
#         seedFrac = trainingExps.get('seedFrac')
#         # monomerFrac = trainingExps.get('monomerFrac')
#         monomerFracTemp = 0.5
#         for i in surfactantConc:
#             monomerFrac.append(monomerFracTemp)
#         trainingExps["surfactantConcMin"] = surfactantConc
#         trainingExps["seedFracMin"] = seedFrac
 

#         particleSize = []
#         for x, y, z, in zip(surfactantConc,
#                             monomerFrac,
#                             seedFrac):
#             particleSizeTemp = round(func((x, y, z), a, b, c, d, e, f, g, h, i, m), 2)
#             particleSize.append(particleSizeTemp)
#         trainingExps["particleSize"] = particleSize

#         particleSizeFunc = []
#         for size in particleSize:
#             particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#             particleSizeFunc.append(particleSizeFuncTemp)
#         trainingExps["particleSizeFunc"] = particleSizeFunc

#         LHCRandomStateValue = trainingExps.pop("Random state value")
#         current_dataSet = trainingExps
#         current_dataSet = pd.DataFrame.from_dict(current_dataSet)

#         iterations = len(particleSize)
#         while iterations < 30:

#             new_conditions = getNextExperiment(inputs, objectives, current_dataSet, iteration=iterations)
#             iterations += 1

#             x = new_conditions.get('surfactantConc')[0]
#             # y = new_conditions.get('monomerFrac')[0]
#             y = 0.5
#             z = new_conditions.get('seedFrac')[0]
#             X = x, y, z
#             newSize = round(func(X, a, b, c, d, e, f, g, h, i, m), 2)
            
#             ### Update data set ###
#             surfactantConc.append(x)
#             monomerFrac.append(monomerFracTemp)
#             seedFrac.append(z)
#             particleSize.append(newSize)
#             particleSizeFunc = []
#             for size in particleSize:
#                 particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#                 particleSizeFunc.append(particleSizeFuncTemp)

#             values = [surfactantConc, seedFrac, particleSizeFunc, surfactantConc, seedFrac]
#             current_dataSet = {}
#             for variable, value in zip(variablesList, values):
#                 current_dataSet[variable] = value
#             current_dataSet = pd.DataFrame.from_dict(current_dataSet)
#             current_dataSet = DataSet(current_dataSet)

#             print(particleSize)
#         print(str(iterations) + ' reactions completed')
            
#         LHC_randomStateExp = []
#         for i in particleSize:
#             LHC_randomStateExp.append(LHCRandomStateValue)

#         optimisationSeries = {
#             "Surfactant concentration" : surfactantConc,
#             "Monomer fraction" : monomerFrac,
#             "Seed fraction" : seedFrac,
#             "Particle sizes" : particleSize,
#             "Particle size function" : particleSizeFunc,
#             "LHC random state" : LHC_randomStateExp
#         }
#         optimisationdf = pd.DataFrame(optimisationSeries)
#         optimisationdf.to_excel((path + r'\TSEMO-STANDARD' + str(exp) + '.xlsx'))

#     minParticleSizeFuncDict = {}
#     surfactantDict = {}
#     monomerDict = {}
#     seedDict = {}
#     particleSizeDict = {}
#     particleSizeFuncDict = {}
#     fileCounter = 0

#     for file in glob.glob(path + r"\*.xlsx"):
#         fileCounter += 1

#         data = pd.read_excel(file)
#         numExps = data.shape[0]

#         minParticleSizeFunc = []
#         minParticleSizeFuncTemp = 99999
#         for exp in range(numExps):
#             if data.iloc[exp, 5] < minParticleSizeFuncTemp:
#                 minParticleSizeFuncTemp = data.iloc[exp, 5]
#             minParticleSizeFunc.append(minParticleSizeFuncTemp)
        
#         minParticleSizeFuncDict[str(fileCounter)] = minParticleSizeFunc
        
#         surfactantConc = []
#         monomerFrac = []
#         seedFrac = []
#         particleSize = []
#         particleSizeFunc = []

#         for exp in range(numExps): # get experimental conditions from each 'reaction' of the experiments
#             surfactantConc.append(data.iloc[exp, 1])
#             monomerFrac.append(data.iloc[exp, 2])
#             seedFrac.append(data.iloc[exp, 3])
#             particleSize.append(data.iloc[exp, 4])
#             particleSizeFunc.append(data.iloc[exp, 5])


#         surfactantDict[str(fileCounter)] = surfactantConc
#         monomerDict[str(fileCounter)] = monomerFrac
#         seedDict[str(fileCounter)] = seedFrac
#         particleSizeDict[str(fileCounter)] = particleSize
#         particleSizeFuncDict[str(fileCounter)] = particleSizeFunc

#     surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
#     monomerSummary = pd.DataFrame.from_dict(monomerDict)
#     seedSummary = pd.DataFrame.from_dict(seedDict)
#     particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
#     particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
#     minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

#     surfactantTot = 0
#     surfactantAveList = []
#     monomerTot = 0
#     monomerAveList = []
#     seedTot = 0
#     seedAveList = []
#     particleSizeTot = 0
#     particleSizeAveList = []
#     particleSizeFuncTot = 0
#     particleSizeFuncAveList = []
#     minParticleSizeFuncTot = 0
#     minParticleSizeFuncAveList = []

#     for exp in range(numExps):
#         surfactantTot = sum(surfactantSummary.iloc[exp, :])
#         surfactantAveList.append(surfactantTot/len(surfactantSummary.columns))
#         monomerTot = sum(monomerSummary.iloc[exp, :])
#         monomerAveList.append(monomerTot/len(monomerSummary.columns))
#         seedTot = sum(seedSummary.iloc[exp, :])
#         seedAveList.append(seedTot/len(seedSummary.columns))
#         particleSizeTot = sum(particleSizeSummary.iloc[exp, :])
#         particleSizeAveList.append(particleSizeTot/len(particleSizeSummary.columns))
#         particleSizeFuncTot = sum(particleSizeFuncSummary.iloc[exp, :])
#         particleSizeFuncAveList.append(particleSizeFuncTot/len(particleSizeFuncSummary.columns))
#         minParticleSizeFuncTot = sum(minParticleSizeFuncSummary.iloc[exp, :])
#         minParticleSizeFuncAveList.append(minParticleSizeFuncTot/len(minParticleSizeFuncSummary.columns))
#     surfactantDict["Average"] = surfactantAveList
#     monomerDict["Average"] = monomerAveList
#     seedDict["Average"] = seedAveList
#     particleSizeDict["Average"] = particleSizeAveList
#     particleSizeFuncDict["Average"] = particleSizeFuncAveList
#     minParticleSizeFuncDict["Average"] = minParticleSizeFuncAveList

#     averagesDict = {
#         "Surfactant concentration" : surfactantAveList,
#         "Monomer fraction" : monomerAveList,
#         "Seed fraction" : seedAveList,
#         "Particle size" : particleSizeAveList,
#         "Particle size function" : particleSizeFuncAveList,
#         "Minimised particle size function" : minParticleSizeFuncAveList
#     }

#     surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
#     monomerSummary = pd.DataFrame.from_dict(monomerDict)
#     seedSummary = pd.DataFrame.from_dict(seedDict)
#     particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
#     particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
#     averagesSummary = pd.DataFrame.from_dict(averagesDict)
#     minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

#     surfactantSummary.to_excel(path + r'\surfactantSummary-' + expString + '.xlsx')
#     monomerSummary.to_excel(path + r'\monomerSummary-' + expString + '.xlsx')
#     seedSummary.to_excel(path + r'\seedSummary-' + expString + '.xlsx')
#     particleSizeSummary.to_excel(path + r'\particleSizeSummary-' + expString + '.xlsx')
#     particleSizeFuncSummary.to_excel(path + r'\particleSizeFuncSummary-' + expString + '.xlsx')
#     averagesSummary.to_excel(path + r'\averagesSummary-' + expString + '.xlsx')
#     minParticleSizeFuncSummary.to_excel(path + r'\minimisationSummary-' + expString + '.xlsx')

# #####################################################################################################################
# '''

# NEXT RUN

# '''
# #####################################################################################################################

# expString = r'standard-targeting-min-seed'
# inputs = {'surfactantConc' : [0.01, 0.05],
#             'seedFrac': [0.02, 0.2086]}
# objectives = {('particleSizeFunc', "MIN") : [0.00001, 50],
#               ('seedFracMin', "MIN") : [0.02, 0.2086]
#             }
# objectivesList = []
# for objective in objectives:
#     objectivesList.append(objective[0])
# variablesList = list(inputs.keys())
# variablesList.extend(objectivesList)

# sizeList = [80, 90, 100, 110, 120]
# n = 5 # number of optimisations to run (repeats)
# LHCs = 10

# for sizeTarget in sizeList:
#     path = str(r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Optimisation\Conventional\Simulated_experiments\CONSOLIDATED\STANDARD\SIZE-TARGETING\MINIMISE-SEED' + r'\_' + str(sizeTarget) + 'nm')

#     for exp in range(n):  
#         monomerFrac = []  
#         parameterNames = list(inputs.keys())
#         parameterBounds = list(inputs.values())
#         trainingExps = latinHypercubeSampling.getHypercubeSamplingParams(LHCs, parameter_names=parameterNames, parameter_bounds=parameterBounds)
#         surfactantConc = trainingExps.get('surfactantConc')
#         seedFrac = trainingExps.get('seedFrac')
#         # monomerFrac = trainingExps.get('monomerFrac')
#         monomerFracTemp = 0.5
#         for i in surfactantConc:
#             monomerFrac.append(monomerFracTemp)
#         trainingExps["surfactantConcMin"] = surfactantConc
#         trainingExps["seedFracMin"] = seedFrac
 

#         particleSize = []
#         for x, y, z, in zip(surfactantConc,
#                             monomerFrac,
#                             seedFrac):
#             particleSizeTemp = round(func((x, y, z), a, b, c, d, e, f, g, h, i, m), 2)
#             particleSize.append(particleSizeTemp)
#         trainingExps["particleSize"] = particleSize

#         particleSizeFunc = []
#         for size in particleSize:
#             particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#             particleSizeFunc.append(particleSizeFuncTemp)
#         trainingExps["particleSizeFunc"] = particleSizeFunc

#         LHCRandomStateValue = trainingExps.pop("Random state value")
#         current_dataSet = trainingExps
#         current_dataSet = pd.DataFrame.from_dict(current_dataSet)

#         iterations = len(particleSize)
#         while iterations < 30:

#             new_conditions = getNextExperiment(inputs, objectives, current_dataSet, iteration=iterations)
#             iterations += 1

#             x = new_conditions.get('surfactantConc')[0]
#             # y = new_conditions.get('monomerFrac')[0]
#             y = 0.5
#             z = new_conditions.get('seedFrac')[0]
#             X = x, y, z
#             newSize = round(func(X, a, b, c, d, e, f, g, h, i, m), 2)
            
#             ### Update data set ###
#             surfactantConc.append(x)
#             monomerFrac.append(monomerFracTemp)
#             seedFrac.append(z)
#             particleSize.append(newSize)
#             particleSizeFunc = []
#             for size in particleSize:
#                 particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#                 particleSizeFunc.append(particleSizeFuncTemp)

#             values = [surfactantConc, seedFrac, particleSizeFunc, surfactantConc, seedFrac]
#             current_dataSet = {}
#             for variable, value in zip(variablesList, values):
#                 current_dataSet[variable] = value
#             current_dataSet = pd.DataFrame.from_dict(current_dataSet)
#             current_dataSet = DataSet(current_dataSet)

#             print(particleSize)
#         print(str(iterations) + ' reactions completed')
            
#         LHC_randomStateExp = []
#         for i in particleSize:
#             LHC_randomStateExp.append(LHCRandomStateValue)

#         optimisationSeries = {
#             "Surfactant concentration" : surfactantConc,
#             "Monomer fraction" : monomerFrac,
#             "Seed fraction" : seedFrac,
#             "Particle sizes" : particleSize,
#             "Particle size function" : particleSizeFunc,
#             "LHC random state" : LHC_randomStateExp
#         }
#         optimisationdf = pd.DataFrame(optimisationSeries)
#         optimisationdf.to_excel((path + r'\TSEMO-STANDARD' + str(exp) + '.xlsx'))

#     minParticleSizeFuncDict = {}
#     surfactantDict = {}
#     monomerDict = {}
#     seedDict = {}
#     particleSizeDict = {}
#     particleSizeFuncDict = {}
#     fileCounter = 0

#     for file in glob.glob(path + r"\*.xlsx"):
#         fileCounter += 1

#         data = pd.read_excel(file)
#         numExps = data.shape[0]

#         minParticleSizeFunc = []
#         minParticleSizeFuncTemp = 99999
#         for exp in range(numExps):
#             if data.iloc[exp, 5] < minParticleSizeFuncTemp:
#                 minParticleSizeFuncTemp = data.iloc[exp, 5]
#             minParticleSizeFunc.append(minParticleSizeFuncTemp)
        
#         minParticleSizeFuncDict[str(fileCounter)] = minParticleSizeFunc
        
#         surfactantConc = []
#         monomerFrac = []
#         seedFrac = []
#         particleSize = []
#         particleSizeFunc = []

#         for exp in range(numExps): # get experimental conditions from each 'reaction' of the experiments
#             surfactantConc.append(data.iloc[exp, 1])
#             monomerFrac.append(data.iloc[exp, 2])
#             seedFrac.append(data.iloc[exp, 3])
#             particleSize.append(data.iloc[exp, 4])
#             particleSizeFunc.append(data.iloc[exp, 5])


#         surfactantDict[str(fileCounter)] = surfactantConc
#         monomerDict[str(fileCounter)] = monomerFrac
#         seedDict[str(fileCounter)] = seedFrac
#         particleSizeDict[str(fileCounter)] = particleSize
#         particleSizeFuncDict[str(fileCounter)] = particleSizeFunc

#     surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
#     monomerSummary = pd.DataFrame.from_dict(monomerDict)
#     seedSummary = pd.DataFrame.from_dict(seedDict)
#     particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
#     particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
#     minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

#     surfactantTot = 0
#     surfactantAveList = []
#     monomerTot = 0
#     monomerAveList = []
#     seedTot = 0
#     seedAveList = []
#     particleSizeTot = 0
#     particleSizeAveList = []
#     particleSizeFuncTot = 0
#     particleSizeFuncAveList = []
#     minParticleSizeFuncTot = 0
#     minParticleSizeFuncAveList = []

#     for exp in range(numExps):
#         surfactantTot = sum(surfactantSummary.iloc[exp, :])
#         surfactantAveList.append(surfactantTot/len(surfactantSummary.columns))
#         monomerTot = sum(monomerSummary.iloc[exp, :])
#         monomerAveList.append(monomerTot/len(monomerSummary.columns))
#         seedTot = sum(seedSummary.iloc[exp, :])
#         seedAveList.append(seedTot/len(seedSummary.columns))
#         particleSizeTot = sum(particleSizeSummary.iloc[exp, :])
#         particleSizeAveList.append(particleSizeTot/len(particleSizeSummary.columns))
#         particleSizeFuncTot = sum(particleSizeFuncSummary.iloc[exp, :])
#         particleSizeFuncAveList.append(particleSizeFuncTot/len(particleSizeFuncSummary.columns))
#         minParticleSizeFuncTot = sum(minParticleSizeFuncSummary.iloc[exp, :])
#         minParticleSizeFuncAveList.append(minParticleSizeFuncTot/len(minParticleSizeFuncSummary.columns))
#     surfactantDict["Average"] = surfactantAveList
#     monomerDict["Average"] = monomerAveList
#     seedDict["Average"] = seedAveList
#     particleSizeDict["Average"] = particleSizeAveList
#     particleSizeFuncDict["Average"] = particleSizeFuncAveList
#     minParticleSizeFuncDict["Average"] = minParticleSizeFuncAveList

#     averagesDict = {
#         "Surfactant concentration" : surfactantAveList,
#         "Monomer fraction" : monomerAveList,
#         "Seed fraction" : seedAveList,
#         "Particle size" : particleSizeAveList,
#         "Particle size function" : particleSizeFuncAveList,
#         "Minimised particle size function" : minParticleSizeFuncAveList
#     }

#     surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
#     monomerSummary = pd.DataFrame.from_dict(monomerDict)
#     seedSummary = pd.DataFrame.from_dict(seedDict)
#     particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
#     particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
#     averagesSummary = pd.DataFrame.from_dict(averagesDict)
#     minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

#     surfactantSummary.to_excel(path + r'\surfactantSummary-' + expString + '.xlsx')
#     monomerSummary.to_excel(path + r'\monomerSummary-' + expString + '.xlsx')
#     seedSummary.to_excel(path + r'\seedSummary-' + expString + '.xlsx')
#     particleSizeSummary.to_excel(path + r'\particleSizeSummary-' + expString + '.xlsx')
#     particleSizeFuncSummary.to_excel(path + r'\particleSizeFuncSummary-' + expString + '.xlsx')
#     averagesSummary.to_excel(path + r'\averagesSummary-' + expString + '.xlsx')
#     minParticleSizeFuncSummary.to_excel(path + r'\minimisationSummary-' + expString + '.xlsx')


# #####################################################################################################################
# '''

# NEXT RUN

# '''
# #####################################################################################################################

# expString = r'flipflop-mapping-min-seed'
# inputs = {'surfactantConc' : [0.01, 0.05],
#             'seedFrac': [0.02, 0.2086]}
# objectives = {('particleSize', "FLIPFLOP") : [10, 150],
#               ('seedFracMin', "MIN") : [0.02, 0.2086]
#             }
# objectivesList = []
# for objective in objectives:
#     objectivesList.append(objective[0])
# variablesList = list(inputs.keys())
# variablesList.extend(objectivesList)

# sizeList = [80, 90, 100, 110, 120]
# n = 5 # number of optimisations to run (repeats)
# LHCs = 10
# path = str(r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Optimisation\Conventional\Simulated_experiments\CONSOLIDATED\FLIP-FLOP\SIZE-MAPPING\MINIMISE-SEED')

# for exp in range(n):  
#     monomerFrac = []  
#     parameterNames = list(inputs.keys())
#     parameterBounds = list(inputs.values())
#     trainingExps = latinHypercubeSampling.getHypercubeSamplingParams(LHCs, parameter_names=parameterNames, parameter_bounds=parameterBounds)
#     surfactantConc = trainingExps.get('surfactantConc')
#     seedFrac = trainingExps.get('seedFrac')
#     # monomerFrac = trainingExps.get('monomerFrac')
#     monomerFracTemp = 0.5
#     for i in surfactantConc:
#         monomerFrac.append(monomerFracTemp)
#     trainingExps["surfactantConcMin"] = surfactantConc
#     trainingExps["seedFracMin"] = seedFrac

#     particleSize = []
#     for x, y, z, in zip(surfactantConc,
#                         monomerFrac,
#                         seedFrac):
#         particleSizeTemp = round(func((x, y, z), a, b, c, d, e, f, g, h, i, m), 2)
#         particleSize.append(particleSizeTemp)
#     trainingExps["particleSize"] = particleSize

#     particleSizeFunc = []
#     for size in particleSize:
#         particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#         particleSizeFunc.append(particleSizeFuncTemp)
#     trainingExps["particleSizeFunc"] = particleSizeFunc

#     LHCRandomStateValue = trainingExps.pop("Random state value")
#     current_dataSet = trainingExps
#     current_dataSet = pd.DataFrame.from_dict(current_dataSet)

#     iterations = len(particleSize)
#     while iterations < 30:

#         new_conditions = getNextExperiment(inputs, objectives, current_dataSet, iteration=iterations)
#         iterations += 1

#         x = new_conditions.get('surfactantConc')[0]
#         # y = new_conditions.get('monomerFrac')[0]
#         y = 0.5
#         z = new_conditions.get('seedFrac')[0]
#         X = x, y, z
#         newSize = round(func(X, a, b, c, d, e, f, g, h, i, m), 2)
        
#         ### Update data set ###
#         surfactantConc.append(x)
#         monomerFrac.append(monomerFracTemp)
#         seedFrac.append(z)
#         particleSize.append(newSize)
#         particleSizeFunc = []
#         for size in particleSize:
#             particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#             particleSizeFunc.append(particleSizeFuncTemp)

#         values = [surfactantConc, seedFrac, particleSizeFunc, surfactantConc, seedFrac]
#         current_dataSet = {}
#         for variable, value in zip(variablesList, values):
#             current_dataSet[variable] = value
#         current_dataSet = pd.DataFrame.from_dict(current_dataSet)
#         current_dataSet = DataSet(current_dataSet)

#         print(particleSize)
#     print(str(iterations) + ' reactions completed')
        
#     LHC_randomStateExp = []
#     for i in particleSize:
#         LHC_randomStateExp.append(LHCRandomStateValue)

#     optimisationSeries = {
#         "Surfactant concentration" : surfactantConc,
#         "Monomer fraction" : monomerFrac,
#         "Seed fraction" : seedFrac,
#         "Particle sizes" : particleSize,
#         "Particle size function" : particleSizeFunc,
#         "LHC random state" : LHC_randomStateExp
#     }
#     optimisationdf = pd.DataFrame(optimisationSeries)
#     optimisationdf.to_excel((path + r'\TSEMO-STANDARD' + str(exp) + '.xlsx'))

# minParticleSizeFuncDict = {}
# surfactantDict = {}
# monomerDict = {}
# seedDict = {}
# particleSizeDict = {}
# particleSizeFuncDict = {}
# fileCounter = 0

# for file in glob.glob(path + r"\*.xlsx"):
#     fileCounter += 1

#     data = pd.read_excel(file)
#     numExps = data.shape[0]

#     minParticleSizeFunc = []
#     minParticleSizeFuncTemp = 99999
#     for exp in range(numExps):
#         if data.iloc[exp, 5] < minParticleSizeFuncTemp:
#             minParticleSizeFuncTemp = data.iloc[exp, 5]
#         minParticleSizeFunc.append(minParticleSizeFuncTemp)
    
#     minParticleSizeFuncDict[str(fileCounter)] = minParticleSizeFunc
    
#     surfactantConc = []
#     monomerFrac = []
#     seedFrac = []
#     particleSize = []
#     particleSizeFunc = []

#     for exp in range(numExps): # get experimental conditions from each 'reaction' of the experiments
#         surfactantConc.append(data.iloc[exp, 1])
#         monomerFrac.append(data.iloc[exp, 2])
#         seedFrac.append(data.iloc[exp, 3])
#         particleSize.append(data.iloc[exp, 4])
#         particleSizeFunc.append(data.iloc[exp, 5])


#     surfactantDict[str(fileCounter)] = surfactantConc
#     monomerDict[str(fileCounter)] = monomerFrac
#     seedDict[str(fileCounter)] = seedFrac
#     particleSizeDict[str(fileCounter)] = particleSize
#     particleSizeFuncDict[str(fileCounter)] = particleSizeFunc

# surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
# monomerSummary = pd.DataFrame.from_dict(monomerDict)
# seedSummary = pd.DataFrame.from_dict(seedDict)
# particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
# particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
# minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

# surfactantTot = 0
# surfactantAveList = []
# monomerTot = 0
# monomerAveList = []
# seedTot = 0
# seedAveList = []
# particleSizeTot = 0
# particleSizeAveList = []
# particleSizeFuncTot = 0
# particleSizeFuncAveList = []
# minParticleSizeFuncTot = 0
# minParticleSizeFuncAveList = []

# for exp in range(numExps):
#     surfactantTot = sum(surfactantSummary.iloc[exp, :])
#     surfactantAveList.append(surfactantTot/len(surfactantSummary.columns))
#     monomerTot = sum(monomerSummary.iloc[exp, :])
#     monomerAveList.append(monomerTot/len(monomerSummary.columns))
#     seedTot = sum(seedSummary.iloc[exp, :])
#     seedAveList.append(seedTot/len(seedSummary.columns))
#     particleSizeTot = sum(particleSizeSummary.iloc[exp, :])
#     particleSizeAveList.append(particleSizeTot/len(particleSizeSummary.columns))
#     particleSizeFuncTot = sum(particleSizeFuncSummary.iloc[exp, :])
#     particleSizeFuncAveList.append(particleSizeFuncTot/len(particleSizeFuncSummary.columns))
#     minParticleSizeFuncTot = sum(minParticleSizeFuncSummary.iloc[exp, :])
#     minParticleSizeFuncAveList.append(minParticleSizeFuncTot/len(minParticleSizeFuncSummary.columns))
# surfactantDict["Average"] = surfactantAveList
# monomerDict["Average"] = monomerAveList
# seedDict["Average"] = seedAveList
# particleSizeDict["Average"] = particleSizeAveList
# particleSizeFuncDict["Average"] = particleSizeFuncAveList
# minParticleSizeFuncDict["Average"] = minParticleSizeFuncAveList

# averagesDict = {
#     "Surfactant concentration" : surfactantAveList,
#     "Monomer fraction" : monomerAveList,
#     "Seed fraction" : seedAveList,
#     "Particle size" : particleSizeAveList,
#     "Particle size function" : particleSizeFuncAveList,
#     "Minimised particle size function" : minParticleSizeFuncAveList
# }

# surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
# monomerSummary = pd.DataFrame.from_dict(monomerDict)
# seedSummary = pd.DataFrame.from_dict(seedDict)
# particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
# particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
# averagesSummary = pd.DataFrame.from_dict(averagesDict)
# minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

# surfactantSummary.to_excel(path + r'\surfactantSummary-' + expString + '.xlsx')
# monomerSummary.to_excel(path + r'\monomerSummary-' + expString + '.xlsx')
# seedSummary.to_excel(path + r'\seedSummary-' + expString + '.xlsx')
# particleSizeSummary.to_excel(path + r'\particleSizeSummary-' + expString + '.xlsx')
# particleSizeFuncSummary.to_excel(path + r'\particleSizeFuncSummary-' + expString + '.xlsx')
# averagesSummary.to_excel(path + r'\averagesSummary-' + expString + '.xlsx')
# minParticleSizeFuncSummary.to_excel(path + r'\minimisationSummary-' + expString + '.xlsx')

# #####################################################################################################################
# '''

# NEXT RUN

# '''
# #####################################################################################################################

# expString = r'flipflop-mapping-min-seed-and-surfactant'
# inputs = {'surfactantConc' : [0.01, 0.05],
#             'seedFrac': [0.02, 0.2086]}
# objectives = {('particleSize', "FLIPFLOP") : [10, 150],
#               ("surfactantConcMin", "MIN") : [0.01, 0.05],
#               ('seedFracMin', "MIN") : [0.02, 0.2086]
#             }
# objectivesList = []
# for objective in objectives:
#     objectivesList.append(objective[0])
# variablesList = list(inputs.keys())
# variablesList.extend(objectivesList)

# n = 5 # number of optimisations to run (repeats)
# LHCs = 10
# path = str(r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Optimisation\Conventional\Simulated_experiments\CONSOLIDATED\FLIP-FLOP\SIZE-MAPPING\MINIMISE-SURFACTANT-AND-SEED')

# for exp in range(n):  
#     monomerFrac = []  
#     parameterNames = list(inputs.keys())
#     parameterBounds = list(inputs.values())
#     trainingExps = latinHypercubeSampling.getHypercubeSamplingParams(LHCs, parameter_names=parameterNames, parameter_bounds=parameterBounds)
#     surfactantConc = trainingExps.get('surfactantConc')
#     seedFrac = trainingExps.get('seedFrac')
#     # monomerFrac = trainingExps.get('monomerFrac')
#     monomerFracTemp = 0.5
#     for i in surfactantConc:
#         monomerFrac.append(monomerFracTemp)
#     trainingExps["surfactantConcMin"] = surfactantConc
#     trainingExps["seedFracMin"] = seedFrac

#     particleSize = []
#     for x, y, z, in zip(surfactantConc,
#                         monomerFrac,
#                         seedFrac):
#         particleSizeTemp = round(func((x, y, z), a, b, c, d, e, f, g, h, i, m), 2)
#         particleSize.append(particleSizeTemp)
#     trainingExps["particleSize"] = particleSize

#     particleSizeFunc = []
#     for size in particleSize:
#         particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#         particleSizeFunc.append(particleSizeFuncTemp)
#     trainingExps["particleSizeFunc"] = particleSizeFunc

#     LHCRandomStateValue = trainingExps.pop("Random state value")
#     current_dataSet = trainingExps
#     current_dataSet = pd.DataFrame.from_dict(current_dataSet)

#     iterations = len(particleSize)
#     while iterations < 30:

#         new_conditions = getNextExperiment(inputs, objectives, current_dataSet, iteration=iterations)
#         iterations += 1

#         x = new_conditions.get('surfactantConc')[0]
#         # y = new_conditions.get('monomerFrac')[0]
#         y = 0.5
#         z = new_conditions.get('seedFrac')[0]
#         X = x, y, z
#         newSize = round(func(X, a, b, c, d, e, f, g, h, i, m), 2)
        
#         ### Update data set ###
#         surfactantConc.append(x)
#         monomerFrac.append(monomerFracTemp)
#         seedFrac.append(z)
#         particleSize.append(newSize)
#         particleSizeFunc = []
#         for size in particleSize:
#             particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#             particleSizeFunc.append(particleSizeFuncTemp)

#         values = [surfactantConc, seedFrac, particleSizeFunc, surfactantConc, seedFrac]
#         current_dataSet = {}
#         for variable, value in zip(variablesList, values):
#             current_dataSet[variable] = value
#         current_dataSet = pd.DataFrame.from_dict(current_dataSet)
#         current_dataSet = DataSet(current_dataSet)

#         print(particleSize)
#     print(str(iterations) + ' reactions completed')
        
#     LHC_randomStateExp = []
#     for i in particleSize:
#         LHC_randomStateExp.append(LHCRandomStateValue)

#     optimisationSeries = {
#         "Surfactant concentration" : surfactantConc,
#         "Monomer fraction" : monomerFrac,
#         "Seed fraction" : seedFrac,
#         "Particle sizes" : particleSize,
#         "Particle size function" : particleSizeFunc,
#         "LHC random state" : LHC_randomStateExp
#     }
#     optimisationdf = pd.DataFrame(optimisationSeries)
#     optimisationdf.to_excel((path + r'\TSEMO-STANDARD' + str(exp) + '.xlsx'))

# minParticleSizeFuncDict = {}
# surfactantDict = {}
# monomerDict = {}
# seedDict = {}
# particleSizeDict = {}
# particleSizeFuncDict = {}
# fileCounter = 0

# for file in glob.glob(path + r"\*.xlsx"):
#     fileCounter += 1

#     data = pd.read_excel(file)
#     numExps = data.shape[0]

#     minParticleSizeFunc = []
#     minParticleSizeFuncTemp = 99999
#     for exp in range(numExps):
#         if data.iloc[exp, 5] < minParticleSizeFuncTemp:
#             minParticleSizeFuncTemp = data.iloc[exp, 5]
#         minParticleSizeFunc.append(minParticleSizeFuncTemp)
    
#     minParticleSizeFuncDict[str(fileCounter)] = minParticleSizeFunc
    
#     surfactantConc = []
#     monomerFrac = []
#     seedFrac = []
#     particleSize = []
#     particleSizeFunc = []

#     for exp in range(numExps): # get experimental conditions from each 'reaction' of the experiments
#         surfactantConc.append(data.iloc[exp, 1])
#         monomerFrac.append(data.iloc[exp, 2])
#         seedFrac.append(data.iloc[exp, 3])
#         particleSize.append(data.iloc[exp, 4])
#         particleSizeFunc.append(data.iloc[exp, 5])


#     surfactantDict[str(fileCounter)] = surfactantConc
#     monomerDict[str(fileCounter)] = monomerFrac
#     seedDict[str(fileCounter)] = seedFrac
#     particleSizeDict[str(fileCounter)] = particleSize
#     particleSizeFuncDict[str(fileCounter)] = particleSizeFunc

# surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
# monomerSummary = pd.DataFrame.from_dict(monomerDict)
# seedSummary = pd.DataFrame.from_dict(seedDict)
# particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
# particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
# minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

# surfactantTot = 0
# surfactantAveList = []
# monomerTot = 0
# monomerAveList = []
# seedTot = 0
# seedAveList = []
# particleSizeTot = 0
# particleSizeAveList = []
# particleSizeFuncTot = 0
# particleSizeFuncAveList = []
# minParticleSizeFuncTot = 0
# minParticleSizeFuncAveList = []

# for exp in range(numExps):
#     surfactantTot = sum(surfactantSummary.iloc[exp, :])
#     surfactantAveList.append(surfactantTot/len(surfactantSummary.columns))
#     monomerTot = sum(monomerSummary.iloc[exp, :])
#     monomerAveList.append(monomerTot/len(monomerSummary.columns))
#     seedTot = sum(seedSummary.iloc[exp, :])
#     seedAveList.append(seedTot/len(seedSummary.columns))
#     particleSizeTot = sum(particleSizeSummary.iloc[exp, :])
#     particleSizeAveList.append(particleSizeTot/len(particleSizeSummary.columns))
#     particleSizeFuncTot = sum(particleSizeFuncSummary.iloc[exp, :])
#     particleSizeFuncAveList.append(particleSizeFuncTot/len(particleSizeFuncSummary.columns))
#     minParticleSizeFuncTot = sum(minParticleSizeFuncSummary.iloc[exp, :])
#     minParticleSizeFuncAveList.append(minParticleSizeFuncTot/len(minParticleSizeFuncSummary.columns))
# surfactantDict["Average"] = surfactantAveList
# monomerDict["Average"] = monomerAveList
# seedDict["Average"] = seedAveList
# particleSizeDict["Average"] = particleSizeAveList
# particleSizeFuncDict["Average"] = particleSizeFuncAveList
# minParticleSizeFuncDict["Average"] = minParticleSizeFuncAveList

# averagesDict = {
#     "Surfactant concentration" : surfactantAveList,
#     "Monomer fraction" : monomerAveList,
#     "Seed fraction" : seedAveList,
#     "Particle size" : particleSizeAveList,
#     "Particle size function" : particleSizeFuncAveList,
#     "Minimised particle size function" : minParticleSizeFuncAveList
# }

# surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
# monomerSummary = pd.DataFrame.from_dict(monomerDict)
# seedSummary = pd.DataFrame.from_dict(seedDict)
# particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
# particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
# averagesSummary = pd.DataFrame.from_dict(averagesDict)
# minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

# surfactantSummary.to_excel(path + r'\surfactantSummary-' + expString + '.xlsx')
# monomerSummary.to_excel(path + r'\monomerSummary-' + expString + '.xlsx')
# seedSummary.to_excel(path + r'\seedSummary-' + expString + '.xlsx')
# particleSizeSummary.to_excel(path + r'\particleSizeSummary-' + expString + '.xlsx')
# particleSizeFuncSummary.to_excel(path + r'\particleSizeFuncSummary-' + expString + '.xlsx')
# averagesSummary.to_excel(path + r'\averagesSummary-' + expString + '.xlsx')
# minParticleSizeFuncSummary.to_excel(path + r'\minimisationSummary-' + expString + '.xlsx')


# #####################################################################################################################
# '''

# NEXT RUN

# '''
# #####################################################################################################################

# expString = r'flipflop-mapping-min-surfactant'
# inputs = {'surfactantConc' : [0.01, 0.05],
#             'seedFrac': [0.02, 0.2086]}
# objectives = {('particleSize', "FLIPFLOP") : [10, 150],
#               ("surfactantConcMin", "MIN") : [0.01, 0.05]
#             }
# objectivesList = []
# for objective in objectives:
#     objectivesList.append(objective[0])
# variablesList = list(inputs.keys())
# variablesList.extend(objectivesList)

# n = 5 # number of optimisations to run (repeats)
# LHCs = 10
# path = str(r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Optimisation\Conventional\Simulated_experiments\CONSOLIDATED\FLIP-FLOP\SIZE-MAPPING\MINIMISE-SURFACTANT')

# for exp in range(n):  
#     monomerFrac = []  
#     parameterNames = list(inputs.keys())
#     parameterBounds = list(inputs.values())
#     trainingExps = latinHypercubeSampling.getHypercubeSamplingParams(LHCs, parameter_names=parameterNames, parameter_bounds=parameterBounds)
#     surfactantConc = trainingExps.get('surfactantConc')
#     seedFrac = trainingExps.get('seedFrac')
#     # monomerFrac = trainingExps.get('monomerFrac')
#     monomerFracTemp = 0.5
#     for i in surfactantConc:
#         monomerFrac.append(monomerFracTemp)
#     trainingExps["surfactantConcMin"] = surfactantConc
#     trainingExps["seedFracMin"] = seedFrac

#     particleSize = []
#     for x, y, z, in zip(surfactantConc,
#                         monomerFrac,
#                         seedFrac):
#         particleSizeTemp = round(func((x, y, z), a, b, c, d, e, f, g, h, i, m), 2)
#         particleSize.append(particleSizeTemp)
#     trainingExps["particleSize"] = particleSize

#     particleSizeFunc = []
#     for size in particleSize:
#         particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#         particleSizeFunc.append(particleSizeFuncTemp)
#     trainingExps["particleSizeFunc"] = particleSizeFunc

#     LHCRandomStateValue = trainingExps.pop("Random state value")
#     current_dataSet = trainingExps
#     current_dataSet = pd.DataFrame.from_dict(current_dataSet)

#     iterations = len(particleSize)
#     while iterations < 30:

#         new_conditions = getNextExperiment(inputs, objectives, current_dataSet, iteration=iterations)
#         iterations += 1

#         x = new_conditions.get('surfactantConc')[0]
#         # y = new_conditions.get('monomerFrac')[0]
#         y = 0.5
#         z = new_conditions.get('seedFrac')[0]
#         X = x, y, z
#         newSize = round(func(X, a, b, c, d, e, f, g, h, i, m), 2)
        
#         ### Update data set ###
#         surfactantConc.append(x)
#         monomerFrac.append(monomerFracTemp)
#         seedFrac.append(z)
#         particleSize.append(newSize)
#         particleSizeFunc = []
#         for size in particleSize:
#             particleSizeFuncTemp = round(((size - sizeTarget)/sizeTarget)**2, 5)
#             particleSizeFunc.append(particleSizeFuncTemp)

#         values = [surfactantConc, seedFrac, particleSizeFunc, surfactantConc, seedFrac]
#         current_dataSet = {}
#         for variable, value in zip(variablesList, values):
#             current_dataSet[variable] = value
#         current_dataSet = pd.DataFrame.from_dict(current_dataSet)
#         current_dataSet = DataSet(current_dataSet)

#         print(particleSize)
#     print(str(iterations) + ' reactions completed')
        
#     LHC_randomStateExp = []
#     for i in particleSize:
#         LHC_randomStateExp.append(LHCRandomStateValue)

#     optimisationSeries = {
#         "Surfactant concentration" : surfactantConc,
#         "Monomer fraction" : monomerFrac,
#         "Seed fraction" : seedFrac,
#         "Particle sizes" : particleSize,
#         "Particle size function" : particleSizeFunc,
#         "LHC random state" : LHC_randomStateExp
#     }
#     optimisationdf = pd.DataFrame(optimisationSeries)
#     optimisationdf.to_excel((path + r'\TSEMO-STANDARD' + str(exp) + '.xlsx'))

# minParticleSizeFuncDict = {}
# surfactantDict = {}
# monomerDict = {}
# seedDict = {}
# particleSizeDict = {}
# particleSizeFuncDict = {}
# fileCounter = 0

# for file in glob.glob(path + r"\*.xlsx"):
#     fileCounter += 1

#     data = pd.read_excel(file)
#     numExps = data.shape[0]

#     minParticleSizeFunc = []
#     minParticleSizeFuncTemp = 99999
#     for exp in range(numExps):
#         if data.iloc[exp, 5] < minParticleSizeFuncTemp:
#             minParticleSizeFuncTemp = data.iloc[exp, 5]
#         minParticleSizeFunc.append(minParticleSizeFuncTemp)
    
#     minParticleSizeFuncDict[str(fileCounter)] = minParticleSizeFunc
    
#     surfactantConc = []
#     monomerFrac = []
#     seedFrac = []
#     particleSize = []
#     particleSizeFunc = []

#     for exp in range(numExps): # get experimental conditions from each 'reaction' of the experiments
#         surfactantConc.append(data.iloc[exp, 1])
#         monomerFrac.append(data.iloc[exp, 2])
#         seedFrac.append(data.iloc[exp, 3])
#         particleSize.append(data.iloc[exp, 4])
#         particleSizeFunc.append(data.iloc[exp, 5])


#     surfactantDict[str(fileCounter)] = surfactantConc
#     monomerDict[str(fileCounter)] = monomerFrac
#     seedDict[str(fileCounter)] = seedFrac
#     particleSizeDict[str(fileCounter)] = particleSize
#     particleSizeFuncDict[str(fileCounter)] = particleSizeFunc

# surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
# monomerSummary = pd.DataFrame.from_dict(monomerDict)
# seedSummary = pd.DataFrame.from_dict(seedDict)
# particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
# particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
# minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

# surfactantTot = 0
# surfactantAveList = []
# monomerTot = 0
# monomerAveList = []
# seedTot = 0
# seedAveList = []
# particleSizeTot = 0
# particleSizeAveList = []
# particleSizeFuncTot = 0
# particleSizeFuncAveList = []
# minParticleSizeFuncTot = 0
# minParticleSizeFuncAveList = []

# for exp in range(numExps):
#     surfactantTot = sum(surfactantSummary.iloc[exp, :])
#     surfactantAveList.append(surfactantTot/len(surfactantSummary.columns))
#     monomerTot = sum(monomerSummary.iloc[exp, :])
#     monomerAveList.append(monomerTot/len(monomerSummary.columns))
#     seedTot = sum(seedSummary.iloc[exp, :])
#     seedAveList.append(seedTot/len(seedSummary.columns))
#     particleSizeTot = sum(particleSizeSummary.iloc[exp, :])
#     particleSizeAveList.append(particleSizeTot/len(particleSizeSummary.columns))
#     particleSizeFuncTot = sum(particleSizeFuncSummary.iloc[exp, :])
#     particleSizeFuncAveList.append(particleSizeFuncTot/len(particleSizeFuncSummary.columns))
#     minParticleSizeFuncTot = sum(minParticleSizeFuncSummary.iloc[exp, :])
#     minParticleSizeFuncAveList.append(minParticleSizeFuncTot/len(minParticleSizeFuncSummary.columns))
# surfactantDict["Average"] = surfactantAveList
# monomerDict["Average"] = monomerAveList
# seedDict["Average"] = seedAveList
# particleSizeDict["Average"] = particleSizeAveList
# particleSizeFuncDict["Average"] = particleSizeFuncAveList
# minParticleSizeFuncDict["Average"] = minParticleSizeFuncAveList

# averagesDict = {
#     "Surfactant concentration" : surfactantAveList,
#     "Monomer fraction" : monomerAveList,
#     "Seed fraction" : seedAveList,
#     "Particle size" : particleSizeAveList,
#     "Particle size function" : particleSizeFuncAveList,
#     "Minimised particle size function" : minParticleSizeFuncAveList
# }

# surfactantSummary = pd.DataFrame.from_dict(surfactantDict)
# monomerSummary = pd.DataFrame.from_dict(monomerDict)
# seedSummary = pd.DataFrame.from_dict(seedDict)
# particleSizeSummary = pd.DataFrame.from_dict(particleSizeDict)
# particleSizeFuncSummary = pd.DataFrame.from_dict(particleSizeFuncDict)
# averagesSummary = pd.DataFrame.from_dict(averagesDict)
# minParticleSizeFuncSummary = pd.DataFrame.from_dict(minParticleSizeFuncDict)

# surfactantSummary.to_excel(path + r'\surfactantSummary-' + expString + '.xlsx')
# monomerSummary.to_excel(path + r'\monomerSummary-' + expString + '.xlsx')
# seedSummary.to_excel(path + r'\seedSummary-' + expString + '.xlsx')
# particleSizeSummary.to_excel(path + r'\particleSizeSummary-' + expString + '.xlsx')
# particleSizeFuncSummary.to_excel(path + r'\particleSizeFuncSummary-' + expString + '.xlsx')
# averagesSummary.to_excel(path + r'\averagesSummary-' + expString + '.xlsx')
# minParticleSizeFuncSummary.to_excel(path + r'\minimisationSummary-' + expString + '.xlsx')
