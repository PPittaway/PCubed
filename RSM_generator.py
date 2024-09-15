import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from types import SimpleNamespace
import matplotlib.pyplot as plt

### Function to fit ###
def func(X, a, b, c, d, e, f, g, h, i, m):
    x, y, z = X
    return a*x + b*y + c*z + d*x**2 + e*y**2 + f*z**2 + g*x*y + h*x*z + i*y*z + m

''' 
Define path to save response surface plot
'''
savePath = r'C:\Users\pm15pmp\OneDrive - University of Leeds\Research\Year 4\Recipe-3\240429_optimisation\Response_surface_modelling\Updated_240911\Syn2-BA_conversion'

'''
Define variable names [x, y, z, response]
'''
varNames = ["init_ratio", "surf_ratio", "seed_frac", "Syn2-BA_conversion"]

'''
Insert experimental input values and response data
'''
xData = [0.45, 0.714, 0.582, 0.494, 0.406, 0.758, 0.362, 0.67, 0.538, 0.626, 0.37369, 0.77806, 0.39281, 0.40633, 0.44177, 0.7796, 0.34, 0.34152, 0.37176, 0.34055, 0.34, 0.34453, 0.34082, 0.48603, 0.59227]
xRange = [0.34, 0.78]
yData = [2.0525, 1.8945, 1.7365, 1.8155, 1.5785, 1.9735, 2.2105, 1.4995, 1.6575, 2.1315, 1.49841, 2.24828, 1.7461, 1.47916, 1.46089, 1.46053, 1.46435, 1.51182, 1.46264, 1.81034, 1.46002, 1.69688, 1.92014, 1.54519, 2.25]
yRange = [1.46, 2.25]
zData = [0.031, 0.087, 0.073, 0.017, 0.129, 0.045, 0.115, 0.059, 0.101, 0.143, 0.08922, 0.07455, 0.02342, 0.07069, 0.02092, 0.14613, 0.01876, 0.05626, 0.01155, 0.08458, 0.01, 0.06405, 0.09039, 0.11445, 0.01093]
zRange = [0.15, 0.15]
response = [0.707605136, 0.797760519, 0.771039559, 0.666715087, 0.735104652, 0.932526199, 0.739772894, 0.647135884, 0.702242519, 0.863055184, 0.675674148, 0.706145165, 0.599420283, 0.557237262, 0.595438406, 0.80394028, 0.553238107, 0.544337868, 0.608538053, 0.72505826, 0.619414949, 0.722443304, 0.674074701, 0.702801874, 0.667100864]

surfaceName = varNames[-1] # Sets the file save name as the response variable
customLabel = ["0.15_seed_slice"] # Add custom label to save file

### Handles file path naming ###
if not customLabel[0]=="":
    surfaceName = surfaceName + f'-{customLabel[0]}'
else: pass

### Provide initial guess for fit parameters ###
p0 = 1, 1, 1, 1, 1, 1, 1, 1, 1, 1

### Pass experimental data to fitting function and store each fitting parameter ###
fit = curve_fit(func, (xData, yData, zData), response, p0, maxfev=100000)

fitVars = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "m"]
fitDict = {}
index = 0
for var in fitVars:
    fitDict[var] = fit[0][index]
    index+=1

### Save fitting parameters to path ###
fit_df = pd.DataFrame(fitDict, index=[0])
fit_df.to_csv(savePath + f'\{surfaceName}-fit_params.csv', index=False)
fit = SimpleNamespace(**fitDict)

### Define resolution of response surface and create mesh bound by the specified variable ranges ###
n = 50
xMesh = np.round(np.linspace(xRange[0], xRange[1], n), 4)
yMesh = np.round(np.linspace(yRange[0], yRange[1], n), 4)
zMesh = np.round(np.linspace(zRange[0], zRange[1], n), 4)

X, Y, Z = np.meshgrid(xMesh, yMesh, zMesh)

### Extract the meshgrid to x, y, z values ###
xValues = []
yValues = []
zValues = []
rows = range(0, (n+1), 1)
for i in range(0, (n+1), 1):
    for row in rows:
        xrow = X[i-1][row-1, :]
        xValues.extend(xrow)

    for row in rows:
        yrow = Y[i-1][row-1, :]
        yValues.extend(yrow)

    for row in rows:
        zrow = Z[i-1][row-1, :]
        zValues.extend(zrow)

### Use response surface function to find a response value for each meshgrid point ###
responseValues = []
for x, y, z in zip(xValues, yValues, zValues):
    responseValueTemp = round(func((x, y, z),
                                   fit.a,
                                   fit.b,
                                   fit.c,
                                   fit.d,
                                   fit.e,
                                   fit.f,
                                   fit.g,
                                   fit.h,
                                   fit.i,
                                   fit.m),
                                   4)
    responseValues.append(responseValueTemp)

### Collect all values and create dict to save to path ###
varValues = [xValues, yValues, zValues, responseValues]

dict = {}
for name, values in zip(varNames, varValues):
    dict[name] = values

df = pd.DataFrame(dict)
df.to_csv(savePath + f'\{surfaceName}.csv')

### Plotting for quick visualisation ###
name_color_map = 'gist_rainbow'
x = df[varNames[0]] 
y = df[varNames[1]] 
z = df[varNames[2]] 
c = df[varNames[3]]

fig = plt.figure(); 
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel(varNames[0]) 
ax.set_ylabel(varNames[1])
ax.set_zlabel(varNames[2])

img = ax.scatter(x, y, z, c = c, cmap = name_color_map)
cbar = fig.colorbar(img, shrink=0.5, aspect=5)
cbar.ax.get_yaxis().labelpad = 15; cbar.ax.set_ylabel(varNames[3], rotation = 270)

# plt.show()


