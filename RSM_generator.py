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
savePath = r'add save path here'

'''
Define variable names [x, y, z, response]
'''
varNames = ["xData", "yData", "zData", "response"]

'''
Insert experimental input values and response data; data = [a, b, c, ..., n]; dataRange = [min, max]
'''
xData = []
xRange = []
yData = []
yRange = []
zData = []
zRange = []
response = []

surfaceName = varNames[-1] # Sets the file save name as the response variable
customLabel = [""] # Add custom label to save file

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


