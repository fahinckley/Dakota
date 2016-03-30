# -*- coding: utf-8 -*-
"""
Analysis file for Cape Thompson data using Dakota

@author: Franklin
"""

import sys
from subprocess import call
import numpy as np


# Set the file names to use for the input file and template input file. The template
# input file is a text file we create that contains the string "{coefficient}". The
# actual input file, 'inputs.txt', will be created automatically by Dakota, which will
# replace "{coefficient}" with a number.
input_file_template = 'input_template.txt'
input_file = 'inputs.txt'

# Run Dakota's "dprepro" (preprocessor) program. This will make a new copy of
# 'my_inputs_template.txt' but replace the string that reads "{coefficient}" with an
# actual numerical value (or rather, a string representing a numeric value). It will then
# save the file with the name 'inputs.txt'. The actual value will be obtained from
# a third, Dakota-created file called 'params.in' (this same name is also passed in as
# an argument, and because it is the first argument, it is in sys.argv[1])
call(['dprepro', sys.argv[1], input_file_template, input_file])

# Get the input value for the model
infile = open(input_file)
tStep = float( infile.readline() )
infile.close()
print 'Running with tstep value', tStep, 'years'

# Define static parameters
k = 2.5; # [W/mK]
kappa = 1e-6; # [m^2/sec]
rhoCp = k/kappa;
Qm = 45.0/1000; # [W/m^2]

# Load and parse Reference Cape Thompson data
cape_thompson = np.loadtxt("cape_thompson.txt")
CT_z = cape_thompson[:,0]
CT_T = cape_thompson[:,1]

# Year of data set
t_CT = 1961;

# Fit step change to Cape Thompson data
# Estimate original surface temperature by extrapolating slope at depth
pd = np.polyfit(CT_z[11:],CT_T[11:],1,rcond=None,full=False,w=None,cov=False)
CT_TS0 = np.polyval(pd,0)

# Estimate original temperature profile for steady state
xp = np.concatenate(([0],CT_z))
CT_T0 = np.polyval(pd,xp);

# Assert current surface temperature is the first data point
CT_TS = CT_T[1]
CT_T = np.concatenate(([CT_T[1]],CT_T))

# Compute magnitude of step
Tstep = (CT_TS - CT_TS0) # [deg C] 

# Determine length step and time step
dz = np.min(np.diff(CT_z))
dt = 0.1*(dz**2/(2*kappa))

# Initialize temperatures
T = CT_T0

# Compute expected temperature profile
secYr = 86400*365.25
t = np.arange(0,tStep*secYr,dt)
for ii in range(0,len(t),1):
    # Get surface temperature at current time
    TsC = Tstep;
    
    # Determine temperature gradient at depth
    dTdzmax = Qm/k # [C/m]

    # Compute temperature gradient
    dTdz = np.diff(T)/dz
    dTdz = np.append(dTdz,dTdzmax)

    # Compute heat flux
    q = -k*dTdz

    # Compute temperature rate
    dqdz = np.diff(q)/dz
    dTdt = (-1.0/(rhoCp))*dqdz

    # Update temperature
    T[1:] = T[1:] + dTdt*dt
    T[0] = TsC


# Find the root-mean-square misfit between model and data
rms_error = np.sqrt( np.sum( (CT_T - T)**2 ))
print 'RMS error with tstep', tStep, 'is', rms_error

# Write the RMS error to a file
outfile = open(sys.argv[2], 'w')
outfile.write(str(rms_error))
outfile.close()
