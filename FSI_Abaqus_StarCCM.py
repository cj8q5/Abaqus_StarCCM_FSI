# This python file runs the other scripts in the FSI building and running process
# 	Written by Casey J. Jesse on June 2014 at the University of Missouri - Columbia
#	Revisions:
#		July 11, 2013 -  Plate mesh parameters were added to create biases in the plate's mesh 
#
import subprocess

def fileReader(fileName):
    parameters = {}
    with open(fileName, 'r') as file:
    	for line in file:
    		if line.startswith('#'):
    			continue
    		elif line.startswith('\n'):
    			continue
    		else:
    			line = line.replace('\t','')
    			name, type, value, comment = line.split(':')
    			if type == 'float':
    				parameters[str(name)] = float(value);
    			elif type == 'integer':
    				parameters[str(name)] = int(value);
    			elif type == 'string':
    				parameters[str(name)] = str(value);
    return parameters

#-----------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
#subprocess.call('mode con cols=150 lines=150', shell = True)
parameters = fileReader('FSI_Input_File.txt')
# Running the Abaqus Python script for building the plate and fluid models
buildAbaqusCall = 'abaqus cae noGUI=FSI_GeometryBuilder.py'
subprocess.call(buildAbaqusCall, shell = True)

# Running the Star-CCM+ Java macro for building the Star-CCM+ simulation
buildStarCCMCall = 'starccm+ -new -np 1 -batch AbaqusMeshingFSI.java'
subprocess.call(buildStarCCMCall)

# Running the FSI simulation
if parameters['runStar'] == 'yes':
	couplingScheme = parameters['couplingScheme']
	numStarProcesses = parameters['starProcesses']
	vel = str(int(parameters['avgChVelocity']))
	plateGeometry = parameters['plateGeometry']
	runStarCall = 'starccm+ -np ' + numStarProcesses + '-time -batch ' + couplingScheme + "_" + vel + '_Star_' + plateGeometry + '_' + str(int(parameters['plateThickness']/0.0254*1000)) + '_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + str(int(parameters['lgChHeight']/0.0254*1000)) + '.sim'
	subprocess.call(runStarCall)
