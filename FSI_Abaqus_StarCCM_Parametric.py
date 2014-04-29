# This python file runs the other scripts in the FSI building and running process
# 	Written by Casey J. Jesse on June 2014 at the University of Missouri - Columbia
#	Revisions:
#		July 11, 2013 -  Plate mesh parameters were added to create biases in the plate's mesh 
#		December 17, 2013 - Added subroutine to change variables in the input file for parametric studies
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
	
def changeParameters(fileName, parameter2Change, newValue):
	with open(fileName, 'r') as file:
		newFileLines = []
		i = 0;
		for line in file:
			if line.startswith('#'):
				newFileLines.append(line)
			elif line.startswith('\n'):
				newFileLines.append(line)
			elif line.startswith(parameter2Change):
				name, type, value, comment = line.split(':')
				value = str(newValue)
				newFileLines.append(name + ':' + type + ': \t\t' + value + ':' + comment)
			else:
				newFileLines.append(line)
			i = i + 1;
		file.close
		file = open(fileName, 'w')
		for line in newFileLines:
			file.write(line)
		

#-----------------------------------------------------------------------------------------------------------
# MAIN FUNCTION
#subprocess.call('mode con cols=150 lines=150', shell = True)
parameters = fileReader('FSI_Input_File.txt')

#changeParameters('FSI_Input_File.txt', 'chBiasDirection', 'Center')
for wallBias in range(750, 1250, 250):
	changeParameters('FSI_Input_File.txt', 'flSmChHeightBias', float(wallBias)/10)
	changeParameters('FSI_Input_File.txt', 'flLgChHeightBias', float(wallBias)/10)
	changeParameters('FSI_Input_File.txt', 'flPlHeightBias', float(wallBias)/10)
	
	# Running the Abaqus Python script for building the plate and fluid models
	buildAbaqusCall = 'abaqus cae noGUI=FSI_GeometryBuilder.py'
	subprocess.call(buildAbaqusCall, shell = True)

	# Running the Star-CCM+ Java macro for building the Star-CCM+ simulation
	buildStarCCMCall = 'starccm+ -new -np 1 -batch AbaqusMeshingFSI.java'
	subprocess.call(buildStarCCMCall)

	parameters = fileReader('FSI_Input_File.txt')
	# Running the FSI simulation
	if parameters['runStar'] == 'yes':
		couplingScheme = parameters['couplingScheme']
		numStarProcesses = parameters['starProcesses']
		vel = str(int(parameters['avgChVelocity']))
		wallHeight = str(int(parameters['flSmChHeightBias']*10))
		plateGeometry = parameters['plateGeometry']
		runStarCall = 'starccm+ -np ' + numStarProcesses + '-time -batch ' + "FSI_" +vel + "_" + wallHeight + "_" + plateGeometry + '_' + str(int(parameters['plateThickness']/0.0254*1000)) + '_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + str(int(parameters['lgChHeight']/0.0254*1000)) + '.sim'
		subprocess.call(runStarCall)