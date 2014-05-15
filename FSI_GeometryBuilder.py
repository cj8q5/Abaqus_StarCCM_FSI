# This python file runs in Abaqus to create the plate and fluid geometry/mesh
# 	Written by Casey J. Jesse on June 2014 at the University of Missouri - Columbia
#	Revisions:
#		July 11, 2013 -  Plate mesh parameters were added to create biases in the plate's mesh 
#
# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import subprocess
import math

#-----------------------------------------------------------------------------------------------------------
'''
	This method creates a flat plate for a FSI simulation
'''
def createFlatPlate(geometryParameters):
	# GRABBING ALL OF THE GEOMETRY PARAMETERS
	plateName = 'Plate'
	modelName = parameters['abaqusModelName']
	plateLength = parameters['plateLength']
	plateWidth = parameters['plateWidth']
	plateThickness = parameters['plateThickness']
	smChHeight = parameters['smChHeight']
	lgChHeight = parameters['lgChHeight']
	inletLength = parameters['inletPlLength']
	outletLength = parameters['outletPlLength']
	numOfPlates = parameters['numOfPlates']

	# Grabbing all of the mesh parameters
	elemType = parameters['elemType']
	plateThickNodes = parameters['plateThickNodes']
	plateLengthNodes = parameters['plateLengthNodes']
	plateLengthBias = parameters['plateLengthBias']
	plateWidthNodes = parameters['plateWidthNodes']
	clampedWidthNodes = parameters['clampedWidthNodes']
	biasDirection = parameters['chBiasDirection']

	# Grabbing all of the material properties of the plate
	plateMaterial = parameters['plateMaterial']
	elasticModulus = parameters['elasticModulus']
	poissonsRatio = parameters['poissonsRatio']
	plateDensity = parameters['plateDensity']
	pinOrCombBC = parameters['pinOrCombBC']

	# Grabbing the FSI coupling parameters
	guessedAbaqusStep = parameters['guessedAbaqusStep']
	couplingScheme = parameters['couplingScheme']
	timeStep = parameters['timeStep']
	minTimeStep = parameters['minTimeStep']
	maxSimTime = parameters['maxSimTime']

	mdb.Model(name=modelName, modelType = STANDARD_EXPLICIT)

	##--------------------------------------------PARTS NODE-------------------------------------------------------
	# Creating part named 'FlatPlate'
	plateSize = [ (-0.0127, 0.0), (plateWidth + 0.0127, plateLength), (plateThickness) ]
	plate = createBox(modelName, plateName, plateSize)

	# Creating three datum planes on the plate part
	plate.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
	plate.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=plateWidth)
	plate.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=plateWidth*0.5)

	# Grabbing the datum planes and the cells of the plate for partitioning the plate
	plateCells = plate.cells
	datumPlanes = plate.datums

	# Creating a set of the entire plate geometry
	entirePlateGeometry = plateCells.getByBoundingBox(
	                          xMin = -1, yMin = -1, zMin = -1,
	                          xMax = 1, yMax = 1, zMax = 1)
	plate.Set(cells = entirePlateGeometry, name = 'EntirePlateGeometry')

	# Partitioning the plate into three sections, wetted region and clamped regions
	pickedCells = plateCells.findAt(((plateWidth*0.25, plateLength*0.5, plateThickness*0.5),))
	plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[2], cells=pickedCells)

	pickedCells = plateCells.findAt(((plateWidth*0.25, plateLength*0.5, plateThickness*0.5),))
	plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[3], cells=pickedCells)

	pickedCells = plateCells.findAt(((plateWidth*0.25, plateLength*0.5, plateThickness*0.5),))
	plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[4], cells=pickedCells)

	# Creating all the sets of edges in the plate part
	plateEdges = plate.edges
	flowLengthEdges = plateEdges.findAt(
	                        ((-0.0127, 				plateLength*0.5,	0.0),),
	                        ((-0.0127, 				plateLength*0.5,	plateThickness),),
	                        ((0.0, 					plateLength*0.5,	0.0),),
	                        ((0.0, 					plateLength*0.5,	plateThickness),),
	                        ((plateWidth*0.5, 		plateLength*0.5,	0.0),),
	                        ((plateWidth*0.5, 		plateLength*0.5,	plateThickness),),
	                        ((plateWidth, 			plateLength*0.5,	0.0),),
	                        ((plateWidth, 			plateLength*0.5,	plateThickness),),
	                        ((plateWidth+0.0127, 	plateLength*0.5,	0.0),),
	                        ((plateWidth+0.0127, 	plateLength*0.5,	plateThickness),))
	plate.Set(edges=flowLengthEdges, name='PlateLength')
	
	flowWidthEdges = plateEdges.findAt(
			                        ((plateWidth*0.25, 		plateLength, 	0.0),),
	                        ((plateWidth*0.25, 		plateLength, 	plateThickness),),
	                        ((plateWidth*0.25, 		0.0, 			0.0),),
	                        ((plateWidth*0.25, 		0.0, 			plateThickness),),
	                        ((plateWidth*0.75, 		plateLength, 	0.0),),
	                        ((plateWidth*0.75, 		plateLength, 	plateThickness),),
	                        ((plateWidth*0.75, 		0.0, 			0.0),),
	                        ((plateWidth*0.75, 		0.0, 			plateThickness),))
	plate.Set(edges=flowWidthEdges, name = 'PlateWidth')

	plateThicknessEdges = plateEdges.findAt(
	                        ((-0.0127, 				0.0, 			plateThickness*0.5),), 
							((0.0, 					0.0, 			plateThickness*0.5),),
	                        ((plateWidth*0.5, 		0.0, 			plateThickness*0.5),), 
							((plateWidth, 			0.0, 			plateThickness*0.5),),
	                        ((plateWidth+0.0127,	0.0, 			plateThickness*0.5),),
							((-0.0127, 				plateLength, 	plateThickness*0.5),), 
	                        ((0.0, 					plateLength, 	plateThickness*0.5),),
							((plateWidth*0.5, 		plateLength, 	plateThickness*0.5),), 
	                        ((plateWidth, 			plateLength, 	plateThickness*0.5),),
							((plateWidth+0.0127, 	plateLength, 	plateThickness*0.5),))
	plate.Set(edges = plateThicknessEdges, name = 'PlateThickness')

	clampedWidthEdges = plateEdges.findAt(
	                        ((-0.0127*0.5,				0.0, 			0.0),),
							((-0.0127*0.5,				0.0, 			plateThickness),),
	                        ((plateWidth+(0.0127*0.5),	0.0, 			0.0),),
							((plateWidth+(0.0127*0.5),	0.0, 			plateThickness),),
	                        ((-0.0127*0.5,				plateLength, 	0.0),),
							((-0.0127*0.5,				plateLength, 	plateThickness),),
	                        ((plateWidth+(0.0127*0.5),	plateLength, 	0.0),),
							((plateWidth+(0.0127*0.5),	plateLength, 	plateThickness),))
	plate.Set(edges = clampedWidthEdges, name = 'ClampedWidth')

	# Creating the set of faces where the plate will be clamped
	plateFaces = plate.faces
	clampedFaces = plateFaces.findAt(
	                        ((-0.0127*0.5,				plateLength*0.5,	0.0),),
							((-0.0127*0.5,				plateLength*0.5,	plateThickness),),
	                        ((plateWidth+(0.0127*0.5),	plateLength*0.5,	0.0),), 
							((plateWidth+(0.0127*0.5),	plateLength*0.5,	plateThickness),))
	plate.Set(faces = clampedFaces, name = 'ClampedFaces')

	# Creating the set of vertices where the plate will be pinned
	if pinOrCombBC == "pin":
		plateVertices = plate.vertices
		pinVertices = plateVertices.findAt(
								((plateWidth*0.5,	plateLength,	0.0),),
								((plateWidth*0.5,	plateLength,	plateThickness),),
								((plateWidth*0.5,	0.0,			0.0),),
								((plateWidth*0.5,	0.0,			plateThickness),))
		plate.Set(vertices = pinVertices, name = 'Pins')

	# Creating the small and larage channel master surfaces
	if pinOrCombBC == "comb":
		smChMasterFaces = plateFaces.findAt( 
								((plateWidth*0.25,		plateLength*0.5,	plateThickness),), 
								((plateWidth*0.75,		plateLength*0.5,	plateThickness),))
		plate.Surface(side2Faces = smChMasterFaces, name = 'SmChMaster')

		lgChMasterFaces = plateFaces.findAt( 
								((plateWidth*0.25,		plateLength*0.5,	0.0),),
								((plateWidth*0.75,		plateLength*0.5,	0.0),))
		plate.Surface(side2Faces = lgChMasterFaces, name = 'LgChMaster')

	if guessedAbaqusStep == "yes":
		guessedPressureFaces = plateFaces.findAt(
							((plateWidth*0.25,		plateLength*0.5,	plateThickness),),
							((plateWidth*0.75,		plateLength*0.5,	plateThickness),))
		plate.Surface(side1Faces = guessedPressureFaces, name = 'Guessed_Pressure_Surface')

	##-----------------------------------------MATERIALS NODE-------------------------------------------------------
	# Creating the plate material aluminum
	mdb.models[modelName].Material(name=plateMaterial)
	mdb.models[modelName].materials[plateMaterial].Density(table=((plateDensity, ), ))
	mdb.models[modelName].materials[plateMaterial].Elastic(table=((elasticModulus, poissonsRatio), ))
	
	##-----------------------------------------SECTIONS NODE-------------------------------------------------------
	# Setting and creating the section assignment for the plate
	plateRegion = plate.sets['EntirePlateGeometry']
	plateCellsRegion = plate.sets['EntirePlateGeometry'].cells
	shellStack = plateFaces.findAt(
							((plateWidth*0.25,		plateLength*0.5,	plateThickness),))

	##-----------------------------------------MESH NODE-------------------------------------------------------
	if elemType == "C3D8I":
		element = C3D8I
		mdb.models[modelName].HomogeneousSolidSection(name='PlateSection_Solid', 
								  material=plateMaterial, thickness=None)
		plate.SectionAssignment(region=plateRegion, sectionName='PlateSection_Solid', offset=0.0, 
			offsetType=MIDDLE_SURFACE, offsetField='', 
			thicknessAssignment=FROM_SECTION)
		elemType = mesh.ElemType(elemCode = element, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
							distortionControl = DEFAULT)
		plate.setElementType(regions = (plateCellsRegion, ), elemTypes = (elemType, ))

	if elemType == "SC8R":
		element = SC8R
		mdb.models[modelName].HomogeneousShellSection(name='PlateSection_Shell', 
			preIntegrate=OFF, material=plateMaterial, thicknessType=UNIFORM, 
			thickness=plateThickness, thicknessField='', idealization=NO_IDEALIZATION, 
			poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
			useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)
		plate.SectionAssignment(region=plateRegion, sectionName='PlateSection_Shell', offset=0.0, 
			offsetType=MIDDLE_SURFACE, offsetField='', 
			thicknessAssignment=FROM_SECTION)
		elemType1 = mesh.ElemType(elemCode = element, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
							distortionControl = DEFAULT)
		elemType2 = mesh.ElemType(elemCode=SC6R, elemLibrary=STANDARD)
		elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)
		plate.setElementType(regions=(plateCellsRegion, ), elemTypes=(elemType1, elemType2, elemType3))
		plate.assignStackDirection(referenceRegion=shellStack[0], cells=plateCellsRegion)

	# Seeding edges of the plate for meshing
	plateLengthEdges = plate.sets['PlateLength'].edges
	plate.seedEdgeByBias(biasMethod=DOUBLE, endEdges=plateLengthEdges, ratio=plateLengthBias,
	    number=plateLengthNodes, constraint=FINER)

	plateWidthEdges = plate.sets['PlateWidth'].edges
	plate.seedEdgeByNumber(edges = plateWidthEdges, number = plateWidthNodes, constraint = FINER)

	plateThicknessEdges = plate.sets['PlateThickness'].edges
	plate.seedEdgeByNumber(edges = plateThicknessEdges, number = plateThickNodes, constraint = FINER)

	clampedWidthEdges = plate.sets['ClampedWidth'].edges
	plate.seedEdgeByNumber(edges = clampedWidthEdges, number = clampedWidthNodes, constraint = FINER)

	# Meshing the part
	plate.generateMesh()

	##------------------------------------------ASSEMBLY NODE-----------------------------------------------------
	# Adding the plates as instances to the assembly
	fsiInterfaceList = []
	plateSpacing = plateThickness + smChHeight
	assembly = mdb.models[modelName].rootAssembly
	for i in range(0, numOfPlates):
		assembly.Instance(name = plateName + "_" + str(i), part = plate, dependent = ON)
		plateFacesTmp = assembly.instances[plateName + '_' + str(i)].faces

		if i >= 1:
			assembly.translate(instanceList=(plateName + "_" + str(i), ), 
								vector=(0.0, 0.0, plateSpacing*i))
		
		tempList = ()
		w = 0.25
		for j in range(0, 2):
			tempList = (((plateWidth*w,		plateLength*0.5,	plateSpacing*i),),) + tempList 
			tempList = (((plateWidth*w,		plateLength*0.5,	plateThickness + plateSpacing*i),),) + tempList
			tempList = (((plateWidth*w,		0.0,				plateThickness*0.5 + plateSpacing*i),),) + tempList
			tempList = (((plateWidth*w,		plateLength,		plateThickness*0.5 + plateSpacing*i),),) + tempList
			w = 0.75
		fsiInterfaceList.append(plateFacesTmp.findAt(tempList[0], tempList[1], tempList[2], tempList[3], 
														tempList[4], tempList[5], tempList[6], tempList[7]))

		# Setting the clamped boundary condition on the clamped regions of the plate
		clampedRegion = assembly.instances[plateName + "_" + str(i)].sets['ClampedFaces']
		mdb.models[modelName].PinnedBC(name = 'Clamps_' + str(i), createStepName = 'Initial', 
								region = clampedRegion, localCsys = None)

		##-----------------------------------------BOUNDARY CONDITIONS -----------------------------------------------
		# Creating either a pinned or combed boundary condition
		if pinOrCombBC == 'pin':
			pinRegion = assembly.instances[plateName + "_" + str(i)].sets['Pins']
			mdb.models[modelName].DisplacementBC(name='Pins', 
				createStepName='Initial', region=pinRegion, u1=UNSET, u2=UNSET, u3=SET, 
				ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, 
				distributionType=UNIFORM, fieldName='', localCsys=None)

		elif pinOrCombBC == 'comb':
			combRadius = 0.086*0.0254*0.5
			combOffset = 0.003*0.0254
			combSize = [ combRadius, 0.01 ]
			combPositions = [ (plateWidth*0.5, plateLength - combRadius, plateThickness + combOffset),
								(plateWidth*0.5, plateLength - combRadius, -(combOffset + combSize[1])),
								(plateWidth*0.5, combRadius, plateThickness + combOffset),
								(plateWidth*0.5, combRadius, -(combOffset + combSize[1])) ]
			mdb.models[modelName].ContactProperty('IntProp')
			for i in range(0,4):
				# Creating the comb part
				comb = createCylinder(modelName,'Comb_' + str(i),combSize)
				combFaces = comb.faces
				combCells = comb.cells

				# Creating the master and slave surface for the interaction and the region set for the rigid body
				if i == 0 or i == 2:
					masterSurf = assembly.instances[plateName + "_" + str(i)].surfaces['SmChMaster']
					slaveSurface = combFaces.findAt( ((0.0, 0.0, 0.0),) )
					comb.Surface(side1Faces = slaveSurface, name = 'Comb_' + str(i) + '_Slave')
				elif i == 1 or i == 3:
					masterSurf = assembly.instances[plateName + "_" + str(i)].surfaces['LgChMaster']
					slaveSurface = combFaces.findAt( ((0.0, 0.0, combSize[1]),) )
					comb.Surface(side1Faces = slaveSurface, name = 'Comb_' + str(i) + '_Slave')

				# Creating the region for the comb
				entireComb = combCells.getByBoundingBox(
										  xMin = -1, yMin = -1, zMin = -1,
										  xMax = 1, yMax = 1, zMax = 1)
				comb.Set(cells = entireComb, name = 'CombCells_' + str(i))
				combRegionCells = comb.Set(cells = entireComb, name = 'CombCells_' + str(i)).cells

				# Meshing the comb
				planesPartition = [YZPLANE, XZPLANE]
				for j in range(0,2):
					comb.DatumPlaneByPrincipalPlane(principalPlane = planesPartition[j], offset = 0.0)
					datumPlanes = comb.datums
					pickedCells = combCells.findAt( ((0.0, combSize[0], combSize[1]*0.5),),
													((-combSize[0], 0.0, combSize[1]*0.5),) )
					comb.PartitionCellByDatumPlane(datumPlane = datumPlanes[(j*2)+5], cells = pickedCells)
				
				comb.seedPart(size = 2*math.pi*combRadius*0.25, deviationFactor = 0.1, minSizeFactor = 0.1)

				if elemType == 'SC8R':
					mdb.models[modelName].HomogeneousSolidSection(name='PlateSection_Solid', 
						material=plateMaterial, thickness=None)
				element = C3D8I
				sectionRegion = comb.sets['CombCells_' + str(i)]
				comb.SectionAssignment(region=sectionRegion, sectionName='PlateSection_Solid', offset=0.0, 
					offsetType=MIDDLE_SURFACE, offsetField='', 
					thicknessAssignment=FROM_SECTION)
				elemType = mesh.ElemType(elemCode = element, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
									distortionControl = DEFAULT)
				comb.setElementType(regions = (combRegionCells, ), elemTypes = (elemType, ))
				comb.generateMesh()

				# Adding the comb to the assembly
				assembly = mdb.models[modelName].rootAssembly
				assembly.Instance(name = 'Comb_' + str(i), part = comb, dependent = ON)

				# Translating the comb
				assembly.translate(instanceList=('Comb_' + str(i), ), vector=(combPositions[i]))

				# Creating the interaction between the plate and the comb
				combSlave = assembly.instances['Comb_' + str(i)].surfaces['Comb_' + str(i) + '_Slave']
				mdb.models[modelName].SurfaceToSurfaceContactStd(
					name='Comb_' + str(i) + 'Int', createStepName='Initial', master=masterSurf, 
					slave=combSlave, sliding=FINITE, thickness=ON, 
					interactionProperty='IntProp', adjustMethod=NONE, 
					initialClearance=OMIT, datumAxis=None, clearanceRegion=None)

			# Creating a set of all four combs
			combList = []
			for i in range(0,4):
				combCells0 = assembly.instances['Comb_' + str(i)].cells
				combCells_0 = combCells0.getByBoundingBox(
												xMin = -1, yMin = -1, zMin = -1,
												xMax = 1, yMax = 1, zMax = 1)
				combList.append(combCells_0)

			combRegions = assembly.Set(cells = combList[0]+combList[1]+combList[2]+combList[3], name = 'Combs')

		
			# Creating a reference point for the rigid body
			assembly.ReferencePoint(point = (0.0, 0.0, 0.0))
			refPts = assembly.referencePoints
			refPt = (refPts[12], )
			refPtRegion = regionToolset.Region(referencePoints=refPt)

			# Setting boundary condition for the reference point
			mdb.models[modelName].EncastreBC(name='RefPfRigid', createStepName='Initial', 
				region=refPtRegion, localCsys=None)

			# Creating the rigid body
			mdb.models[modelName].RigidBody(name='RigidComb', 
				refPointRegion=refPtRegion, bodyRegion=combRegions, refPointAtCOM=ON)

	# Creating the FSI interface for all plates in the stack
	assembly.Surface(side1Faces = fsiInterfaceList, name = 'FSI_INTERFACE')


	##-------------------------------------------JOBS NODE-------------------------------------------------------
	# Creating and writing the input file to the Abaqus work directory
	if pinOrCombBC == 'pin':
		BC = str(int(plateThickness/0.0254*1000)) + '_Pinned'

	elif pinOrCombBC == 'comb':
		BC = str(int(plateThickness/0.0254*1000)) + '_Combed'

	elif pinOrCombBC == 'none':
		BC = str(int(plateThickness/0.0254*1000)) + '_Free'

	# Creating the job for creating and appending the input file
	fileName = couplingScheme + '_' + BC + 'Plate'
	createInputFile(fileName, modelName)

	# Appending the written input file for explicit FSI coupling
	appendInputFile(couplingScheme, BC, timeStep, maxSimTime, minTimeStep)


'''
	This module creates the fluid geometry around a flat plate
'''
def createFlatFluid(parameters):
	modelName = parameters['abaqusModelName']
	plateName = 'Plate'
	pinOrCombBC = parameters['pinOrCombBC']
	numOfPlates = parameters['numOfPlates']

	# Grabbing all of the geometry variables for the fluid model
	plateLength = parameters['plateLength']
	plateWidth = parameters['plateWidth']
	plateThickness = parameters['plateThickness']
	inletPlLength = parameters['inletPlLength']
	outletPlLength = parameters['outletPlLength']
	smChHeight = parameters['smChHeight']
	lgChHeight = parameters['lgChHeight']

	# Grabbing all of the fluid mesh variables for the fluid model
	flPlLenNodes = parameters['flPlLenNodes']
	flPlLenBias = parameters['flPlLenBias']
	flPlWidthNodes = parameters['flPlWidthNodes']
	flPlWidthBias = parameters['flPlWidthBias']
	flInletNodes = parameters['flInletNodes']
	flInletBias = parameters['flInletBias']
	flOutletNodes = parameters['flOutletNodes']
	flOutletBias = parameters['flOutletBias']
	flSmChHeightNodes = parameters['flSmChHeightNodes']
	flSmChHeightBias = parameters['flSmChHeightBias']
	flLgChHeightNodes = parameters['flLgChHeightNodes']
	flLgChHeightBias = parameters['flLgChHeightBias']
	flPlHeightNodes = parameters['flPlHeightNodes']
	flPlHeightBias = parameters['flPlHeightBias']
	biasDirection = parameters['chBiasDirection']

	# Creating part named 'BulkFluid'
	smplateSpacing = plateThickness + smChHeight
	lgplateSpacing = plateThickness + lgChHeight
	bulkFluidThickness = (plateThickness + smChHeight + lgChHeight) + smplateSpacing*(numOfPlates - 1)
	size = [ (0.0, -outletPlLength), (plateWidth, plateLength + inletPlLength), (bulkFluidThickness) ]
	bulkFluid = createBox(modelName, 'BulkFluid', size)

	# Adding the bulk fluid as an instance to the assembly
	assembly = mdb.models[modelName].rootAssembly
	assembly.Instance(name = 'BulkFluid', part = bulkFluid, dependent = ON)

	# Suppressing the combs if they exist
	if pinOrCombBC == "comb":
		for i in range(0,4):
			assembly.features['Comb_' + str(i)].suppress()

	# Translating the bulkFluid part the height of the large channel 
	assembly.translate(instanceList=('BulkFluid', ), vector=(0.0, 0.0, -lgChHeight))

	# Cutting the plate instance from the bulkFluid instance to create the Fluid part
	plates = []
	for i in range(0, numOfPlates):
		plate = assembly.instances[plateName + "_" + str(i)]
		plates.append(plate)

	assembly.InstanceFromBooleanCut(name='Fluid', 
		instanceToBeCut=assembly.instances['BulkFluid'], 
		cuttingInstances=plates, 
		originalInstances=SUPPRESS)
	assembly.features.changeKey( fromName = 'Fluid-1', toName = 'Fluid' )

	# Creating the datum planes on the fluid part for partitioning
	fluid = mdb.models[modelName].parts['Fluid']
	fluid.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=plateLength)
	fluid.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
	for i in range(0, numOfPlates):
		fluid.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=smplateSpacing*i)
		fluid.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=(smplateSpacing*i) + plateThickness)

	# Grabbing the datum planes and the cells of the fluid for partitioning the fluid
	fluidCells = fluid.cells
	datumPlanes = fluid.datums

	# Creating a set of the entire plate geometry
	entireFluidGeometry = fluidCells.getByBoundingBox(
								xMin = -1, yMin = -1, zMin = -1,
								xMax = 1, yMax = 1, zMax = 1)
	fluid.Set(cells = entireFluidGeometry, name = 'EntireFluidGeometry')

	# Partitioning the fluid for meshing
	for i in range(2, (numOfPlates*2)+4):
		pickedCells = fluid.sets['EntireFluidGeometry'].cells
		fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[i], cells=pickedCells)

	# Creating all the sets of edges in the fluid part
	fluidEdges = fluid.edges
	pos = [ 'FluidPlateLength', plateLength*0.5, 
			'FluidOutletLength', -outletPlLength*0.5,
			'FluidInletLength', plateLength + (inletPlLength*0.5) ]
	for k in range(0, 3):
		flowPlateEdges = []
		for i in range(0, numOfPlates+1):
			w = 0.0
			for j in range(0, 2):
				flowPlateEdges.append(fluidEdges.findAt( 
										((w,	pos[(k*2)+1],	-lgChHeight + lgplateSpacing*i),) ))
				flowPlateEdges.append(fluidEdges.findAt( 
										((w,	pos[(k*2)+1],	smplateSpacing*i),) ))
				w = plateWidth
		fluid.Set(edges=flowPlateEdges, name=pos[k*2])


	zPos = ['LargeChHeight', -lgChHeight*0.5, 
			'SmallChHeight', plateThickness + smChHeight*0.5,
			'PlateHeight', plateThickness*0.5]
	yPos = [ -outletPlLength, 0.0, plateLength, plateLength + inletPlLength ] 
	for k in range(0, len(zPos)/2):
		flowPlateEdges = []
		w = 0.0
		for i in range(0, numOfPlates):
			for l in range(0,len(yPos)):
				w = 0.0
				for j in range(0, 2):
					if zPos[k*2] == 'PlateHeight':
						flowPlateEdges.append(fluidEdges.findAt( 
												((w, yPos[l], zPos[(k*2)+1] + smplateSpacing*i),) ))
					else:
						flowPlateEdges.append(fluidEdges.findAt( 
												((w, yPos[l], zPos[(k*2)+1] + smplateSpacing*i*2),) ))
					w = plateWidth
		fluid.Set(edges=flowPlateEdges, name=zPos[k*2])

	halfPlateWidth = plateWidth*0.5
	yPos = [ -outletPlLength,	0.0,			plateLength,	plateLength + inletPlLength ]
	zPos = [ 0.0,				plateThickness,	-lgChHeight,	smChHeight + plateThickness ]
	flowPlateEdges = []
	for i in range(0, numOfPlates+1):
		for j in range(0, len(yPos)):
			flowPlateEdges.append(fluidEdges.findAt( 
									((halfPlateWidth, yPos[j], -lgChHeight + lgplateSpacing*i),) ))
			flowPlateEdges.append(fluidEdges.findAt( 
									((halfPlateWidth, yPos[j], smplateSpacing*i),) ))
		fluid.Set(edges=flowPlateEdges, name='FluidWidth')

	# Creating the set of faces for the inlet
	pos = [ 'Inlet', plateLength + inletPlLength, 
			'Outlet', -outletPlLength ]
	fluidFaces = fluid.faces
	for i in range(0, len(pos)/2):
		inletFaces = []
		inletFaces.append(fluidFaces.findAt( ((plateWidth*0.5, pos[(i*2)+1], -lgChHeight*0.5),) ))
		for j in range(0, numOfPlates+1):
			inletFaces.append(fluidFaces.findAt( 
								((halfPlateWidth, pos[(i*2)+1], plateThickness*0.5 + smplateSpacing*j),) ))
			inletFaces.append(fluidFaces.findAt( 
								((halfPlateWidth, pos[(i*2)+1], plateThickness + lgChHeight*0.5 + smplateSpacing*j),) ))
		fluid.Surface(side2Faces = inletFaces, name = pos[(i*2)])

	# Creating the FSI surfaces
	yPos = [ plateLength*0.5, plateLength*0.5, plateLength, 0.0 ]
	zPos = [ 0.0, plateThickness, plateThickness*0.5, plateThickness*0.5 ]
	fsiSurfNames = ['FSI_Back', 'FSI_Front', 'FSI_Top', 'FSI_Bottom']
	for i in range(0, numOfPlates):
		for j in range(0, len(yPos)):
			fsiBack = fluidFaces.findAt( 
							   ((plateWidth*0.5, yPos[j], zPos[j] + smplateSpacing*i),) )
			fluid.Surface(side1Faces = fsiBack, name = fsiSurfNames[j] + '_' + str(i))

	# Seeding edges of the fluid for meshing
	fluidPlateLengthEdges = fluid.sets['FluidPlateLength'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidPlateLengthEdges, ratio=flPlLenBias,
        number=flPlLenNodes, constraint=FINER)

	fluidWidthEdges = fluid.sets['FluidWidth'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidWidthEdges, ratio=flPlWidthBias,
        number=flPlWidthNodes, constraint=FINER)

	fluidOutletEdges = fluid.sets['FluidOutletLength'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidOutletEdges, ratio=flOutletBias,
        number=flOutletNodes, constraint=FINER)

	fluidInletEdges = fluid.sets['FluidInletLength'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidInletEdges, ratio=flInletBias,
        number=flInletNodes, constraint=FINER)

	fluidSmallChEdges = fluid.sets['SmallChHeight'].edges
	fluidLargeChEdges = fluid.sets['LargeChHeight'].edges
	fluidPlateEdges = fluid.sets['PlateHeight'].edges
	if biasDirection == 'Center':
		fluid.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=fluidSmallChEdges, ratio=flSmChHeightBias,
    		number=flSmChHeightNodes, constraint=FINER)

		fluid.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=fluidLargeChEdges, ratio=flLgChHeightBias,
    		number=flLgChHeightNodes, constraint=FINER)

		fluid.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=fluidPlateEdges, ratio=flPlHeightBias,
			number=flPlHeightNodes, constraint=FINER)
	elif biasDirection == 'End':
		fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidSmallChEdges, ratio=flSmChHeightBias,
    		number=flSmChHeightNodes, constraint=FINER)

		fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidLargeChEdges, ratio=flLgChHeightBias,
    		number=flLgChHeightNodes, constraint=FINER)

		fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidPlateEdges, ratio=flPlHeightBias,
    		number=flPlHeightNodes, constraint=FINER)

	# Setting the mesh element type and meshing the part
	elemType = mesh.ElemType(elemCode = C3D8I, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
                              distortionControl = DEFAULT)
	fluidCellsRegion = fluid.sets['EntireFluidGeometry'].cells
	fluid.setElementType(regions = (fluidCellsRegion, ), elemTypes = (elemType, ))
	fluid.generateMesh()

	# Creating a job for input file creation
	if numOfPlates == 1:
		fileName = ( 'Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
				 '_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
				 '_' + str(int(parameters['lgChHeight']/0.0254*1000)) )
	else:
		fileName = ( 'Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
				'_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
				'_' + str(numOfPlates) + '_' + parameters['plateGeometry'] + '_Plate_Stack' )

	inputFileName = createInputFile(fileName,modelName)


def createCurvedPlate(geometryParameters):
	""" Creates input files for a curved plate model
	    Inputs:
	    model_name:    Abaqus model name
	    parameters:    model parameters
	    Outputs:
	"""
	# GRABBING ALL OF THE GEOMETRY PARAMETERS
	plateName = 'Plate'
	modelName = parameters['abaqusModelName']
	plateLength = parameters['plateLength']
	plateWidth = parameters['plateWidth']
	plateThickness = parameters['plateThickness']
	smChHeight = parameters['smChHeight']
	lgChHeight = parameters['lgChHeight']
	inletLength = parameters['inletPlLength']
	outletLength = parameters['outletPlLength']
	
	# Grabbing all of the mesh parameters
	plateThickNodes = parameters['plateThickNodes']
	plateLengthNodes = parameters['plateLengthNodes']
	plateWidthNodes = parameters['plateWidthNodes']
	clampedWidthNodes = parameters['clampedWidthNodes']

	# Grabbing all of the material properties of the plate
	numOfPlates = parameters['numOfPlates']
	plateMaterial = parameters['plateMaterial']
	elasticModulus = parameters['elasticModulus']
	poissonsRatio = parameters['poissonsRatio']
	plateDensity = parameters['plateDensity']
	pinOrCombBC = parameters['pinOrCombBC']

	# Grabbing the FSI coupling parameters
	couplingScheme = parameters['couplingScheme']
	timeStep = parameters['timeStep']
	minTimeStep = parameters['minTimeStep']
	maxSimTime = parameters['maxSimTime']

	plateSpacing = plateThickness + smChHeight
	mdb.Model(name=modelName, modelType = STANDARD_EXPLICIT)

	##------------------------------------------MATERIALS NODE-----------------------------------------------------
	# Creating the plate material aluminum
	mdb.models[modelName].Material(name=plateMaterial)
	mdb.models[modelName].materials[plateMaterial].Density(table=((plateDensity, ), ))
	mdb.models[modelName].materials[plateMaterial].Elastic(table=((elasticModulus, poissonsRatio), ))

	##------------------------------------------SECTIONS NODE-----------------------------------------------------
	# Setting and creating the section assignment for the plate
	mdb.models[modelName].HomogeneousSolidSection(name='PlateSection', 
								material=plateMaterial, thickness=None)

	##--------------------------------------------PARTS NODE-------------------------------------------------------
	# Building the curved plate stack
	# Creating a array of the plate's radii
	r_in = []
	r_out = []
	plateArcLen = []
	r_in.append((4*plateWidth/math.pi)-plateThickness*0.5)
	r_out.append(r_in[0] + plateThickness)
	plateArcLen.append(plateWidth)
	wettedTheta = math.pi/4
	plateTheta = (plateWidth + 0.0254)/r_in[0]
	clampedTheta = 0.5*(plateTheta - wettedTheta)
	emptyTheta = (math.pi/2 - plateTheta)/2

	# Sines and cosines for the clamped, wetted, and axial midpoint edges and the midpoints between the three edges
	cos_0C = math.cos(math.pi/2 - emptyTheta)
	cos_0_C = math.cos(math.pi/2 - emptyTheta - clampedTheta/2)
	cos_0 = math.cos(math.pi/2 - clampedTheta - emptyTheta)
	cos_mid_0 = math.cos(emptyTheta + clampedTheta + wettedTheta*0.75)
	cos_mid = math.cos(math.pi/4)
	cos_mid_1 = math.cos(emptyTheta + clampedTheta + wettedTheta*0.25)
	cos_1 = math.cos(emptyTheta + clampedTheta)
	cos_1_C = math.cos(emptyTheta + clampedTheta/2)
	cos_1C = math.cos(emptyTheta)

	sin_0C = math.sin(math.pi/2 - emptyTheta)
	sin_0_C = math.sin(math.pi/2 - emptyTheta - clampedTheta/2)
	sin_0 = math.sin(math.pi/2 - clampedTheta - emptyTheta)
	sin_mid_0 = math.sin(emptyTheta + clampedTheta + wettedTheta*0.75)
	sin_mid = math.sin(math.pi/4)
	sin_mid_1 = math.sin(emptyTheta + clampedTheta + wettedTheta*0.25)
	sin_1 = math.sin(emptyTheta + clampedTheta)
	sin_1_C = math.sin(emptyTheta + clampedTheta/2)
	sin_1C = math.sin(emptyTheta)

	halfPlateLength = plateLength*0.5
	fsiInterfaceList = []
	for i in range(0, numOfPlates):
		if i > 0:
			r_in.append(r_in[i-1] - plateSpacing)
			r_out.append(r_out[i-1] - plateSpacing)
			plateArcLen.append((r_in[i] + r_out[i])/2*plateTheta)
		plate = mdb.models[modelName].Part(name=plateName + '_' + str(i), dimensionality=THREE_D, type=DEFORMABLE_BODY )
		plateSketch = curved_sketch(plate, modelName, r_in[i], emptyTheta, plateThickness, 'Plate')

		plate.BaseSolidExtrude(sketch=plateSketch, depth=plateLength)
		plateSketch.unsetPrimaryObject()
		del mdb.models[modelName].sketches['__profile__']
		plate.regenerate
		mdb.models[modelName].parts[plateName + '_' + str(i)].setValues(geometryRefinement=EXTRA_FINE)
	
		# Creating datum points for the three datum planes
		plate.DatumPointByCoordinate(coords=(r_in[i]*cos_0,		plateLength, r_in[i]*sin_0))
		plate.DatumPointByCoordinate(coords=(r_in[i]*cos_1,		plateLength, r_in[i]*sin_1))
		plate.DatumPointByCoordinate(coords=(r_out[i]*cos_0,	plateLength, r_out[i]*sin_0))
		plate.DatumPointByCoordinate(coords=(r_out[i]*cos_1,	plateLength, r_out[i]*sin_1))

		plate.DatumPointByCoordinate(coords=(r_in[i]*cos_0, 0.0, r_in[i]*sin_0))
		plate.DatumPointByCoordinate(coords=(r_in[i]*cos_1, 0.0, r_in[i]*sin_1))
	
		plate.DatumPointByCoordinate(coords=(r_in[i]*cos_mid,	plateLength,	r_in[i]*sin_mid))
		plate.DatumPointByCoordinate(coords=(r_out[i]*cos_mid,	plateLength,	r_out[i]*sin_mid))
		plate.DatumPointByCoordinate(coords=(r_in[i]*cos_mid,	0.0,			r_in[i]*sin_mid))
		datumPoints = plate.datums

		# Creating three datum planes on the plate part
		plate.DatumPlaneByThreePoints(point1=datumPoints[10], point2=datumPoints[11], point3=datumPoints[12])
		plate.DatumPlaneByThreePoints(point1=datumPoints[4], point2=datumPoints[6], point3=datumPoints[8])
		plate.DatumPlaneByThreePoints(point1=datumPoints[5], point2=datumPoints[7], point3=datumPoints[9])

		# Grabbing the datum planes and the cells of the plate for partitioning the plate
		plateCells = plate.cells
		datumPlanes = plate.datums

		# Creating a set of the entire plate geometry
		entirePlateGeometry = plateCells.getByBoundingBox(
									xMin = -1, yMin = -1, zMin = -1,
									xMax = 1, yMax = 1, zMax = 1)
		plate.Set(cells = entirePlateGeometry, name = 'EntirePlateGeometry')

		# Partitioning the plate into three sections, wetted region and clamped regions
		pickedCells = plateCells.getByBoundingBox(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0)
		plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[15], cells=pickedCells)
			    
		pickedCells = plateCells.findAt(
				((r_in[i]*math.cos(plateTheta/2), halfPlateLength, (r_in[i] + (plateThickness*0.5))*math.sin(plateTheta/2),)))
		plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[13], cells=pickedCells)
			    
		pickedCells = plateCells.findAt(
				((r_in[i]*math.cos(plateTheta), plateLength*0.25, (r_in[i] + (plateThickness*0.5))*math.sin(plateTheta),)))
		plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[14], cells=pickedCells)
			    
		# Creating all the sets of edges in the plate part
		plateEdges = plate.edges
		flowLengthEdges = plateEdges.findAt(
			((r_in[i]*cos_0C,	halfPlateLength, 	r_in[i]*sin_0C),),
			((r_out[i]*cos_0C,	halfPlateLength,	r_out[i]*sin_0C),),
			((r_in[i]*cos_0,	halfPlateLength,	r_in[i]*sin_0),),
			((r_out[i]*cos_0,	halfPlateLength,	r_out[i]*sin_0),),
			((r_in[i]*cos_mid,	halfPlateLength,	r_in[i]*sin_mid),),
			((r_out[i]*cos_mid, halfPlateLength,	r_out[i]*sin_mid),),
			((r_in[i]*cos_1, 	halfPlateLength,	r_in[i]*sin_1),),
			((r_out[i]*cos_1, 	halfPlateLength,	r_out[i]*sin_1),),
			((r_in[i]*cos_1C, 	halfPlateLength,	r_in[i]*sin_1C),),
			((r_out[i]*cos_1C,	halfPlateLength,	r_out[i]*sin_1C),),)
		plate.Set(edges=flowLengthEdges, name='PlateLength')
			    
		flowWidthEdges = plateEdges.findAt(
			((r_in[i]*cos_mid_0, 	plateLength,	r_in[i]*sin_mid_0),),
			((r_out[i]*cos_mid_0, 	plateLength, 	r_out[i]*sin_mid_0),),
			((r_in[i]*cos_mid_1, 	plateLength, 	r_in[i]*sin_mid_1),),
			((r_out[i]*cos_mid_1, 	plateLength, 	r_out[i]*sin_mid_1),),
			((r_in[i]*cos_mid_0, 	0.0,			r_in[i]*sin_mid_0),),
			((r_out[i]*cos_mid_0, 	0.0, 			r_out[i]*sin_mid_0),),
			((r_in[i]*cos_mid_1, 	0.0, 			r_in[i]*sin_mid_1),),
			((r_out[i]*cos_mid_1, 	0.0, 			r_out[i]*sin_mid_1),))		
		plate.Set(edges=flowWidthEdges, name = 'PlateWidth')
			    
		r_plateThick = r_in[i] + plateThickness/2
		plateThicknessEdges = plateEdges.findAt(
			((r_plateThick*cos_0C,	0.0,			r_plateThick*sin_0C),),
			((r_plateThick*cos_0,	0.0,			r_plateThick*sin_0),),
			((r_plateThick*cos_mid,	0.0,			r_plateThick*sin_mid),),
			((r_plateThick*cos_1, 	0.0,			r_plateThick*sin_1),),
			((r_plateThick*cos_1C, 	0.0,			r_plateThick*sin_1C),),
			((r_plateThick*cos_0C,	plateLength,	r_plateThick*sin_0C),),
			((r_plateThick*cos_0,	plateLength,	r_plateThick*sin_0),),
			((r_plateThick*cos_mid,	plateLength,	r_plateThick*sin_mid),),
			((r_plateThick*cos_1, 	plateLength,	r_plateThick*sin_1),),
			((r_plateThick*cos_1C, 	plateLength,	r_plateThick*sin_1C),),)
		plate.Set(edges = plateThicknessEdges, name = 'PlateThickness')
			    
		clampedWidthEdges = plateEdges.findAt(
			((r_in[i]*cos_0_C, 	plateLength,	r_in[i]*sin_0_C),),
			((r_out[i]*cos_0_C, plateLength, 	r_out[i]*sin_0_C),),
			((r_in[i]*cos_1_C, 	plateLength, 	r_in[i]*sin_1_C),),
			((r_out[i]*cos_1_C,	plateLength, 	r_out[i]*sin_1_C),),
			((r_in[i]*cos_0_C, 	0.0,			r_in[i]*sin_0_C),),
			((r_out[i]*cos_0_C, 0.0, 			r_out[i]*sin_0_C),),
			((r_in[i]*cos_1_C, 	0.0, 			r_in[i]*sin_1_C),),
			((r_out[i]*cos_1_C,	0.0, 			r_out[i]*sin_1_C),))
		plate.Set(edges = clampedWidthEdges, name = 'ClampedWidth')
			    
		# Creating the set of faces where the plate will be clamped
		plateFaces = plate.faces
		clampedFaces = plateFaces.findAt(
			((r_in[i]*cos_0_C,	halfPlateLength,	r_in[i]*sin_0_C),),
			((r_out[i]*cos_0_C,	halfPlateLength,	r_out[i]*sin_0_C),),
			((r_in[i]*cos_1_C,	halfPlateLength,	r_in[i]*sin_1_C),),
			((r_out[i]*cos_1_C,	halfPlateLength,	r_out[i]*sin_1_C),))
		plate.Set(faces = clampedFaces, name = 'ClampedFaces')
			    
		# Creating the set of vertices where the plate will be pinned
		plateVertices = plate.vertices
		pinVertices = plateVertices.findAt(
			((r_in[i]*cos_mid,		plateLength,	r_in[i]*sin_mid),),
			((r_in[i]*cos_mid,		0.0,			r_in[i]*sin_mid),),)
		plate.Set(vertices = pinVertices, name = 'Pins')

		    
		# Setting the section and materal of the plate
		plateRegion = plate.sets['EntirePlateGeometry']
		plate.SectionAssignment(region=plateRegion, sectionName='PlateSection', offset=0.0, 
			offsetType=MIDDLE_SURFACE, offsetField='', 
			thicknessAssignment=FROM_SECTION)

		##--------------------------------------------MESH NODE-------------------------------------------------------
		# Seeding edges of the plate for meshing
		plateLengthEdges = plate.sets['PlateLength'].edges

		plate.seedEdgeByNumber(edges = plateLengthEdges, number = plateLengthNodes, constraint = FINER)
		plateWidthEdges = plate.sets['PlateWidth'].edges
		plate.seedEdgeByNumber(edges = plateWidthEdges, number = plateWidthNodes, constraint = FINER)

		plateThicknessEdges = plate.sets['PlateThickness'].edges
		plate.seedEdgeByNumber(edges = plateThicknessEdges, number = plateThickNodes, constraint = FINER)

		clampedWidthEdges = plate.sets['ClampedWidth'].edges
		plate.seedEdgeByNumber(edges = clampedWidthEdges, number = clampedWidthNodes, constraint = FINER)

		# Setting the mesh element type and meshing the part
		elemType = mesh.ElemType(elemCode = C3D8I, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
									distortionControl = DEFAULT)
		plateCellsRegion = plate.sets['EntirePlateGeometry'].cells
		plate.setElementType(regions = (plateCellsRegion, ), elemTypes = (elemType, ))
		plate.generateMesh()

		##-----------------------------------------ASSEMBLY NODE------------------------------------------------------
		# Adding the plate as an instance to the assembly
		assembly = mdb.models[modelName].rootAssembly
		assembly.Instance(name = plateName + '_' + str(i), part = plate, dependent = ON)

		# Creating the FSI_INTERFACE surface
		r_fsi = (r_in[i] + plateThickness/2)
		plateFacesTmp = assembly.instances['Plate_' + str(i)].faces
		tempList = ()
		tempList = (((r_in[i]*cos_mid_0,	halfPlateLength,	r_in[i]*sin_mid_0),),) + tempList 
		tempList = (((r_out[i]*cos_mid_0,	halfPlateLength, 	r_out[i]*sin_mid_0),),) + tempList
		tempList = (((r_fsi*cos_mid_0,		0.0, 				r_fsi*sin_mid_0),),) + tempList
		tempList = (((r_fsi*cos_mid_0,		plateLength, 		r_fsi*sin_mid_0),),) + tempList
		tempList = (((r_in[i]*cos_mid_1,	halfPlateLength,	r_in[i]*sin_mid_1),),) + tempList 
		tempList = (((r_out[i]*cos_mid_1,	halfPlateLength, 	r_out[i]*sin_mid_1),),) + tempList
		tempList = (((r_fsi*cos_mid_1,		0.0, 				r_fsi*sin_mid_1),),) + tempList
		tempList = (((r_fsi*cos_mid_1,		plateLength, 		r_fsi*sin_mid_1),),) + tempList

		fsiInterfaceList.append(plateFacesTmp.findAt(tempList[0], tempList[1], tempList[2], tempList[3], 
														tempList[4], tempList[5], tempList[6], tempList[7]))

		# Setting the clamped boundary condition on the clamped regions of the plate
		clampedRegion = assembly.instances[plateName + '_' + str(i)].sets['ClampedFaces']
		mdb.models[modelName].PinnedBC(name = 'Clamps', createStepName = 'Initial', 
									region = clampedRegion, localCsys = None)

		# Setting the pins boundary condition on the 'Pins' set
		if pinOrCombBC == 'pin':
			pinRegion = assembly.instances[plateName + '_' + str(i)].sets['Pins']
			mdb.models[modelName].DisplacementBC(name='Pins', 
    			createStepName='Initial', region=pinRegion, u1=SET, u2=UNSET, u3=SET, 
    			ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, 
    			distributionType=UNIFORM, fieldName='', localCsys=None)

	# Creating the FSI interface for all plates in the stack
	assembly.Surface(side1Faces = fsiInterfaceList, name = 'FSI_INTERFACE')

	##-------------------------------------------JOBS NODE-------------------------------------------------------
	# Creating and writing the input file to the Abaqus work directory
	if pinOrCombBC == 'pin':
		BC = str(int(plateThickness/0.0254*1000)) + '_Pinned'

	elif pinOrCombBC == 'comb':
		BC = str(int(plateThickness/0.0254*1000)) + '_Combed'

	elif pinOrCombBC == 'none':
		BC = str(int(plateThickness/0.0254*1000)) + '_Free'

	# Creating the job for creating and appending the input file
	fileName = couplingScheme + '_' + BC + 'Plate'
	createInputFile(fileName, modelName)

	# Appending the written input file for explicit FSI coupling
	appendInputFile(couplingScheme, BC, timeStep, maxSimTime, minTimeStep)

def createCurvedFluid(parameters):
	modelName = parameters['abaqusModelName']
	plateName = 'Plate'
	fluidName = 'Fluid'
	# Grabbing all of the geometry variables for the fluid model
	numOfPlates = parameters['numOfPlates']
	plateLength = parameters['plateLength']
	plateWidth = parameters['plateWidth']
	plateThickness = parameters['plateThickness']
	inletPlLength = parameters['inletPlLength']
	outletPlLength = parameters['outletPlLength']
	smChHeight = parameters['smChHeight']
	lgChHeight = parameters['lgChHeight']

	# Grabbing all of the fluid mesh variables for the fluid model
	flPlLenNodes = parameters['flPlLenNodes']
	flPlLenBias = parameters['flPlLenBias']
	flPlWidthNodes = parameters['flPlWidthNodes']
	flPlWidthBias = parameters['flPlWidthBias']
	flInletNodes = parameters['flInletNodes']
	flInletBias = parameters['flInletBias']
	flOutletNodes = parameters['flOutletNodes']
	flOutletBias = parameters['flOutletBias']
	flSmChHeightNodes = parameters['flSmChHeightNodes']
	flSmChHeightBias = parameters['flSmChHeightBias']
	flLgChHeightNodes = parameters['flLgChHeightNodes']
	flLgChHeightBias = parameters['flLgChHeightBias']
	flPlHeightNodes = parameters['flPlHeightNodes']
	flPlHeightBias = parameters['flPlHeightBias']
	biasDirection = parameters['chBiasDirection']

	##----------------------------------------PARTS/ASSEMBLY NODES----------------------------------------------------
	plateSpacing = smChHeight + plateThickness
	fluidThickness = plateSpacing*numOfPlates + smChHeight
	r_in = plateWidth/(math.pi/4) - plateThickness/2 - smChHeight - plateSpacing*(numOfPlates-1)
	r_out = r_in + fluidThickness
	plateTheta = (plateWidth + 0.0254)/r_in
	wettedTheta = math.pi/4
	clampedTheta = 0.5*(plateTheta - wettedTheta)
	emptyTheta = (math.pi/2 - plateTheta)/2
	cos_0 = math.cos(math.pi/2 - emptyTheta - clampedTheta)
	cos_1 = math.cos(emptyTheta + clampedTheta)
	sin_0 = math.sin(math.pi/2 - emptyTheta - clampedTheta)
	sin_1 = math.sin(emptyTheta + clampedTheta)
	cos_mid = math.cos(math.pi/4)
	sin_mid = math.sin(math.pi/4)
	    
	# Creating part named 'BulkFluid'
	bulkFluid = mdb.models[modelName].Part(name='BulkFluid', dimensionality=THREE_D, type=DEFORMABLE_BODY)
	bulkFluid = mdb.models[modelName].parts['BulkFluid']

	bulkFluid.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0)
	bulkFluid.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
	datumPlanes = bulkFluid.datums
	t = bulkFluid.MakeSketchTransform(sketchPlane=datumPlanes[1], origin=(0.0, 0.0, 0.0), 
								sketchOrientation=RIGHT, sketchPlaneSide=SIDE1, sketchUpEdge=datumPlanes[2])

	bulkFluidSketch = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=2.0, 
								transform=t)
	
	bulkFluidSketch = curved_sketch(bulkFluid, modelName, r_in, emptyTheta, fluidThickness, 'Fluid')  

	bulkFluid.BaseSolidExtrude(sketch=bulkFluidSketch, depth=(inletPlLength + plateLength + outletPlLength))
	bulkFluidSketch.unsetPrimaryObject()
	del mdb.models[modelName].sketches['__profile__']
	bulkFluid.regenerate
	mdb.models[modelName].parts['BulkFluid'].setValues(geometryRefinement=EXTRA_FINE)

	# Adding the bulk fluid as an instance to the assembly
	assembly = mdb.models[modelName].rootAssembly
	assembly.Instance(name = 'BulkFluid', part = bulkFluid, dependent = ON)

	# Translating the bulkFluid part the height of the large channel 
	assembly.translate(instanceList=('BulkFluid', ), vector=(0.0, -outletPlLength, 0.0))

	# Cutting the plate instance from the bulkFluid instance to create the Fluid part
	plates = []
	for i in range(0, numOfPlates):
		plates.append(assembly.instances[plateName + '_' + str(i)])

	assembly.InstanceFromBooleanCut(name=fluidName, 
		instanceToBeCut=assembly.instances['BulkFluid'], 
		cuttingInstances=(plates),
		originalInstances=SUPPRESS)
	assembly.features.changeKey( fromName = 'Fluid-1', toName = fluidName )
	mdb.models[modelName].parts[fluidName].setValues(geometryRefinement=EXTRA_FINE)
	    
	# Creating the datum points for creating the two datum planes
	fluid = mdb.models[modelName].parts[fluidName]
	fluid.DatumPointByCoordinate(coords=(r_in*cos_0,	plateLength, r_in*sin_0))
	fluid.DatumPointByCoordinate(coords=(r_in*cos_1, 	plateLength, r_in*sin_1))
	fluid.DatumPointByCoordinate(coords=(r_out*cos_0, 	plateLength, r_out*sin_0))
	fluid.DatumPointByCoordinate(coords=(r_out*cos_1,	plateLength, r_out*sin_1))
	    
	fluid.DatumPointByCoordinate(coords=(r_in*cos_0,	0.0, 		r_in*sin_0))
	fluid.DatumPointByCoordinate(coords=(r_in*cos_1, 	0.0, 		r_in*sin_1))
	fluid.DatumPointByCoordinate(coords=(r_out*cos_0, 	0.0, 		r_out*sin_0))
	fluid.DatumPointByCoordinate(coords=(r_out*cos_1,	0.0, 		r_out*sin_1))
	datumPoints = fluid.datums
	    
	# Creating three datum planes on the fluid part
	fluid.DatumPlaneByThreePoints(point1=datumPoints[2], point2=datumPoints[3], point3=datumPoints[4])
	fluid.DatumPlaneByThreePoints(point1=datumPoints[6], point2=datumPoints[7], point3=datumPoints[8])

	# Grabbing the datum planes and the cells of the fluid for partitioning the fluid
	fluidCells = fluid.cells
	datumPlanes = fluid.datums

	# Creating a set of the entire plate geometry
	entireFluidGeometry = fluidCells.getByBoundingBox(
								xMin = -1, yMin = -1, zMin = -1,
								xMax = 1, yMax = 1, zMax = 1)
	fluid.Set(cells = entireFluidGeometry, name = 'EntireFluidGeometry')

	# Partitioning the plate into three sections, wetted region and clamped regions
	pickedCells = fluid.sets['EntireFluidGeometry'].cells
	fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[10], cells=pickedCells)

	pickedCells = fluid.sets['EntireFluidGeometry'].cells
	fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[11], cells=pickedCells)
	    
	# Creating all the sets of edges in the fluid part
	halfPlateLength = plateLength*0.5
	r_i = []
	r_o = []
	flowPlateLengthEdges = []
	plateEdges = []
	smChannelEdges = []
	lgChannelEdges = []
	flowWidthEdges = []
	fluidEdges = fluid.edges
	fluidFaces = fluid.faces
	for i in range(0, numOfPlates+1):
		if i == 0:
			r_i.append(plateWidth/(math.pi/4) - plateThickness/2 - smChHeight - plateSpacing*(numOfPlates-1))
			r_o.append(r_i[i] + plateThickness + smChHeight + lgChHeight)
		else:
			r_i.append(r_i[i-1] + plateSpacing)
			r_o.append(r_o[i-1] + plateSpacing)

		r_plFlow = (r_i[i] + lgChHeight)
		flowPlateLengthEdges.append( fluidEdges.findAt(
			((r_i[i]*cos_0,		halfPlateLength,	r_i[i]*sin_0),),
			((r_i[i]*cos_1,		halfPlateLength,	r_i[i]*sin_1),),
			((r_plFlow*cos_0,	halfPlateLength,	r_plFlow*sin_0),),
			((r_plFlow*cos_1,	halfPlateLength,	r_plFlow*sin_1),),) )

		r_plateEdges = r_i[i] + plateThickness/2 + lgChHeight
		plateEdges.append( fluidEdges.findAt(
			((r_plateEdges*cos_0,	plateLength,	r_plateEdges*sin_0),),
			((r_plateEdges*cos_1,	plateLength,	r_plateEdges*sin_1),),
			((r_plateEdges*cos_0,	0.0,			r_plateEdges*sin_0),),
			((r_plateEdges*cos_1,	0.0,			r_plateEdges*sin_1),),) )

		if i % 2 == 0:
			r_lgChEdges = r_i[i] + lgChHeight/2 + plateThickness
			lgChannelEdges.append( fluidEdges.findAt(
				((r_lgChEdges*cos_0,	plateLength,	r_lgChEdges*sin_0),),
				((r_lgChEdges*cos_1,	plateLength,	r_lgChEdges*sin_1),),
				((r_lgChEdges*cos_0,	0.0,			r_lgChEdges*sin_0),),
				((r_lgChEdges*cos_1,	0.0,			r_lgChEdges*sin_1),),) )

		else:
			r_smChEdges = r_i[i] + smChHeight/2
			smChannelEdges.append( fluidEdges.findAt(
				((r_smChEdges*cos_0,	plateLength,	r_smChEdges*sin_0),),
				((r_smChEdges*cos_1,	plateLength,	r_smChEdges*sin_1),),
				((r_smChEdges*cos_0,	0.0,			r_smChEdges*sin_0),),
				((r_smChEdges*cos_1,	0.0,			r_smChEdges*sin_1),),) )

		r_flowWidth = r_i[i] + lgChHeight
		if i == 0:
			flowWidthEdges.append( fluidEdges.findAt( 
									((r_i[i]*cos_mid,		-outletPlLength,			r_i[i]*sin_mid),),
									((r_i[i]*cos_mid,		plateLength+inletPlLength,	r_i[i]*sin_mid),),
									((r_i[i]*cos_mid,		0.0,						r_i[i]*sin_mid),),
									((r_flowWidth*cos_mid,	0.0,						r_flowWidth*sin_mid),),
									((r_i[i]*cos_mid,		plateLength,				r_i[i]*sin_mid),),
									((r_flowWidth*cos_mid,	plateLength,				r_flowWidth*sin_mid),), ) )
		elif i == numOfPlates:
			flowWidthEdges.append( fluidEdges.findAt( 
									((r_flowWidth*cos_mid,	-outletPlLength,			r_flowWidth*sin_mid),),
									((r_flowWidth*cos_mid,	plateLength+inletPlLength,	r_flowWidth*sin_mid),),
									((r_i[i]*cos_mid,		0.0,						r_i[i]*sin_mid),),
									((r_flowWidth*cos_mid,	0.0,						r_flowWidth*sin_mid),),
									((r_i[i]*cos_mid,		plateLength,				r_i[i]*sin_mid),),
									((r_flowWidth*cos_mid,	plateLength,				r_flowWidth*sin_mid),), ) )
		else:
			flowWidthEdges.append( fluidEdges.findAt(
									((r_i[i]*cos_mid,		0.0,			r_i[i]*sin_mid),),
									((r_flowWidth*cos_mid,	0.0,			r_flowWidth*sin_mid),),
									((r_i[i]*cos_mid,		plateLength,	r_i[i]*sin_mid),),
									((r_flowWidth*cos_mid,	plateLength,	r_flowWidth*sin_mid),), ) )

		if i == numOfPlates:
			continue

		r_fsiBack = r_i[i] + lgChHeight
		fsiBack = fluidFaces.findAt(
									((r_fsiBack*cos_mid,	plateLength*0.5, r_fsiBack*cos_mid),),)
		fluid.Surface(side1Faces = fsiBack, name = 'FSI_Back_' + str(i))

		r_fsiFront = r_i[i] + lgChHeight + plateThickness
		fsiFront = fluidFaces.findAt(
									((r_fsiFront*cos_mid,	plateLength*0.5, r_fsiFront*cos_mid),),)
		fluid.Surface(side1Faces = fsiFront, name = 'FSI_Front_' + str(i))

		r_fsiTopBot = r_i[i] + lgChHeight + plateThickness/2
		fsiTop = fluidFaces.findAt(
									((r_fsiTopBot*cos_mid,	plateLength,	r_fsiTopBot*cos_mid),),)
		fluid.Surface(side1Faces = fsiTop, name = 'FSI_Top_' + str(i))

		fsiBot = fluidFaces.findAt(
									((r_fsiTopBot*cos_mid,	0.0,			r_fsiTopBot*cos_mid),),)
		fluid.Surface(side1Faces = fsiBot, name = 'FSI_Bottom_' + str(i))

	fluid.Set(edges=flowPlateLengthEdges, name='FluidPlateLength')
	fluid.Set(edges=lgChannelEdges, name='LargeChHeight')
	fluid.Set(edges=smChannelEdges, name='SmallChHeight')
	fluid.Set(edges=plateEdges, name='PlateHeight')
	fluid.Set(edges=flowWidthEdges, name='FluidWidth')
			   
	halfOutletLength = outletPlLength*0.5 
	flowOutletLengthEdges = fluidEdges.findAt(
		((r_in*cos_0,	-halfOutletLength,	r_in*sin_0),),
		((r_in*cos_1,	-halfOutletLength,	r_in*sin_1),),
		((r_out*cos_0,	-halfOutletLength,	r_out*sin_0),),
		((r_out*cos_1,	-halfOutletLength,	r_out*sin_1),),)
	fluid.Set(edges=flowOutletLengthEdges, name='FluidOutletLength')
	
	halfInletPlateLength = plateLength + inletPlLength*0.5
	flowInletLengthEdges = fluidEdges.findAt(
		((r_in*cos_0,	halfInletPlateLength,	r_in*sin_0),),
		((r_in*cos_1,	halfInletPlateLength,	r_in*sin_1),),
		((r_out*cos_0,	halfInletPlateLength,	r_out*sin_0),),
		((r_out*cos_1,	halfInletPlateLength,	r_out*sin_1),),)
	fluid.Set(edges=flowInletLengthEdges, name='FluidInletLength')

	r_inOut = r_in + plateThickness/2
	inOutEdges = fluidEdges.findAt(
		((r_inOut*cos_0,	plateLength + inletPlLength,	r_inOut*sin_0),),
		((r_inOut*cos_1,	plateLength + inletPlLength,	r_inOut*sin_1),),
		((r_inOut*cos_0,	-outletPlLength,				r_inOut*sin_0),),
		((r_inOut*cos_1,	-outletPlLength,				r_inOut*sin_1),),)
	fluid.Set(edges=inOutEdges, name='InletOutletHeight')
	        
	# Creating the set of faces for the inlet
	inletFaces = fluidFaces.findAt(
		((r_i[numOfPlates/2]*cos_mid,	plateLength + inletPlLength, r_i[numOfPlates/2]*sin_mid),),)
	fluid.Surface(side1Faces = inletFaces, name = 'Inlet')
	    
	outletFaces = fluidFaces.findAt(
		((r_i[numOfPlates/2]*cos_mid,	-outletPlLength,			r_i[numOfPlates/2]*sin_mid),),)
	fluid.Surface(side1Faces = outletFaces, name = 'Outlet')

	##--------------------------------------------MESH NODE-------------------------------------------------------
	# Seeding edges of the fluid for meshing
	fluidPlateLengthEdges = fluid.sets['FluidPlateLength'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidPlateLengthEdges, ratio=flPlLenBias,
		number=flPlLenNodes, constraint=FINER)

	fluidWidthEdges = fluid.sets['FluidWidth'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidWidthEdges, ratio=flPlWidthBias,
		number=flPlWidthNodes, constraint=FINER)

	fluidOutletEdges = fluid.sets['FluidOutletLength'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidOutletEdges, ratio=flOutletBias,
		number=flOutletNodes, constraint=FINER)

	fluidInletEdges = fluid.sets['FluidInletLength'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidInletEdges, ratio=flInletBias,
		number=flInletNodes, constraint=FINER)

	fluidSmallChEdges = fluid.sets['SmallChHeight'].edges
	if biasDirection == 'Center':
		fluid.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=fluidSmallChEdges, ratio=flSmChHeightBias,
    		number=flSmChHeightNodes, constraint=FINER)
	elif biasDirection == 'End':
		fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidSmallChEdges, ratio=flSmChHeightBias,
    		number=flSmChHeightNodes, constraint=FINER)

	fluidLargeChEdges = fluid.sets['LargeChHeight'].edges
	if biasDirection == 'Center':
		fluid.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=fluidLargeChEdges, ratio=flLgChHeightBias,
    		number=flLgChHeightNodes, constraint=FINER)
	elif biasDirection == 'End':
		fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidLargeChEdges, ratio=flLgChHeightBias,
    		number=flLgChHeightNodes, constraint=FINER)

	fluidPlateEdges = fluid.sets['PlateHeight'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidPlateEdges, ratio=flPlHeightBias,
		number=flPlHeightNodes, constraint=FINER)

	flPlHeightNodes = flPlHeightNodes*numOfPlates
	if numOfPlates+1 % 2 == 0:
		flLgChHeightNodes = flLgChHeightNodes*(numOfPlates+1)/2
		flSmChHeightNodes = flSmChHeightNodes*(numOfPlates+1)/2
	else:
		flLgChHeightNodes = int(flLgChHeightNodes*math.ceil((numOfPlates+1)/2))
		flSmChHeightNodes = int(flSmChHeightNodes*math.floor((numOfPlates+1)/2))
	
	fluidInletOutletEdges = fluid.sets['InletOutletHeight'].edges
	fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidInletOutletEdges, ratio=flPlHeightBias,
		number=flPlHeightNodes + flLgChHeightNodes + flSmChHeightNodes, constraint=FINER)

	# Setting the mesh element type and meshing the part
	elemType = mesh.ElemType(elemCode = C3D8I, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
								distortionControl = DEFAULT)
	fluidCellsRegion = fluid.sets['EntireFluidGeometry'].cells
	fluid.setElementType(regions = (fluidCellsRegion, ), elemTypes = (elemType, ))
	fluid.generateMesh()

	##--------------------------------------------JOBS NODE-------------------------------------------------------
	# Creating a job for input file creation
	if numOfPlates == 1:
		fileName = ( 'Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
				 '_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
				 '_' + str(int(parameters['lgChHeight']/0.0254*1000)) )
	else:
		fileName = ( 'Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
				'_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
				'_' + str(numOfPlates) + '_' + parameters['plateGeometry'] + '_Plate_Stack' )

	inputFileName = createInputFile(fileName,modelName)

def createBox(modelName, boxName, size):
    pt1 = size[0]
    pt2 = size[1]
    extDepth = size[2]

    box = mdb.models[modelName].Part(name=boxName, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    boxSketch = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=2.0)
    boxSketch.setPrimaryObject(option=STANDALONE)

    boxSketch.rectangle(point1=pt1, point2=pt2)

    box.BaseSolidExtrude(sketch=boxSketch, depth=extDepth)
    boxSketch.unsetPrimaryObject()
    del mdb.models[modelName].sketches['__profile__']
    return box

def createCylinder(modelName,cylinderName,size):
	r = size[0]
	h = size[1]

	cylinder = mdb.models[modelName].Part(name=cylinderName, dimensionality=THREE_D, type=DEFORMABLE_BODY)

	cylSketch = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=1.0)
	cylSketch.setPrimaryObject(option=STANDALONE)

	cylSketch.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, r))

	cylinder.BaseSolidExtrude(sketch=cylSketch, depth=h)
	cylSketch.unsetPrimaryObject()
	del mdb.models[modelName].sketches['__profile__']
	return cylinder

def createInputFile(jobName, modelName):
	mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS, 
		atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
		memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
		explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
		modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
		scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
		activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
	mdb.jobs[jobName].writeInput(consistencyChecking=OFF)

def appendInputFile(couplingScheme, BC, timeStep, maxSimTime, minTimeStep):
	# Appending the written input file for explicit FSI coupling
	with open(couplingScheme + "_" + BC + "Plate.inp", "a") as inputFile:
		if couplingScheme == 'Explicit':
			couplingScheme = 'GAUSS-SEIDEL'
		elif couplingScheme == 'Implicit':
			couplingScheme = 'ITERATIVE'

		inputFileLines = [
    		"** ----------------------------------------------------------------\n",
    		"**\n",
    		"** STEP: FSI\n",
    		"**\n",
    		"*Step, name=FSI, nlgeom=YES, inc=1000000\n",
    		"*Dynamic,application=QUASI-STATIC \n",
    		"" + str(timeStep) + "," + str(maxSimTime) + "," + str(minTimeStep) + "," + str(timeStep) + "\n",
    		"**\n",
			"** OUTPUT REQUESTS\n",
			"**\n",
    		"** FIELD OUTPUT: F-Output-1\n",
    		"**\n",
    		"*Output, field, variable=PRESELECT, frequency=1\n",
    		"**\n",
    		"** HISTORY OUTPUT: H-Output-1\n",
    		"**\n",
    		"*Output, history, variable=PRESELECT\n",
    		"**\n",
    		"*CO-SIMULATION, NAME=FSI_Mech, PROGRAM=MULTIPHYSICS, CONTROLS=Control-1\n",
    		"*CO-SIMULATION REGION, TYPE=SURFACE, EXPORT\n",
    		"ASSEMBLY_FSI_INTERFACE, U\n",
    		"ASSEMBLY_FSI_INTERFACE, V\n",
    		"*CO-SIMULATION REGION, TYPE=SURFACE, IMPORT\n",
    		"ASSEMBLY_FSI_INTERFACE, CF\n",
    		"*CO-SIMULATION CONTROLS, NAME=Control-1, COUPLING SCHEME=" + couplingScheme +", SCHEME MODIFIER=LAG,\n",
    		"STEP SIZE=" + str(timeStep) + ", TIME INCREMENTATION=SUBCYCLE, TIME MARKS=YES\n",
    		"**\n",
    		"*End Step\n"]			
		inputFile.writelines(inputFileLines)
	inputFile.close()

def curved_sketch(part, modelName, r_in, emptyTheta, plateThickness, plateOrFluid):
	part.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0)
	part.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
	datumPlanes = part.datums
	t = part.MakeSketchTransform(sketchPlane=datumPlanes[1], origin=(0.0, 0.0, 0.0), 
								sketchOrientation=RIGHT, sketchPlaneSide=SIDE1, sketchUpEdge=datumPlanes[2])
	    
	plateSketch = mdb.models[modelName].ConstrainedSketch(
								name='__profile__', sheetSize=3.0*r_in*(math.pi/2 - emptyTheta), transform=t)
							 
	#    Inner radius:
	ri = r_in
	#    Outer radius:
	ro = ri + plateThickness

	#    Point 0: Center of radius
	point0 = (0.0, 0.0)
	if plateOrFluid == 'Plate':
		#    Point 1: Inside, left
		point1 = (-ri*math.sin(math.pi/2 - emptyTheta), ri*math.cos(math.pi/2 - emptyTheta))
		#    Point 2: Inside, right
		point2 = (-ri*math.sin(emptyTheta), ri*math.cos(emptyTheta))
		#    Point 3: Outside, left
		point3 = (-ro*math.sin(math.pi/2 - emptyTheta), ro*math.cos(math.pi/2 - emptyTheta))
		#    Point 4: Outside, right
		point4 = (-ro*math.sin(emptyTheta), ro*math.cos(emptyTheta))
	elif plateOrFluid == 'Fluid':
		plateTheta = math.pi/2 - 2*emptyTheta
		clampedTheta = 0.5*(plateTheta - math.pi/4)
		#    Point 1: Inside, left
		point1 = (-ri*sin(math.pi/2 - emptyTheta - clampedTheta), ri*cos(math.pi/2 - emptyTheta - clampedTheta))
		#    Point 2: Inside, right
		point2 = (-ri*sin(emptyTheta + clampedTheta), ri*cos(emptyTheta + clampedTheta))
		#    Point 3: Outside, left
		point3 = (-ro*sin(math.pi/2 - emptyTheta - clampedTheta), ro*cos(math.pi/2 - emptyTheta - clampedTheta))
		#    Point 4: Outside, right
		point4 = (-ro*sin(emptyTheta + clampedTheta), ro*cos(emptyTheta + clampedTheta))
	#    Draw curved profile
	g, v, d, c = plateSketch.geometry, plateSketch.vertices, plateSketch.dimensions, plateSketch.constraints
	plateSketch.setPrimaryObject(option=STANDALONE )
	plateSketch.ArcByCenterEnds(center=point0, point1=point1, point2=point2, direction=CLOCKWISE)
	plateSketch.ArcByCenterEnds(center=point0, point1=point3, point2=point4, direction=CLOCKWISE)
	plateSketch.Line(point1=point1, point2=point3)
	plateSketch.Line(point1=point2, point2=point4)
	return plateSketch

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
print('\n\n Script is running \n')

# Reading the .txt file for the FSI model's parameters
parameters = fileReader('FSI_Input_File.txt')
print('\n\n The plates geometry and mesh parameters have been read from external file')

# Running the plate and fluid subroutines for creating the plate and fluid models
if parameters['plateGeometry'] == 'Flat':
	createFlatPlate(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		if parameters['numOfPlates'] == 1:
			mdb.saveAs('FSI_SolidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' +
				    str(int(parameters['lgChHeight']/0.0254*1000)))
		else:
			mdb.saveAs('FSI_SolidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + 
					str(parameters['numOfPlates']) + '_' + parameters['plateGeometry'] + '_' + 'Plate_Stack')
	print('\n\n The Abaqus script has completed building the plate model\n')

	createFlatFluid(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		if parameters['numOfPlates'] == 1:
			mdb.saveAs('FSI_FluidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + 
					str(int(parameters['lgChHeight']/0.0254*1000)))
		else:
			mdb.saveAs('FSI_FluidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + 
					str(parameters['numOfPlates']) + '_' + parameters['plateGeometry'] + '_' + 'Plate_Stack')
	print('\n\n The Abaqus script has completed building the fluid model\n')

if parameters['plateGeometry'] == 'Curved':
	createCurvedPlate(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		if parameters['numOfPlates'] == 1:
			mdb.saveAs('FSI_SolidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' +
				    str(int(parameters['lgChHeight']/0.0254*1000)))
		else:
			mdb.saveAs('FSI_SolidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + 
					str(parameters['numOfPlates']) + '_' + parameters['plateGeometry'] + '_' + 'Plate_Stack')

	print('\n\n The Abaqus script has completed building the plate model\n')

	createCurvedFluid(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		if parameters['numOfPlates'] == 1:
			mdb.saveAs('FSI_FluidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + 
					str(int(parameters['lgChHeight']/0.0254*1000)))
		else:
			mdb.saveAs('FSI_FluidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
					'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + 
					str(parameters['numOfPlates']) + '_' + parameters['plateGeometry'] + '_' + 'Plate_Stack')

	print('\n\n The Abaqus script has completed building the fluid model\n')
	