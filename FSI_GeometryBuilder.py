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
	plateVertices = plate.vertices
	pinVertices = plateVertices.findAt(
							((plateWidth*0.5,	plateLength,	0.0),),
							((plateWidth*0.5,	plateLength,	plateThickness),),
							((plateWidth*0.5,	0.0,			0.0),),
							((plateWidth*0.5,	0.0,			plateThickness),))
	plate.Set(vertices = pinVertices, name = 'Pins')
	
	# Creating the FSI_INTERFACE surface
	fsiFaces = plateFaces.findAt(
	                        ((plateWidth*0.25,		plateLength*0.5,	0.0),), 
							((plateWidth*0.25,		plateLength*0.5,	plateThickness),),
	                        ((plateWidth*0.25,		0.0,				plateThickness*0.5),),
							((plateWidth*0.25,		plateLength,		plateThickness*0.5),),
	                        ((plateWidth*0.75,		plateLength*0.5,	0.0),),
							((plateWidth*0.75,		plateLength*0.5,	plateThickness),),
	                        ((plateWidth*0.75,		0.0,				plateThickness*0.5),),
							((plateWidth*0.75,		plateLength,		plateThickness*0.5),))
	plate.Surface(side1Faces = fsiFaces, name = 'FSI_INTERFACE')

	# Creating the small channel master surface
	smChMasterFaces = plateFaces.findAt( 
							((plateWidth*0.25,		plateLength*0.5,	plateThickness),), 
							((plateWidth*0.75,		plateLength*0.5,	plateThickness),))
	plate.Surface(side2Faces = smChMasterFaces, name = 'SmChMaster')

	# Creating the large channel master surface
	lgChMasterFaces = plateFaces.findAt( 
	                        ((plateWidth*0.25,		plateLength*0.5,	0.0),),
							((plateWidth*0.75,		plateLength*0.5,	0.0),))
	plate.Surface(side2Faces = lgChMasterFaces, name = 'LgChMaster')

	if guessedAbaqusStep == "yes":
		guessedPressureFaces = plateFaces.findAt(
							((plateWidth*0.25,		plateLength*0.5,	plateThickness),),
							((plateWidth*0.75,		plateLength*0.5,	plateThickness),))
		plate.Surface(side1Faces = guessedPressureFaces, name = 'Guessed_Pressure_Surface')

	# Creating the plate material aluminum
	mdb.models[modelName].Material(name=plateMaterial)
	mdb.models[modelName].materials[plateMaterial].Density(table=((plateDensity, ), ))
	mdb.models[modelName].materials[plateMaterial].Elastic(table=((elasticModulus, poissonsRatio), ))

	# Setting and creating the section assignment for the plate
	plateRegion = plate.sets['EntirePlateGeometry']
	plateCellsRegion = plate.sets['EntirePlateGeometry'].cells
	shellStack = plateFaces.findAt(
							((plateWidth*0.25,		plateLength*0.5,	plateThickness),))
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

	# Adding the plate as an instance to the assembly
	assembly = mdb.models[modelName].rootAssembly
	assembly.Instance(name = plateName, part = plate, dependent = ON)

	# Setting the clamped boundary condition on the clamped regions of the plate
	clampedRegion = assembly.instances[plateName].sets['ClampedFaces']
	mdb.models[modelName].PinnedBC(name = 'Clamps', createStepName = 'Initial', 
                              region = clampedRegion, localCsys = None)

	# Creating either a pinned or combed boundary condition
	if pinOrCombBC == 'pin':
		pinRegion = assembly.instances[plateName].sets['Pins']
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
				masterSurf = assembly.instances[plateName].surfaces['SmChMaster']
				slaveSurface = combFaces.findAt( ((0.0, 0.0, 0.0),) )
				comb.Surface(side1Faces = slaveSurface, name = 'Comb_' + str(i) + '_Slave')
			elif i == 1 or i == 3:
				masterSurf = assembly.instances[plateName].surfaces['LgChMaster']
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


	# Creating and writing the input file to the Abaqus work directory
	if pinOrCombBC == 'pin':
		BC = str(int(plateThickness/0.0254*1000)) + '_Pinned'

	elif pinOrCombBC == 'comb':
		BC = str(int(plateThickness/0.0254*1000)) + '_Combed'

	elif pinOrCombBC == 'no':
		BC = str(int(plateThickness/0.0254*1000)) + '_Free'

	# Creating the job for creating and appending the input file
	fileName = couplingScheme + '_' + BC + 'Plate'
	createInputFile(fileName, modelName)

	# Appending the written input file for explicit FSI coupling
	appendInputFile(couplingScheme, BC, timeStep, maxSimTime, minTimeStep)

def createFlatFluid(parameters):
	modelName = parameters['abaqusModelName']
	plateName = 'Plate'
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
	size = [ (0.0, -outletPlLength), 
			(plateWidth, plateLength + inletPlLength), 
			(plateThickness + smChHeight + lgChHeight) ]
	bulkFluid = createBox(modelName, 'BulkFluid', size)
	bulkFluid = mdb.models[modelName].Part(name='BulkFluid', dimensionality=THREE_D, type=DEFORMABLE_BODY)
	bulkFluid = mdb.models[modelName].parts['BulkFluid']

	plateSketch = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=2.0)
	plateSketch.setPrimaryObject(option=STANDALONE)

	plateSketch.rectangle(point1=(0.0, -outletPlLength), point2=(plateWidth, plateLength + inletPlLength))

	bulkFluid.BaseSolidExtrude(sketch=plateSketch, depth=plateThickness + smChHeight + lgChHeight)
	plateSketch.unsetPrimaryObject()

	# Adding the bulk fluid as an instance to the assembly
	assembly = mdb.models[modelName].rootAssembly
	assembly.Instance(name = 'BulkFluid', part = bulkFluid, dependent = ON)

	# Suppressing the combs if they exist
	if pinOrCombBC.equals("comb"):
		for i in range(0,4):
			assembly.features['Comb_' + str(i)].suppress()

	# Translating the bulkFluid part the height of the large channel 
	assembly.translate(instanceList=('BulkFluid', ), vector=(0.0, 0.0, -lgChHeight))

	# Cutting the plate instance from the bulkFluid instance to create the Fluid part
	assembly.InstanceFromBooleanCut(name='Fluid', 
		instanceToBeCut=assembly.instances['BulkFluid'], 
		cuttingInstances=(assembly.instances[plateName], ), 
		originalInstances=SUPPRESS)
	assembly.features.changeKey( fromName = 'Fluid-1', toName = 'Fluid' )

	# Creating three datum planes on the fluid part
	fluid = mdb.models[modelName].parts['Fluid']
	fluid.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
	fluid.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=plateLength)
	fluid.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
	fluid.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=plateThickness)

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
	fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[2], cells=pickedCells)

	pickedCells = fluid.sets['EntireFluidGeometry'].cells
	fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[3], cells=pickedCells)

	pickedCells = fluid.sets['EntireFluidGeometry'].cells
	fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[4], cells=pickedCells)

	pickedCells = fluid.sets['EntireFluidGeometry'].cells
	fluid.PartitionCellByDatumPlane(datumPlane=datumPlanes[5], cells=pickedCells)

	# Creating all the sets of edges in the fluid part
	fluidEdges = fluid.edges
	flowPlateLengthEdges = fluidEdges.findAt(
                            ((0.0, 			plateLength*0.5,	0.0),),
                            ((0.0, 			plateLength*0.5,	plateThickness),),
                            ((0.0, 			plateLength*0.5,	-lgChHeight),),
                            ((0.0, 			plateLength*0.5,	smChHeight + plateThickness),),
                            ((plateWidth, 	plateLength*0.5,	0.0),),
                            ((plateWidth, 	plateLength*0.5,	plateThickness),),
                            ((plateWidth, 	plateLength*0.5,	-lgChHeight),),
                            ((plateWidth, 	plateLength*0.5,	smChHeight + plateThickness),))
	fluid.Set(edges=flowPlateLengthEdges, name='FluidPlateLength')

	flowOutletLengthEdges = fluidEdges.findAt(
                            ((0.0, 			-outletPlLength*0.5,	0.0),),
                            ((0.0, 			-outletPlLength*0.5,	plateThickness),),
                            ((0.0, 			-outletPlLength*0.5,	-lgChHeight),),
                            ((0.0, 			-outletPlLength*0.5,	smChHeight + plateThickness),),
                            ((plateWidth, 	-outletPlLength*0.5,	0.0),),
                            ((plateWidth, 	-outletPlLength*0.5,	plateThickness),),
                            ((plateWidth, 	-outletPlLength*0.5,	-lgChHeight),),
                            ((plateWidth, 	-outletPlLength*0.5,	smChHeight + plateThickness),))
	fluid.Set(edges=flowOutletLengthEdges, name='FluidOutletLength')

	flowInletLengthEdges = fluidEdges.findAt(
                            ((0.0, 			plateLength + (inletPlLength*0.5),	0.0),),
                            ((0.0, 			plateLength + (inletPlLength*0.5),	plateThickness),),
                            ((0.0, 			plateLength + (inletPlLength*0.5),	-lgChHeight),),
                            ((0.0, 			plateLength + (inletPlLength*0.5),	smChHeight + plateThickness),),
                            ((plateWidth, 	plateLength + (inletPlLength*0.5),	0.0),),
                            ((plateWidth, 	plateLength + (inletPlLength*0.5),	plateThickness),),
                            ((plateWidth, 	plateLength + (inletPlLength*0.5),	-lgChHeight),),
                            ((plateWidth, 	plateLength + (inletPlLength*0.5),	smChHeight + plateThickness),))
	fluid.Set(edges=flowInletLengthEdges, name='FluidInletLength')

	lgChannelEdges = fluidEdges.findAt(
                            ((0.0, 			-outletPlLength,				-lgChHeight*0.5),),
                            ((plateWidth,	-outletPlLength,				-lgChHeight*0.5),),
                            ((0.0, 			0.0,							-lgChHeight*0.5),),
                            ((plateWidth,	0.0,							-lgChHeight*0.5),),
                            ((0.0, 			plateLength,					-lgChHeight*0.5),),
                            ((plateWidth, 	plateLength,					-lgChHeight*0.5),),
                            ((0.0, 			plateLength + inletPlLength,	-lgChHeight*0.5),),
                            ((plateWidth, 	plateLength + inletPlLength,	-lgChHeight*0.5),))
	fluid.Set(edges=lgChannelEdges, name='LargeChHeight')

	smChannelEdges = fluidEdges.findAt(
                            ((0.0, 			-outletPlLength,				plateThickness + smChHeight*0.5),),
                            ((plateWidth,	-outletPlLength,				plateThickness + smChHeight*0.5),),
                            ((0.0, 			0.0,							plateThickness + smChHeight*0.5),),
                            ((plateWidth,	0.0,							plateThickness + smChHeight*0.5),),
                            ((0.0, 			plateLength,					plateThickness + smChHeight*0.5),),
                            ((plateWidth, 	plateLength,					plateThickness + smChHeight*0.5),),
                            ((0.0, 			plateLength + inletPlLength,	plateThickness + smChHeight*0.5),),
                            ((plateWidth, 	plateLength + inletPlLength,	plateThickness + smChHeight*0.5),))
	fluid.Set(edges=smChannelEdges, name='SmallChHeight')

	plateEdges = fluidEdges.findAt(
                            ((0.0, 			-outletPlLength,				plateThickness*0.5),),
                            ((plateWidth,	-outletPlLength,				plateThickness*0.5),),
                            ((0.0, 			0.0,							plateThickness*0.5),),
                            ((plateWidth,	0.0,							plateThickness*0.5),),
                            ((0.0, 			plateLength,					plateThickness*0.5),),
                            ((plateWidth, 	plateLength,					plateThickness*0.5),),
                            ((0.0, 			plateLength + inletPlLength,	plateThickness*0.5),),
                            ((plateWidth, 	plateLength + inletPlLength,	plateThickness*0.5),))
	fluid.Set(edges=plateEdges, name='PlateHeight')

	flowWidthEdges = fluidEdges.findAt(
                        ((plateWidth*0.5, 	-outletPlLength,				0.0),),
                        ((plateWidth*0.5, 	-outletPlLength,				plateThickness),),
                        ((plateWidth*0.5, 	-outletPlLength,				-lgChHeight),),
                        ((plateWidth*0.5, 	-outletPlLength,				smChHeight + plateThickness),),
                        ((plateWidth*0.5, 	0.0,							0.0),),
                        ((plateWidth*0.5, 	0.0,							plateThickness),),
                        ((plateWidth*0.5, 	0.0,							-lgChHeight),),
                        ((plateWidth*0.5, 	0.0,							smChHeight + plateThickness),),
                        ((plateWidth*0.5, 	plateLength,					0.0),),
                        ((plateWidth*0.5, 	plateLength,					plateThickness),),
                        ((plateWidth*0.5, 	plateLength,					-lgChHeight),),
                        ((plateWidth*0.5, 	plateLength,					smChHeight + plateThickness),),
                        ((plateWidth*0.5, 	plateLength + inletPlLength,	0.0),),
                        ((plateWidth*0.5, 	plateLength + inletPlLength,	plateThickness),),
                        ((plateWidth*0.5, 	plateLength + inletPlLength,	-lgChHeight),),
                        ((plateWidth*0.5, 	plateLength + inletPlLength,	smChHeight + plateThickness),))
	fluid.Set(edges=flowWidthEdges, name='FluidWidth')

	# Creating the set of faces for the inlet
	fluidFaces = fluid.faces
	inletFaces = fluidFaces.findAt(
                            ((plateWidth*0.5,	plateLength + inletPlLength,	plateThickness*0.5),),
    						((plateWidth*0.5,	plateLength + inletPlLength,	-smChHeight*0.5),),
                            ((plateWidth*0.5,	plateLength + inletPlLength,	plateThickness + lgChHeight*0.5),))
	fluid.Surface(side1Faces = inletFaces, name = 'Inlet')

	outletFaces = fluidFaces.findAt(
                            ((plateWidth*0.5,	-outletPlLength,	plateThickness*0.5),),
    						((plateWidth*0.5,	-outletPlLength,	-smChHeight*0.5),),
                            ((plateWidth*0.5,	-outletPlLength,	plateThickness + lgChHeight*0.5),))
	fluid.Surface(side1Faces = outletFaces, name = 'Outlet')

	# Creating the FSI surfaces
	fsiBack = fluidFaces.findAt(
                            ((plateWidth*0.5,		plateLength*0.5,	0.0),))
	fluid.Surface(side1Faces = fsiBack, name = 'FSI_Back')

	fsiFront = fluidFaces.findAt(
    						((plateWidth*0.5,		plateLength*0.5,	plateThickness),))
	fluid.Surface(side1Faces = fsiFront, name = 'FSI_Front')

	fsiTop = fluidFaces.findAt(
                            ((plateWidth*0.5,		plateLength,		plateThickness*0.5),),)
	fluid.Surface(side1Faces = fsiTop, name = 'FSI_Top')

	fsiBottom = fluidFaces.findAt(
    						((plateWidth*0.5,		0.0,				plateThickness*0.5),),)
	fluid.Surface(side1Faces = fsiBottom, name = 'FSI_Bottom')

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
	fileName = ( 'Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
			 '_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
			 '_' + str(int(parameters['lgChHeight']/0.0254*1000)) )
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
    plateMaterial = parameters['plateMaterial']
    elasticModulus = parameters['elasticModulus']
    poissonsRatio = parameters['poissonsRatio']
    plateDensity = parameters['plateDensity']
    pinnedBC = parameters['pinnedBC']

    # Grabbing the FSI coupling parameters
    couplingScheme = parameters['couplingScheme']
    timeStep = parameters['timeStep']
    minTimeStep = parameters['minTimeStep']
    maxSimTime = parameters['maxSimTime']
    r_in = 4*plateWidth/math.pi
    r_out = r_in + plateThickness
    plateTheta = (plateWidth + 0.0254)/r_in
    wettedTheta = math.pi/4
    clampedTheta = 0.5*(plateTheta - wettedTheta)
    emptyTheta = (math.pi/2 - plateTheta)/2

    #    Create Plate Part
    mdb.Model(name=modelName, modelType = STANDARD_EXPLICIT)
    plate = mdb.models[modelName].Part(name=plateName, dimensionality=THREE_D, type=DEFORMABLE_BODY )

    plate.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0)
    plate.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
    datumPlanes = plate.datums
    t = plate.MakeSketchTransform(sketchPlane=datumPlanes[1], origin=(0.0, 0.0, 0.0), 
                             sketchOrientation=RIGHT, sketchPlaneSide=SIDE1, sketchUpEdge=datumPlanes[2])
    
    plateSketch = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=3.0*plateWidth, 
                             transform=t)
							 
    #    Inner radius:
    ri = r_in
    #    Outer radius:
    ro = ri + plateThickness

    plateOrFluid = 'plate'
    if plateOrFluid == 'plate':
		theta = (plateWidth + 0.0254)/ri
         
    if plateOrFluid == 'fluid':
    	theta = math.pi/4

    #    Point 0: Center of radius
    point0 = (0.0, 0.0)
    #    Point 1: Inside, left
    point1 = (-ri*math.sin(math.pi/2 - emptyTheta), ri*math.cos(math.pi/2 - emptyTheta))
    #    Point 2: Inside, right
    point2 = (-ri*math.sin(emptyTheta), ri*math.cos(emptyTheta))
    #    Point 3: Outside, left
    point3 = (-ro*math.sin(math.pi/2 - emptyTheta), ro*math.cos(math.pi/2 - emptyTheta))
    #    Point 4: Outside, right
    point4 = (-ro*math.sin(emptyTheta), ro*math.cos(emptyTheta))
    #    Draw curved profile
    g, v, d, c = plateSketch.geometry, plateSketch.vertices, plateSketch.dimensions, plateSketch.constraints
    plateSketch.setPrimaryObject(option=STANDALONE )
    plateSketch.ArcByCenterEnds(center=point0, point1=point1, point2=point2, direction=CLOCKWISE)
    plateSketch.ArcByCenterEnds(center=point0, point1=point3, point2=point4, direction=CLOCKWISE)
    plateSketch.Line(point1=point1, point2=point3)
    plateSketch.Line(point1=point2, point2=point4)
    	
    plate.BaseSolidExtrude(sketch=plateSketch, depth=plateLength)
    plateSketch.unsetPrimaryObject()
    del mdb.models[modelName].sketches['__profile__']
    plate.regenerate
    mdb.models[modelName].parts[plateName].setValues(geometryRefinement=EXTRA_FINE)
    
    # Creating datum points for the three datum planes
    plate.DatumPointByCoordinate(coords=(r_in*math.cos(math.pi/2 - clampedTheta - emptyTheta), plateLength, r_in*math.sin(math.pi/2 - clampedTheta - emptyTheta)))
    plate.DatumPointByCoordinate(coords=(r_in*math.cos(emptyTheta + clampedTheta), plateLength, r_in*math.sin(emptyTheta + clampedTheta)))
    plate.DatumPointByCoordinate(coords=(r_out*math.cos(math.pi/2 - clampedTheta - emptyTheta), plateLength, r_out*math.sin(math.pi/2 - clampedTheta - emptyTheta)))
    plate.DatumPointByCoordinate(coords=(r_out*math.cos(emptyTheta + clampedTheta), plateLength, r_out*math.sin(emptyTheta + clampedTheta)))

    plate.DatumPointByCoordinate(coords=(r_in*math.cos(math.pi/2 - clampedTheta - emptyTheta), 0.0, r_in*math.sin(math.pi/2 - clampedTheta - emptyTheta)))
    plate.DatumPointByCoordinate(coords=(r_in*math.cos(emptyTheta + clampedTheta), 0.0, r_in*math.sin(emptyTheta + clampedTheta)))
	
    plate.DatumPointByCoordinate(coords=(r_in*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)), plateLength, r_in*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))))
    plate.DatumPointByCoordinate(coords=(r_out*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)), plateLength, r_out*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))))
    plate.DatumPointByCoordinate(coords=(r_in*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)), 0.0, r_in*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))))
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
    
    pickedCells = plateCells.findAt(((r_in*math.cos(plateTheta/2), plateLength*0.5, (r_in + (plateThickness*0.5))*math.sin(plateTheta/2),)))
    plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[13], cells=pickedCells)
    
    pickedCells = plateCells.findAt(((r_in*math.cos(plateTheta), plateLength*0.25, (r_in + (plateThickness*0.5))*math.sin(plateTheta),)))
    plate.PartitionCellByDatumPlane(datumPlane=datumPlanes[14], cells=pickedCells)
    
    # Creating all the sets of edges in the plate part
    plateEdges = plate.edges
    flowLengthEdges = plateEdges.findAt(
        ((r_in*math.cos(math.pi/2 - emptyTheta),						plateLength*0.5, 	r_in*math.sin(math.pi/2 - emptyTheta)),),
        ((r_out*math.cos(math.pi/2 - emptyTheta),						plateLength*0.5,	r_out*math.sin(math.pi/2 - emptyTheta)),),
        ((r_in*math.cos(math.pi/2 - emptyTheta - clampedTheta),			plateLength*0.5,	r_in*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        ((r_out*math.cos(math.pi/2 - emptyTheta - clampedTheta),		plateLength*0.5,	r_out*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        ((r_in*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)),	plateLength*0.5,	r_in*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))),),
        ((r_out*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)), 	plateLength*0.5,	r_out*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))),),
        ((r_in*math.cos(emptyTheta + clampedTheta), 					plateLength*0.5,	r_in*math.sin(emptyTheta + clampedTheta)),),
        ((r_out*math.cos(emptyTheta + clampedTheta), 					plateLength*0.5,	r_out*math.sin(emptyTheta + clampedTheta)),),
        ((r_in*math.cos(emptyTheta), 									plateLength*0.5,	r_in*math.sin(emptyTheta)),),
        ((r_out*math.cos(emptyTheta),									plateLength*0.5,	r_out*math.sin(emptyTheta)),),)
    plate.Set(edges=flowLengthEdges, name='PlateLength')
    
    flowWidthEdges = plateEdges.findAt(
        ((r_in*math.cos(emptyTheta + clampedTheta + wettedTheta*0.75), 		plateLength,	r_in*math.sin(emptyTheta + clampedTheta + wettedTheta*0.75)),),
        ((r_out*math.cos(emptyTheta + clampedTheta + wettedTheta*0.75), 	plateLength, 	r_out*math.sin(emptyTheta + clampedTheta + wettedTheta*0.75)),),
        ((r_in*math.cos(emptyTheta + clampedTheta + wettedTheta*0.25), 		plateLength, 	r_in*math.sin(emptyTheta + clampedTheta + wettedTheta*0.25)),),
        ((r_out*math.cos(emptyTheta + clampedTheta + wettedTheta*0.25), 	plateLength, 	r_out*math.sin(emptyTheta + clampedTheta + wettedTheta*0.25)),),
        ((r_in*math.cos(emptyTheta + clampedTheta + wettedTheta*0.75), 		0.0,			r_in*math.sin(emptyTheta + clampedTheta + wettedTheta*0.75)),),
        ((r_out*math.cos(emptyTheta + clampedTheta + wettedTheta*0.75), 	0.0, 			r_out*math.sin(emptyTheta + clampedTheta + wettedTheta*0.75)),),
        ((r_in*math.cos(emptyTheta + clampedTheta + wettedTheta*0.25), 		0.0, 			r_in*math.sin(emptyTheta + clampedTheta + wettedTheta*0.25)),),
        ((r_out*math.cos(emptyTheta + clampedTheta + wettedTheta*0.25), 	0.0, 			r_out*math.sin(emptyTheta + clampedTheta + wettedTheta*0.25)),))
    plate.Set(edges=flowWidthEdges, name = 'PlateWidth')
    
    plateThicknessEdges = plateEdges.findAt(
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta),						0.0,			(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta)),),
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		0.0,			(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)),	0.0,			(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta), 					0.0,			(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta), 									0.0,			(r_in + plateThickness/2)*math.sin(emptyTheta)),),
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta),						plateLength,	(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta)),),
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		plateLength,	(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta + (wettedTheta/2)),	plateLength,	(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta + (wettedTheta/2))),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta), 					plateLength,	(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta), 									plateLength,	(r_in + plateThickness/2)*math.sin(emptyTheta)),))
    plate.Set(edges = plateThicknessEdges, name = 'PlateThickness')
    
    clampedWidthEdges = plateEdges.findAt(
        ((r_in*math.cos(math.pi/2 - emptyTheta - (clampedTheta/2)), 	plateLength,	r_in*math.sin(math.pi/2 - emptyTheta - (clampedTheta/2))),),
        ((r_out*math.cos(math.pi/2 - emptyTheta - (clampedTheta/2)), 	plateLength, 	r_out*math.sin(math.pi/2 - emptyTheta - (clampedTheta/2))),),
        ((r_in*math.cos(emptyTheta + (clampedTheta/2)), 				plateLength, 	r_in*math.sin(emptyTheta + (clampedTheta/2))),),
        ((r_out*math.cos(emptyTheta + (clampedTheta/2)),				plateLength, 	r_out*math.sin(emptyTheta + (clampedTheta/2))),),
        ((r_in*math.cos(math.pi/2 - emptyTheta - (clampedTheta/2)), 	0.0,			r_in*math.sin(math.pi/2 - emptyTheta - (clampedTheta/2))),),
        ((r_out*math.cos(math.pi/2 - emptyTheta - (clampedTheta/2)), 	0.0, 			r_out*math.sin(math.pi/2 - emptyTheta - (clampedTheta/2))),),
        ((r_in*math.cos(emptyTheta + (clampedTheta/2)), 				0.0, 			r_in*math.sin(emptyTheta + (clampedTheta/2))),),
        ((r_out*math.cos(emptyTheta + (clampedTheta/2)),				0.0, 			r_out*math.sin(emptyTheta + (clampedTheta/2))),))
    plate.Set(edges = clampedWidthEdges, name = 'ClampedWidth')
    
    # Creating the set of faces where the plate will be clamped
    plateFaces = plate.faces
    clampedFaces = plateFaces.findAt(
        ((r_in*math.cos(math.pi/2 - emptyTheta - (clampedTheta/2)),		plateLength*0.5,	r_in*math.sin(math.pi/2 - emptyTheta - (clampedTheta/2))),),
		((r_out*math.cos(math.pi/2 - emptyTheta - (clampedTheta/2)),	plateLength*0.5,	r_out*math.sin(math.pi/2 - emptyTheta - (clampedTheta/2))),),
        ((r_in*math.cos(emptyTheta + (clampedTheta/2)),					plateLength*0.5,	r_in*math.sin(emptyTheta + (clampedTheta/2))),),
		((r_out*math.cos(emptyTheta + (clampedTheta/2)),				plateLength*0.5,	r_out*math.sin(emptyTheta + (clampedTheta/2))),))
    plate.Set(faces = clampedFaces, name = 'ClampedFaces')
    
    # Creating the set of vertices where the plate will be pinned
    plateVertices = plate.vertices
    pinVertices = plateVertices.findAt(
		((r_in*math.cos(math.pi/4),		plateLength,	r_in*math.sin(math.pi/4)),),
		((r_in*math.cos(math.pi/4),		0.0,			r_in*math.sin(math.pi/4)),),)
    plate.Set(vertices = pinVertices, name = 'Pins')
    
    # Creating the FSI_INTERFACE surface
    fsiFaces = plateFaces.findAt(
        ((r_in*math.cos(emptyTheta + clampedTheta + wettedTheta/4),							plateLength*0.5,	r_in*math.sin(emptyTheta + clampedTheta + wettedTheta/4)),),
        ((r_out*math.cos(emptyTheta + clampedTheta + wettedTheta/4), 						plateLength*0.5, 	r_out*math.sin(emptyTheta + clampedTheta + wettedTheta/4)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta + wettedTheta/4),	0.0, 				(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta + wettedTheta/4)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta + wettedTheta/4),	plateLength, 		(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta + wettedTheta/4)),),
        ((r_in*math.cos(emptyTheta + clampedTheta + wettedTheta*3/4), 						plateLength*0.5,	r_in*math.sin(emptyTheta + clampedTheta + wettedTheta*3/4)),),
        ((r_out*math.cos(emptyTheta + clampedTheta + wettedTheta*3/4),						plateLength*0.5, 	r_out*math.sin(emptyTheta + clampedTheta + wettedTheta*3/4)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta + wettedTheta*3/4), 	0.0, 				(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta + wettedTheta*3/4)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta + wettedTheta*3/4), 	plateLength, 		(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta + wettedTheta*3/4)),),)
    plate.Surface(side1Faces = fsiFaces, name = 'FSI_INTERFACE')
    
    # Creating the plate material aluminum
    mdb.models[modelName].Material(name=plateMaterial)
    mdb.models[modelName].materials[plateMaterial].Density(table=((plateDensity, ), ))
    mdb.models[modelName].materials[plateMaterial].Elastic(table=((elasticModulus, poissonsRatio), ))

    # Setting and creating the section assignment for the plate
    mdb.models[modelName].HomogeneousSolidSection(name='PlateSection', 
                              material=plateMaterial, thickness=None)
    plateRegion = plate.sets['EntirePlateGeometry']
    plate.SectionAssignment(region=plateRegion, sectionName='PlateSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    
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

    # Adding the plate as an instance to the assembly
    assembly = mdb.models[modelName].rootAssembly
    assembly.Instance(name = plateName, part = plate, dependent = ON)

    # Setting the clamped boundary condition on the clamped regions of the plate
    clampedRegion = assembly.instances[plateName].sets['ClampedFaces']
    mdb.models[modelName].PinnedBC(name = 'Clamps', createStepName = 'Initial', 
                              region = clampedRegion, localCsys = None)

    # Setting the pins boundary condition on the 'Pins' set
    if pinnedBC == 'yes':
    	pinRegion = assembly.instances[plateName].sets['Pins']
    	mdb.models[modelName].DisplacementBC(name='Pins', 
    		createStepName='Initial', region=pinRegion, u1=SET, u2=UNSET, u3=SET, 
    		ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, 
    		distributionType=UNIFORM, fieldName='', localCsys=None)

    # Creating and writing the input file to the Abaqus work directory
    if pinnedBC == 'yes':
		x = str(int(plateThickness/0.0254*1000)) + '_Pinned'

    if pinnedBC == 'no':
		x = str(int(plateThickness/0.0254*1000)) + '_Free'
		
    mdb.Job(name=couplingScheme + '_' + x + 'Plate', model=modelName, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
    mdb.jobs[couplingScheme + '_' + x + 'Plate'].writeInput(consistencyChecking=OFF)

    # Appending the written input file for explicit FSI coupling
    with open(couplingScheme + "_" + x + "Plate.inp", "a") as inputFile:
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
    		"*Dynamic,application=QUASI-STATIC\n",
    		"" + str(timeStep) + "," + str(maxSimTime) + "," + str(minTimeStep) + "," + str(timeStep) + "\n",
    		"**\n",
    		"** OUTPUT REQUESTS\n",
    		"**\n",
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
    		"ASSEMBLY_PLATE_FSI_INTERFACE, U\n",
    		"ASSEMBLY_PLATE_FSI_INTERFACE, V\n",
    		"*CO-SIMULATION REGION, TYPE=SURFACE, IMPORT\n",
    		"ASSEMBLY_PLATE_FSI_INTERFACE, CF\n",
    		"*CO-SIMULATION CONTROLS, NAME=Control-1, COUPLING SCHEME=" + couplingScheme +", SCHEME MODIFIER=LAG,\n",
    		"STEP SIZE=" + str(timeStep) + ", TIME INCREMENTATION=SUBCYCLE, TIME MARKS=YES\n",
    		"**\n",
    		"*End Step\n"]			
    	inputFile.writelines(inputFileLines)
    inputFile.close()

def createCurvedFluid(parameters):
    modelName = parameters['abaqusModelName']
    plateName = 'Plate'
    fluidName = 'Fluid'
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

    r_in = plateWidth/(math.pi/4)
    r_out = r_in + plateThickness
    plateTheta = (plateWidth + 0.0254)/r_in
    wettedTheta = math.pi/4
    clampedTheta = 0.5*(plateTheta - wettedTheta)
    emptyTheta = (math.pi/2 - plateTheta)/2

    r_lgCh = r_in - lgChHeight
    
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
    
	#    Inner radius:
    ri = r_lgCh
    #    Outer radius:
    ro = ri + plateThickness + smChHeight + lgChHeight

    plateOrFluid = 'fluid'
    if plateOrFluid == 'plate':
    	theta = (plateWidth + 0.0254)/ri
         
    if plateOrFluid == 'fluid':
    	theta = math.pi/4
    	plateTheta = (plateWidth + 0.0254)/(plateWidth/(math.pi/4))
    	clampedTheta = 0.5*(plateTheta - theta)
    	#    sine and cosine factors
    	#cos = math.cos(theta)
    	#sin = math.sin(theta)
    	#    Point 0: Center of radius
    	point0 = (0.0, 0.0)
    	#    Point 1: Inside, left
    	point1 = (-ri*sin(math.pi/2 - emptyTheta - clampedTheta), ri*cos(math.pi/2 - emptyTheta - clampedTheta))
    	#    Point 2: Inside, right
    	point2 = (-ri*sin(emptyTheta + clampedTheta), ri*cos(emptyTheta + clampedTheta))
    	#    Point 3: Outside, left
    	point3 = (-ro*sin(math.pi/2 - emptyTheta - clampedTheta), ro*cos(math.pi/2 - emptyTheta - clampedTheta))
    	#    Point 4: Outside, right
    	point4 = (-ro*sin(emptyTheta + clampedTheta), ro*cos(emptyTheta + clampedTheta))

    #    Draw curved profile
    g, v, d, c = bulkFluidSketch.geometry, bulkFluidSketch.vertices, bulkFluidSketch.dimensions, bulkFluidSketch.constraints
    bulkFluidSketch.setPrimaryObject(option=STANDALONE )
    bulkFluidSketch.ArcByCenterEnds(center=point0, point1=point1, point2=point2, direction=CLOCKWISE)
    bulkFluidSketch.ArcByCenterEnds(center=point0, point1=point3, point2=point4, direction=CLOCKWISE)
    bulkFluidSketch.Line(point1=point1, point2=point3)
    bulkFluidSketch.Line(point1=point2, point2=point4)
    
	#curved_sketch(bulkFluidSketch, r_lgCh, plateWidth, plateThickness + smChHeight + lgChHeight, 'fluid')

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
    assembly.InstanceFromBooleanCut(name=fluidName, 
        instanceToBeCut=assembly.instances['BulkFluid'], 
        cuttingInstances=(assembly.instances[plateName], ),
        originalInstances=SUPPRESS)
    assembly.features.changeKey( fromName = 'Fluid-1', toName = 'Fluid' )
    mdb.models[modelName].parts[fluidName].setValues(geometryRefinement=EXTRA_FINE)
    
    # Creating the datum points for creating the four datum planes
    fluid = mdb.models[modelName].parts[fluidName]
    fluid.DatumPointByCoordinate(coords=(r_in*math.cos(math.pi/2 - emptyTheta - clampedTheta),		plateLength, r_in*math.sin(math.pi/2 - emptyTheta - clampedTheta)))
    fluid.DatumPointByCoordinate(coords=(r_in*math.cos(emptyTheta + clampedTheta), 					plateLength, r_in*math.sin(emptyTheta + clampedTheta)))
    fluid.DatumPointByCoordinate(coords=(r_out*math.cos(math.pi/2 - emptyTheta - clampedTheta), 	plateLength, r_out*math.sin(math.pi/2 - emptyTheta - clampedTheta)))
    fluid.DatumPointByCoordinate(coords=(r_out*math.cos(emptyTheta + clampedTheta),					plateLength, r_out*math.sin(emptyTheta + clampedTheta)))
    
    fluid.DatumPointByCoordinate(coords=(r_in*math.cos(math.pi/2 - emptyTheta - clampedTheta),		0.0, 		r_in*math.sin(math.pi/2 - emptyTheta - clampedTheta)))
    fluid.DatumPointByCoordinate(coords=(r_in*math.cos(emptyTheta + clampedTheta), 					0.0, 		r_in*math.sin(emptyTheta + clampedTheta)))
    fluid.DatumPointByCoordinate(coords=(r_out*math.cos(math.pi/2 - emptyTheta - clampedTheta), 	0.0, 		r_out*math.sin(math.pi/2 - emptyTheta - clampedTheta)))
    fluid.DatumPointByCoordinate(coords=(r_out*math.cos(emptyTheta + clampedTheta),					0.0, 		r_out*math.sin(emptyTheta + clampedTheta)))
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
    fluidEdges = fluid.edges
    flowPlateLengthEdges = fluidEdges.findAt(
        ((r_in*math.cos(math.pi/2 - emptyTheta - clampedTheta),										plateLength*0.5,	r_in*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        ((r_out*math.cos(math.pi/2 - emptyTheta - clampedTheta), 									plateLength*0.5,	r_out*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        ((r_in*math.cos(emptyTheta + clampedTheta), 												plateLength*0.5,	r_in*math.sin(emptyTheta + clampedTheta)),),
        ((r_out*math.cos(emptyTheta + clampedTheta),												plateLength*0.5,	r_out*math.sin(emptyTheta + clampedTheta)),),
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		plateLength*0.5,	(r_in + smChHeight + plateThickness)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + smChHeight + plateThickness)*math.cos(emptyTheta + clampedTheta),					plateLength*0.5,	(r_in + smChHeight + plateThickness)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in - lgChHeight)*math.cos(math.pi/2 - emptyTheta - clampedTheta),						plateLength*0.5,	(r_in - lgChHeight)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in - lgChHeight)*math.cos(emptyTheta + clampedTheta),									plateLength*0.5,	(r_in - lgChHeight)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=flowPlateLengthEdges, name='FluidPlateLength')
    
    flowOutletLengthEdges = fluidEdges.findAt(
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/2 - emptyTheta - clampedTheta),	-outletPlLength*0.5,	(r_in + smChHeight + plateThickness)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + smChHeight + plateThickness)*math.cos(emptyTheta + clampedTheta),				-outletPlLength*0.5,	(r_in + smChHeight + plateThickness)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in - lgChHeight)*math.cos(math.pi/2 - emptyTheta - clampedTheta),					-outletPlLength*0.5,	(r_in - lgChHeight)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in - lgChHeight)*math.cos(emptyTheta + clampedTheta),								-outletPlLength*0.5,	(r_in - lgChHeight)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=flowOutletLengthEdges, name='FluidOutletLength')
    
    flowInletLengthEdges = fluidEdges.findAt(
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/2 - emptyTheta - clampedTheta),	plateLength + inletPlLength*0.5,	(r_in + smChHeight + plateThickness)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + smChHeight + plateThickness)*math.cos(emptyTheta + clampedTheta),				plateLength + inletPlLength*0.5,	(r_in + smChHeight + plateThickness)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in - lgChHeight)*math.cos(math.pi/2 - emptyTheta - clampedTheta),					plateLength + inletPlLength*0.5,	(r_in - lgChHeight)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in - lgChHeight)*math.cos(emptyTheta + clampedTheta),								plateLength + inletPlLength*0.5,	(r_in - lgChHeight)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=flowInletLengthEdges, name='FluidInletLength')
    
    smChannelEdges = fluidEdges.findAt(
        (((r_in + smChHeight/2 + plateThickness)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		plateLength,	(r_in + smChHeight/2 + plateThickness)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + smChHeight/2 + plateThickness)*math.cos(emptyTheta + clampedTheta),	plateLength,	(r_in + smChHeight/2 + plateThickness)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in + smChHeight/2 + plateThickness)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		0.0,			(r_in + smChHeight/2 + plateThickness)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + smChHeight/2 + plateThickness)*math.cos(emptyTheta + clampedTheta),					0.0,			(r_in + smChHeight/2 + plateThickness)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=smChannelEdges, name='SmallChHeight')
    
    lgChannelEdges = fluidEdges.findAt(
        (((r_in - lgChHeight/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		plateLength,	(r_in - lgChHeight/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in - lgChHeight/2)*math.cos(emptyTheta + clampedTheta),					plateLength,	(r_in - lgChHeight/2)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in - lgChHeight/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),		0.0,			(r_in - lgChHeight/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in - lgChHeight/2)*math.cos(emptyTheta + clampedTheta),					0.0,			(r_in - lgChHeight/2)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=lgChannelEdges, name='LargeChHeight')

    plateEdges = fluidEdges.findAt(
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),	plateLength,	(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta),				plateLength,	(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),	0.0,			(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta),				0.0,			(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=plateEdges, name='PlateHeight')

    inOutEdges = fluidEdges.findAt(
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),	plateLength + inletPlLength,	(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta),				plateLength + inletPlLength,	(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(math.pi/2 - emptyTheta - clampedTheta),	-outletPlLength,				(r_in + plateThickness/2)*math.sin(math.pi/2 - emptyTheta - clampedTheta)),),
        (((r_in + plateThickness/2)*math.cos(emptyTheta + clampedTheta),				-outletPlLength,				(r_in + plateThickness/2)*math.sin(emptyTheta + clampedTheta)),),)
    fluid.Set(edges=inOutEdges, name='InletOutletHeight')
    
    flowWidthEdges = fluidEdges.findAt(
        (((r_in - lgChHeight)*math.cos(math.pi/4), 					-outletPlLength,				(r_in - lgChHeight)*math.sin(math.pi/4)),),
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/4), -outletPlLength,				(r_in + smChHeight + plateThickness)*math.sin(math.pi/4)),),
        ((r_in*math.cos(math.pi/4), 								0.0,							r_in*math.sin(math.pi/4)),),
        (((r_in + plateThickness)*math.cos(math.pi/4), 				0.0,							(r_in + plateThickness)*math.sin(math.pi/4)),),
        (((r_in - lgChHeight)*math.cos(math.pi/4), 					0.0,							(r_in - lgChHeight)*math.sin(math.pi/4)),),
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/4), 0.0,							(r_in + smChHeight + plateThickness)*math.sin(math.pi/4)),),
        ((r_in*math.cos(math.pi/4), 								plateLength,					r_in*math.sin(math.pi/4)),),
        (((r_in + plateThickness)*math.cos(math.pi/4), 				plateLength,					(r_in + plateThickness)*math.sin(math.pi/4)),),
        (((r_in - lgChHeight)*math.cos(math.pi/4), 					plateLength,					(r_in - lgChHeight)*math.sin(math.pi/4)),),
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/4), plateLength,					(r_in + smChHeight + plateThickness)*math.sin(math.pi/4)),),
        (((r_in - lgChHeight)*math.cos(math.pi/4), 					plateLength + inletPlLength,	(r_in - lgChHeight)*math.sin(math.pi/4)),),
        (((r_in + smChHeight + plateThickness)*math.cos(math.pi/4), plateLength + inletPlLength,	(r_in + smChHeight + plateThickness)*math.sin(math.pi/4)),))
    fluid.Set(edges=flowWidthEdges, name='FluidWidth')
    
    # Creating the set of faces for the inlet
    fluidFaces = fluid.faces
    inletFaces = fluidFaces.findAt(
        (((r_in + (plateThickness)*0.5)*math.cos(math.pi/4),		plateLength + inletPlLength,	(r_in + (plateThickness + smChHeight + lgChHeight)*0.5)*math.sin(math.pi/4)),),)
    fluid.Surface(side1Faces = inletFaces, name = 'Inlet')
    
    outletFaces = fluidFaces.findAt(
        (((r_in + (plateThickness)*0.5)*math.cos(math.pi/4),		-outletPlLength,				(r_in + (plateThickness + smChHeight + lgChHeight)*0.5)*math.sin(math.pi/4)),),)
    fluid.Surface(side1Faces = outletFaces, name = 'Outlet')

    # Creating the FSI surfaces
    fsiBack = fluidFaces.findAt(
        ((r_in*math.cos(math.pi/4),									plateLength*0.5,	r_in*math.sin(math.pi/4)),),)
    fluid.Surface(side1Faces = fsiBack, name = 'FSI_Back')
    
    fsiFront = fluidFaces.findAt(
        (((r_in + plateThickness)*math.cos(math.pi/4),				plateLength*0.5,	(r_in + plateThickness)*math.sin(math.pi/4)),),)
    fluid.Surface(side1Faces = fsiFront, name = 'FSI_Front')

    fsiTop = fluidFaces.findAt(
        (((r_in + (plateThickness*0.5))*math.cos(math.pi/4),		plateLength,		(r_in + (plateThickness*0.5))*math.sin(math.pi/4)),),)
    fluid.Surface(side1Faces = fsiTop, name = 'FSI_Top')

    fsiBottom = fluidFaces.findAt(
        (((r_in + (plateThickness*0.5))*math.cos(math.pi/4),		0.0,				(r_in + (plateThickness*0.5))*math.sin(math.pi/4)),),)
    fluid.Surface(side1Faces = fsiBottom, name = 'FSI_Bottom')

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

    fluidInletOutletEdges = fluid.sets['InletOutletHeight'].edges
    fluid.seedEdgeByBias(biasMethod=DOUBLE, endEdges=fluidInletOutletEdges, ratio=flPlHeightBias,
        number=flPlHeightNodes + flLgChHeightNodes + flSmChHeightNodes, constraint=FINER)

    # Setting the mesh element type and meshing the part
    elemType = mesh.ElemType(elemCode = C3D8I, elemLibrary = STANDARD, secondOrderAccuracy = OFF,
                              distortionControl = DEFAULT)
    fluidCellsRegion = fluid.sets['EntireFluidGeometry'].cells
    fluid.setElementType(regions = (fluidCellsRegion, ), elemTypes = (elemType, ))
    fluid.generateMesh()

    # Creating a job for input file creation
    mdb.Job(name='Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
        '_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
        '_' + str(int(parameters['lgChHeight']/0.0254*1000)),
        model=modelName, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)

    # Writing the input file to the Abaqus work directory
    mdb.jobs['Star_Fluid_' + str(int(parameters['plateThickness']/0.0254*1000)) + 
        '_' + str(int(parameters['smChHeight']/0.0254*1000)) + 
        '_' + str(int(parameters['lgChHeight']/0.0254*1000))].writeInput(consistencyChecking=OFF)

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
    		"ASSEMBLY_PLATE_FSI_INTERFACE, U\n",
    		"ASSEMBLY_PLATE_FSI_INTERFACE, V\n",
    		"*CO-SIMULATION REGION, TYPE=SURFACE, IMPORT\n",
    		"ASSEMBLY_PLATE_FSI_INTERFACE, CF\n",
    		"*CO-SIMULATION CONTROLS, NAME=Control-1, COUPLING SCHEME=" + couplingScheme +", SCHEME MODIFIER=LAG,\n",
    		"STEP SIZE=" + str(timeStep) + ", TIME INCREMENTATION=SUBCYCLE, TIME MARKS=YES\n",
    		"**\n",
    		"*End Step\n"]			
		inputFile.writelines(inputFileLines)
	inputFile.close()

def curved_sketch(s, innerRadius, plateWidth, thickness, plateOrFluid):
    """ Creates a sketch of the plate cross section
        Inputs:
        s:            Constrained sketch part
        width:      Plate width on the inner radius
        thickness:      Plate thickness
        plateOrFluid:	'plate' for plate part and 'fluid' for fluid part
        Outputs:
        s:            Completed constrained sketch
    """
    #    Inner radius:
    ri = innerRadius
    #    Outer radius:
    ro = ri + thickness

    if plateOrFluid == 'plate':
		theta = (plateWidth + 0.0254)/ri
		#    sine and cosine factors
		cos = math.cos(theta)
		sin = math.sin(theta)
		#    Point 0: Center of radius
		point0 = (0.0, 0.0)
		#    Point 1: Inside, left
		point1 = (0, ri)
		#    Point 2: Inside, right
		point2 = (-ri*sin, ri*cos)
		#    Point 3: Outside, left
		point3 = (0, ro)
		#    Point 4: Outside, right
		point4 = (-ro*sin, ro*cos)
         
    if plateOrFluid == 'fluid':
    	theta = math.pi/4
    	plateTheta = (plateWidth + 0.0254)/(plateWidth/(math.pi/4))
    	clampedTheta = 0.5*(plateTheta - theta)
    	#    Point 0: Center of radius
    	point0 = (0.0, 0.0)
    	#    Point 1: Inside, left
    	point1 = (-ri*math.sin(clampedTheta), ri*math.cos(clampedTheta))
    	#    Point 2: Inside, right
    	point2 = (-ri*math.sin(plateTheta - clampedTheta), ri*math.cos(plateTheta - clampedTheta))
    	#    Point 3: Outside, left
    	point3 = (-ro*math.sin(clampedTheta), ro*math.cos(clampedTheta))
    	#    Point 4: Outside, right
    	point4 = (-ro*math.sin(plateTheta - clampedTheta), ro*math.cos(plateTheta - clampedTheta))

    #    Draw curved profile
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE )
    s.ArcByCenterEnds(center=point0, point1=point2, point2=point1, direction=CLOCKWISE)
    s.ArcByCenterEnds(center=point0, point1=point4, point2=point3, direction=CLOCKWISE)
    s.Line(point1=point1, point2=point3)
    s.Line(point1=point2, point2=point4)
    return s

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
		mdb.saveAs('FSI_PlateGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
				'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + str(int(parameters['lgChHeight']/0.0254*1000)))
	print('\n\n The Abaqus script has completed building the plate model\n')

	createFlatFluid(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		mdb.saveAs('FSI_FluidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
				'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + str(int(parameters['lgChHeight']/0.0254*1000)))
	print('\n\n The Abaqus script has completed building the fluid model\n')

if parameters['plateGeometry'] == 'Curved':
	createCurvedPlate(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		mdb.saveAs('FSI_PlateGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
				'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + str(int(parameters['lgChHeight']/0.0254*1000)))
	print('\n\n The Abaqus script has completed building the plate model\n')

	createCurvedFluid(parameters)
	if parameters['saveCAEFiles'] == 'yes':
		mdb.saveAs('FSI_FluidGeometry_' + str(int(parameters['plateThickness']/0.0254*1000)) +
				'_' + str(int(parameters['smChHeight']/0.0254*1000)) + '_' + str(int(parameters['lgChHeight']/0.0254*1000)))
	print('\n\n The Abaqus script has completed building the fluid model\n')
	