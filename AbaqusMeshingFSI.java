/** This java script is to be run using the FSI_Abaqus_StarCCM.py python script within Star-CCM+
 * 		Written by Casey J. Jesse on June 2014 at the University of Missouri - Columbia	
 * 		Revisions:
 *			July 17, 2013 - Changed the boundaries included in the wall y+ scene and plot
 *			July 17, 2013 - Added a section for creating a solution history file
 *			
 */
package myStarJavaMacros;

import java.io.File;
import java.io.IOException;

import star.base.report.AreaAverageReport;
import star.base.report.MaxReport;
import star.base.report.MinReport;
import star.common.ImplicitUnsteadyModel;
import star.common.PhysicsContinuum;
import star.common.PrimitiveFieldFunction;
import star.common.Simulation;
import star.common.StarMacro;
import star.common.VectorMagnitudeFieldFunction;
import star.common.XYPlot;
import star.cosimulation.abaqus.AbaqusCoSimulation;
import star.cosimulation.abaqus.AbaqusCoSimulationModel;
import star.cosimulation.common.*;
import star.flow.ConstantDensityModel;
import star.keturb.KEpsilonTurbulence;
import star.keturb.KeTwoLayerAllYplusWallTreatment;
import star.keturb.RkeTwoLayerTurbModel;
import star.material.SingleComponentLiquidModel;
import star.metrics.ThreeDimensionalModel;
import star.segregatedflow.SegregatedFlowModel;
import star.coupledflow.CoupledFlowModel;
import star.turbulence.RansTurbulenceModel;
import star.turbulence.TurbulentModel;
import star.vis.LinePart;
import star.vis.Scene;
import starClasses.CoSimulationAbaqus;
import starClasses.ContiuumBuilder;
import starClasses.DerivedParts;
import starClasses.FieldFunctions;
import starClasses.ImportCAE;
import starClasses.MeshMorpher;
import starClasses.NewDataReader;
import starClasses.RegionBuilder;
import starClasses.ReportsMonitorsPlots;
import starClasses.Scenes;
import starClasses.SolutionHistoryCreator;
import starClasses.SolversNode;
import starClasses.StoppingCriteria;
import starClasses.Tools;

public class AbaqusMeshingFSI extends StarMacro 
{
	public void execute() 
	{
		String currentDirectory = System.getProperty("user.dir");
		String abqExecutableFileName = "abq6122.bat";
		
		// Star-CCM+ settings and variables
		Simulation activeSim = getActiveSimulation();
		
		// Reading in the geometry parameters from the external file
		NewDataReader reader = new NewDataReader();
		try 
		{
			reader.readGeometryData(currentDirectory + File.separator + "FSI_Input_File.txt");
		}
		catch (NumberFormatException e1) 
		{
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		catch (IOException e1) 
		{
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		// Grabbing the Co-Simulation settings and variables
		String plateGeometry = reader.getStringData("plateGeometry");
		String couplingScheme = reader.getStringData("couplingScheme");
		double couplingTimeStep = reader.getDoubleData("timeStep");
		int numAbaqusCPUs = reader.getIntData("abaqusCPUs");
		int numExchanges = reader.getIntData("numImplicitExch");
		int iterationsPerExchange = reader.getIntData("iterPerExch");
		int iterationsPerTS = reader.getIntData("iterPerTS");
		double deflectionUnderRelax = reader.getDoubleData("plateUnderRelax");
		String morphAtInnIter = reader.getStringData("morphAtInnIter");
		double maxSimTime = reader.getDoubleData("maxSimTime");
		double courantNumber = reader.getDoubleData("courantNumber");
		boolean fluidStopCriteria = Boolean.parseBoolean(reader.getStringData("fluidStopCriteria"));
		String pinOrCombBC = reader.getStringData("pinOrCombBC");
		String SSorFSI = reader.getStringData("steadyStateOrFSI");
		double wallHeight = reader.getDoubleData("flSmChHeightBias");
		
		// Grabbing the geometry parameters
		double plateLength = reader.getDoubleData("plateLength");
		double plateHeight = reader.getDoubleData("plateThickness");
		double wettedPlateWidth = reader.getDoubleData("plateWidth");
		double smallChannelHeight = reader.getDoubleData("smChHeight");
		double largeChannelHeight = reader.getDoubleData("lgChHeight");
		double inletLength = reader.getDoubleData("inletPlLength");
		double outletLength = reader.getDoubleData("outletPlLength");
		double avgChVel = reader.getDoubleData("avgChVelocity");
		
		double r_in = wettedPlateWidth/(Math.PI/4);
		
		double[] initialVel = {0.0, -avgChVel, 0.0};
		double inletVel = Math.abs((smallChannelHeight + largeChannelHeight)/(smallChannelHeight + largeChannelHeight + plateHeight)*avgChVel);		
		
		String abaqusInputFilePath = null;
		//String abaqusInputFilePath = currentDirectory + File.separator + couplingScheme + "_PinnedPlate.inp";;
		if(pinOrCombBC.equals("pin"))
		{
			abaqusInputFilePath = currentDirectory + File.separator + couplingScheme + "_" + (int)(plateHeight/0.0254*1000) + "_PinnedPlate.inp";
		}
		if(pinOrCombBC.equals("comb"))
		{
			abaqusInputFilePath = currentDirectory + File.separator + couplingScheme + "_" + (int)(plateHeight/0.0254*1000) + "_CombedPlate.inp";
		}
		else if(pinOrCombBC.equals("none"))
		{
			abaqusInputFilePath = currentDirectory + File.separator + couplingScheme + "_" + (int)(plateHeight/0.0254*1000) + "_FreePlate.inp";
		}	
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			IMPORTED PARTS NODE */
		// Importing Abaqus fluid part
		String fileLocation = currentDirectory + File.separator;
		String abqFileName = "Star_Fluid_" + (int)(plateHeight/0.0254*1000) + 
			"_" + (int)(smallChannelHeight/0.0254*1000) + 
			"_" + (int)(largeChannelHeight/0.0254*1000);
		//String inputFileLocation = resolvePath(fileLocation + abqFileName);
		ImportCAE cae = new ImportCAE(activeSim, fileLocation, "Fluid");
		cae.importAbaqusInputFile(abqFileName, true, false);
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			REGION NODE */
		// Building the fluid regions
		RegionBuilder fluidRegion = new RegionBuilder(activeSim, "Fluid");
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			PHYSICS NODE */
		// Creating the physics continuum
		ContiuumBuilder physics = new ContiuumBuilder(activeSim);
		PhysicsContinuum fluidPhysics = physics.setPhysicsName("Physics 1", "Water");
		//PhysicsContinuum fluidPhysics = physics.createPhysicsContinua("Water");
		
		fluidPhysics.enable(ThreeDimensionalModel.class);
		fluidPhysics.enable(ImplicitUnsteadyModel.class);
		fluidPhysics.enable(SingleComponentLiquidModel.class);
		
		fluidPhysics.enable(SegregatedFlowModel.class);
		//fluidPhysics.enable(CoupledFlowModel.class);
		
		fluidPhysics.enable(ConstantDensityModel.class);
		fluidPhysics.enable(TurbulentModel.class);
		fluidPhysics.enable(RansTurbulenceModel.class);
		fluidPhysics.enable(KEpsilonTurbulence.class);
		fluidPhysics.enable(RkeTwoLayerTurbModel.class);
		fluidPhysics.enable(KeTwoLayerAllYplusWallTreatment.class);
		fluidPhysics.enable(CoSimulationModel.class);
		fluidPhysics.enable(AbaqusCoSimulationModel.class);
		
		// Setting the initial velocity in all cells
		physics.setInitialConditionsVel(fluidPhysics, initialVel);
		
		// Setting the inlet surfaces/boundaries to velocity inlets and setting the inlet velocity
		fluidRegion.setBoundaryCondition("Fluid.Inlet", "Velocity Inlet", new double[] {0, -1, 0}, inletVel);
		
		// Setting the outlet surfaces/boundaries to pressure outlets
		fluidRegion.setBoundaryCondition("Fluid.Outlet", "Pressure Outlet", new double[] {0, 1, 0}, 0);
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			CO-SIMULATIONS NODE */
		// Creating the Abaqus Co-Simulation
		CoSimulationAbaqus abaqus = new CoSimulationAbaqus(activeSim, "Fluid");
		abaqus.setCouplingBoundaries(new String[] {"Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"});
		abaqus.setAbaqusExecutionSettings(
				"FSI_" + (int)Math.abs(initialVel[1]) + "_" + (int)(wallHeight*10) + "_Abaqus_" + 
					(int)(plateHeight/0.0254*1000) + "_" + (int)(smallChannelHeight/0.0254*1000) + "_" + (int)(largeChannelHeight/0.0254*1000), 
				abaqusInputFilePath, abqExecutableFileName, numAbaqusCPUs);
		abaqus.abaqusCouplingAlgorithm(couplingScheme, "Star Leads", couplingTimeStep);
		abaqus.setFieldExchangeControls(numExchanges, iterationsPerExchange, deflectionUnderRelax);
		abaqus.setMapperTolSettings(0.01,0.01);

		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			SOLVERS NODE */
		SolversNode solvers = new SolversNode(activeSim);
		solvers.setKepsilonRelax(0.6);
		solvers.setUnsteadyTimeStep(couplingTimeStep, 1);
		//solvers.setCourantNumber(courantNumber);
		
		// Setting up the mesh morpher solver
		MeshMorpher morpher =  new MeshMorpher(activeSim, "Fluid");
		morpher.addRegionBoundary("Fluid.FSI_Back", "FSI");
		morpher.addRegionBoundary("Fluid.FSI_Front", "FSI");
		morpher.addRegionBoundary("Fluid.FSI_Bottom", "FSI");
		morpher.addRegionBoundary("Fluid.FSI_Top", "FSI");
		if(morphAtInnIter.equals("yes"))
		{
			morpher.innerIterationMorphing(true);
		}
		else if(morphAtInnIter.equals("no"))
		{
			morpher.innerIterationMorphing(false);
		}
		
		if(SSorFSI.equals("SS"))
		{
			solvers.setSSorFSI(true);
		}
		else if(SSorFSI.equals("FSI"))
		{
			solvers.setSSorFSI(false);
		}
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			STOPPING CRITERIA NODE */
		StoppingCriteria stoppingCriteria = new StoppingCriteria(activeSim);
		stoppingCriteria.innerIterationStoppingCriteriaController(iterationsPerTS, "OR", true);
		stoppingCriteria.maxPhysicalTime(maxSimTime, "OR", true);
		
		// Creating stopping criteria for the static pressure on the plate's wall
		ReportsMonitorsPlots reports = new ReportsMonitorsPlots(activeSim);
		AreaAverageReport avgPressure = reports.createAverageReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "AvgPressure");
		MaxReport maxPressure = reports.createMaxReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "MaxPressure");
		MinReport minPressure = reports.createMinReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "MinPressure");
		
		FieldFunctions fieldFunctions = new FieldFunctions(activeSim);
		PrimitiveFieldFunction staticPressure = fieldFunctions.getStaticPressureFunction();
		avgPressure.setScalar(staticPressure);
		maxPressure.setScalar(staticPressure);
		minPressure.setScalar(staticPressure);
		
		reports.createMonitorPlot2(new String[] {"MinPressure", "MaxPressure", "AvgPressure"}, 
				"Avg, Max, and Min Plate Pressure", new String[] {"Iteration", "Static Pressure (Pa)"}, true);
		
		stoppingCriteria.createAsymStoppingCriteria("AvgPressure Monitor", "Report", 1.0, 25, "AND", "Inner", fluidStopCriteria);
		stoppingCriteria.createAsymStoppingCriteria("MaxPressure Monitor", "Report", 1.0, 25, "AND", "Inner", fluidStopCriteria);
		stoppingCriteria.createAsymStoppingCriteria("MinPressure Monitor", "Report", 1.0, 25, "AND", "Inner", fluidStopCriteria);
		
		// Creating stopping criteria for the wall shear stres on the plate's wall
		AreaAverageReport avgStress = reports.createAverageReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "AvgWallShearStress");
		MaxReport maxStress = reports.createMaxReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "MaxWallShearStress");
		MinReport minStress = reports.createMinReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "MinWallShearStress");
		
		PrimitiveFieldFunction wallShearStress = fieldFunctions.getWallShearStress();
		VectorMagnitudeFieldFunction vectorMagnitudeFieldFunction = ((VectorMagnitudeFieldFunction) wallShearStress.getMagnitudeFunction());
		avgStress.setScalar(vectorMagnitudeFieldFunction);
		maxStress.setScalar(vectorMagnitudeFieldFunction);
		minStress.setScalar(vectorMagnitudeFieldFunction);
		
		reports.createMonitorPlot2(new String[] {"MinWallShearStress", "MaxWallShearStress", "AvgWallShearStress"}, 
				"Avg, Max, and Min Plate Wall Shear Stress", new String[] {"Iteration", "Static Pressure (Pa)"}, true);
		
		stoppingCriteria.createAsymStoppingCriteria("AvgWallShearStress Monitor", "Report", 1.0, 25, "AND", "Inner", fluidStopCriteria);
		stoppingCriteria.createAsymStoppingCriteria("MaxWallShearStress Monitor", "Report", 1.0, 25, "AND", "Inner", fluidStopCriteria);
		stoppingCriteria.createAsymStoppingCriteria("MinWallShearStress Monitor", "Report", 1.0, 25, "AND", "Inner", fluidStopCriteria);
		
		// Creating stopping criteria for the plate deflection 
		AreaAverageReport avgDeflection = reports.createAverageReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "AvgDeflection");
		MaxReport maxDeflection = reports.createMaxReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "MaxDeflection");
		MinReport minDeflection = reports.createMinReport(
				new String[] {"Fluid", "Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"}, "MinDeflection");
		
		PrimitiveFieldFunction plateDeflection = fieldFunctions.getNodalDisplacement();
		VectorMagnitudeFieldFunction plateDeflection_Mag = ((VectorMagnitudeFieldFunction) plateDeflection.getMagnitudeFunction());
		avgDeflection.setScalar(plateDeflection_Mag);
		maxDeflection.setScalar(plateDeflection_Mag);
		minDeflection.setScalar(plateDeflection_Mag);
		
		reports.createMonitorPlot2(new String[] {"MinDeflection", "MaxDeflection", "AvgDeflection"}, "Avg, Max, and Min PlateDeflection",
				new String[] {"Iteration", "Plate Deflection (m)"}, true);
		
		stoppingCriteria.createAsymStoppingCriteria("AvgDeflection Monitor", "Report", 1.27e-6, 5, "AND", "Outer", fluidStopCriteria);
		stoppingCriteria.createAsymStoppingCriteria("MaxDeflection Monitor", "Report", 1.27e-6, 5, "AND", "Outer", fluidStopCriteria);
		stoppingCriteria.createAsymStoppingCriteria("MinDeflection Monitor", "Report", 1.27e-6, 5, "AND", "Outer", fluidStopCriteria);
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			DERIVED PARTS NODE */
		DerivedParts centerPlane = new DerivedParts(activeSim, new String[] {"Fluid"});
		if(plateGeometry.equals("Flat"))
		{
			centerPlane.createSectionPlane(new double[] {1, 0, 0}, new double[] {wettedPlateWidth * 0.5, 0, 0}, "CenterPlane");
		}
		else if(plateGeometry.equals("Curved"))
		{		
			centerPlane.createSectionPlane(new double[] {1, 0, -1}, new double[] {0.0, 0, 0.0}, "CenterPlane");
		}
		
		// Creating line probes throughout the model for plotting the pressure profile through the entire model
		String[] lineProbeRegions = {"Fluid"};
		DerivedParts smChLineProbe = new DerivedParts(activeSim, lineProbeRegions);
		LinePart smChLinePart = null;
		LinePart lgChLinePart = null;
		DerivedParts lgChLineProbe = new DerivedParts(activeSim, lineProbeRegions);
		if(plateGeometry.equals("Flat"))
		{		
			double[] lgChLineProbeCoord_0 = {wettedPlateWidth*0.5, -outletLength, -largeChannelHeight*0.5};
			double[] lgChLineProbeCoord_1 = {wettedPlateWidth*0.5, plateLength + inletLength, -(largeChannelHeight*0.5)};
			smChLinePart = smChLineProbe.createLineProbe(lgChLineProbeCoord_0, lgChLineProbeCoord_1, 255, "LargeChannelLineProbe");
			
			double[] smChLineProbeCoord_0 = {wettedPlateWidth*0.5, -outletLength, smallChannelHeight*0.5 + plateHeight};
			double[] smChLineProbeCoord_1 = {wettedPlateWidth*0.5, plateLength + inletLength, smallChannelHeight*0.5 + plateHeight};
			lgChLinePart = lgChLineProbe.createLineProbe(smChLineProbeCoord_0, smChLineProbeCoord_1, 255, "SmallChannelLineProbe");
		}
		else if(plateGeometry.equals("Curved"))
		{		
			double[] lgChLineProbeCoord_0 = {(r_in - largeChannelHeight*0.5)*Math.cos(Math.PI/4), -outletLength, (r_in - largeChannelHeight*0.5)*Math.cos(Math.PI/4)};
			double[] lgChLineProbeCoord_1 = {r_in*Math.cos(Math.PI/4), plateLength + inletLength, r_in*Math.cos(Math.PI/4)};
			smChLinePart = smChLineProbe.createLineProbe(lgChLineProbeCoord_0, lgChLineProbeCoord_1, 255, "LargeChannelLineProbe");
			
			double[] smChLineProbeCoord_0 = {(r_in + plateHeight + smallChannelHeight*0.5)*Math.cos(Math.PI/4), -outletLength, (r_in + plateHeight + smallChannelHeight*0.5)*Math.cos(Math.PI/4)};
			double[] smChLineProbeCoord_1 = {(r_in + plateHeight + smallChannelHeight*0.5)*Math.cos(Math.PI/4), plateLength + inletLength, (r_in + plateHeight + smallChannelHeight*0.5)*Math.cos(Math.PI/4)};
			lgChLinePart = lgChLineProbe.createLineProbe(smChLineProbeCoord_0, smChLineProbeCoord_1, 255, "SmallChannelLineProbe");
		}
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
		 	PLOTS NODE */
		FieldFunctions fieldFunction = new FieldFunctions(activeSim);
		
		// Turning off the "Auto" normalization option for all the residual monitors
		ReportsMonitorsPlots reportsMonitorsPlots = new ReportsMonitorsPlots(activeSim);
		reportsMonitorsPlots.residualNormalization(new String[] {"Continuity", "Tdr", "Tke", "X-momentum", "Y-momentum", "Z-momentum"});
		
		// Creating an XY plot of the pressure profiles throughout the model
		ReportsMonitorsPlots pressureProfilePlot = new ReportsMonitorsPlots(activeSim);
		XYPlot pressureProfile_XYPlot = pressureProfilePlot.createXYPlot(new double[] {0, 1, 0}, "PressureProfiles", "Static Pressure (Pa)");
		fieldFunction.setXYPlotFieldFunction(pressureProfile_XYPlot, "StaticPressure", "0");
		pressureProfile_XYPlot.getParts().addObjects(smChLinePart);
		pressureProfile_XYPlot.getParts().addObjects(lgChLinePart);
		
		// Creating an XY plot of the plate's wall y+ values
		ReportsMonitorsPlots wallYplusPlot = new ReportsMonitorsPlots(activeSim);
		XYPlot wallYplus_XYPlot = wallYplusPlot.createXYPlot(new double[] {0, 1, 0}, "Plate Wall y+ Values", "Wall y+ Values");
		fieldFunction.setXYPlotFieldFunction(wallYplus_XYPlot, "WallYplus", "0");
		wallYplusPlot.addObjects2XYPlot(wallYplus_XYPlot, "Fluid", 
				new String[] {"Fluid.FSI_Back", "Fluid.FSI_Front"});
				
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
 			SCENES NODE */
		// Creating a scene of pressure
		Scenes deflectionScene = new Scenes(activeSim, "Deflection");
		Scene deflection_Scene = deflectionScene.createScalarScene();
		//fieldFunction.setSceneFieldFunction(pressure_Scene, "Morpher Displacement", "Magnitude");
		PrimitiveFieldFunction nodalDisplacement = fieldFunction.getNodalDisplacement();
		deflectionScene.setSceneFieldFunction(nodalDisplacement, 4);
		deflectionScene.addObject2Scene(deflection_Scene, "Fluid", new String[] {"Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"});
		deflectionScene.addDerivedPart2Scene(deflection_Scene, new String[] {"CenterPlane"});
		
		// Creating a scene of velocity on the "CenterPlane"
		Scenes velocityScene = new Scenes(activeSim, "Velocity");
		Scene velocity_Scene = velocityScene.createScalarScene();
		fieldFunction.setSceneFieldFunction(velocity_Scene, "Velocity", "1");
		velocityScene.addDerivedPart2Scene(velocity_Scene, new String[] {"CenterPlane"});
		
		// Creating a scene of wall y+ values on the region "Fluid"
		Scenes wallYScene = new Scenes(activeSim, "WallY+");
		Scene wallY_Scene = wallYScene.createScalarScene();
		fieldFunction.setSceneFieldFunction(wallY_Scene, "WallYplus", "0");
		wallYScene.addObject2Scene(wallY_Scene, "Fluid", new String[] {"Fluid.FSI_Back", "Fluid.FSI_Front", "Fluid.FSI_Bottom", "Fluid.FSI_Top"});
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			SOLUTION HISTORY NODE */
			/*
		String simhFileLocation =  currentDirectory + File.separator + couplingScheme + "_FSI_" + (int)Math.abs(initialVel[1]) + "_StarCCM_" + 
				(int)(plateHeight/0.0254*1000) + "_" + (int)(smallChannelHeight/0.0254*1000) + "_" + (int)(largeChannelHeight/0.0254*1000) + 
				"SolutionHistory.simh";
		String[] scalarFieldFunctions = {"StaticPressure", "Volume"};
		String[] vectorFieldFunctions = {"NodalDisplacement", "Morpher Displacement"};
		SolutionHistoryCreator solutionHistory = new SolutionHistoryCreator(activeSim, simhFileLocation);
		solutionHistory.addScalarFieldFunction(scalarFieldFunctions);
		solutionHistory.addVectorFieldFunction(vectorFieldFunctions);
		solutionHistory.setUpdateSettings("Time Step", 1);
		
		*/
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			TOOLS NODE */
		//Tools tools = new Tools(activeSim);
		//tools.createXYZInternalTable("Fluid", new String[] {"FSI_Back", "FSI_Front"});
		//tools.setXYZInternalTableFieldFunction("StaticPressure");
		//tools.extractAndExportXYZInternalTableData(currentDirectory + File.separator + "XYZ_Pressure.csv");
		
		/**-----------------------------------------------------------------------------------------------------------------------------------------------------
			SAVING AND RUNNING NODE */		
		// Saving and running the simulation
		String saveLocation = currentDirectory + File.separator + 
			"FSI_" + (int)Math.abs(initialVel[1]) + "_" +  (int)(wallHeight*10) + "_" +
				plateGeometry + '_' + (int)(plateHeight/0.0254*1000) + "_" + (int)(smallChannelHeight/0.0254*1000) + "_" + (int)(largeChannelHeight/0.0254*1000) +".sim";
		activeSim.saveState(saveLocation);
		AbaqusCoSimulation coSim = (AbaqusCoSimulation) activeSim.get(CoSimulationManager.class).getCoSimulation("Abaqus Co-Simulation 1");
		//activeSim.getSimulationIterator().run();
		try {
			coSim.terminate();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//activeSim.saveState(saveLocation);
		activeSim.close();
	}

}
