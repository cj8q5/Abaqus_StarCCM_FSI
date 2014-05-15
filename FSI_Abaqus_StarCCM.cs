using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;

/// <summary>
/// This class is used for creating and running FSI models using Star-CCM+ and Abaqus in either parametric mode 
/// or without
/// </summary>
public class FSI_Abaqus_StarCCM
{
    private static Hashtable m_doubleParameters = new Hashtable();
    private static Hashtable m_stringParameters = new Hashtable();
    private static Hashtable m_intParameters = new Hashtable();

    static void Main()
    {
        string currentDirectory = Directory.GetCurrentDirectory();
        readGeometryData(Path.Combine(currentDirectory, "FSI_Input_File.txt"));
        
        // Running the Abaqus Python script for building the plate and fluid models
        if (getDoubleData("smChHeight") != getDoubleData("lgChHeight") && getIntData("numOfPlates") > 1)
        {
            Console.WriteLine("\n Plate stack must have equal small and large channel heights! \n");
        }
        else
        {
            Console.WriteLine("\n Abaqus is building the fluid and solid domain geomtries...\n");
            string buildAbaqusCall = "/C abaqus cae noGUI=FSI_GeometryBuilder.py";
            var abaqusProcess = Process.Start("cmd.exe", buildAbaqusCall);
            abaqusProcess.WaitForExit();
            Console.WriteLine("\n Abaqus has finished building the fluid and solid domains!\n");

            // Running Star-CCM+ to build the fluid model and setup the FSI problem
            if (getStringData("createStarFile").Equals("yes"))
            {
                Console.WriteLine("\n Star-CCM+ is now building the fluid model and setting up the FSI problem...\n");
                string buildStarCall = "/C starccm+ -new -np 1 -batch AbaqusMeshingFSI.java";
                var starProcess = Process.Start("cmd.exe", buildStarCall);
                starProcess.WaitForExit();
                Console.WriteLine("\n Star-CCM+ has finished building the fluid model and the FSI problem!\n");
            }

            // Running the FSI simulation
            if (getStringData("runStar").Equals("yes"))
            {
                string couplingScheme = getStringData("couplingScheme");
                string numStarProcesses = getStringData("starProcesses");
                double vel = getDoubleData("avgChVelocity");
                string plateGeometry = getStringData("plateGeometry");
                int intPlateThickness = (int)(getDoubleData("plateThickness") / 0.0254 * 1000);
                int intSmChHeight = (int)(getDoubleData("smChHeight") / 0.0254 * 1000);
                int intLgHeight = (int)(getDoubleData("lgChHeight") / 0.0254 * 1000);

                string runStarCall = "/C starccm+ -np " + numStarProcesses + "-time -batch " + "FSI_" +
                    vel + "_" + plateGeometry + "_" + intPlateThickness + "_" + intSmChHeight + "_" +
                    intLgHeight + ".sim";
                var runStarProcess = Process.Start("cmd.exe", runStarCall);
                runStarProcess.WaitForExit();
            }
        }
    }

    /**
     * Method for reading data from an input file
     * 
     */
	public static void readGeometryData(string file2Read)
	{
        try
        {
            StreamReader file = new StreamReader(file2Read);
            string line;
            while((line = file.ReadLine()) != null)
            {
                if(line.StartsWith("#"))
                {
                    continue;
                }
                else if(line.Length == 0)
                {
                    continue;
                }
                else
                {
                    string[] things2Replace = new string[] {"\t", " "};
                    for (int i = 0; i <= 1; i++)
                    {
                        line = line.Replace(things2Replace[i], String.Empty);
                    }

                    string[] stringSeparators = new string[] { ":" };
                    string[] lineParameters = line.Split(stringSeparators, StringSplitOptions.RemoveEmptyEntries);

                    if(lineParameters[1].Equals("float"))
                    {
                        m_doubleParameters.Add(lineParameters[0], Double.Parse(lineParameters[2]));
                    }
                    else if(lineParameters[1].Equals("string"))
                    {
                        m_stringParameters.Add(lineParameters[0], lineParameters[2]);
                    }
                    else if(lineParameters[1].Equals("integer"))
                    {
                        m_intParameters.Add(lineParameters[0], int.Parse(lineParameters[2]));
                    }
                }
            }
            file.Close();
        }
        catch(FileNotFoundException e)
        {
            Console.WriteLine(e);
        }
	}

    /**
     * Getter methods for getting the double, string, and integer data
     */
    public static double getDoubleData(string variableName)
    {
        double variable = (double)m_doubleParameters[variableName];
        return variable;
    }

    public static string getStringData(string variableName)
    {
        string variable = (string)m_stringParameters[variableName];
        return variable;
    }

    public static int getIntData(string variableName)
    {
        int variable = (int)m_intParameters[variableName];
        return variable;
    }

    /*
     * This method changes a variable in the input file
     */
    public static void changeInputParameters(string fileName, string parameter2Change, string newValue)
    {
        List<string> inputFileLines = new List<string>();
        using (StreamReader file = new StreamReader(fileName))
        {
            string line;
            while ((line = file.ReadLine()) != null)
            {
                if (line.StartsWith("#") || line.Length == 0)
                {
                    inputFileLines.Add(line);
                }
                else if (line.StartsWith(parameter2Change))
                {
                    string[] stringSeparators = new string[] { ":" };
                    string[] lineParameters = line.Split(stringSeparators, StringSplitOptions.RemoveEmptyEntries);
                    string newLine = lineParameters[0] + ":" + lineParameters[1] + ": \t" + newValue + ":" + lineParameters[3];
                    inputFileLines.Add(newLine);
                }
                else
                {
                    inputFileLines.Add(line);
                }
            }
            file.Close();
        }
        
        // Rewriting the input file
        using (StreamWriter file = new StreamWriter(fileName))
        {
            foreach (string line in inputFileLines)
            {
                file.Write(line +"\r\n");
            }
        }
    }
}