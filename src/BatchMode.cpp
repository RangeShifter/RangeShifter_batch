/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell
 *
 *	This file is part of RangeShifter.
 *
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *
 --------------------------------------------------------------------------*/


 //---------------------------------------------------------------------------

#include "BatchMode.h"
//---------------------------------------------------------------------------

// Note - all batch files are prefixed 'b' here for reasons concerned with RS v1.0
ifstream ifsSimFile, ifsParamFile, ifsLandFile, ifsSpLandFile, ifsDynLandFile;
ifstream ifsSpDistFile, ifsStageStructFile, ifsTransMatrix;
ifstream ifsStageWeightsFile;
ifstream ifsEmigrationFile, ifsTransferFile, ifsSettlementFile;
ifstream ifsTraitsFile, ifsGeneticsFile;
ifstream ifsInitFile, ifsInitIndsFile;
ifstream ifsFecDens, ifsDevDens, ifsSurvDens;

ofstream batchLogOfs;

// global variables passed between parsing functions...
// should be removed eventually, maybe share variables through members of a class
int gUsesPatches, gUsesStageStruct, gResol;
int gTransferType, gLandType, gMaxNbHab;
bool gAnyUsesGenetics;
int gNbLandscapes = 0;

int gEnvStochType;
bool gStochInK;
map<int, int> gNbReplicates;
set<species_id> gSpeciesNames;
map<species_id, bool> gUseSpeciesDist;
map<species_id, bool> gUseSMSCosts;

// sim x species grid of parameters to check coherency between input files
map<int, map<species_id, spInputOptions>> gSpInputOpt;

rasterdata landRaster;
// ...including names of the input files
string gSimFile, gParametersFile;
string landFile;
string gHabMapName, gDynLandFileName;
string gSpLandName;
string stageStructFile, transMatrix;
string emigrationFile, transferFile, settleFile, geneticsFile, traitsFile, initialFile;
string prevInitialIndsFile = " ";

const string gNbLinesStr = "No. of lines for final Simulation ";
const string gShouldBeStr = " should be ";
const string gResolOfStr = "*** Resolution of ";
const string gResolNotMatchStr = " does not match Resolution in Control file ";
const string gHeadersOfStr = "*** Headers of ";
const string gHeadersNotMatchStr = " do not match headers of LandscapeFile";
const string gPatchReqdStr = " is required for patch-based model";
const string gSpecMustMatchStr = " must match the specification exactly";
const string gCaseSensitiveStr = " case-sensitive parameter names";

float** gMatrix = nullptr;	// temporary matrix used in batch mode
int gMatrixSize = 0; 		// size of temporary matrix

//---------------------------------------------------------------------------
// Fractal dimensions must be some power of two, plus one
bool isValidFractalDim(int x) {
	if (x < 2) return false;
	int r = x % 2;
	while (r == 0) {
		x /= 2; 
		r = x % 2;
	}
	return x == 1;
}

//---------------------------------------------------------------------------
bool checkInputFiles(string pathToControlFile, string inputDir, string outputDir)
{
	int lines, nSimuls;
	int nbErrors = 0;
	string paramName, filename, pathToFile, batchLogPath, header;
	string whichInputFile = "Control file";
	bool anyFormatError = false;
	
	// Open batch log
	batchLogPath = outputDir + "BatchLog.txt";
	batchLogOfs.open(batchLogPath.c_str());
	if (!batchLogOfs.is_open()) {
		cout << "Error opening batch output log file " << batchLogPath << endl;
		return false;
	}

	// Open control file
	ifstream controlIfs{ pathToControlFile.c_str() };
	if (!controlIfs.is_open()) {
		cout << "Error opening Control file: " << pathToControlFile << endl;
		batchLogOfs << "Error opening Control file: " << pathToControlFile << endl;
		if (batchLogOfs.is_open()) { 
			batchLogOfs.close(); 
			batchLogOfs.clear(); 
		}
		return false;
	}
	else batchLogOfs << "Checking Control file " << pathToControlFile << endl;

	// Check batch parameters
	int batchNb;
	controlIfs >> paramName >> batchNb;
	if (paramName == "BatchNum") {
		if (batchNb < 0) {
			BatchError(whichInputFile, -999, 19, "BatchNum");
			nbErrors++;
		}
		else paramsSim->setBatchNum(batchNb);
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gUsesPatches;
	if (paramName == "PatchModel") {
		if (gUsesPatches != 0 && gUsesPatches != 1) {
			BatchError(whichInputFile, -999, 1, "PatchModel"); 
			nbErrors++;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gResol;
	if (paramName == "Resolution") {
		if (gResol < 1) {
			BatchError(whichInputFile, -999, 11, "Resolution");
			nbErrors++;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gLandType;
	if (paramName == "LandType") {
		if (gLandType != 0 && gLandType != 2 && gLandType != 9) {
			BatchError(whichInputFile, -999, 0, "LandType");
			batchLogOfs << "LandType must be 0, 2 or 9" << endl;
			nbErrors++;
		}
		else {
			if (gLandType == 9 && gUsesPatches) {
				BatchError(whichInputFile, -999, 0, "LandType");
				batchLogOfs << "LandType may not be 9 for a patch-based model" << endl;
				nbErrors++;
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gMaxNbHab;
	if (paramName == "MaxHabitats") {
		if (gLandType == 0) { // raster with unique habitat codes
			if (gMaxNbHab < 2) {
				BatchError(whichInputFile, -999, 12, "MaxHabitats"); nbErrors++;
			}
		}
		else if (gMaxNbHab != 1) { // habitat quality or artificial landscape
			BatchError(whichInputFile, -999, 0, " "); nbErrors++;
			batchLogOfs << "MaxHabitats must be 1 for LandType = " << gLandType << endl;
		}
	}
	else anyFormatError = true;

	controlIfs >> paramName >> gUsesStageStruct;
	if (paramName == "StageStruct") {
		if (gUsesStageStruct != 0 && gUsesStageStruct != 1) {
			BatchError(whichInputFile, -999, 1, "StageStruct"); 
			nbErrors++;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gTransferType;
	if (paramName == "Transfer") {
		if (gTransferType < 0 || gTransferType > 2) {
			BatchError(whichInputFile, -999, 2, "Transfer"); 
			nbErrors++;
		}
	}
	else anyFormatError = true; // wrong control file format

	if (anyFormatError || nbErrors > 0) { // terminate batch error checking
		if (anyFormatError) printControlFormatError();
		batchLogOfs << endl
			<< "*** Model parameters in Control file must be corrected before further input file checks are conducted"
			<< endl;
		batchLogOfs.close(); 
		batchLogOfs.clear();
		controlIfs.close(); 
		controlIfs.clear();
		return false;
	}

	bool areInputFilesOk = true;

	// Check simulation file
	controlIfs >> paramName >> filename;
	if (paramName == "SimFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsSimFile.open(pathToFile.c_str());
		if (ifsSimFile.is_open()) {
			if (!CheckSimFile())
				areInputFilesOk = false;
			else {
				FileOK(paramName, gSpInputOpt.size(), 0);
				gSimFile = pathToFile;
			}
			ifsSimFile.close();
		}
		else {
			OpenError(paramName, pathToFile);
			areInputFilesOk = false;
			cout << "Unable to open SimFile" << endl;
		}
		ifsSimFile.clear();
		if (!areInputFilesOk) {
			batchLogOfs << endl
				<< "*** SimFile must be corrected before further input file checks are conducted"
				<< endl;
			batchLogOfs.close();
			batchLogOfs.clear();
			controlIfs.close();
			controlIfs.clear();
			return false;
		}
	}
	else anyFormatError = true; // wrong control file format
	if (ifsSimFile.is_open())
		ifsSimFile.close();
	ifsSimFile.clear();

	// Check parameter file
	controlIfs >> paramName >> filename;
	if (paramName == "ParameterFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsParamFile.open(pathToFile.c_str());
		if (ifsParamFile.is_open()) {
			if (!CheckParameterFile())
				areInputFilesOk = false;
			else {
				FileOK(paramName, gSpInputOpt.size(), 0);
				gParametersFile = pathToFile;
			}
			ifsParamFile.close();
		}
		else {
			OpenError(paramName, pathToFile); 
			areInputFilesOk = false;
			cout << "Unable to open ParameterFile" << endl;
		}
		ifsParamFile.clear();
		if (!areInputFilesOk) {
			batchLogOfs << endl
				<< "*** ParameterFile must be corrected before further input file checks are conducted"
				<< endl;
			batchLogOfs.close(); 
			batchLogOfs.clear();
			controlIfs.close(); 
			controlIfs.clear();
			return false;
		}
	}
	else anyFormatError = true; // wrong control file format
	if (ifsParamFile.is_open()) 
		ifsParamFile.close();
	ifsParamFile.clear();

	// Check land file
	controlIfs >> paramName >> filename;
	if (paramName == "LandFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsLandFile.open(pathToFile.c_str());
		if (ifsLandFile.is_open()) {
			if (!CheckLandFile(gLandType, inputDir)) {
				areInputFilesOk = false;
				batchLogOfs << "*** Format error in " << paramName << endl;
			}
			else {
				FileOK(paramName, gNbLandscapes, 1);
				landFile = pathToFile;
			}
			ifsLandFile.close();
		}
		else {
			OpenError(paramName, pathToFile); 
			areInputFilesOk = false;
		}
		ifsLandFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check stage structure file if required file
	controlIfs >> paramName >> filename;
	batchLogOfs << endl;
	if (paramName == "StageStructFile" && !anyFormatError) {
		if (filename == "NULL") {
			if (gUsesStageStruct) {
				batchLogOfs << "*** File name is required for " << paramName << endl;
				areInputFilesOk = false;
			}
		}
		else { // filename is not NULL
			if (gUsesStageStruct) { // check file only if it is required
				pathToFile = inputDir + filename;
				batchLogOfs << "Checking " << paramName << " " << pathToFile << endl;
				ifsStageStructFile.open(pathToFile.c_str());
				if (ifsStageStructFile.is_open()) {
					nSimuls = CheckStageFile(inputDir);
					if (nSimuls < 0) areInputFilesOk = false;
					else {
						FileOK(paramName, nSimuls, 0);
						if (nSimuls != gSpInputOpt.size()) {
							SimulnCountError(filename); 
							areInputFilesOk = false;
						}
						else stageStructFile = pathToFile;
					}
					ifsStageStructFile.close();
				}
				else {
					OpenError(paramName, pathToFile);
					areInputFilesOk = false;
				}
				ifsStageStructFile.clear();
			} // end of required
			else { // file is not required, and filename should be NULL
				if (filename != "NULL") {
					batchLogOfs << "*** File name for stageStructFile should be NULL as StageStruct = "
						<< gUsesStageStruct << endl;
					areInputFilesOk = false;
				}
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	// Check emigration file
	controlIfs >> paramName >> filename;
	if (paramName == "EmigrationFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsEmigrationFile.open(pathToFile.c_str());
		if (ifsEmigrationFile.is_open()) {
			nSimuls = CheckEmigFile();
			if (nSimuls < 0) 
				areInputFilesOk = false;
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != gSpInputOpt.size()) {
					SimulnCountError(filename); 
					areInputFilesOk = false;
				}
				else emigrationFile = pathToFile;
			}
			ifsEmigrationFile.close();
		}
		else {
			OpenError(paramName, pathToFile); 
			areInputFilesOk = false;
		}
		ifsEmigrationFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check transfer file
	controlIfs >> paramName >> filename;
	if (paramName == "TransferFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsTransferFile.open(pathToFile.c_str());
		if (ifsTransferFile.is_open()) {
			nSimuls = CheckTransferFile(inputDir);
			if (nSimuls < 0) {
				areInputFilesOk = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != gSpInputOpt.size()) {
					SimulnCountError(filename); 
					areInputFilesOk = false;
				}
				else transferFile = pathToFile;
			}
			ifsTransferFile.close(); 
			ifsTransferFile.clear();
		}
		else {
			OpenError(paramName, pathToFile); 
			areInputFilesOk = false;
		}
		ifsTransferFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check settlement file
	controlIfs >> paramName >> filename;
	if (paramName == "SettlementFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsSettlementFile.open(pathToFile.c_str());
		if (ifsSettlementFile.is_open()) {
			nSimuls = CheckSettleFile();
			if (nSimuls < 0) {
				areInputFilesOk = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != gSpInputOpt.size()) {
					SimulnCountError(filename); 
					areInputFilesOk = false;
				}
				else settleFile = pathToFile;
			}
			ifsSettlementFile.close();
		}
		else {
			OpenError(paramName, pathToFile); 
			areInputFilesOk = false;
		}
		ifsSettlementFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check genetics file if required file
	controlIfs >> paramName >> filename;
	batchLogOfs << endl;
	if (paramName == "GeneticsFile" && !anyFormatError) {
		if (filename == "NULL") {
			gAnyUsesGenetics = false;
			for (auto const& [simNb, simOpt] : gSpInputOpt) {
				for (auto const& [sp, spOpt] : simOpt) {
					if (spOpt.isEmigIndVar) {
						batchLogOfs << "Error: GeneticsFile is NULL but Emigration is set to IndVar for species " 
							<< to_string(sp) << " in simulation " << to_string(simNb) << endl;
						areInputFilesOk = false;
					}
					if (spOpt.isSettIndVar) {
						batchLogOfs << "Error: GeneticsFile is NULL but Settlement is set to IndVar for species "
							<< to_string(sp) << " in simulation " << to_string(simNb) << endl;
						areInputFilesOk = false;
					}
					if (spOpt.isKernTransfIndVar) {
						batchLogOfs << "Error: GeneticsFile is NULL but Transfer is set to IndVar for species "
							<< to_string(sp) << " in simulation " << to_string(simNb) << endl;
						areInputFilesOk = false;
					}
					if (spOpt.isSMSTransfIndVar) {
						batchLogOfs << "Error: GeneticsFile is NULL but Transfer is set to IndVar for species "
							<< to_string(sp) << " in simulation " << to_string(simNb) << endl;
						areInputFilesOk = false;
					}
				}
			}
		}
		else {
			gAnyUsesGenetics = true;
			pathToFile = inputDir + filename;
			batchLogOfs << "Checking " << paramName << " " << pathToFile << endl;
			ifsGeneticsFile.open(pathToFile.c_str());
			if (ifsGeneticsFile.is_open()) {
				nSimuls = CheckGeneticsFile(inputDir);
				if (nSimuls < 0) {
					areInputFilesOk = false;
				}
				else {
					FileOK(paramName, nSimuls, 0);
					geneticsFile = pathToFile;
				}
				ifsGeneticsFile.close();
			}
			else {
				OpenError(paramName, pathToFile); 
				areInputFilesOk = false;
			}
			ifsGeneticsFile.clear();
		}
	}
	else anyFormatError = true; // wrong control file format

	// Check TraitsFile
	controlIfs >> paramName >> filename;
	batchLogOfs << endl;

	if (paramName == "TraitsFile" && !anyFormatError) {
		if (filename == "NULL") {
			if (gAnyUsesGenetics) {
				batchLogOfs << "Error: Genetics are enabled but no TraitsFile is provided." << endl;
				areInputFilesOk = false;
			}
		}
		else {
			pathToFile = inputDir + filename;
			batchLogOfs << "Checking " << paramName << " " << pathToFile << endl;
			ifsTraitsFile.open(pathToFile.c_str());
			if (ifsTraitsFile.is_open()) {
				nSimuls = CheckTraitsFile(inputDir);
				if (nSimuls < 0) {
					areInputFilesOk = false;
				}
				else {
					FileOK(paramName, nSimuls, 0);
					traitsFile = pathToFile;
				}
				ifsTraitsFile.close();
			}
			else {
				OpenError(paramName, filename);
				areInputFilesOk = false;
			}
			if (ifsTraitsFile.is_open()) ifsTraitsFile.close();
			ifsTraitsFile.clear();
		}
	}

	// Check initialisation file
	controlIfs >> paramName >> filename;
	if (paramName == "InitialisationFile" && !anyFormatError) {
		pathToFile = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << pathToFile << endl;
		ifsInitFile.open(pathToFile.c_str());
		if (ifsInitFile.is_open()) {
			nSimuls = CheckInitFile(inputDir);
			if (nSimuls < 0) {
				areInputFilesOk = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != gSpInputOpt.size()) {
					SimulnCountError(filename); 
					areInputFilesOk = false;
				}
				else initialFile = pathToFile;
			}
			ifsInitFile.close();
		}
		else {
			OpenError(paramName, pathToFile);
			areInputFilesOk = false;
		}
		ifsInitFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	if (anyFormatError) {
		printControlFormatError();
		areInputFilesOk = false;
	}

	if (controlIfs.is_open()) { controlIfs.close(); controlIfs.clear(); }
	if (batchLogOfs.is_open()) { batchLogOfs.close(); batchLogOfs.clear(); }

	return areInputFilesOk;
}

bool CheckSimFile() {
	int nbErrors = 0;
	string header;
	ifsSimFile >> header; if (header != "Simulation") nbErrors++;
	ifsSimFile >> header; if (header != "Replicates") nbErrors++;
	ifsSimFile >> header; if (header != "Years") nbErrors++;
	ifsSimFile >> header; if (header != "Absorbing") nbErrors++;
	ifsSimFile >> header; if (header != "FixReplicateSeed") nbErrors++;
	ifsSimFile >> header; if (header != "EnvStoch") nbErrors++;
	ifsSimFile >> header; if (header != "EnvStochType") nbErrors++;
	ifsSimFile >> header; if (header != "ac") nbErrors++;
	ifsSimFile >> header; if (header != "std") nbErrors++;

	string whichFile = "SimFile";
	if (nbErrors > 0) {
		FormatError(whichFile, nbErrors);
		batchLogOfs << "*** SimFile column headers are incorrect." << endl;
		return false;
	}

	// Parse data lines
	int whichLine = 1;
	int nbSims = 0, prevSim;
	const int errSimNb = -98765;
	int simNb = errSimNb;
	ifsSimFile >> simNb; // first simulation number
	if (simNb == errSimNb) {
		batchLogOfs << "*** Error in SimFile - first simulation number could not be read." << endl;
		nbErrors++;
	}
	else if (simNb < 0) {
		batchLogOfs << "*** Error in SimFile - first simulation number must be >= 0" << endl;
		nbErrors++;
	}
	else {
		prevSim = simNb;
		nbSims++;
	}

	int inReplicates, inYears, inAbsorb, inFixReplicateSeed, inStochInK;
	float inStochAC, inStochStD;

	while (simNb != -98765) {

		// Initialise input option map with simulation numbers
		gSpInputOpt.emplace(simNb, map<species_id, spInputOptions>());

		ifsSimFile >> inReplicates;
		if (inReplicates <= 0) {
			BatchError(whichFile, whichLine, 11, "Replicates");
			nbErrors++;
		}
		else {
			gNbReplicates.emplace(simNb, inReplicates);
		}
		ifsSimFile >> inYears;
		if (inYears <= 0) {
			BatchError(whichFile, whichLine, 11, "Years");
			nbErrors++;
		}
		ifsSimFile >> inAbsorb;
		if (inAbsorb != 0 && inAbsorb != 1) {
			BatchError(whichFile, whichLine, 1, "Absorbing");
			nbErrors++;
		}
		ifsSimFile >> inFixReplicateSeed;
		if (inFixReplicateSeed != 0 && inFixReplicateSeed != 1) {
			BatchError(whichFile, whichLine, 1, "FixReplicateSeed");
			nbErrors++;
		}
		ifsSimFile >> gEnvStochType;
		if (gUsesPatches == 0) { // cell-based model
			if (gEnvStochType != 0 && gEnvStochType != 1 && gEnvStochType != 2) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "EnvStoch must be 0, 1 or 2 for cell-based model" << endl;
				nbErrors++;
			}
		}
		else { // patch-based model
			if (gEnvStochType != 0 && gEnvStochType != 1) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "EnvStoch must be 0 or 1 for patch-based model" << endl;
				nbErrors++;
			}
		}
		ifsSimFile >> inStochInK;
		if (gEnvStochType) {
			if (inStochInK < 0 || inStochInK > 1) {
				BatchError(whichFile, whichLine, 1, "EnvStochType");
				nbErrors++;
			}
			else gStochInK = inStochInK == 1;
		}
		ifsSimFile >> inStochAC;
		if (gEnvStochType && (inStochAC < 0.0 || inStochAC >= 1.0)) {
			BatchError(whichFile, whichLine, 20, "ac");
			nbErrors++;
		}
		ifsSimFile >> inStochStD;
		if (gEnvStochType && (inStochStD <= 0.0 || inStochStD > 1.0)) {
			BatchError(whichFile, whichLine, 20, "std");
			nbErrors++;
		}

		whichLine++;
		// read next simulation number
		simNb = -98765;
		ifsSimFile >> simNb;
		if (ifsSimFile.eof()) {
			simNb = -98765;
		}
		else { // check for valid simulation number
			if (simNb != prevSim + 1) {
				BatchError(whichFile, whichLine, 222, " ");
				nbErrors++;
			}
			prevSim = simNb;
			nbSims++;
		}
	} // end of while loop
	if (!ifsSimFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}
	return nbErrors == 0;
}

//---------------------------------------------------------------------------
bool CheckParameterFile()
{
	string header, Kheader, intext;
	int i, simNb, inYears;
	species_id inSp;
	int inAbsorb, inGradient, inShifting, inShiftStart, inShiftEnd;
	int inOptimum, inSaveMaps;
	int prevsimul = 0;
	int inRepro, inNbStages, inRepSeasons;
	float inMinR, inMaxR, inMinK, inMaxK, sum_K, min_K, max_K;
	float inGradSteep, inGradScalingFactor, inLocalExtOpt, inShiftRate;
	float inBc, inRmax, inK, inLocalExtProb, inPropMales, inHarem;
	int inOutStartPop, inOutStartInd, inOutStartTraitCell, inOutStartTraitRow;
	int inOutStartConn, inOutIntRange, inOutIntOcc, inOutIntPop, inOutIntInd;
	int inOutIntTraitCell, inOutIntTraitRow, inOutIntConn, inMapsInterval;
	int inSMSHeatMap, inDrawLoadedSp, inFixReplicateSeed;
	int nbErrors = 0;
	int nbKerrors = 0;
	string whichFile = "ParameterFile";

	// Parse header line;
	ifsParamFile >> header; if (header != "Simulation") nbErrors++;
	ifsParamFile >> header; if (header != "Species") nbErrors++;
	ifsParamFile >> header; if (header != "Gradient") nbErrors++;
	ifsParamFile >> header; if (header != "GradSteep") nbErrors++;
	ifsParamFile >> header; if (header != "Optimum") nbErrors++;
	ifsParamFile >> header; if (header != "f") nbErrors++;
	ifsParamFile >> header; if (header != "LocalExtOpt") nbErrors++;
	ifsParamFile >> header; if (header != "Shifting") nbErrors++;
	ifsParamFile >> header; if (header != "ShiftRate") nbErrors++;
	ifsParamFile >> header; if (header != "ShiftStart") nbErrors++;
	ifsParamFile >> header; if (header != "ShiftEnd") nbErrors++;
	ifsParamFile >> header; if (header != "minR") nbErrors++;
	ifsParamFile >> header; if (header != "maxR") nbErrors++;
	ifsParamFile >> header; if (header != "minK") nbErrors++;
	ifsParamFile >> header; if (header != "maxK") nbErrors++;
	ifsParamFile >> header; if (header != "LocalExtProb") nbErrors++;
	ifsParamFile >> header; if (header != "NbStages") nbErrors++;
	ifsParamFile >> header; if (header != "Reproduction") nbErrors++;
	ifsParamFile >> header; if (header != "RepSeasons") nbErrors++;
	ifsParamFile >> header; if (header != "PropMales") nbErrors++;
	ifsParamFile >> header; if (header != "Harem") nbErrors++;
	ifsParamFile >> header; if (header != "bc") nbErrors++;
	ifsParamFile >> header; if (header != "Rmax") nbErrors++;
	for (i = 0; i < gMaxNbHab; i++) {
		Kheader = "K" + to_string(i + 1);
		ifsParamFile >> header; 
		if (header != Kheader) nbKerrors++;
	}
	ifsParamFile >> header; if (header != "OutStartPop") nbErrors++;
	ifsParamFile >> header; if (header != "OutStartInd") nbErrors++;
	ifsParamFile >> header; if (header != "OutStartTraitCell") nbErrors++;
	ifsParamFile >> header; if (header != "OutStartTraitRow") nbErrors++;
	ifsParamFile >> header; if (header != "OutStartConn") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntRange") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntOcc") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntPop") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntInd") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntTraitCell") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntTraitRow") nbErrors++;
	ifsParamFile >> header; if (header != "OutIntConn") nbErrors++;
	ifsParamFile >> header; if (header != "SMSHeatMap") nbErrors++;

	if (nbErrors > 0 || nbKerrors > 0) {
		FormatError(whichFile, nbErrors);
		batchLogOfs << "*** ParameterFile column headers are incorrect." << endl;
		if (nbKerrors > 0) {
			BatchError(whichFile, -999, 333, "K");
		}
		return false;
	}

	// Parse data lines
	int whichLine = 1;
	const int errSimNb = -98765;
	simNb = errSimNb;
	ifsParamFile >> simNb; // first simulation number
	if (simNb == errSimNb) {
		batchLogOfs << "*** Error in ParameterFile - first simulation number could not be read." << endl;
		nbErrors++;
	}
	else if (simNb < 0) {
		batchLogOfs << "*** Error in ParameterFile - first simulation number must be >= 0" << endl;
		nbErrors++;
	}
	while (simNb != -98765) {

		if (!gSpInputOpt.contains(simNb)) {
			batchLogOfs << "Error in ParameterFile - simulation number " << simNb << " is inconsistent with SimFile" << endl;
			nbErrors++;
			break;
		}

		ifsParamFile >> inSp;

		// Initialise input option map for each species
		gSpInputOpt.at(simNb).emplace(inSp, spInputOptions());

		// Record which species exist for landscape checks
		if (!gSpeciesNames.contains(inSp))
			gSpeciesNames.insert(inSp);
		if (!gUseSpeciesDist.contains(inSp)) {
			gUseSpeciesDist.emplace(inSp, false); // false unless set true in SpLandFile
		}
		if (!gUseSMSCosts.contains(inSp)) {
			gUseSMSCosts.emplace(inSp, false);
		}

		ifsParamFile >> inGradient;
		if (gUsesPatches) {
			if (inGradient != 0) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "Gradient must be 0 for patch-based model" << endl;
				nbErrors++;
				inGradient = 0; // to prevent checking of subsequent fields
			}
			inGradient = 0; // to prevent unnecessary checking of subsequent fields
		}
		else { // cell-based model
			if (inGradient < 0 || inGradient > 3) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "Gradient must be between 0 and 3 for cell-based model" << endl;
				nbErrors++;
			}
		}
		ifsParamFile >> inGradSteep;
		if (inGradient && inGradSteep < 0.0) {
			BatchError(whichFile, whichLine, 19, "GradSteep"); 
			nbErrors++; 
		}
		ifsParamFile >> inOptimum;
		if (inGradient && inOptimum < 0) {
			BatchError(whichFile, whichLine, 19, "Optimum"); 
			nbErrors++; 
		}
		ifsParamFile >> inGradScalingFactor;
		if (inGradient && inGradScalingFactor < 0.0) {
			BatchError(whichFile, whichLine, 19, "f"); 
			nbErrors++; 
		}
		ifsParamFile >> inLocalExtOpt;
		if (inGradient == 4 
			&& (inLocalExtOpt < 0.0 || inLocalExtOpt >= 1.0)) {
			BatchError(whichFile, whichLine, 20, "LocalExtOpt"); 
			nbErrors++;
		}
		ifsParamFile >> inShifting;
		if (inGradient && (inShifting != 0 && inShifting != 1)) { 
			BatchError(whichFile, whichLine, 1, "Shifting");
			nbErrors++; 
		}
		ifsParamFile >> inShiftRate;
		if (inGradient && inShifting && inShiftRate <= 0.0) {
			BatchError(whichFile, whichLine, 10, "ShiftRate"); 
			nbErrors++; 
		}
		ifsParamFile >> inShiftStart;
		if (inGradient && inShifting && inShiftStart <= 0) {
			BatchError(whichFile, whichLine, 10, "ShiftStart");
			nbErrors++;
		}
		ifsParamFile >> inShiftEnd;
		if (inGradient && inShifting && inShiftEnd <= inShiftStart) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "ShiftEnd must be greater than ShiftStart" << endl;
			nbErrors++;
		}
		ifsParamFile >> inMinR;
		if (gEnvStochType > 0 && !gStochInK && inMinR <= 0.0) {
			BatchError(whichFile, whichLine, 10, "minR"); 
			nbErrors++; 
		}
		ifsParamFile >> inMaxR;
		if (gEnvStochType > 0 && !gStochInK && inMaxR <= inMinR) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "maxR must be greater than minR" << endl;
			nbErrors++;
		}
		ifsParamFile >> inMinK >> inMaxK;
		if (gEnvStochType > 0 && gStochInK) {
			if (inMinK <= 0.0) { 
				BatchError(whichFile, whichLine, 10, "minK"); 
				nbErrors++; 
			}
			if (inMaxK <= inMinK) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "maxK must be greater than minK" << endl;
				nbErrors++;
			}
		}

		ifsParamFile >> inLocalExtProb;
		if (gUsesPatches && inLocalExtProb != 0.0) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "LocalExtProb must be zero for patch-based models." << endl;
			nbErrors++;
		}
		else if (inLocalExtProb < 0.0 || inLocalExtProb > 1.0) {
			BatchError(whichFile, whichLine, 20, "LocalExtProb");
			nbErrors++;
		}
		
		ifsParamFile >> inNbStages;
		if (gUsesStageStruct) {
			if (inNbStages < 1) {
				BatchError(whichFile, whichLine, 21, "NbStages");
				nbErrors++;
			}
			else if (inNbStages < 2 || inNbStages > 10) {
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << "NbStages must be between 2 and 10." << endl;
			}
			else gSpInputOpt.at(simNb).at(inSp).nbStages = inNbStages;
		}
		else if (inNbStages != gEmptyVal) {
			BatchError(whichFile, whichLine, 0, " ");
			nbErrors++;
			batchLogOfs << "NbStages should left empty (-9) if stage-structure is disabled." << endl;
		}

		ifsParamFile >> inRepro;
		if (!(inRepro == 0 || inRepro == 1 || inRepro == 2)) {
			BatchError(whichFile, whichLine, 2, "ReproType");
			nbErrors++;
		}
		else gSpInputOpt.at(simNb).at(inSp).reproType = inRepro;
		
		ifsParamFile >> inRepSeasons;
		if (inRepSeasons < 1) {
			BatchError(whichFile, whichLine, 21, "RepSeasons");
			nbErrors++;
		}

		ifsParamFile >> inPropMales;
		if (inRepro > 0 
			&& (inPropMales <= 0.0 || inPropMales >= 1.0)) {
			BatchError(whichFile, whichLine, 0, "");
			batchLogOfs << "PropMales should be above 0 and below 1 for sexual models" << endl;
			nbErrors++;
		}
		ifsParamFile >> inHarem;
		if (inRepro == 2 && inHarem <= 0.0) {
			BatchError(whichFile, whichLine, 10, "Harem"); 
			nbErrors++; 
		}
		ifsParamFile >> inBc;
		if (gUsesStageStruct == 0 && inBc <= 0.0) {
			BatchError(whichFile, whichLine, 10, "bc"); 
			nbErrors++; 
		}
		ifsParamFile >> inRmax;
		if (gUsesStageStruct == 0 && inRmax <= 0.0) {
			BatchError(whichFile, whichLine, 10, "Rmax");
			nbErrors++; 
		}
		sum_K = 0.0; 
		min_K = 9999999.0; 
		max_K = 0.0;
		for (i = 0; i < gMaxNbHab; i++) {
			ifsParamFile >> inK;
			if (inK < 0.0) {
				Kheader = "K" + to_string(i + 1);
				BatchError(whichFile, whichLine, 19, Kheader); 
				nbErrors++;
			}
			else {
				sum_K += inK;
				if (inK > 0.0) {
					if (inK < min_K) min_K = inK;
					if (inK > max_K) max_K = inK;
				}
			}
		}
		if (sum_K <= 0.0) {
			BatchError(whichFile, whichLine, 0, " "); 
			nbErrors++;
			batchLogOfs << "At least one K column must be non-zero" << endl;
		}
		else {
			if (gEnvStochType > 0 && gStochInK) { // environmental stochasticity in K
				if (min_K < inMinK || max_K > inMaxK) {
					BatchError(whichFile, whichLine, 0, " "); 
					nbErrors++;
					batchLogOfs << "Non-zero K values must lie between minK and maxK" << endl;
				}
			}
		}

		ifsParamFile >> inOutStartPop;
		if (inOutStartPop < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartPop"); 
			nbErrors++; 
		}
		ifsParamFile >> inOutStartInd;
		if (inOutStartInd < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartInd"); 
			nbErrors++; 
		}
		ifsParamFile >> inOutStartTraitCell;
		if (inOutStartTraitCell < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartTraitCell"); 
			nbErrors++; 
		}
		ifsParamFile >> inOutStartTraitRow;
		if (inOutStartTraitRow < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartTraitRow"); 
			nbErrors++; 
		}
		ifsParamFile >> inOutStartConn;
		if (inOutStartConn < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartConn"); 
			nbErrors++; 
		}
		ifsParamFile >> inOutIntRange;
		if (inOutIntRange < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntRange"); 
			nbErrors++;
		}
		ifsParamFile >> inOutIntOcc;
		if (inOutIntOcc < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntOcc"); 
			nbErrors++; 
		}
		else {
			if (gLandType == 9) {
				if (inOutIntOcc > 0) {
					BatchError(whichFile, whichLine, 0, " "); 
					nbErrors++;
					batchLogOfs << "OutIntOcc must be zero for a generated landscape" << endl;
				}
			}
			else {
				if (gNbReplicates.at(simNb) < 2 && inOutIntOcc > 0) {
					BatchError(whichFile, whichLine, 0, " "); 
					nbErrors++;
					batchLogOfs << "OutIntOcc may be non-zero only if Replicates >= 2" << endl;
				}
			}
		}
		ifsParamFile >> inOutIntPop;
		if (inOutIntPop < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntPop");
			nbErrors++; 
		}
		ifsParamFile >> inOutIntInd;
		if (inOutIntInd < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntInd");
			nbErrors++;
		}
		ifsParamFile >> inOutIntTraitCell;
		if (inOutIntTraitCell < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntTraitCell");
			nbErrors++; 
		}
		ifsParamFile >> inOutIntTraitRow;
		if (inOutIntTraitRow < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntTraitRow"); 
			nbErrors++;
		}
		ifsParamFile >> inOutIntConn;
		if (inOutIntConn < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntConn"); 
			nbErrors++; 
		}
		else if (gUsesPatches != 1 && inOutIntConn > 0) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "OutIntConn may be >0 only if PatchModel is 1" << endl;
			nbErrors++;
		}
		
		ifsParamFile >> inSMSHeatMap;
		if (inSMSHeatMap != 0 && inSMSHeatMap != 1) {
			BatchError(whichFile, whichLine, 1, "SMSHeatMap");
			nbErrors++;
		}

		whichLine++;
		// read next simulation number
		simNb = -98765;
		ifsParamFile >> simNb;
		if (ifsParamFile.eof()) {
			simNb = -98765;
		}
		else { // check for valid simulation number
			if (simNb != prevsimul + 1) {
				BatchError(whichFile, whichLine, 222, " ");
				nbErrors++;
			}
			prevsimul = simNb; 
		}
	} // end of while loop
	if (!ifsParamFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}
	return nbErrors == 0;
}


bool CheckLandFile(int landtype, string inputDir)
{
	string header, inSpLand, inDynLand, inLandscape, whichInputFile;
	int landNb, inNbHab, whichLine;
	float infloat;
	rasterdata patchRaster, spdistraster, costraster;
	int nbErrors = 0;
	vector <int> landlist;
	string whichFile = "LandFile";

	if (landtype == 0 || landtype == 2) { // imported landscape
		// Parse header line;
		ifsLandFile >> header; if (header != "LandNum") nbErrors++;
		ifsLandFile >> header; if (header != "Nhabitats") nbErrors++;
		ifsLandFile >> header; if (header != "LandscapeFile") nbErrors++;
		ifsLandFile >> header; if (header != "SpeciesLandFile") nbErrors++;
		ifsLandFile >> header; if (header != "DynLandFile") nbErrors++;
		if (nbErrors > 0) {
			FormatError(whichFile, 0);
			batchLogOfs << "*** Ensure format is correct for imported landscape" << endl;
			return false;
		}
		// Parse data lines
		whichLine = 1;
		landNb = -98765;
		ifsLandFile >> landNb;
		while (landNb != -98765) {
			if (landNb < 1) {
				BatchError(whichFile, whichLine, 11, "LandNum"); 
				nbErrors++;
			}
			else {
				// landscape number must be unique - retain in list to check
				for (int j = 0; j < (int)landlist.size(); j++) {
					if (landNb == landlist[j]) {
						BatchError(whichFile, whichLine, 666, "LandNum"); 
						j = (int)landlist.size() + 1; 
						nbErrors++;
					}
				}
				landlist.push_back(landNb);
			}
			ifsLandFile >> inNbHab;
			if (landtype == 0) { // raster map with unique habitat codes
				if (inNbHab < 0) {
					BatchError(whichFile, whichLine, 10, "Nhabitats"); 
					nbErrors++;
				}
				if (inNbHab > gMaxNbHab) {
					BatchError(whichFile, whichLine, 0, " ");
					batchLogOfs << "Nhabitats may not exceed MaxHabitats in Control file" << endl;
					nbErrors++;
				}
			}

			// check landscape filename
			whichInputFile = "LandscapeFile";
			ifsLandFile >> inLandscape;
			string pathToLandscape = inputDir + inLandscape;
			landRaster = CheckRasterFile(pathToLandscape);
			if (landRaster.ok) {
				if (landRaster.cellsize == gResol)
					batchLogOfs << whichInputFile << " headers OK: " << pathToLandscape << endl;
				else {
					nbErrors++;
					batchLogOfs << gResolOfStr << whichInputFile << " " << pathToLandscape
						<< gResolNotMatchStr << endl;
				}
			}
			else {
				nbErrors++;
				if (landRaster.errors == -111)
					OpenError(whichInputFile, pathToLandscape);
				else FormatError(pathToLandscape, landRaster.errors);
			}

			// Check species-specific landscape parameters
			whichInputFile = "SpeciesLandFile";
			ifsLandFile >> inSpLand;
			if (inSpLand == "NULL") {
				if (gUsesPatches) {
					BatchError(whichFile, whichLine, 0, " ");
					nbErrors++;
					batchLogOfs << whichInputFile << "SpLandFile is required for patch-based models." << endl;
				}
				if (gTransferType == 1 && gLandType == 2) { // SMS
					BatchError(whichFile, whichLine, 0, " ");
					nbErrors++;
					batchLogOfs << whichInputFile << "SpLandFile is required to specify SMS costs for habitat quality landscapes." << endl;

				}
			}
			else {
				string pathToSpLand = inputDir + inSpLand;
				batchLogOfs << "Checking " << whichInputFile << " " << pathToSpLand << endl;
				ifsSpLandFile.open(pathToSpLand.c_str());
				if (ifsSpLandFile.is_open()) {
					if (!CheckSpLandFile(inputDir, true))
						nbErrors++;
					ifsSpLandFile.close();
					ifsSpLandFile.clear();
				}
				else {
					ifsSpLandFile.clear();
					nbErrors++;
					OpenError(whichInputFile, pathToSpLand);
				}
			}
			
			// check dynamic landscape filename
			whichInputFile = "DynLandFile";
			ifsLandFile >> inDynLand;
			if (inDynLand != "NULL") { // landscape is dynamic
				string pathToDyn = inputDir + inDynLand;
				batchLogOfs << "Checking " << whichInputFile << " " << pathToDyn << endl;
				ifsDynLandFile.open(pathToDyn.c_str());
				if (ifsDynLandFile.is_open()) {
					int errCode = CheckDynamicFile(inputDir);
					if (errCode < 0) {
						nbErrors++;
					}
					ifsDynLandFile.close(); 
					ifsDynLandFile.clear();
				}
				else {
					ifsDynLandFile.clear();
					nbErrors++;
					OpenError(whichInputFile, pathToDyn);
				}
			}

			gNbLandscapes++; 
			whichLine++;
			// read first field on next line
			landNb = -98765;
			ifsLandFile >> landNb;
		} // end of while loop
		landlist.clear();
	} // end of imported landscape
	else {
		if (landtype == 9) { // artificial landscape
			int isFractal, type, Xdim, Ydim;
			float minhab, maxhab;
			// Parse header line;
			ifsLandFile >> header; if (header != "LandNum") nbErrors++;
			ifsLandFile >> header; if (header != "Fractal") nbErrors++;
			ifsLandFile >> header; if (header != "Type") nbErrors++;
			ifsLandFile >> header; if (header != "Xdim") nbErrors++;
			ifsLandFile >> header; if (header != "Ydim") nbErrors++;
			ifsLandFile >> header; if (header != "MinHab") nbErrors++;
			ifsLandFile >> header; if (header != "MaxHab") nbErrors++;
			ifsLandFile >> header; if (header != "Psuit") nbErrors++;
			ifsLandFile >> header; if (header != "H") nbErrors++;
			if (nbErrors > 0) {
				FormatError(whichFile, 0);
				batchLogOfs << "*** Ensure format is correct for artificial landscape" << endl;
				return -111;
			}
			// Parse data lines
			whichLine = 1;
			landNb = -98765;
			ifsLandFile >> landNb;
			while (landNb != -98765) {
				for (int j = 0; j < (int)landlist.size(); j++) {
					if (landNb < 1 || landNb == landlist[j]) {
						BatchError(whichFile, whichLine, 666, "LandNum"); j = (int)landlist.size() + 1; nbErrors++;
					}
				}
				landlist.push_back(landNb);
				ifsLandFile >> isFractal;
				if (isFractal < 0 || isFractal > 1) {
					BatchError(whichFile, whichLine, 1, "Fractal"); nbErrors++;
				}
				ifsLandFile >> type;
				if (type < 0 || type > 1) {
					BatchError(whichFile, whichLine, 1, "Type"); nbErrors++;
				}
				ifsLandFile >> Xdim >> Ydim;
				if (isFractal == 1) {
					if (Xdim < 3) {
						BatchError(whichFile, whichLine, 13, "Xdim"); nbErrors++;
					}
					if (Ydim < 3) {
						BatchError(whichFile, whichLine, 13, "Ydim"); nbErrors++;
					}
				}
				else {
					if (Xdim < 1) {
						BatchError(whichFile, whichLine, 11, "Xdim"); nbErrors++;
					}
					if (Ydim < 1) {
						BatchError(whichFile, whichLine, 11, "Ydim"); nbErrors++;
					}
				}
				if (isFractal == 1) {
					if (Ydim < Xdim) {
						BatchError(whichFile, whichLine, 0, " ");
						batchLogOfs << "Y dimension may not be less than X dimension" << endl; nbErrors++;
					}
					if ((Xdim > 2 && !isValidFractalDim(Xdim - 1))
						|| (Ydim > 2 && !isValidFractalDim(Ydim - 1))) {
						BatchError(whichFile, whichLine, 0, " ");
						batchLogOfs << "X and Y dimensions must be a power of 2 plus 1" << endl; nbErrors++;
					}
				}
				ifsLandFile >> minhab >> maxhab;
				if (type == 1) { // continuous landscape
					if (minhab <= 0.0 || minhab >= 100.0) {
						BatchError(whichFile, whichLine, 100, "MinHab"); nbErrors++;
					}
					if (maxhab <= 0.0 || maxhab > 100.0) {
						BatchError(whichFile, whichLine, 100, "MaxHab"); nbErrors++;
					}
					if (maxhab <= minhab) {
						BatchError(whichFile, whichLine, 0, " ");
						batchLogOfs << "MaxHab must exceed MinHab" << endl; nbErrors++;
					}
				}
				ifsLandFile >> infloat;
				if (infloat < 0.0 || infloat > 1.0) {
					BatchError(whichFile, whichLine, 20, "Psuit"); nbErrors++;
				}
				ifsLandFile >> infloat;
				if (isFractal == 1) {
					if (infloat <= 0.0 || infloat >= 1.0) {
						BatchError(whichFile, whichLine, 20, "H"); nbErrors++;
					}
				}
				gNbLandscapes++; 
				whichLine++;
				// read first field on next line
				landNb = -98765;
				ifsLandFile >> landNb;
			} // end of while loop
		} // end of artificial landscape
		else { // ERROR condition which should not occur
			batchLogOfs << "*** Critical error in land file. "
				<< "Invalid value of landscape type passed to function ParseLandFile()" << endl;
			nbErrors++;
		}
	}
	if (!ifsLandFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}
	return nbErrors > 0;
}

bool CheckSpLandFile(string inputDir, bool isInitial) {

	string header;
	int nbErrors = 0, whichLine = 1;
	string whichFile = "SpLandFile";
	string inPatchFile, inCostFile, inSpDistFile;

	ifsSpLandFile >> header; if (header != "Species") nbErrors++;
	ifsSpLandFile >> header; if (header != "PatchFile") nbErrors++;
	ifsSpLandFile >> header; if (header != "CostMapFile") nbErrors++;
	if (isInitial) { // column absent for dynamic input
		ifsSpLandFile >> header; if (header != "SpDistFile") nbErrors++;
	}
	const int errSpNb = -978;
	int spNb = errSpNb;
	ifsSpLandFile >> spNb;
	while (spNb != errSpNb) {

		if (!gSpeciesNames.contains(spNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Species number " << to_string(spNb) << " doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}

		// check patch map filename
		string whichInputFile = "PatchFile";
		ifsSpLandFile >> inPatchFile;
		if (inPatchFile == "NULL") {
			if (gUsesPatches) {
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << whichInputFile << gPatchReqdStr << endl;
			}
		}
		else {
			string pathToPatchFile = inputDir + inPatchFile;
			rasterdata patchRaster = CheckRasterFile(pathToPatchFile);
			if (patchRaster.ok) {
				if (patchRaster.cellsize == gResol) {
					if (patchRaster.ncols == landRaster.ncols
						&& patchRaster.nrows == landRaster.nrows
						&& patchRaster.cellsize == landRaster.cellsize
						&& (int)patchRaster.xllcorner == (int)landRaster.xllcorner
						&& (int)patchRaster.yllcorner == (int)landRaster.yllcorner) {
						batchLogOfs << whichInputFile << " headers OK: " << pathToPatchFile << endl;
					}
					else {
						batchLogOfs << gHeadersOfStr << whichInputFile << " " << pathToPatchFile
							<< gHeadersNotMatchStr << endl;
						nbErrors++;
					}
				}
				else {
					batchLogOfs << gResolOfStr << whichInputFile << " " << pathToPatchFile
						<< gResolNotMatchStr << endl;
					nbErrors++;
				}
			}
			else {
				nbErrors++;
				if (patchRaster.errors == -111)
					OpenError(whichInputFile, pathToPatchFile);
				else FormatError(pathToPatchFile, patchRaster.errors);
			}

		}

		// check cost map filename
		whichInputFile = "CostMapFile";
		ifsSpLandFile >> inCostFile;
		if (inCostFile == "NULL") {
			if (gTransferType == 1 && gLandType == 2) { // SMS
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << whichInputFile << " is required for a habitat quality landscape" << endl;
			}
		}
		else {
			if (gTransferType != 1) { // SMS
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << whichInputFile << " must be NULL if transfer model is not SMS" << endl;
			}
			else {
				gUseSMSCosts.at(spNb) = true;
				string pathToCosts = inputDir + inCostFile;
				rasterdata costRaster = CheckRasterFile(pathToCosts);
				if (costRaster.ok) {
					if (costRaster.cellsize == gResol) {
						if (costRaster.ncols == landRaster.ncols
							&& costRaster.nrows == landRaster.nrows
							&& costRaster.cellsize == landRaster.cellsize
							&& (int)costRaster.xllcorner == (int)landRaster.xllcorner
							&& (int)costRaster.yllcorner == (int)landRaster.yllcorner) {
							batchLogOfs << whichInputFile << " headers OK: " << pathToCosts << endl;
						}
						else {
							batchLogOfs << gHeadersOfStr << whichInputFile << " " << pathToCosts
								<< gHeadersNotMatchStr << endl;
							nbErrors++;
						}
					}
					else {
						batchLogOfs << gResolOfStr << whichInputFile << " " << pathToCosts
							<< gResolNotMatchStr << endl;
						nbErrors++;
					}
				}
				else {
					nbErrors++;
					if (costRaster.errors == -111)
						OpenError(whichInputFile, pathToCosts);
					else FormatError(pathToCosts, costRaster.errors);
				}
			}
		}

		if (isInitial) { // column absent if called from DynLandFile

			// check initial distribution map filename
			whichInputFile = "SpDistFile";
			ifsSpLandFile >> inSpDistFile;
			if (inSpDistFile != "NULL") {
				gUseSpeciesDist.at(spNb) = true;
				string pathToSpDist = inputDir + inSpDistFile;
				rasterdata spDistRaster = CheckRasterFile(pathToSpDist);
				if (spDistRaster.ok) {
					if (spDistRaster.cellsize < gResol) {
						batchLogOfs << "*** Resolution of " << whichInputFile << " " << pathToSpDist
							<< " must not be smaller than landscape resolution" << endl;
						nbErrors++;
					}
					if (gResol % spDistRaster.cellsize != 0) {
						batchLogOfs << "*** Resolution of " << whichInputFile << " " << pathToSpDist
							<< " must be an integer multiple of landscape resolution" << endl;
						nbErrors++;
					}
					if (spDistRaster.cellsize == landRaster.cellsize) {
						// check that extent matches landscape extent
						if (spDistRaster.ncols != landRaster.ncols
							|| spDistRaster.nrows != landRaster.nrows) {
							batchLogOfs << "*** Extent of " << whichInputFile
								<< " does not match extent of LandscapeFile" << endl;
							nbErrors++;
						}
					}
					// check origins match
					if ((int)spDistRaster.xllcorner == (int)landRaster.xllcorner
						&& (int)spDistRaster.yllcorner == (int)landRaster.yllcorner) {
						batchLogOfs << whichInputFile << " headers OK: " << pathToSpDist << endl;
					}
					else {
						batchLogOfs << "*** Origin co-ordinates of " << whichInputFile
							<< " do not match those of LandscapeFile" << endl;
						nbErrors++;
					}
				}
				else {
					nbErrors++;
					if (spDistRaster.errors == -111)
						OpenError(whichInputFile, pathToSpDist);
					else FormatError(pathToSpDist, spDistRaster.errors);
				}
			}
		}

		spNb = errSpNb;
		ifsSpLandFile >> spNb;
		whichLine++;
	};
	return nbErrors == 0;
}

int CheckDynamicFile(string inputDir) {

	string header, inLandChgFile, inDynSpLand;
	int change, prevChg, year, prevYr = 0;
	rasterdata landChgRaster;
	int nbErrors = 0;
	string whichFile = "DynLandFile";

	ifsDynLandFile >> header; if (header != "Change") nbErrors++;
	ifsDynLandFile >> header; if (header != "Year") nbErrors++;
	ifsDynLandFile >> header; if (header != "LandChangeFile") nbErrors++;
	ifsDynLandFile >> header; if (header != "SpLandChangeFile") nbErrors++;

	if (nbErrors > 0) {
		FormatError(whichFile, nbErrors);
		return -111;
	}

	// Parse data lines
	int whichLine = 1;
	change = -98765;
	ifsDynLandFile >> change; // first change number
	if (change != 1) {
		batchLogOfs << "*** Error in DynLandFile - first change number must be 1" << endl;
		nbErrors++;
	}
	else {
		prevChg = change;
	}
	while (change != -98765) {

		ifsDynLandFile >> year; 
		if (year <= 0) { 
			BatchError(whichFile, whichLine, 10, "Year"); 
			nbErrors++; 
		}
		if (whichLine > 1) {
			if (year <= prevYr) {
				BatchError(whichFile, whichLine, 1, "Year", "previous Year"); 
				nbErrors++;
			}
		}
		prevYr = year;

		// check landscape filename
		string strLandChg = "LandChangeFile";
		ifsDynLandFile >> inLandChgFile;
		string pathToLandChg = inputDir + inLandChgFile;
		landChgRaster = CheckRasterFile(pathToLandChg);
		if (landChgRaster.ok) {
			if (landChgRaster.cellsize == gResol)
				if (landChgRaster.ncols == landRaster.ncols
					&& landChgRaster.nrows == landRaster.nrows
					&& landChgRaster.cellsize == landRaster.cellsize
					&& (int)landChgRaster.xllcorner == (int)landRaster.xllcorner
					&& (int)landChgRaster.yllcorner == (int)landRaster.yllcorner) {
					batchLogOfs << strLandChg << " headers OK: " << pathToLandChg << endl;
				}
				else {
					batchLogOfs << gHeadersOfStr << strLandChg << " " << pathToLandChg
						<< gHeadersNotMatchStr << endl;
					nbErrors++;
				}
			else {
				nbErrors++;
				batchLogOfs << gResolOfStr << strLandChg << " " << pathToLandChg << gResolNotMatchStr << endl;
			}
		}
		else {
			nbErrors++;
			if (landChgRaster.errors == -111)
				OpenError(strLandChg, pathToLandChg);
			else FormatError(pathToLandChg, landChgRaster.errors);
		}

		// Check changes in species-specific parameters
		ifsDynLandFile >> inDynSpLand;
		if (inDynSpLand == "NULL") {
			if (gUsesPatches) {
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << "SpLandChangeFile is required for patch-based models." << endl;
			}
			if (gTransferType == 1 && gLandType == 2) { // SMS
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << "v is required to specify SMS costs for habitat quality landscapes." << endl;

			}
		}
		else {
			string pathToSpLand = inputDir + inDynSpLand;
			batchLogOfs << "Checking SpLandChangeFile " << pathToSpLand << endl;
			ifsSpLandFile.open(pathToSpLand.c_str());
			if (ifsSpLandFile.is_open()) {
				if (!CheckSpLandFile(inputDir, false))
					nbErrors++;
				ifsSpLandFile.close();
				ifsSpLandFile.clear();
			}
			else {
				ifsSpLandFile.clear();
				nbErrors++;
				OpenError("SpLandChangeFile", pathToSpLand);
			}
		}

		whichLine++;
		// read first field on next line
		change = -98765;
		ifsDynLandFile >> change;
		if (ifsDynLandFile.eof()) {
			change = -98765;
		}
		else { // check for valid change number
			if (change != prevChg + 1) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "Change numbers must be sequential integers" << endl;
				nbErrors++;
			}
			prevChg = change;
		}
	}

	if (nbErrors > 0) return -111;
	else return 0;
}

//---------------------------------------------------------------------------
int CheckStageFile(string indir)
{
	string header, fname;
	int simNb, inSp, i, inFecDensDep;
	int usesFecStgWeights, usesDevDensDep, usesDevStgWts, usesSurvDensDep, usesSurvStgWts;
	int inPostDest, inPRep, inRepInt, inMaxAge, inSurvSched;
	float inDevDensCoeff, inSurvDensDepCoeff;
	float infloat;
	int nbErrors = 0;
	int nbSims = 0;
	int prevSim = 0;
	string inTrMatrixFile, fecStgWtFile, inDevStgWtsFile, inSurvWtsFile;
	vector <string> transfiles, wtsfiles;
	const string strStageFile = "StageStructFile";
	const string strTrMatrix = "TransMatrixFile";

	// Parse header line;
	ifsStageStructFile >> header; if (header != "Simulation") nbErrors++;
	ifsStageStructFile >> header; if (header != "Species") nbErrors++;
	ifsStageStructFile >> header; if (header != "PostDestructn") nbErrors++;
	ifsStageStructFile >> header; if (header != "PRep") nbErrors++;
	ifsStageStructFile >> header; if (header != "RepInterval") nbErrors++;
	ifsStageStructFile >> header; if (header != "MaxAge") nbErrors++;
	ifsStageStructFile >> header; if (header != "TransMatrixFile") nbErrors++;
	ifsStageStructFile >> header; if (header != "SurvSched") nbErrors++;
	ifsStageStructFile >> header; if (header != "FecDensDep") nbErrors++;
	ifsStageStructFile >> header; if (header != "FecStageWts") nbErrors++;
	ifsStageStructFile >> header; if (header != "FecStageWtsFile") nbErrors++;
	ifsStageStructFile >> header; if (header != "DevDensDep") nbErrors++;
	ifsStageStructFile >> header; if (header != "DevDensCoeff") nbErrors++;
	ifsStageStructFile >> header; if (header != "DevStageWts") nbErrors++;
	ifsStageStructFile >> header; if (header != "DevStageWtsFile") nbErrors++;
	ifsStageStructFile >> header; if (header != "SurvDensDep") nbErrors++;
	ifsStageStructFile >> header; if (header != "SurvDensCoeff") nbErrors++;
	ifsStageStructFile >> header; if (header != "SurvStageWts") nbErrors++;
	ifsStageStructFile >> header; if (header != "SurvStageWtsFile") nbErrors++;
	if (nbErrors > 0) {
		FormatError(strStageFile, nbErrors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	simNb = -98765;
	ifsStageStructFile >> simNb;
	
	while (simNb != -98765) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(strStageFile, line, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in SimFile" << endl;
			nbErrors++;
		}
		nbSims++;

		ifsStageStructFile >> inSp;
		if (!gSpInputOpt.at(simNb).contains(inSp)) {
			BatchError(strStageFile, line, 0, " ");
			batchLogOfs << "Species number " << to_string(inSp) << " doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(inSp);
		int nbStg = inputOpt.nbStages;
		int nbSexDem = inputOpt.reproType == 2 ? 2 : 1;

		ifsStageStructFile >> inPostDest;
		if (inPostDest < 0 || inPostDest > 1) {
			BatchError(strStageFile, line, 1, "PostDestructn"); 
			nbErrors++; 
		}
		ifsStageStructFile >> inPRep;
		if (inPRep <= 0 || inPRep > 1.0) {
			BatchError(strStageFile, line, 20, "PRep"); 
			nbErrors++;
		}
		
		ifsStageStructFile >> inRepInt;
		if (inRepInt < 0) {
			BatchError(strStageFile, line, 19, "RepInterval"); 
			nbErrors++;
		}
		ifsStageStructFile >> inMaxAge;
		if (inMaxAge < 2) {
			BatchError(strStageFile, line, 12, "MaxAge"); 
			nbErrors++; 
		}
		ifsStageStructFile >> inTrMatrixFile;
		// transition matrix file - compulsory
		bool mustCheck = true;
		for (i = 0; i < (int)transfiles.size(); i++) {
			if (inTrMatrixFile == transfiles[i]) { // file has already been checked
				mustCheck = false;
			}
		}
		if (mustCheck) {
			if (inTrMatrixFile == "NULL") {
				batchLogOfs << "*** " << strTrMatrix << " is compulsory for stage-structured model" << endl;
				nbErrors++;
			}
			else {
				fname = indir + inTrMatrixFile;
				batchLogOfs << "Checking " << strTrMatrix << " " << fname << endl;
				ifsTransMatrix.open(fname.c_str());
				if (ifsTransMatrix.is_open()) {
					if (CheckTransitionFile(nbStg, nbSexDem)) 
						FileHeadersOK(strTrMatrix);
					else nbErrors++;
					ifsTransMatrix.close();
				}
				else {
					OpenError(strTrMatrix, fname); 
					nbErrors++;
				}
				if (ifsTransMatrix.is_open()) 
					ifsTransMatrix.close();
				ifsTransMatrix.clear();
			}
		}
		transfiles.push_back(inTrMatrixFile);

		ifsStageStructFile >> inSurvSched;
		if (inSurvSched < 0 || inSurvSched > 2) {
			BatchError(strStageFile, line, 2, "SurvSched"); 
			nbErrors++; 
		}

		ifsStageStructFile >> inFecDensDep;
		if (inFecDensDep < 0 || inFecDensDep > 1) {
			BatchError(strStageFile, line, 1, "FecDensDep"); 
			nbErrors++; 
			inFecDensDep = 1;
		}
		ifsStageStructFile >> usesFecStgWeights;
		if (inFecDensDep) {
			if (usesFecStgWeights < 0 || usesFecStgWeights > 1) {
				BatchError(strStageFile, line, 1, "FecStageWts");
				nbErrors++; 
				usesFecStgWeights = 1;
			}
		}
		else if (usesFecStgWeights != 0) {
			BatchError(strStageFile, line, 0, " ");
			batchLogOfs << "FecStageWts must be 0 if FecDensDep is 0" << endl; 
			nbErrors++; 
			usesFecStgWeights = 1;
		}
		// fecundity stage weights file - optional
		const string strFecStgWt = "FecStageWtsFile";
		ifsStageStructFile >> fecStgWtFile;
		if (fecStgWtFile == "NULL") {
			if (usesFecStgWeights) {
				BatchError(strStageFile, line, 0, " ");
				batchLogOfs << strFecStgWt << " is compulsory unless FecStageWts is 0" << endl;
				nbErrors++;
			}
		}
		else {
			mustCheck = true;
			for (i = 0; i < (int)wtsfiles.size(); i++) {
				if (fecStgWtFile == wtsfiles[i])
					mustCheck = false; // file has already been checked
			}
			if (mustCheck) {
				fname = indir + fecStgWtFile;
				batchLogOfs << "Checking " << strTrMatrix << " " << fname << endl;
				ifsStageWeightsFile.open(fname.c_str());
				if (ifsStageWeightsFile.is_open()) {
					if (CheckWeightsFile(strTrMatrix, nbStg, nbSexDem))
						FileHeadersOK(strTrMatrix);
					else nbErrors++;
					ifsStageWeightsFile.close();
				}
				else {
					OpenError(strTrMatrix, fname); 
					nbErrors++;
				}
				if (ifsStageWeightsFile.is_open()) 
					ifsStageWeightsFile.close();
				ifsStageWeightsFile.clear();
			}
			wtsfiles.push_back(fecStgWtFile);
		}

		ifsStageStructFile >> usesDevDensDep;
		if (usesDevDensDep < 0 || usesDevDensDep > 1) {
			BatchError(strStageFile, line, 1, "DevDensDep"); 
			nbErrors++; 
			usesDevDensDep = 1;
		}
		ifsStageStructFile >> inDevDensCoeff >> usesDevStgWts;
		if (usesDevDensDep) {
			if (inDevDensCoeff <= 0.0) {
				BatchError(strStageFile, line, 10, "DevDensCoeff"); 
				nbErrors++;
			}
			if (usesDevStgWts < 0 || usesDevStgWts > 1) {
				BatchError(strStageFile, line, 1, "DevStageWts");
				nbErrors++; 
				usesDevStgWts = 1;
			}
		}
		else if (usesDevStgWts != 0) {
			BatchError(strStageFile, line, 0, " ");
			batchLogOfs << "DevStageWts must be 0 if DevDensDep is 0" << endl;
			nbErrors++;
			usesDevStgWts = 1;
		}

		// development stage weights file - optional
		const string strDevStgWts = "DevStageWtsFile";
		ifsStageStructFile >> inDevStgWtsFile;
		if (inDevStgWtsFile == "NULL") {
			if (usesDevStgWts) {
				BatchError(strStageFile, line, 0, " ");
				batchLogOfs << strDevStgWts << " is compulsory unless DevStageWts is 0" << endl;
				nbErrors++;
			}
		}
		else {
			mustCheck = true;
			for (i = 0; i < (int)wtsfiles.size(); i++) {
				if (inDevStgWtsFile == wtsfiles[i])
					mustCheck = false; // file has already been checked
			}
			if (mustCheck) {
				fname = indir + inDevStgWtsFile;
				batchLogOfs << "Checking " << strDevStgWts << " " << fname << endl;
				ifsStageWeightsFile.open(fname.c_str());
				if (ifsStageWeightsFile.is_open()) {
					if (CheckWeightsFile(strDevStgWts, nbStg, nbSexDem))
						FileHeadersOK(strDevStgWts);
					else nbErrors++;
					ifsStageWeightsFile.close();
				}
				else {
					OpenError(strDevStgWts, fname); 
					nbErrors++;
				}
				if (ifsStageWeightsFile.is_open()) 
					ifsStageWeightsFile.close();
				ifsStageWeightsFile.clear();
			}
			wtsfiles.push_back(inDevStgWtsFile);
		}

		ifsStageStructFile >> usesSurvDensDep;
		if (usesSurvDensDep < 0 || usesSurvDensDep > 1) {
			BatchError(strStageFile, line, 1, "SurvDensDep"); 
			nbErrors++; 
			usesSurvDensDep = 1;
		}

		ifsStageStructFile >> inSurvDensDepCoeff >> usesSurvStgWts;
		if (usesSurvDensDep) {
			if (inSurvDensDepCoeff <= 0.0) {
				BatchError(strStageFile, line, 10, "SurvDensCoeff"); 
				nbErrors++;
			}
			if (usesSurvStgWts < 0 || usesSurvStgWts > 1) {
				BatchError(strStageFile, line, 1, "SurvStageWts");
				nbErrors++; 
				usesSurvStgWts = 1;
			}
		}
		else if (usesSurvStgWts != 0) {
			BatchError(strStageFile, line, 0, " ");
			batchLogOfs << "SurvStageWts must be 0 if SurvDensDep is 0" << endl; nbErrors++;
			nbErrors++; 
			usesSurvStgWts = 1;
		}
		// survival stage weights file - optional
		const string strSurvStgWts = "SurvStageWtsFile";
		ifsStageStructFile >> inSurvWtsFile;
		if (inSurvWtsFile == "NULL") {
			if (usesSurvStgWts) {
				BatchError(strStageFile, line, 0, " ");
				batchLogOfs << strSurvStgWts << " is compulsory unless SurvStageWts is 0" << endl;
				nbErrors++;
			}
		}
		else {
			mustCheck = true;
			for (i = 0; i < (int)wtsfiles.size(); i++) {
				if (inSurvWtsFile == wtsfiles[i])
					mustCheck = false; // file has already been checked
			}
			if (mustCheck) {
				fname = indir + inSurvWtsFile;
				batchLogOfs << "Checking " << strSurvStgWts << " " << fname << endl;
				ifsStageWeightsFile.open(fname.c_str());
				if (ifsStageWeightsFile.is_open()) {
					if (CheckWeightsFile(strSurvStgWts, nbStg, nbSexDem))
						FileHeadersOK(strSurvStgWts);
					else nbErrors++;
					ifsStageWeightsFile.close();
				}
				else {
					OpenError(strSurvStgWts, fname);
					nbErrors++;
				}
				if (ifsStageWeightsFile.is_open()) 
					ifsStageWeightsFile.close();
				ifsStageWeightsFile.clear();
			}
			wtsfiles.push_back(inSurvWtsFile);
		}

		// read next simulation
		line++;
		simNb = -98765;
		ifsStageStructFile >> simNb;
		if (ifsStageStructFile.eof()) {
			simNb = -98765;
		}
		else { // check for valid simulation number
			if (simNb != prevSim + 1) {
				BatchError(strStageFile, line, 222, " ");
				nbErrors++;
			}
			prevSim = simNb;
		}
	}
	if (!ifsStageStructFile.eof()) {
		EOFerror(strStageFile);
		nbErrors++;
	}

	transfiles.clear();
	wtsfiles.clear();

	if (nbErrors > 0) return -111;
	else return nbSims;
}

//---------------------------------------------------------------------------
// Check transition matrix file
bool CheckTransitionFile(short nstages, short nsexesDem)
{
	string header, expectedHeader;
	int iStage, iSex, stage, sex, line, minage;
	//int prevminage;
	float infloat;
	int errors = 0;
	string filetype = "TransMatrixFile";

	// check header records
	ifsTransMatrix >> header; 
	if (header != "Transition") errors++;

	for (iStage = 0; iStage < nstages; iStage++) {
		for (iSex = 0; iSex < nsexesDem; iSex++) {
			ifsTransMatrix >> header;
			expectedHeader = to_string(iStage);
			if (nsexesDem != 1) expectedHeader += iSex == 0 ? "m" : "f";
			if (header != expectedHeader) errors++;
		}
	}
	ifsTransMatrix >> header; 
	if (header != "MinAge") errors++;

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// check matrix, including row headers

	// single row for juveniles
	line = 1;
	ifsTransMatrix >> header;
	if (header != "0") {
		BatchError(filetype, line, 0, " ");
		batchLogOfs << "Invalid row header" << endl; errors++;
	}
	float totfecundity = 0.0;
	for (iStage = 0; iStage < nstages; iStage++) {
		for (iSex = 0; iSex < nsexesDem; iSex++) {
			ifsTransMatrix >> infloat;
			if (iStage > 0) {
				if (infloat < 0.0) {
					BatchError(filetype, line, 19, "Fecundity"); errors++;
				}
				totfecundity += infloat;
			}
		}
	}
	if (totfecundity <= 0.0) {
		BatchError(filetype, line, 10, "Total fecundity"); errors++;
	}
	ifsTransMatrix >> minage;
	if (minage != 0) {
		BatchError(filetype, line, 0, " ");
		batchLogOfs << "MinAge must be zero for juvenile stage" << endl; errors++;
	}

	// one row for each stage/sex combination
	for (stage = 1; stage < nstages; stage++) {
		for (sex = 0; sex < nsexesDem; sex++) {
			line++;
			// row header
			ifsTransMatrix >> header;
			expectedHeader = to_string(stage);
			if (nsexesDem == 2) expectedHeader += sex == 0 ? "m" : "f";
			if (header != expectedHeader) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "Invalid row header" << endl; errors++;
			}
			for (iStage = 0; iStage < nstages; iStage++) {
				for (iSex = 0; iSex < nsexesDem; iSex++) {
					ifsTransMatrix >> infloat;
					if (infloat < 0.0 || infloat > 1) {
						BatchError(filetype, line, 20, "Transition probability"); errors++;
					}
				}
			}
			ifsTransMatrix >> minage;
			if (stage == 1 && minage != 0) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "MinAge must be zero for stage 1" << endl; 
				errors++;
			}
			if (stage > 1) {
				if (minage < 0) {
					BatchError(filetype, line, 19, "MinAge"); 
					errors++;
				}
			}
		}
	}
	// final read should hit EOF
	ifsTransMatrix >> header;

	if (!ifsTransMatrix.eof()) {
		EOFerror(filetype);
		errors++;
	}

	return errors == 0;

}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Check stage weights matrix file
bool CheckWeightsFile(string filetype, int nbStages, int nbSexes)
{
	string header, hhh;
	int i, j, stage, sex, line;
	float infloat;
	int errors = 0;

	// check header records
	ifsStageWeightsFile >> header; if (header != "StageWts") errors++;
	for (i = 0; i < nbStages; i++) {
		for (j = 0; j < nbSexes; j++) {
			ifsStageWeightsFile >> header;
			if (nbSexes == 1) hhh = to_string(i);
			else {
				if (j == 0) hhh = to_string(i) + "m"; else hhh = to_string(i) + "f";
			}
			if (header != hhh) errors++;
		}
	}

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// check matrix, including row headers
	// one row for each stage/sex combination
	line = 0;
	for (stage = 0; stage < nbStages; stage++) {
		for (sex = 0; sex < nbSexes; sex++) {
			line++;
			// row header
			ifsStageWeightsFile >> header;
			if (nbSexes == 1) hhh = to_string(stage);
			else {
				if (sex == 0) hhh = to_string(stage) + "m"; else hhh = to_string(stage) + "f";
			}
			if (header != hhh) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "Invalid row header" << endl; errors++;
			}
			for (i = 0; i < nbStages; i++) {
				for (j = 0; j < nbSexes; j++) {
					ifsStageWeightsFile >> infloat;
					// NOTE - any real number is acceptable - no check required
				}
			}
		}
	}
	// final read should hit EOF
	ifsStageWeightsFile >> header;

	if (!ifsStageWeightsFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	return errors == 0;

}

//---------------------------------------------------------------------------
int CheckEmigFile()
{
	string header;
	int inDensDep, inUseFullKern, inStgDep, inSexDep, inIndVar, inEmigStg, inStage, inSex;
	bool isDensDep, isIndVar;
	float inSp, inEP, inD0, inAlpha, inBeta;
	int nbErrors = 0;
	int nbSims = 0;
	string whichInputFile = "EmigrationFile";

	isDensDep = false;
	isIndVar = false;
	inEP = 0.0;

	// Parse header line;
	ifsEmigrationFile >> header; if (header != "Simulation") nbErrors++;
	ifsEmigrationFile >> header; if (header != "Species") nbErrors++;
	ifsEmigrationFile >> header; if (header != "DensDep") nbErrors++;
	ifsEmigrationFile >> header; if (header != "UseFullKern") nbErrors++;
	ifsEmigrationFile >> header; if (header != "StageDep") nbErrors++;
	ifsEmigrationFile >> header; if (header != "SexDep") nbErrors++;
	ifsEmigrationFile >> header; if (header != "IndVar") nbErrors++;
	ifsEmigrationFile >> header; if (header != "EmigStage") nbErrors++;
	ifsEmigrationFile >> header; if (header != "Stage") nbErrors++;
	ifsEmigrationFile >> header; if (header != "Sex") nbErrors++;
	ifsEmigrationFile >> header; if (header != "EP") nbErrors++;
	ifsEmigrationFile >> header; if (header != "D0") nbErrors++;
	ifsEmigrationFile >> header; if (header != "alpha") nbErrors++;
	ifsEmigrationFile >> header; if (header != "beta") nbErrors++;

	if (nbErrors > 0) {
		FormatError(whichInputFile, nbErrors);
		return -111;
	}

	// Parse data lines
	bool readNextLine = true;
	int lineNb = 1;
	simCheck currentLine, prevLine;
	int simNb = prevLine.simNb = -999;
	prevLine.simLines = prevLine.reqdSimLines = 0;
	ifsEmigrationFile >> simNb;

	while (readNextLine) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in SimFile" << endl;
			nbErrors++;
		}

		ifsEmigrationFile >> inSp;
		if (!gSpInputOpt.at(simNb).contains(inSp)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Species number " << to_string(inSp) << " doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(inSp);

		// read and validate columns relating to stage and sex-dependency and to IIV
		ifsEmigrationFile >> inDensDep >> inUseFullKern >> inStgDep >> inSexDep;
		ifsEmigrationFile >> inIndVar >> inEmigStg >> inStage >> inSex;
		currentLine = CheckStageSex(whichInputFile, lineNb, simNb, inSp, prevLine, inStgDep, inSexDep, inStage, inSex, inIndVar, true, false);
		if (currentLine.isNewSim) nbSims++;
		nbErrors += currentLine.errors;
		prevLine = currentLine;

		// validate density dependency
		if (inDensDep != 0 && inDensDep != 1) {
			BatchError(whichInputFile, lineNb, 1, "DensDep"); 
			nbErrors++;
		}
		else {
			inputOpt.isEmigDensDep = (inDensDep == 1);
		}

		// validate individual variation
		if (inIndVar != 0 && inIndVar != 1) {
			BatchError(whichInputFile, lineNb, 1, "IndVar");
			nbErrors++;
		}
		else {
			inputOpt.isEmigIndVar = (inIndVar == 1);
		}

		// validate use full kernel
		if (inUseFullKern != 0 && inUseFullKern != 1) {
			BatchError(whichInputFile, lineNb, 1, "UseFullKern"); 
			nbErrors++;
		}
		if (inDensDep == 1 && inUseFullKern != 0) {
				BatchError(whichInputFile, lineNb, 0, "UseFullKern"); 
				nbErrors++;
				batchLogOfs << "UseFullKern must be 0 if there is density-dependent emigration" << endl;
		}
		// validate emigration stage
		if (gUsesStageStruct && !inStgDep && inIndVar == 1
			&& inStage == 0 && inSex == 0
			&& (inEmigStg < 0 || inEmigStg >= inputOpt.nbStages)) {
			BatchError(whichInputFile, lineNb, 0, "EmigStage");
			nbErrors++;
			batchLogOfs << "EmigStage must be from 0 to " << to_string(inputOpt.nbStages - 1) << endl;
		}
		if (inSexDep != 0 && inSexDep != 1) {
			BatchError(whichInputFile, lineNb, 1, "SexDep");
			nbErrors++;
		} 
		else {
			inputOpt.isEmigSexDep = (inSexDep == 1);
		}

		if (inStage == 0 && inSex == 0) { // first line of a simulation
			// record whether density dependence and individual variability are applied
			isDensDep = (inDensDep == 1);
			isIndVar = (inIndVar == 1);
		}

		// read remaining columns of the current record
		ifsEmigrationFile >> inEP >> inD0 >> inAlpha >> inBeta;

		if (inputOpt.isEmigIndVar) {
			if (inEP != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if individual variability is enabled EP must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inD0 != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if individual variability is enabled D0 must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inAlpha != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if individual variability is enabled alpha must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inBeta != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if individual variability is enabled beta must be " << gEmptyVal << endl;
				nbErrors++;
			}
		}
		else if (isDensDep) {
			if (inEP != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is enabled EP must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inD0 < 0.0 || inD0 > 1.0) {
				BatchError(whichInputFile, lineNb, 20, "D0"); 
				nbErrors++;
			}
		}
		else { // !densdepset
			if (inEP < 0.0 || inEP > 1.0) {
				BatchError(whichInputFile, lineNb, 20, "EP"); 
				nbErrors++;
			}
			if (inD0 != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled D0 must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inAlpha != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled alpha must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inBeta != gEmptyVal) {
				batchLogOfs << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled beta must be " << gEmptyVal << endl;
				nbErrors++;
			}
		}

		// read next simulation
		lineNb++;
		const int errSimNb = -98765;
		simNb = errSimNb;
		ifsEmigrationFile >> simNb;
		if (simNb == errSimNb || ifsEmigrationFile.eof()) 
			readNextLine = false;

	} // end of while loop

	// check for correct number of lines for previous simulation
	if (currentLine.simLines != currentLine.reqdSimLines) {
		BatchError(whichInputFile, lineNb, 0, " "); 
		nbErrors++;
		batchLogOfs << gNbLinesStr << currentLine.simNb
			<< gShouldBeStr << currentLine.reqdSimLines << endl;
	}
	if (!ifsEmigrationFile.eof()) {
		EOFerror(whichInputFile);
		nbErrors++;
	}
	if (nbErrors > 0) return -111;
	else return nbSims;
}

//---------------------------------------------------------------------------
int CheckTransferFile(string indir)
{
	string header, colheader, intext, fname, ftype;
	int i, simNb, inSp, inStageDep, inSexDep, inKernelType, inDistMort, inIndVar, inStage, inSex;
	int	inPercRangeMethod, inSMType, inStraightenPath;
	float inPerceptualRange, inDirPersistence, inSMConst;
	int inGoalType, inMemSize, inBetaDispBias; float inGoalBias, inAlphaDispBias;
	float meanDistI, meanDistII, ProbKernelI;
	float mortProb, slope, inflPoint;
	float morthab, mortmatrix;
	int costhab, costmatrix;
	float inStepLength, inStepCorr;

	vector <string> costsfiles;

	int errors = 0; 
	int morthaberrors = 0; 
	int costerrors = 0; 
	int hrerrors = 0;
	int simuls = 0;
	string whichFile = "TransferFile";

	// Parse header line;
	ifsTransferFile >> header; 
	if (header != "Simulation") errors++;
	ifsTransferFile >> header; 
	if (header != "Species") errors++;

	switch (gTransferType) {

	case 0: { // negative exponential dispersal kernel
		batchLogOfs << "Checking dispersal kernel format file" << endl;
		ifsTransferFile >> header; if (header != "StageDep") errors++;
		ifsTransferFile >> header; if (header != "SexDep") errors++;
		ifsTransferFile >> header; if (header != "KernelType") errors++;
		ifsTransferFile >> header; if (header != "DistMort") errors++;
		ifsTransferFile >> header; if (header != "IndVar") errors++;
		ifsTransferFile >> header; if (header != "Stage") errors++;
		ifsTransferFile >> header; if (header != "Sex") errors++;
		ifsTransferFile >> header; if (header != "meanDistI") errors++;
		ifsTransferFile >> header; if (header != "meanDistII") errors++;
		ifsTransferFile >> header; if (header != "ProbKernelI") errors++;
		ifsTransferFile >> header; if (header != "MortProb") errors++;
		ifsTransferFile >> header; if (header != "Slope") errors++;
		ifsTransferFile >> header; if (header != "InflPoint") errors++;
		break;
	} // end of negative exponential dispersal kernel

	case 1: { // SMS
		batchLogOfs << "Checking SMS format file ";
		ifsTransferFile >> header; if (header != "IndVar") errors++;
		ifsTransferFile >> header; if (header != "PR") errors++;
		ifsTransferFile >> header; if (header != "PRMethod") errors++;
		ifsTransferFile >> header; if (header != "DP") errors++;
		ifsTransferFile >> header; if (header != "MemSize") errors++;
		ifsTransferFile >> header; if (header != "GB") errors++;
		ifsTransferFile >> header; if (header != "GoalType") errors++;
		ifsTransferFile >> header; if (header != "AlphaDB") errors++;
		ifsTransferFile >> header; if (header != "BetaDB") errors++;
		ifsTransferFile >> header; if (header != "StraightenPath") errors++;
		ifsTransferFile >> header; if (header != "SMtype") errors++;
		ifsTransferFile >> header; if (header != "SMconst") errors++;
		switch (gLandType) {
		case 0: { // raster map with unique habitat codes
			batchLogOfs << "for LandType = 0" << endl;
			for (i = 0; i < gMaxNbHab; i++) {
				colheader = "MortHab" + to_string(i + 1);
				ifsTransferFile >> header; if (header != colheader) morthaberrors++;
			}
			for (i = 0; i < gMaxNbHab; i++) {
				colheader = "CostHab" + to_string(i + 1);
				ifsTransferFile >> header; if (header != colheader) costerrors++;
			}
			break;
		} // end of raster map with unique habitat codes
		case 2: { // raster map with habitat quality
			batchLogOfs << "for LandType = 2" << endl;
			break;
		} // end of raster map with habitat quality
		case 9: { // artificial landscape
			batchLogOfs << "for LandType = 9" << endl;
			ifsTransferFile >> header; if (header != "MortHabitat") errors++;
			ifsTransferFile >> header; if (header != "MortMatrix") errors++;
			ifsTransferFile >> header; if (header != "CostHabitat") errors++;
			ifsTransferFile >> header; if (header != "CostMatrix") errors++;
			break;
		} // end of artificial landscape
		} // end of switch (landtype)
		break;
	} // end of SMS

	case 2: { // CRW
		batchLogOfs << "Checking CRW format file" << endl;
		ifsTransferFile >> header; if (header != "IndVar") errors++;
		ifsTransferFile >> header; if (header != "SL") errors++;
		ifsTransferFile >> header; if (header != "Rho") errors++;
		ifsTransferFile >> header; if (header != "StraightenPath") errors++;
		ifsTransferFile >> header; if (header != "SMtype") errors++;
		ifsTransferFile >> header; if (header != "SMconst") errors++;
		if (gLandType == 0) {
			for (i = 0; i < gMaxNbHab; i++) {
				colheader = "MortHab" + to_string(i + 1);
				ifsTransferFile >> header; if (header != colheader) morthaberrors++;
			}
		}
		break;
	} // end of CRW

	} // end of switch (transfer)
	// report any errors in headers, and if so, terminate validation
	if (errors > 0 || morthaberrors > 0 || costerrors > 0 || hrerrors > 0) {
		FormatError(whichFile, errors + morthaberrors + costerrors);
		if (morthaberrors > 0) BatchError(whichFile, -999, 333, "MortHab");
		if (costerrors > 0) BatchError(whichFile, -999, 333, "CostHab");
		if (hrerrors > 0) BatchError(whichFile, -999, 444, "Hr");
		return -111;
	}

	// Parse data lines
	int whichLine = 1;
	simCheck current, prev;
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	ifsTransferFile >> simNb;
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in SimFile." << endl;
			errors++;
		}

		ifsTransferFile >> inSp;
		if (!gSpInputOpt.at(simNb).contains(inSp)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Species number " << to_string(inSp) << " doesn't match those in ParametersFile" << endl;
			errors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(inSp);

		switch (gTransferType) {

		case 0: { // negative exponential dispersal kernel
			// read and validate columns relating to stage and sex-dependency and to IIV
			ifsTransferFile >> inStageDep >> inSexDep >> inKernelType >> inDistMort;
			ifsTransferFile >> inIndVar >> inStage >> inSex;
			current = CheckStageSex(whichFile, whichLine, simNb, inSp, prev, inStageDep, inSexDep, inStage, inSex, inIndVar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			// validate kernel type
			if (inKernelType != 0 && inKernelType != 1) {
				BatchError(whichFile, whichLine, 1, "KernelType"); errors++;
			}
			else {
				inputOpt.usesTwoKernels = (inKernelType == 1);
			}
			// validate mortality
			if (inDistMort != 0 && inDistMort != 1) {
				BatchError(whichFile, whichLine, 1, "DistMort"); 
				errors++;
			}
			// read remaining columns of the current record
			ifsTransferFile >> meanDistI >> meanDistII >> ProbKernelI;
			ifsTransferFile >> mortProb >> slope >> inflPoint;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); 
				errors++;
			}
			else {
				inputOpt.isKernTransfIndVar = (inIndVar == 1);
			}

			if (inSexDep != 0 && inSexDep != 1) {
				BatchError(whichFile, whichLine, 1, "SexDep"); 
				errors++;
			}
			else {
				inputOpt.isKernTransfSexDep = (inSexDep == 1);
			}

			// validate mortality
			if (inDistMort != 0 && inDistMort != 1) {
				BatchError(whichFile, whichLine, 1, "DistMort"); 
				errors++;
			}

			if (inputOpt.isKernTransfIndVar) {
				if (meanDistI != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled meanDistI must be " << gEmptyVal << endl;
					errors++;
				}
				if (meanDistII != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled meanDistII must be " << gEmptyVal << endl;
					errors++;
				}
				if (ProbKernelI != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled ProbKernelI must be " << gEmptyVal << endl;
					errors++;
				}
			}
			else {
				if (meanDistI < gResol) {
					// NOTE - should also check whether emigration prob is constant and equal to 1
					//but checks across diffferent input files are not yet implemented
					BatchError(whichFile, whichLine, 2, "meanDistI", "Resolution"); errors++;
				}
				if (inKernelType != 0) {
					if (meanDistII < gResol) {
						// NOTE - DITTO
						BatchError(whichFile, whichLine, 2, "meanDistII", "Resolution"); errors++;
					}
					if (ProbKernelI <= 0.0 || ProbKernelI >= 1.0) {
						BatchError(whichFile, whichLine, 20, "ProbKernelI"); errors++;
					}
				}
			}

			if (inStage == 0 && inSex == 0) {
				if (inDistMort) { // distance-dependent mortality
					// WHAT CONDITIONS APPLY TO MORTALITY SLOPE AND INFLECTION POINT?
				}
				else { // constant mortality
					if (mortProb < 0.0 || mortProb >= 1.0) {
						BatchError(whichFile, whichLine, 20, "MortProb"); errors++;
					}
				}
			}

			break;
		} // end of negative exponential dispersal kernel

		case 1: { // SMS
			ifsTransferFile >> inIndVar;
			ifsTransferFile >> inPerceptualRange >> inPercRangeMethod >> inDirPersistence;
			ifsTransferFile >> inMemSize >> inGoalBias >> inGoalType >> inAlphaDispBias >> inBetaDispBias;
			current = CheckStageSex(whichFile, whichLine, simNb, inSp, prev, 0, 0, 0, 0, 0, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				inputOpt.isSMSTransfIndVar = (inIndVar == 1);
			}

			// validate SMS movement parameters
			if (inPerceptualRange < 1) {
				BatchError(whichFile, whichLine, 11, "PR"); errors++;
			}
			if (inPercRangeMethod < 1 || inPercRangeMethod > 3) {
				BatchError(whichFile, whichLine, 33, "PRmethod"); errors++;
			}
			if (inputOpt.isSMSTransfIndVar) {
				if (inGoalBias != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled GB must be " << gEmptyVal << endl;
					errors++;
				}
				if (inDirPersistence != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled DP must be " << gEmptyVal << endl;
					errors++;
				}
				if (inAlphaDispBias != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled AlphaDB must be " << gEmptyVal << endl;
					errors++;
				}
				if (inBetaDispBias != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled BetaDB must be " << gEmptyVal << endl;
					errors++;
				}
			}

			if (!inIndVar && inDirPersistence < 1.0) {
				BatchError(whichFile, whichLine, 11, "DP"); errors++;
			}
			if (inMemSize < 1 || inMemSize > 14) {
				BatchError(whichFile, whichLine, 0, "MemSize"); errors++;
				batchLogOfs << "MemSize must be from 1 to 14" << endl;
			}
			if (!inIndVar && inGoalBias < 1.0) {
				BatchError(whichFile, whichLine, 11, "GB"); errors++;
			}
			if (inGoalType != 0 && inGoalType != 2) {
				BatchError(whichFile, whichLine, 2, "GoalType"); errors++;
			}
			else {
				inputOpt.usesSMSGoalBias = (inGoalType == 2);
			}
			ifsTransferFile >> inStraightenPath >> inSMType >> inSMConst;
			if (inStraightenPath != 0 && inStraightenPath != 1) {
				BatchError(whichFile, whichLine, 1, "StraightenPath"); errors++;
			}
			if (gLandType == 2) // habitat quality landscape 
			{ // must have constant mortality
				if (inSMType != 0) {
					BatchError(whichFile, whichLine, 0, " "); errors++;
					batchLogOfs << "SMtype must be 0 for LandType 2" << endl;
				}
			}
			else if (inSMType != 0 && inSMType != 1) {
				BatchError(whichFile, whichLine, 1, "SMtype"); errors++;
			}
			if (inSMType == 0) {
				if (inSMConst < 0.0 || inSMConst >= 1.0) {
					BatchError(whichFile, whichLine, 20, "SMconst"); errors++;
				}
			}
			switch (gLandType) {

			case 0: { // raster map with unique habitat codes
				for (i = 0; i < gMaxNbHab; i++) {
					ifsTransferFile >> morthab;
					if (inSMType == 1) {
						if (morthab < 0.0 || morthab >= 1.0) {
							colheader = "MortHab" + to_string(i + 1);
							BatchError(whichFile, whichLine, 20, colheader); errors++;
						}
					}
				}
				for (i = 0; i < gMaxNbHab; i++) {
					ifsTransferFile >> costhab;
					if (!gUseSMSCosts.at(inSp) && costhab < 1) {
						colheader = "CostHab" + to_string(i + 1);
						BatchError(whichFile, whichLine, 11, colheader);
						errors++;
					}
				}
				break;
			} // end of raster map with unique habitat codes

			case 2: { // raster map with habitat quality
				break;
			}

			case 9: { // artificial landscape
				ifsTransferFile >> morthab >> mortmatrix;
				ifsTransferFile >> costhab >> costmatrix;
				if (inSMType) { // validate habitat-dependent mortality
					if (morthab < 0.0 || morthab >= 1.0) {
						BatchError(whichFile, whichLine, 20, "MortHabitat"); errors++;
					}
					if (mortmatrix < 0.0 || mortmatrix >= 1.0) {
						BatchError(whichFile, whichLine, 20, "MortMatrix"); errors++;
					}
				}
				if (costhab < 1) {
					BatchError(whichFile, whichLine, 11, "CostHabitat"); errors++;
				}
				if (costmatrix < 1) {
					BatchError(whichFile, whichLine, 11, "CostMatrix"); errors++;
				}
				break;
			} // end of artificial landscape

			} // end of switch (landtype)

			break;

		} // end of SMS

		case 2: { // CRW
			ifsTransferFile >> inIndVar >> inStepLength >> inStepCorr >> inStraightenPath >> inSMType >> inSMConst;
			current = CheckStageSex(whichFile, whichLine, simNb, inSp, prev, 0, 0, 0, 0, inIndVar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				inputOpt.isCRWTransfIndVar = (inIndVar == 1);
			}

			if (inputOpt.isCRWTransfIndVar) {
				if (inStepLength != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled SL must be " << gEmptyVal << endl;
					errors++;
				}
				if (inStepCorr != gEmptyVal) {
					batchLogOfs << "*** Error in " << whichFile << ": "
						<< "if individual variability is enabled Rho must be " << gEmptyVal << endl;
					errors++;
				}
			}
			else {
				if (inStepLength <= 0.0) {
					BatchError(whichFile, whichLine, 10, "SL"); errors++;
				}
				if (inStepCorr <= 0.0 || inStepCorr >= 1.0) {
					BatchError(whichFile, whichLine, 20, "Rho"); errors++;
				}
			}
			
			if (inStraightenPath != 0 && inStraightenPath != 1) {
				BatchError(whichFile, whichLine, 1, "StraightenPath"); errors++;
			}
			if (gLandType == 0) { // imported landscape with habitat types
				if (inSMType != 0 && inSMType != 1) {
					BatchError(whichFile, whichLine, 1, "SMtype"); errors++;
				}
				if (!inSMType) {
					if (inSMConst < 0.0 || inSMConst >= 1.0) {
						BatchError(whichFile, whichLine, 20, "SMconst"); errors++;
					}
				}
				for (int i = 0; i < gMaxNbHab; i++) {
					ifsTransferFile >> morthab;
					if (inSMType) {
						if (morthab < 0.0 || morthab >= 1.0) {
							colheader = "MortHab" + to_string(i + 1);
							BatchError(whichFile, whichLine, 20, colheader); errors++;
						}
					}
				}
			}
			else { // imported landscape with quality OR artificial landscape
				if (inSMType != 0) {
					BatchError(whichFile, whichLine, 0, " "); errors++;
					batchLogOfs << "SMtype must be 0 for LandType 2 or 9" << endl;
				}
				if (inSMConst < 0.0 || inSMType >= 1.0) {
					BatchError(whichFile, whichLine, 20, "SMconst"); errors++;
				}
			}
			break;
		} // end of CRW

		} // end of switch (transfer)

		// read next simulation
		whichLine++;
		simNb = -98765;
		ifsTransferFile >> simNb;
		if (ifsTransferFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (gTransferType == 0 // no. of lines checked for dispersal kernel transfer method only
		&& current.simLines != current.reqdSimLines) {
		BatchError(whichFile, whichLine, 0, " "); 
		errors++;
		batchLogOfs << gNbLinesStr << current.simNb
			<< gShouldBeStr << current.reqdSimLines << endl;
	}
	if (!ifsTransferFile.eof()) {
		EOFerror(whichFile);
		errors++;
	}
	costsfiles.clear();

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int CheckSettleFile()
{
	string header;
	int simNb, inSp, inStageDep, inSexDep, inStage, inSex, inSettleType;
	int inDensDep, inIndVar, inFindMate, inMinSteps, inMaxSteps, inMaxStepsYear;
	float inS0, inAlphaS, inBetaS;
	int nbErrors = 0;
	int nbSims = 0;
	string whichFile = "SettlementFile";

	// Parse header line;
	ifsSettlementFile >> header; if (header != "Simulation") nbErrors++;
	ifsSettlementFile >> header; if (header != "Species") nbErrors++;
	ifsSettlementFile >> header; if (header != "StageDep") nbErrors++;
	ifsSettlementFile >> header; if (header != "SexDep") nbErrors++;
	ifsSettlementFile >> header; if (header != "Stage") nbErrors++;
	ifsSettlementFile >> header; if (header != "Sex") nbErrors++;
	if (gTransferType == 0) { 
		// dispersal kernel
		ifsSettlementFile >> header; if (header != "SettleType") nbErrors++;
		ifsSettlementFile >> header; if (header != "FindMate") nbErrors++;
	}
	else { // movement method
		ifsSettlementFile >> header; if (header != "DensDep") nbErrors++;
		ifsSettlementFile >> header; if (header != "IndVar") nbErrors++;
		ifsSettlementFile >> header; if (header != "FindMate") nbErrors++;
		ifsSettlementFile >> header; if (header != "MinSteps") nbErrors++;
		ifsSettlementFile >> header; if (header != "MaxSteps") nbErrors++;
		ifsSettlementFile >> header; if (header != "MaxStepsYear") nbErrors++;
		ifsSettlementFile >> header; if (header != "S0") nbErrors++;
		ifsSettlementFile >> header; if (header != "AlphaS") nbErrors++;
		ifsSettlementFile >> header; if (header != "BetaS") nbErrors++;
	}
	if (nbErrors > 0) {
		FormatError(whichFile, nbErrors);
		return -111;
	}

	// Parse data lines
	int whichLine = 1;
	simCheck current, prev;
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	ifsSettlementFile >> simNb;
	
	current.simNb = 0;
	while (simNb != -98765) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in SimFile" << endl;
			nbErrors++;
		}

		ifsSettlementFile >> inSp;
		if (!gSpInputOpt.at(simNb).contains(inSp)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Species number " << to_string(inSp) << " doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(inSp);
		int nbSexDisp = inputOpt.reproType == 0 ? 1 : 2;

		if (gTransferType == 0) { 
			// dispersal kernel
			// read and validate columns relating to stage and sex-dependency (NB no IIV here)
			ifsSettlementFile >> inStageDep >> inSexDep >> inStage >> inSex >> inSettleType >> inFindMate;
			current = CheckStageSex(whichFile, whichLine, simNb, inSp, prev, inStageDep, inSexDep, inStage, inSex, 0, true, false);
			if (current.isNewSim) nbSims++;
			nbErrors += current.errors;
			prev = current;
			if (inSettleType < 0 || inSettleType > 3) {
				BatchError(whichFile, whichLine, 3, "SettleType"); nbErrors++;
			}
			if (!gUsesStageStruct && (inSettleType == 1 || inSettleType == 3)) {
				BatchError(whichFile, whichLine, 0, " "); nbErrors++;
				batchLogOfs << "Invalid SettleType for a non-stage-structured population" << endl;
			}
			if (nbSexDisp > 1) {
				if (inFindMate < 0 || inFindMate > 1) {
					BatchError(whichFile, whichLine, 1, "FindMate"); 
					nbErrors++;
				}
			}
		}
		else { // movement method
			// read and validate columns relating to stage and sex-dependency (IIV psossible)
			ifsSettlementFile >> inStageDep >> inSexDep >> inStage >> inSex >> inDensDep >> inIndVar >> inFindMate;
			current = CheckStageSex(whichFile, whichLine, simNb, inSp, prev, inStageDep, inSexDep, inStage, inSex, inIndVar, true, false);
			if (current.isNewSim) nbSims++;
			nbErrors += current.errors;
			prev = current;

			if (inDensDep != 0 && inDensDep != 1) {
				BatchError(whichFile, whichLine, 1, "DensDep");
				nbErrors++;
			}
			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar");
				nbErrors++;
			}

			if (inDensDep == 0 && inIndVar != 0) {
				BatchError(whichFile, whichLine, 0, " ");
				nbErrors++;
				batchLogOfs << "IndVar must be 0 if DensDep is 0" << endl;
			}
			else {
				inputOpt.isSettIndVar = inIndVar == 1;
			}

			if (inSexDep != 0 && inSexDep != 1) {
				BatchError(whichFile, whichLine, 1, "SexDep");
				nbErrors++;
			}
			else {
				inputOpt.isSettSexDep = inSexDep == 1;
			}

			if (inputOpt.reproType != 0 && nbSexDisp > 1) {
				if (inFindMate != 0 && inFindMate != 1) {
					BatchError(whichFile, whichLine, 1, "FindMate"); 
					nbErrors++;
				}
			}
			ifsSettlementFile >> inMinSteps >> inMaxSteps >> inMaxStepsYear;
			if (inStage == 0 && inSex == 0) {
				if (inMinSteps < 0) {
					BatchError(whichFile, whichLine, 19, "MinSteps"); 
					nbErrors++;
				}
				if (inMaxSteps < 0) {
					BatchError(whichFile, whichLine, 19, "MaxSteps");
					nbErrors++;
				}
			}
			if (inMaxStepsYear < 0) {
				BatchError(whichFile, whichLine, 19, "MaxStepsYear");
				nbErrors++;
			}
			ifsSettlementFile >> inS0 >> inAlphaS >> inBetaS;

			if (inputOpt.isSettIndVar) {
					if (inS0 != gEmptyVal) {
						batchLogOfs << "*** Error in " << whichFile << ": "
							<< "if individual variability is enabled S0 must be " << gEmptyVal << endl;
						nbErrors++;
					}
					if (inAlphaS != gEmptyVal) {
						batchLogOfs << "*** Error in " << whichFile << ": "
							<< "if individual variability is enabled AlphaS must be " << gEmptyVal << endl;
						nbErrors++;
					}
					if (inBetaS != gEmptyVal) {
						batchLogOfs << "*** Error in " << whichFile << ": "
							<< "if individual variability is enabled BetaS must be " << gEmptyVal << endl;
						nbErrors++;
					}
			}
			else if (inDensDep == 1) {

				if (inS0 <= 0.0 || inS0 > 1.0) {
					BatchError(whichFile, whichLine, 20, "S0"); 
					nbErrors++;
				}
				// alphaS and betaS can take any value
			}
		}
		// read next simulation
		whichLine++;
		simNb = -98765;
		ifsSettlementFile >> simNb;
		if (ifsSettlementFile.eof())
			simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation

	if (current.simLines != current.reqdSimLines) {
		BatchError(whichFile, whichLine, 0, " "); 
		nbErrors++;
		batchLogOfs << gNbLinesStr << current.simNb
			<< gShouldBeStr << current.reqdSimLines << endl;
	}
	if (!ifsSettlementFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}

	if (nbErrors > 0) return -111;
	else return nbSims;
}

//---------------------------------------------------------------------------
int CheckTraitsFile(string indir)
{
	string header, colheader;
	int simNb, inSp, nextLineSimNb, nextLineSp;
	string filename, inTraitType, inSex, inInitDist, inInitParams, inInitDomDist, 
		inInitDomParams, inDominanceDist, inDominanceParams, inIsInherited, inMutationDist, 
		inMutationParams, inPositions, inNbPositions, inExpressionType, inMutationRate, inIsOutput;
	int nbErrors = 0;
	int nbSims = 0;
	int nbGenLoadTraits = 0;
	const string whichInputFile = "TraitsFile";
	vector<TraitType> allReadTraits;

	// Parse header line
	ifsTraitsFile >> header; if (header != "Simulation") nbErrors++;
	ifsTraitsFile >> header; if (header != "Species") nbErrors++;
	ifsTraitsFile >> header; if (header != "TraitType") nbErrors++;
	ifsTraitsFile >> header; if (header != "ExprSex") nbErrors++;
	ifsTraitsFile >> header; if (header != "Positions") nbErrors++;
	ifsTraitsFile >> header; if (header != "NbrOfPositions") nbErrors++;
	ifsTraitsFile >> header; if (header != "ExpressionType") nbErrors++;
	ifsTraitsFile >> header; if (header != "InitialAlleleDist") nbErrors++;
	ifsTraitsFile >> header; if (header != "InitialAlleleParams") nbErrors++;
	ifsTraitsFile >> header; if (header != "InitialDomDist") nbErrors++;
	ifsTraitsFile >> header; if (header != "InitialDomParams") nbErrors++;
	ifsTraitsFile >> header; if (header != "IsInherited") nbErrors++;
	ifsTraitsFile >> header; if (header != "MutationDistribution") nbErrors++;
	ifsTraitsFile >> header; if (header != "MutationParameters") nbErrors++;
	ifsTraitsFile >> header; if (header != "DominanceDistribution") nbErrors++;
	ifsTraitsFile >> header; if (header != "DominanceParameters") nbErrors++;
	ifsTraitsFile >> header; if (header != "MutationRate") nbErrors++;
	ifsTraitsFile >> header; if (header != "OutputValues") nbErrors++;

	if (nbErrors > 0) {
		FormatError(whichInputFile, nbErrors);
		return -111;
	}

	// Parse data lines
	int lineNb = 1;
	simCheck current, prev;		
	constexpr int simNbNotRead = -98765;
	simNb = simNbNotRead;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;

	ifsTraitsFile >> simNb >> inSp;

	bool stopReading = (simNb == simNbNotRead);
	int nbRowsToRead = 0;

	while (!stopReading) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in SimFile" << endl;
			nbErrors++;
		}

		if (!gSpInputOpt.at(simNb).contains(inSp)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Species number " << to_string(inSp) << " doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(inSp);

		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		ifsTraitsFile >> inTraitType >> inSex >> inPositions >> inNbPositions 
			>> inExpressionType >> inInitDist >> inInitParams >> inInitDomDist >> inInitDomParams
			>> inIsInherited >> inMutationDist >> inMutationParams >> inDominanceDist >> inDominanceParams
			>> inMutationRate >> inIsOutput;

		current = CheckStageSex(whichInputFile, lineNb, simNb, inSp, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) nbSims++;
		nbErrors += current.errors;
		prev = current;
		nbRowsToRead++;

		////  Validate parameters

		// Check sex is valid
		sex_t sex = stringToSex(inSex);
		if (sex == sex_t::INVALID_SEX) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << inSex << " is invalid: ExprSex must be either female, male, or # (if not applicable)." << endl;
			nbErrors++;
		}

		// Check trait type is legal
		TraitType tr = stringToTraitType(inTraitType);
		if (tr == TraitType::INVALID_TRAIT) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << inTraitType << " is not a valid TraitType." << endl;
			nbErrors++;
		}
		// Can trait be sex-dependent?
		const bool canBeSexDep = tr == E_D0 || tr == E_ALPHA || tr == E_BETA
			|| tr == S_S0 || tr == S_ALPHA || tr == S_BETA
			|| tr == KERNEL_MEANDIST_1 || tr == KERNEL_MEANDIST_2
			|| tr == KERNEL_PROBABILITY;
		if (!canBeSexDep && (sex == FEM || sex == MAL)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << inTraitType << " cannot be sex-dependent so ExprSex must be left blank (#)." << endl;
			nbErrors++;
		}
		if (sex != NA) // add sex to trait if present
			tr = addSexDepToTrait(tr, sex);

		// There can be up to 5 genetic load traits
		if (tr == GENETIC_LOAD) {
			nbGenLoadTraits++;
			if (nbGenLoadTraits > 5) {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "There cannot be more than 5 genetic load traits." << endl;
				nbErrors++;
			}
		}
		else if (traitExists(tr, allReadTraits)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Trait " << to_string(tr) << " is supplied multiple times." << endl;
			nbErrors++;
		}
		allReadTraits.push_back(tr);

		// Check Positions and NbrOfPositions
		const regex patternPositions{ "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$" };
		bool isMatch = regex_search(inPositions, patternPositions);
		if (!isMatch && inPositions != "random") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Positions must be either a comma-separated list of integer ranges, or random." << endl;
			nbErrors++;
		}
		if (inPositions == "random") {
			if (stoi(inNbPositions) <= 0) {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "NbrOfPositions must be a strictly positive integrer." << endl;
				nbErrors++;
			}
		}
		else if (inNbPositions != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "If Positions is not random NbrOfPositions must be blank (#)." << endl;
			nbErrors++;
		}

		// Check ExpressionType
		if (tr == NEUTRAL && inExpressionType != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "ExpressionType must be left blank (#) for the neutral trait." << endl;
			nbErrors++;
		}
		if (tr == GENETIC_LOAD && inExpressionType != "multiplicative") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "ExpressionType must be \"multiplicative\" for genetic load traits." << endl;
			nbErrors++;
		}
		const bool isDisp = tr != NEUTRAL && tr != GENETIC_LOAD && tr != INVALID_TRAIT;
		if (isDisp && inExpressionType != "additive" && inExpressionType != "average") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "ExpressionType must be \"additive\" or \"average\" for dispersal traits." << endl;
			nbErrors++;
		}

		// Check InitialAlleleDist
		if (tr == NEUTRAL && inInitDist != "uniform") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "InitialAlleleDist must be uniform for the neutral trait." << endl;
			nbErrors++;
		}
		if (isDisp && inInitDist != "normal" && inInitDist != "uniform") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "InitialAlleleDist must be either normal or uniform for dispersal traits." << endl;
			nbErrors++;
		}

		// Check InitialAlleleParams
		const regex patternParamsUnif{ "^\"?min=[-]?([0-9]*[.])?[0-9]+,max=[-]?([0-9]*[.])?[0-9]+\"?$" };
		const regex patternParamsNormal{ "^\"?mean=[-]?([0-9]*[.])?[0-9]+,sd=[-]?([0-9]*[.])?[0-9]+\"?$" };
		const regex patternParamsGamma{ "^\"?shape=[-]?([0-9]*[.])?[0-9]+,scale=[-]?([0-9]*[.])?[0-9]+\"?$" };
		const regex patternParamsMean{ "^\"?mean=[-]?([0-9]*[.])?[0-9]+\"?$" };
		const regex patternParamsNeutral{ "^\"?max=[0-9]+\"?$" };

		if (tr == NEUTRAL) {
			if (inInitDist == "uniform") {
				isMatch = regex_search(inInitParams, patternParamsNeutral);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For neutral trait with uniform initialisation, InitialAlleleParams must have form max=int" << endl;
					nbErrors++;
				}
				else {
					const int maxVal = stoi(inInitParams.substr(4));
					if (maxVal > 255) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLogOfs << "For neutral trait with uniform initialisation, max parameter must be between 0 and 255." << endl;
						nbErrors++;
					}
				}
			}
			// if not uniform then initDist must be blank, no params
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "For neutral trait with uniform initialisation, InitialAlleleParams must have form max=int" << endl;
				nbErrors++;
			}
		}

		if (isDisp) {
			if (inInitDist == "uniform") {
				isMatch = regex_search(inInitParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For dispersal trait uniform initialisation, InitialAlleleParams must have form min=float,max=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDist == "normal") {
				isMatch = regex_search(inInitParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For normal initialisation, InitialAlleleParams must have form mean=float,sd=float" << endl;
					nbErrors++;
				}
			}
		}
		if (tr == GENETIC_LOAD) {
			if (inInitDist == "uniform") {
				isMatch = regex_search(inInitParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a uniform distribution, InitialAlleleParams must have form min=float,max=float." << endl;
					nbErrors++;
				}
			}
			else if (inInitDist == "normal") {
				isMatch = regex_search(inInitParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a normal distribution, InitialAlleleParams must have form mean=float,sd=float." << endl;
					nbErrors++;
				}
			}
			else if (inInitDist == "gamma") {
				isMatch = regex_search(inInitParams, patternParamsGamma);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a Gamma distribution, InitialAlleleParams must have form shape=float,scale=float." << endl;
					nbErrors++;
				}
			}
			else if (inInitDist == "negExp") {
				isMatch = regex_search(inInitParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a negative exponential distribution, InitialAlleleParams must have form mean=float." << endl;
					nbErrors++;
				}
			}
			else if (inInitDist == "#") {
				if (inInitParams != "#") {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "If InitialAlleleDist is left blank, InitialAlleleParams must also be blank." << endl;
					nbErrors++;
				}
				// otherwise fine!
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "For genetic load traits, InitialAlleleDist must be either blank (#), uniform, gamma, negExp or normal" << endl;
				nbErrors++;
			}
		}

		// Check InitialDomDist and InitialDomParams
		if ((isDisp || tr == NEUTRAL)
			&& (inInitDomDist != "#" || inInitDomParams != "#")){
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "InitialDomDist and InitialDomParams must be blank (#) for dispersal and neutral traits." << endl;
			nbErrors++;
		}
		else if (tr == GENETIC_LOAD) {
			if (inInitDomDist == "normal") {
				isMatch = regex_search(inInitDomParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a normal dominance distribution, InitialDomParams must have form mean=float,sd=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDomDist == "gamma") {
				isMatch = regex_search(inInitDomParams, patternParamsGamma);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a Gamma dominance distribution, InitialDomParams must have form shape=float,scale=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDomDist == "uniform") {
				isMatch = regex_search(inInitDomParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a uniform dominance distribution, InitialDomParams must have form min=float,max=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDomDist == "negExp") {
				isMatch = regex_search(inInitDomParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a negative exponential dominance distribution, InitialDomParams must have form mean=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDomDist == "scaled") {
				if (inInitDist == "#") {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "Initial scaled dominance distribution requires InitialAlleleDist to be non-blank." << endl;
					nbErrors++;
				}
				isMatch = regex_search(inInitDomParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a scaled dominance distribution, InitialDomParams must have form mean=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDomDist == "#") {
				if (inInitDomParams != "#") {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "If InitialDomDist is left blank, InitialDomParams must also be blank." << endl;
					nbErrors++;
				}
				// otherwise fine
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "InitialDomDist must be either blank (#), normal, gamma, uniform, negExp or scaled for genetic load traits." << endl;
				nbErrors++;
			}
		}

		// Check isInherited and MutationRate
		if ((tr == NEUTRAL || tr == GENETIC_LOAD) && inIsInherited != "TRUE") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "isInherited must always be TRUE for neutral and genetic load traits." << endl;
			nbErrors++;
		}
		else if (isDisp) {
			if (inIsInherited != "TRUE" && inIsInherited != "FALSE") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "IsInherited must be either TRUE or FALSE for dispersal traits." << endl;
				nbErrors++;
			}
		}
		if ((inIsInherited == "TRUE") 
			&& (stof(inMutationRate) < 0.0 || stof(inMutationRate) > 1.0)) {
			BatchError(whichInputFile, lineNb, 20, "mutationRate"); 
			nbErrors++;
		}
		else if (inIsInherited == "FALSE" && inMutationRate != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "If isInherited if off, mutationRate must be blank (#)." << endl;
			nbErrors++;
		}

		// Check MutationDistribution and MutationParameters
		if (tr == NEUTRAL) {
			if (inMutationDist == "KAM" || inMutationDist == "SSM") {
				isMatch = regex_search(inMutationParams, patternParamsNeutral);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a neutral trait, mutationParams must have form max=int." << endl;
					nbErrors++;
				}
				else {
					const int maxVal = stoi(inMutationParams.substr(4));
					if (maxVal > 255) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLogOfs << "For the neutral trait mutation max parameter must be between 0 and 255." << endl;
						nbErrors++;
					}
				}
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "For a neutral trait, mutationDistribution must be either KAM or SSM." << endl;
				nbErrors++;
			}
		}
		if (isDisp) {
			if (inIsInherited == "TRUE") {
				if (inMutationDist == "uniform") {
					isMatch = regex_search(inMutationParams, patternParamsUnif);
					if (!isMatch) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLogOfs << "For a uniform distribution, mutationParams must have form min=float,max=float." << endl;
						nbErrors++;
					}
				}
				else if (inMutationDist == "normal") {
					isMatch = regex_search(inMutationParams, patternParamsNormal);
					if (!isMatch) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLogOfs << "For a normal distribution, mutationParams must have form mean=float,sd=float." << endl;
						nbErrors++;
					}
				}
				else {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For dispersal traits, mutationDistribution must be either uniform or normal" << endl;
					nbErrors++;
				}
			}
			else { // not inherited
				if (inMutationDist != "#" || inMutationParams != "#") {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "If isInherited is turned off, mutationDistribution and mutationParameters must be left blank (#)." << endl;
					nbErrors++;
				}
			}
		}
		if (tr == GENETIC_LOAD) {
			if (inMutationDist == "uniform") {
				isMatch = regex_search(inMutationParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a uniform distribution, mutationParams must have form min=float,max=float." << endl;
					nbErrors++;
				}
			}
			else if (inMutationDist == "normal") {
				isMatch = regex_search(inMutationParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a normal distribution, mutationParams must have form mean=float,sd=float." << endl;
					nbErrors++;
				}
			}
			else if (inMutationDist == "gamma") {
				isMatch = regex_search(inMutationParams, patternParamsGamma);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a Gamma distribution, mutationParams must have form shape=float,scale=float." << endl;
					nbErrors++;
				}
			}
			else if (inMutationDist == "negExp") {
				isMatch = regex_search(inMutationParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a negative exponential distribution, mutationParams must have form mean=float." << endl;
					nbErrors++;
				}
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "For genetic load traits, mutationDistribution must be either uniform, gamma, negExp or normal" << endl;
				nbErrors++;
			}
		}

		// Check DominanceDistribution and DominanceParameters
		if (tr == NEUTRAL) {
			if (inDominanceDist != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "DominanceDistribution must be left blank (#) for the neutral trait." << endl;
				nbErrors++;
			}
			if (inDominanceParams != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "DominanceParameters must be left blank (#) for the neutral trait." << endl;
				nbErrors++;
			}
		}
		if (isDisp) {
			if (inDominanceDist != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "DominanceDistribution must be left blank (#) for dispersal traits." << endl;
				nbErrors++;
			}
			if (inDominanceParams != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "DominanceParameters must be left blank (#) for dispersal traits." << endl;
				nbErrors++;
			}
		}
		if (tr == GENETIC_LOAD) {
			if (inDominanceDist == "normal") {
				isMatch = regex_search(inDominanceParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a normal dominance distribution, DominanceParams must have form mean=float,sd=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "gamma") {
				isMatch = regex_search(inDominanceParams, patternParamsGamma);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a Gamma dominance distribution, DominanceParams must have form shape=float,scale=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "uniform") {
				isMatch = regex_search(inDominanceParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a uniform dominance distribution, DominanceParams must have form min=float,max=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "negExp") {
				isMatch = regex_search(inDominanceParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a negative exponential dominance distribution, DominanceParams must have form mean=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "scaled") {
				isMatch = regex_search(inDominanceParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLogOfs << "For a scaled dominance distribution, DominanceParams must have form mean=float" << endl;
					nbErrors++;
				}
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLogOfs << "DominanceDistribution must be either normal, gamma, uniform, negExp or scaled for genetic load traits." << endl;
				nbErrors++;
			}
		}

		if (inIsOutput != "TRUE" && inIsOutput != "FALSE") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "OutputValues must be either TRUE or FALSE." << endl;
			nbErrors++;
		}

		// Preview next line
		nextLineSimNb = simNbNotRead;
		ifsTraitsFile >> nextLineSimNb >> nextLineSp;

		if (nextLineSimNb == simNbNotRead
			|| ifsTraitsFile.eof()) {
			// Exit loop
			stopReading = true;
			nbErrors += checkTraitSetCoherency(allReadTraits, simNb, inSp);
			inputOpt.nbTraitFileRows = nbRowsToRead;
		}
		else if (nextLineSimNb != simNb || nextLineSp != inSp) {
			// About to change sim or species, conduct checks of all read traits
			nbErrors += checkTraitSetCoherency(allReadTraits, simNb, inSp);
			// Store nb of rows to help reading file later on
			inputOpt.nbTraitFileRows = nbRowsToRead;
			nbRowsToRead = 0; // reset for next sim or species
			nbGenLoadTraits = 0;
			allReadTraits.clear();
			simNb = nextLineSimNb;
			inSp = nextLineSp;
		} // else continue reading traits for same sim
		lineNb++; 
	} // end of while loop

	if (!ifsTraitsFile.eof()) {
		EOFerror(whichInputFile);
		nbErrors++;
	}

	if (nbErrors > 0) 
		return -111;
	else return 0;
}

int checkTraitSetCoherency(const vector <TraitType>& allReadTraits, const int& simNb, const species_id& sp) {
	int nbErrors = 0;
	const string whichInputFile = "TraitsFile";

	const spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(sp);

	if (inputOpt.anyNeutral && !traitExists(NEUTRAL, allReadTraits)) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLogOfs << "Neutral statistics enabled but neutral trait is missing." << endl;
		nbErrors++;
	}

	//// Check dispersal traits and sex-dependencies are complete 
	// and consistent with parameters in dispersal input files

	// Emigration traits
	bool hasD0 = traitExists(E_D0, allReadTraits) || traitExists(E_D0_F, allReadTraits) || traitExists(E_D0_M, allReadTraits);
	bool hasEmigAlpha = (traitExists(E_ALPHA, allReadTraits) || traitExists(E_ALPHA_F, allReadTraits) || traitExists(E_ALPHA_M, allReadTraits));
	bool hasEmigBeta = (traitExists(E_BETA, allReadTraits) || traitExists(E_BETA_F, allReadTraits) || traitExists(E_BETA_M, allReadTraits));

	bool anyEmigNeitherSex = traitExists(E_D0, allReadTraits) || traitExists(E_ALPHA, allReadTraits) || traitExists(E_BETA, allReadTraits);
	bool eitherSexD0 = traitExists(E_D0_F, allReadTraits) || traitExists(E_D0_M, allReadTraits);
	bool bothSexesD0 = traitExists(E_D0_F, allReadTraits) && traitExists(E_D0_M, allReadTraits);
	bool eitherSexEmigAlpha = traitExists(E_ALPHA_F, allReadTraits) || traitExists(E_ALPHA_M, allReadTraits);
	bool bothSexesEmigAlpha = traitExists(E_ALPHA_F, allReadTraits) && traitExists(E_ALPHA_M, allReadTraits);
	bool eitherSexEmigBeta = traitExists(E_BETA_F, allReadTraits) || traitExists(E_BETA_M, allReadTraits);
	bool bothSexesEmigBeta = traitExists(E_BETA_F, allReadTraits) && traitExists(E_BETA_M, allReadTraits);
	bool anyEmigSexDep = eitherSexD0 || eitherSexEmigAlpha || eitherSexEmigBeta;

	if (inputOpt.isEmigIndVar) {
		if (!hasD0) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "EP or d0 is missing." << endl;
			nbErrors++;
		}
		if (inputOpt.isEmigSexDep) {
			if (anyEmigNeitherSex) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Emigration SexDep is on but a trait has been supplied without a sex." << endl;
				nbErrors++;
			}
			if (!bothSexesD0) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Either sex is missing for D0 trait." << endl;
				nbErrors++;
			}
		}
		else if (anyEmigSexDep) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Emigration SexDep is off but a trait has been supplied with a sex." << endl;
			nbErrors++;
		}

		if (inputOpt.isEmigDensDep) {
			if (!hasEmigAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Emigration alpha is missing." << endl;
				nbErrors++;
			}
			if (!hasEmigBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Emigration beta is missing." << endl;
				nbErrors++;
			}
			if (inputOpt.isEmigSexDep) {
				if (!bothSexesEmigAlpha) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLogOfs << "Either sex is missing for emigration alpha trait." << endl;
					nbErrors++;
				}
				if (!bothSexesEmigBeta) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLogOfs << "Either sex is missing for emigration beta trait." << endl;
					nbErrors++;
				}
			}
		}
		else {
			if (hasEmigAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Specified emigration alpha, but emigration is not density-dependent." << endl;
				nbErrors++;
			}
			if (hasEmigBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Specified emigration beta, but emigration is not density-dependent." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasD0 || hasEmigAlpha || hasEmigBeta) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLogOfs << "Specified emigration trait, but emigration is not variable." << endl;
		nbErrors++;
	}

	// Transfer traits
	/// Kernels
	bool hasKern1 = traitExists(KERNEL_MEANDIST_1, allReadTraits) || traitExists(KERNEL_MEANDIST_1_F, allReadTraits) || traitExists(KERNEL_MEANDIST_1_M, allReadTraits);
	bool hasKern2 = traitExists(KERNEL_MEANDIST_2, allReadTraits) || traitExists(KERNEL_MEANDIST_2_F, allReadTraits) || traitExists(KERNEL_MEANDIST_2_M, allReadTraits);
	bool hasKernProb = traitExists(KERNEL_PROBABILITY, allReadTraits) || traitExists(KERNEL_PROBABILITY_F, allReadTraits) || traitExists(KERNEL_PROBABILITY_M, allReadTraits);

	bool anyKernelNeitherSex = traitExists(KERNEL_MEANDIST_1, allReadTraits) || traitExists(KERNEL_MEANDIST_2, allReadTraits) || traitExists(KERNEL_PROBABILITY, allReadTraits);
	bool eitherSexMeanDist1 = traitExists(KERNEL_MEANDIST_1_F, allReadTraits) || traitExists(KERNEL_MEANDIST_1_M, allReadTraits);
	bool bothSexesMeanDist1 = traitExists(KERNEL_MEANDIST_1_F, allReadTraits) && traitExists(KERNEL_MEANDIST_1_M, allReadTraits);
	bool eitherSexMeanDist2 = traitExists(KERNEL_MEANDIST_2_F, allReadTraits) || traitExists(KERNEL_MEANDIST_2_M, allReadTraits);
	bool bothSexesMeanDist2 = traitExists(KERNEL_MEANDIST_2_F, allReadTraits) && traitExists(KERNEL_MEANDIST_2_F, allReadTraits);
	bool eitherSexKernProb = traitExists(KERNEL_PROBABILITY_F, allReadTraits) || traitExists(KERNEL_PROBABILITY_M, allReadTraits);
	bool bothSexesKernProb = traitExists(KERNEL_PROBABILITY_F, allReadTraits) && traitExists(KERNEL_PROBABILITY_M, allReadTraits);
	bool anyKernelSexDep = eitherSexMeanDist1 || eitherSexMeanDist2 || eitherSexKernProb;

	if (inputOpt.isKernTransfIndVar) {
		if (!hasKern1) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "(First) kernel mean is missing." << endl;
			nbErrors++;
		}
		if (inputOpt.isKernTransfSexDep) {
			if (anyKernelNeitherSex) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Kernel SexDep is on but a trait has been supplied without a sex." << endl;
				nbErrors++;
			}
			if (!bothSexesMeanDist1) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Either sex is missing for first kernel mean trait." << endl;
				nbErrors++;
			}
		}
		else if (anyKernelSexDep) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Kernel SexDep is off but a trait has been supplied with a sex." << endl;
			nbErrors++;
		}
		if (inputOpt.usesTwoKernels) {
			if (!hasKern2) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Second kernel mean is missing." << endl;
				nbErrors++;
			}
			if (!hasKernProb) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Kernel probability is missing." << endl;
				nbErrors++;
			}
			if (inputOpt.isKernTransfSexDep) {
				if (!bothSexesMeanDist2) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLogOfs << "Either sex is missing for second kernel mean trait." << endl;
					nbErrors++;
				}
				if (!bothSexesKernProb) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLogOfs << "Either sex is missing for kernel probability trait." << endl;
					nbErrors++;
				}
			}
		}
		else {
			if (hasKern2) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Specified second kernel, but only one kernel is used." << endl;
				nbErrors++;
			}
			if (hasKernProb) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Specified kernel probability, but only one kernel is used." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasKern1 || hasKern2 || hasKernProb) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLogOfs << "Specified kernel transfer trait, but kernel transfer is not variable." << endl;
		nbErrors++;
	}

	/// SMS
	bool hasDP = traitExists(SMS_DP, allReadTraits);
	bool hasGB = traitExists(SMS_GB, allReadTraits);
	bool hasSMSAlpha = traitExists(SMS_ALPHADB, allReadTraits);
	bool hasSMSBeta = traitExists(SMS_BETADB, allReadTraits);
	if (inputOpt.isSMSTransfIndVar) {
		if (!hasDP) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "SMS directional persistence trait is missing." << endl;
			nbErrors++;
		}
		if (inputOpt.usesSMSGoalBias) {
			if (!hasGB) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "SMS goal bias trait is missing." << endl;
				nbErrors++;
			}
			if (!hasSMSAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "SMS alpha direction bias trait is missing." << endl;
				nbErrors++;
			}
			if (!hasSMSBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "SMS beta direction bias trait is missing." << endl;
				nbErrors++;
			}
		}
		else {
			if (hasGB) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "SMS goal bias trait supplied, but SMS GoalType not set to option 2." << endl;
				nbErrors++;
			}
			if (hasSMSAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "SMS alpha direction bias trait supplied, but SMS GoalType not set to option 2." << endl;
				nbErrors++;
			}
			if (hasSMSBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "SMS beta direction bias trait supplied, but SMS GoalType not set to option 2." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasDP || hasGB || hasSMSAlpha || hasSMSBeta) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLogOfs << "Specified SMS trait, but SMS not set to be variable." << endl;
		nbErrors++;
	}

	/// CRW
	bool hasStepLen = traitExists(CRW_STEPLENGTH, allReadTraits);
	bool hasRho = traitExists(CRW_STEPCORRELATION, allReadTraits);
	if (inputOpt.isCRWTransfIndVar) {
		if (!hasStepLen) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "CRW step length trait is missing." << endl;
			nbErrors++;
		}
		if (!hasRho) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "CRW step correlation trait is missing." << endl;
			nbErrors++;
		}
	}
	else if (hasStepLen || hasRho) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLogOfs << "Specified CRW trait, but CRW not set to be variable." << endl;
		nbErrors++;
	}

	// Settlement traits
	bool hasS0 = traitExists(S_S0, allReadTraits) || traitExists(S_S0_F, allReadTraits) || traitExists(S_S0_M, allReadTraits);
	bool hasSettAlpha = traitExists(S_ALPHA, allReadTraits) || traitExists(S_ALPHA_F, allReadTraits) || traitExists(S_ALPHA_M, allReadTraits);
	bool hasSettBeta = traitExists(S_BETA, allReadTraits) || traitExists(S_BETA_F, allReadTraits) || traitExists(S_BETA_M, allReadTraits);

	bool anySettNeitherSex = traitExists(S_S0, allReadTraits) || traitExists(S_ALPHA, allReadTraits) || traitExists(S_BETA, allReadTraits);
	bool eitherSexS0 = traitExists(S_S0_F, allReadTraits) || traitExists(S_S0_M, allReadTraits);
	bool bothSexesS0 = traitExists(S_S0_F, allReadTraits) && traitExists(S_S0_M, allReadTraits);
	bool eitherSexSettAlpha = traitExists(S_ALPHA_F, allReadTraits) || traitExists(S_ALPHA_M, allReadTraits);
	bool bothSexesSettAlpha = traitExists(S_ALPHA_F, allReadTraits) && traitExists(S_ALPHA_M, allReadTraits);
	bool eitherSexSettBeta = traitExists(S_BETA_F, allReadTraits) || traitExists(S_BETA_M, allReadTraits);
	bool bothSexesSettBeta = traitExists(S_BETA_F, allReadTraits) && traitExists(S_BETA_M, allReadTraits);
	bool anySettSexDep = eitherSexS0 || eitherSexSettAlpha || eitherSexSettBeta;

	if (inputOpt.isSettIndVar) {
		if (!hasS0) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Settlement probability trait is missing." << endl;
			nbErrors++;
		}
		if (inputOpt.isSettSexDep) {
			if (anySettNeitherSex) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Settlement SexDep is on but a trait has been supplied without a sex." << endl;
				nbErrors++;
			}
			if (!bothSexesS0) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Either sex is missing for settlement probabibility trait." << endl;
				nbErrors++;
			}
		}
		else if (anySettSexDep) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Settlement SexDep is off but a trait has been supplied with a sex." << endl;
			nbErrors++;
		}
		// if settlement is IndVar, it is always density-dependent
		if (!hasSettAlpha) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Settlement alpha trait is missing." << endl;
			nbErrors++;
		}
		if (!hasSettBeta) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Settlement beta trait is missing." << endl;
			nbErrors++;
		}
		if (inputOpt.isSettSexDep) {
			if (!bothSexesSettAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Either sex is missing for settlement alpha trait." << endl;
				nbErrors++;
			}
			if (!bothSexesSettBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLogOfs << "Either sex is missing for settlement beta trait." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasS0 || hasSettAlpha || hasSettBeta) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLogOfs << "Specified settlement trait, but settlement not set to be variable." << endl;
		nbErrors++;
	}
	return nbErrors;
}

bool traitExists(const TraitType& tr, const vector<TraitType>& existingTraits) {
	return std::find(existingTraits.begin(), existingTraits.end(), tr) != existingTraits.end();
}

TraitType addSexDepToTrait(const TraitType& t, const sex_t& sex) {
	if (sex == FEM) {
		if (t == E_D0) return E_D0_F; // EP uses d0 for trait data
		else if (t == E_ALPHA) return E_ALPHA_F;
		else if (t == E_BETA) return E_BETA_F;
		else if (t == S_S0) return S_S0_F;
		else if (t == S_ALPHA) return S_ALPHA_F;
		else if (t == S_BETA) return S_BETA_F;
		else if (t == KERNEL_MEANDIST_1) return KERNEL_MEANDIST_1_F;
		else if (t == KERNEL_MEANDIST_2) return KERNEL_MEANDIST_2_F;
		else if (t == KERNEL_PROBABILITY) return KERNEL_PROBABILITY_F;
		else return INVALID_TRAIT;
	}
	else if (sex == MAL) {
		if (t == E_D0) return E_D0_M; // EP uses d0 for trait data
		else if (t == E_ALPHA) return E_ALPHA_M;
		else if (t == E_BETA) return E_BETA_M;
		else if (t == S_S0) return S_S0_M;
		else if (t == S_ALPHA) return S_ALPHA_M;
		else if (t == S_BETA) return S_BETA_M;
		else if (t == KERNEL_MEANDIST_1) return KERNEL_MEANDIST_1_M;
		else if (t == KERNEL_MEANDIST_2) return KERNEL_MEANDIST_2_M;
		else if (t == KERNEL_PROBABILITY) return KERNEL_PROBABILITY_M;
		else return INVALID_TRAIT;
	}
	else return INVALID_TRAIT;
}

//---------------------------------------------------------------------------

int CheckGeneticsFile(string inputDirectory) {

	string header;
	int simNb, inSp, prevSimNb, errCode;
	string inChromosomeEnds, inRecombinationRate, inTraitsFile, inPatchList, inStages,
		inOutGeneValues, inOutWeirCockerham, inOutWeirHill,
		inOutStartGenetics, inOutputInterval, inNbrPatchesToSample, inNIndsToSample;
	int inGenomeSize;
	int nbErrors = 0;
	int nbSims = 0;
	string whichFile = "GeneticsFile";

	const regex patternIntList{ "^\"?([0-9]+,)*[0-9]+\"?$" }; // comma-separated integer list
	bool isMatch = false;

	// Parse header line;
	ifsGeneticsFile >> header; if (header != "Simulation") nbErrors++;
	ifsGeneticsFile >> header; if (header != "Species") nbErrors++;
	ifsGeneticsFile >> header; if (header != "GenomeSize") nbErrors++;
	ifsGeneticsFile >> header; if (header != "ChromosomeEnds") nbErrors++;
	ifsGeneticsFile >> header; if (header != "RecombinationRate") nbErrors++;
	ifsGeneticsFile >> header; if (header != "OutputGeneValues") nbErrors++;
	ifsGeneticsFile >> header; if (header != "OutputFstatsWeirCockerham") nbErrors++;
	ifsGeneticsFile >> header; if (header != "OutputFstatsWeirHill") nbErrors++;
	ifsGeneticsFile >> header; if (header != "OutputStartGenetics") nbErrors++;
	ifsGeneticsFile >> header; if (header != "OutputInterval") nbErrors++;
	ifsGeneticsFile >> header; if (header != "PatchList") nbErrors++;
	ifsGeneticsFile >> header; if (header != "NbrPatchesToSample") nbErrors++;
	ifsGeneticsFile >> header; if (header != "nIndividualsToSample") nbErrors++;
	ifsGeneticsFile >> header; if (header != "Stages") nbErrors++;

	if (nbErrors > 0) {
		FormatError(whichFile, nbErrors);
		return -111;
	}

	// Parse data lines
	int whichLine = 1;
	simNb = -98765;
	ifsGeneticsFile >> simNb;
	while (simNb != -98765) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}

		ifsGeneticsFile >> inSp;
		if (!gSpInputOpt.at(simNb).contains(inSp)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Species number " << to_string(inSp) << " doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(inSp);

		ifsGeneticsFile >> inGenomeSize >> inChromosomeEnds >> inRecombinationRate >> inOutGeneValues >> inOutWeirCockerham >>
			inOutWeirHill >> inOutStartGenetics >> inOutputInterval >> inPatchList >> inNbrPatchesToSample
			>> inNIndsToSample >> inStages;

		//// Validate parameters
		
		// Check GenomeSize
		if (inGenomeSize <= 0) {
			BatchError(whichFile, whichLine, 10, "GenomeSize");
			nbErrors++;
		}
		
		// Check ChromosomeEnds
		isMatch = regex_search(inChromosomeEnds, patternIntList);
		if (!isMatch && inChromosomeEnds != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "ChromosomeEnds must be either a comma-separated list of integers, or blank (#)." << endl;
			nbErrors++;
		}
		set<int> chrEnds = stringToChromosomeEnds(inChromosomeEnds, inGenomeSize);
		const int maxVal = *chrEnds.rbegin();
		if (maxVal >= inGenomeSize) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Positions for ChromosomeEnds cannot exceed GenomeSize." << endl;
			nbErrors++;
		}

		// Check RecombinationRate
		if (inputOpt.reproType == 0 && inRecombinationRate != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Do not specify a recombination rate for haploid/asexual systems." << endl;
			nbErrors++;
		}
		else if (inRecombinationRate != "#") {
			float recombinationRate = stof(inRecombinationRate);
			if (recombinationRate < 0.0 || recombinationRate > 0.5) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "RecombinationRate must be positive and not exceed 0.5." << endl;
				nbErrors++;
			}
		}

		// Check genetic output fields
		if (inOutGeneValues != "TRUE" && inOutGeneValues != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "OutGeneValues must be either TRUE or FALSE" << endl;
			nbErrors++;
		}
		if (inOutWeirCockerham != "TRUE" && inOutWeirCockerham != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "OutputFstatsWeirCockerham must be either TRUE or FALSE" << endl;
			nbErrors++;
		}
		if (inOutWeirHill != "TRUE" && inOutWeirHill != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "OutputFstatsWeirHill must be either TRUE or FALSE" << endl;
			nbErrors++;
		}
		inputOpt.anyNeutral = inOutWeirCockerham == "TRUE"
			|| inOutWeirHill == "TRUE";
		bool anyGeneticsOutput = inOutGeneValues == "TRUE" 
			|| inputOpt.anyNeutral;

		if (anyGeneticsOutput) {
			if (inOutStartGenetics == "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "OutStartGenetics cannot be left blank (#) if any genetic output option is TRUE." << endl;
				nbErrors++;
			}
			else {
				int outStartGenetics = stoi(inOutStartGenetics);
				if (outStartGenetics < 0) {
					BatchError(whichFile, whichLine, 10, "OutStartGenetics");
					nbErrors++;
				}
			}
			if (inOutputInterval == "#" || inOutputInterval == "0") {
				// Minimum interval is 1, not 0
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "OutputInterval cannot be left blank (#) or 0 if any genetic output option is TRUE." << endl;
				nbErrors++;
			}
			else {
				int outputInterval = stoi(inOutputInterval);
				if (outputInterval < 0) {
					BatchError(whichFile, whichLine, 10, "OutputInterval");
					nbErrors++;
				}
			}
		} // no genetics output
		else {
			if (inOutStartGenetics != "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "OutStartGenetics should be blank (#) if all genetic output options are FALSE." << endl;
				nbErrors++;
			}
			if (inOutputInterval != "#" && inOutputInterval != "0") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "OutputInterval should be blank (#) or 0 if all genetic output options are FALSE." << endl;
				nbErrors++;
			}
		}

		// Check PatchList
		if (anyGeneticsOutput) {
			if (inPatchList == "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "PatchList cannot be left blank (#) if any genetic output option is TRUE." << endl;
				nbErrors++;
			}
			else {
				isMatch = regex_search(inPatchList, patternIntList);
				if (!isMatch && inPatchList != "random" && inPatchList != "all" && inPatchList != "random_occupied") {
					BatchError(whichFile, whichLine, 0, " ");
					batchLogOfs << "PatchList must be either a comma-separated list of integers, random, random_occupied or all." << endl;
					nbErrors++;
				}
			}
		}
		else if (inPatchList != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "PatchList should be blank (#) if all genetic output options are FALSE." << endl;
			nbErrors++;
		}

		// Check NbrPatchesToSample
		if (inPatchList == "random" || inPatchList == "random_occupied") {
			if (inNbrPatchesToSample == "#" || inNbrPatchesToSample == "0") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "NbrPatchesToSample cannot be blank (#) or 0 if PatchList is random or random_occupied." << endl;
				nbErrors++;
			}
			else {
				int nbPatches = stoi(inNbrPatchesToSample);
				if (nbPatches <= 0) {
					BatchError(whichFile, whichLine, 10, "NbrPatchesToSample");
					nbErrors++;
				}
			}
		}
		else if (inNbrPatchesToSample != "#" && inNbrPatchesToSample != "0") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "NbrPatchesToSample must be blank (#) or zero if PatchList is not random or random_occupied." << endl;
			nbErrors++;
		}

		// Check IndividualsToSample
		if (anyGeneticsOutput) {
			if (inNIndsToSample == "#" || inNIndsToSample == "0") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "NIndsToSample cannot be blank (#) or zero if any genetics output option is TRUE." << endl;
				nbErrors++;
			}
			else if (inNIndsToSample != "all") {
				int nIndsToSample = stoi(inNIndsToSample);
				if (nIndsToSample <= 0) {
					BatchError(whichFile, whichLine, 10, "nIndsToSample");
					nbErrors++;
				}
			}
		}
		else if (inNIndsToSample != "#" && inNIndsToSample != "0") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "NIndsToSample must be blank (#) or zero if all genetics output options are FALSE." << endl;
			nbErrors++;
		}

		// Check Stages
		if (anyGeneticsOutput) {
			if (inStages == "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "Stages cannot be blank (#) if any genetic output option is TRUE." << endl;
				nbErrors++;
			}
			else {
				isMatch = regex_search(inStages, patternIntList);
				if (!isMatch && inStages != "all") {
					BatchError(whichFile, whichLine, 0, " ");
					batchLogOfs << "Stages must be either a comma-separated list of integers, or \"all\"." << endl;
					nbErrors++;
				}
			}
		}
		else if (inStages != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Stages must be blank (#) if all genetic output options are FALSE." << endl;
			nbErrors++;
		}

		// read next simulation
		whichLine++;
		simNb = -98765;
		ifsGeneticsFile >> simNb;
		if (ifsGeneticsFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation

	if (!ifsGeneticsFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}

	if (nbErrors > 0) return -111;
	else return nbSims;

}

//---------------------------------------------------------------------------
int CheckInitFile(string indir)
{
	string header, colheader;
	int i, simNb, spNb;
	int seedtype, freetype, sptype, initdens, indscell = 0, minX, maxX, minY, maxY;
	int nCells, nSpCells, initAge;
	int initFreezeYear, restrictRows, restrictFreq, finalFreezeYear;
	float inds_per_ha;

	int errors = 0; int propnerrors = 0;
	int simuls = 0;
	string filetype = "InitialisationFile";

	int maxNbStages = 0;

	for (auto& [sim, inputOptSp] : gSpInputOpt)
		for (auto [sp, inputOpt] : inputOptSp)
			maxNbStages = max(inputOpt.nbStages, maxNbStages);

	// Parse header line;
	ifsInitFile >> header; if (header != "Simulation") errors++;
	ifsInitFile >> header; if (header != "Species") errors++;
	ifsInitFile >> header; if (header != "SeedType") errors++;
	ifsInitFile >> header; if (header != "FreeType") errors++;
	ifsInitFile >> header; if (header != "SpType") errors++;
	ifsInitFile >> header; if (header != "InitDens") errors++;
	ifsInitFile >> header;
	if (gUsesPatches) { if (header != "IndsHa") errors++; }
	else { if (header != "IndsCell") errors++; }
	ifsInitFile >> header; if (header != "minX") errors++;
	ifsInitFile >> header; if (header != "maxX") errors++;
	ifsInitFile >> header; if (header != "minY") errors++;
	ifsInitFile >> header; if (header != "maxY") errors++;
	ifsInitFile >> header; if (header != "NCells") errors++;
	ifsInitFile >> header; if (header != "NSpCells") errors++;
	ifsInitFile >> header; if (header != "InitFreezeYear") errors++;
	ifsInitFile >> header; if (header != "RestrictRows") errors++;
	ifsInitFile >> header; if (header != "RestrictFreq") errors++;
	ifsInitFile >> header; if (header != "FinalFreezeYear") errors++;
	ifsInitFile >> header; if (header != "InitIndsFile") errors++;
	ifsInitFile >> header; if (header != "InitAge") errors++;
	for (i = 1; i < maxNbStages; i++) {
		colheader = "PropStage" + to_string(i);
		ifsInitFile >> header; if (header != colheader) propnerrors++;
	}
	// report any errors in headers, and if so, terminate validation
	if (errors > 0 || propnerrors > 0) {
		FormatError(filetype, errors + propnerrors);
		if (propnerrors > 0) BatchError(filetype, -999, 444, "PropStage");
		return -111;
	}

	// Parse data lines
	int line = 1;
	int err;
	simCheck current, prev;
	bool checkfile;
	string filename, ftype2, fname;
	vector <string> indsfiles;
	ftype2 = "InitIndsFile";
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	ifsInitFile >> simNb;
	
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		if (!gSpInputOpt.contains(simNb)) {
			BatchError(filetype, line, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in SimFile" << endl;
			errors++;
		}

		ifsInitFile >> spNb;
		if (!gSpInputOpt.at(simNb).contains(spNb)) {
			BatchError(filetype, line, 0, " ");
			batchLogOfs << "Species number " << to_string(spNb) << " doesn't match those in ParametersFile" << endl;
			errors++;
		}
		spInputOptions& inputOpt = gSpInputOpt.at(simNb).at(spNb);

		current = CheckStageSex(filetype, line, simNb, spNb, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) simuls++;
		errors += current.errors;
		prev = current;

		ifsInitFile >> seedtype >> freetype >> sptype >> initdens >> inds_per_ha;
		if (!gUsesPatches) indscell = (int)inds_per_ha;
		if (seedtype < 0 || seedtype > 2) {
			BatchError(filetype, line, 2, "SeedType"); errors++;
		}
		if (gLandType == 9 && seedtype != 0) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLogOfs << "SeedType must be 0 for an artificial landscape"
				<< endl;
		}
		if (!gUseSpeciesDist.at(spNb)) {
			if (seedtype == 1) {
				BatchError(filetype, line, 0, " ");
				errors++;
				batchLogOfs << "SeedType is 1 but there is no species distribution map in the SpeciesLandFile"
					<< endl;
			}
		}
		else if (seedtype != 1) {
			BatchError(filetype, line, 0, " ");
			errors++;
			batchLogOfs << "Species distribution map specified in SpeciesLandFile, but SeedType is not 1."
				<< endl;
		}
		if (seedtype == 0) {
			if (freetype < 0 || freetype > 1) {
				BatchError(filetype, line, 1, "FreeType"); errors++;
			}
		}
		if (seedtype == 1) {
			if (sptype < 0 || sptype > 1) {
				BatchError(filetype, line, 1, "SpType"); errors++;
			}
		}
		if (initdens < 0 || initdens > 2) {
			BatchError(filetype, line, 2, "initDens"); errors++;
		}
		if (seedtype < 2) {
			if (initdens == 2) { // specified density
				if (gUsesPatches) {
					if (inds_per_ha <= 0.0) {
						BatchError(filetype, line, 10, "IndsHa"); errors++;
					}
				}
				else {
					if (indscell < 1) {
						BatchError(filetype, line, 11, "IndsCell"); errors++;
					}
				}
			}
		}

		ifsInitFile >> minX >> maxX >> minY >> maxY >> nCells >> nSpCells;
		if (seedtype == 0) {
			if (maxX < minX) {
				BatchError(filetype, line, 2, "maxX", "minX"); errors++;
			}
			if (maxY < minY) {
				BatchError(filetype, line, 2, "maxY", "minY"); errors++;
			}
		}
		if (seedtype == 0 && freetype == 0) {
			if (nCells < 1) {
				BatchError(filetype, line, 11, "NCells"); errors++;
			}
			int range_cells;
			range_cells = (maxX - minX) * (maxY - minY);
			if (nCells > range_cells) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << "NCells may not be greater than the area specified (i.e. "
					<< range_cells << " cells)" << endl;
			}
		}
		if (seedtype == 1 && sptype == 1 && nSpCells < 1) {
			BatchError(filetype, line, 11, "NSpCells"); errors++;
		}

		ifsInitFile >> initFreezeYear >> restrictRows >> restrictFreq >> finalFreezeYear;
		if (seedtype == 0) {
			if (initFreezeYear < 0) {
				BatchError(filetype, line, 19, "InitFreezeYear"); errors++;
			}
			if (restrictRows < 0) {
				BatchError(filetype, line, 19, "RestrictRows"); errors++;
			}
			if (restrictRows > 0 && restrictFreq <= 0) {
				BatchError(filetype, line, 10, "RestrictFreq"); errors++;
			}
			if (finalFreezeYear < 0) {
				BatchError(filetype, line, 19, "FinalFreezeYear"); errors++;
			}
			else {
				if (finalFreezeYear > 0 && finalFreezeYear <= initFreezeYear) {
					BatchError(filetype, line, 1, "FinalFreezeYear", "InitFreezeYear"); errors++;
				}
			}
		}

		ifsInitFile >> filename;
		if (filename == "NULL") {
			if (seedtype == 2) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << ftype2 << " is compulsory for SeedType 2" << endl;
			}
		}
		else {
			if (seedtype == 2) {
				checkfile = true;
				for (i = 0; i < (int)indsfiles.size(); i++) {
					if (filename == indsfiles[i]) { // file has already been checked
						checkfile = false;
					}
				}
				if (checkfile) {
					fname = indir + filename;
					batchLogOfs << "Checking " << ftype2 << " " << fname << endl;
					ifsInitIndsFile.open(fname.c_str());
					if (ifsInitIndsFile.is_open()) {
						err = CheckInitIndsFile(simNb, spNb);
						if (err == 0) FileHeadersOK(ftype2); 
						else errors++;
						ifsInitIndsFile.close();
					}
					else {
						OpenError(ftype2, fname); errors++;
					}
					if (ifsInitIndsFile.is_open()) 
						ifsInitIndsFile.close();
					ifsInitIndsFile.clear();
					indsfiles.push_back(filename);
				}
			}
			else {
				BatchError(filetype, line, 0, " "); 
				errors++;
				batchLogOfs << ftype2 << " must be NULL for SeedType "
					<< seedtype << endl;
			}
		}

		if (gUsesStageStruct && maxNbStages > 1) {
			ifsInitFile >> initAge;
			if (seedtype != 2 && (initAge < 0 || initAge > 2)) {
				BatchError(filetype, line, 2, "initAge"); 
				errors++;
			}
			float propstage;
			float cumprop = 0.0;
			for (i = 1; i < maxNbStages; i++) {
				ifsInitFile >> propstage;
				cumprop += propstage;
				if (seedtype != 2 && (propstage < 0.0 || propstage > 1.0)) {
					colheader = "PropStage" + to_string(i);
					BatchError(filetype, line, 20, colheader); 
					errors++;
				}
			}
			if (seedtype != 2 && (cumprop < 0.99999 || cumprop > 1.00001)) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << "Initial proportions must sum to 1.0" << endl;
			}
		}

		// read next simulation
		line++;
		simNb = -98765;
		ifsInitFile >> simNb;
		if (ifsInitFile.eof()) simNb = -98765;

	} // end of while loop
	// check for correct number of lines for previous simulation
	if (current.simLines != current.reqdSimLines) {
		BatchError(filetype, line, 0, " "); errors++;
		batchLogOfs << gNbLinesStr << current.simNb
			<< gShouldBeStr << current.reqdSimLines << endl;
	}
	if (!ifsInitFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int CheckInitIndsFile(int simNb, species_id sp) {
	string header;
	int year, species, patchID, x, y, ninds, sex, age, stage, prevyear;

	int nbErrors = 0;
	string filetype = "InitIndsFile";
	spInputOptions inputOpt = gSpInputOpt.at(simNb).at(sp);

	// Parse header line
	ifsInitIndsFile >> header; 
	if (header != "Year") nbErrors++;
	ifsInitIndsFile >> header; 
	if (gUsesPatches) {
		ifsInitIndsFile >> header; 
		if (header != "PatchID") nbErrors++;
	}
	else {
		ifsInitIndsFile >> header; 
		if (header != "X") nbErrors++;
		ifsInitIndsFile >> header;
		if (header != "Y") nbErrors++;
	}
	ifsInitIndsFile >> header; 
	if (header != "Ninds") nbErrors++;
	ifsInitIndsFile >> header;
	if (header != "Sex") nbErrors++;
	if (gUsesStageStruct) {
		ifsInitIndsFile >> header; 
		if (header != "Age") nbErrors++;
		ifsInitIndsFile >> header; 
		if (header != "Stage") nbErrors++;
	}

	// Report any errors in headers, and if so, terminate validation
	if (nbErrors > 0) {
		FormatError(filetype, nbErrors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	string filename, ftype2, fname;
	year = prevyear = -98765;
	ifsInitIndsFile >> year;
	while (year != -98765) {
		if (year < 0) {
			BatchError(filetype, line, 19, "Year"); 
			nbErrors++;
		}
		else {
			if (year < prevyear) {
				BatchError(filetype, line, 2, "Year", "previous Year"); 
				nbErrors++;
			}
		}
		prevyear = year;
		ifsInitIndsFile >> species;
		if (species != 0) {
			BatchError(filetype, line, 0, " "); 
			nbErrors++;
			batchLogOfs << "Species must be 0" << endl;
		}
		if (gUsesPatches) {
			ifsInitIndsFile >> patchID;
			if (patchID < 1) {
				BatchError(filetype, line, 11, "PatchID"); 
				nbErrors++;
			}
		}
		else {
			ifsInitIndsFile >> x >> y;
			if (x < 0 || y < 0) {
				BatchError(filetype, line, 19, "X and Y"); 
				nbErrors++;
			}
		}
		ifsInitIndsFile >> ninds;
		if (ninds < 1) {
			BatchError(filetype, line, 11, "Ninds"); 
			nbErrors++;
		}
		if (inputOpt.reproType > 0) {
			ifsInitIndsFile >> sex;
			if (sex < 0 || sex > 1) {
				BatchError(filetype, line, 1, "Sex"); 
				nbErrors++;
			}
		}
		if (gUsesStageStruct) {
			ifsInitIndsFile >> age >> stage;
			if (age < 1) {
				BatchError(filetype, line, 11, "Age"); 
				nbErrors++;
			}
			if (stage < 1) {
				BatchError(filetype, line, 11, "Stage");
				nbErrors++;
			}
			if (stage >= inputOpt.nbStages) {
				BatchError(filetype, line, 4, "Stage", "no. of stages"); 
				nbErrors++;
			}
		}
		line++;
		year = -98765;
		ifsInitIndsFile >> year;
		if (ifsInitIndsFile.eof()) 
			year = -98765;

	} // end of while loop
	if (!ifsInitIndsFile.eof()) {
		EOFerror(filetype);
		nbErrors++;
	}
	return nbErrors;
}

//---------------------------------------------------------------------------
/*
Check stage- and sex-dependency fields for any of the dispersal files.
Check that the number of records for a simulation matches the stage-
and sex-dependency settings (unless checklines is false).
Validate the IIV field (if present).
*/
simCheck CheckStageSex(string whichInputFile, int whichLine, int simNb, species_id sp, 
	simCheck prev, int isStageDep, int isSexDep, int stage, int sex, 
	int isIndVar, bool mustCheckLines, bool mustCheckStgDepWithIndVar)
{
	simCheck current;
	current.errors = 0;
	int expectedStage;

	int nbSexDisp = gSpInputOpt.at(simNb).at(sp).reproType == 0 ? 1 : 2;
	int nbStg = gSpInputOpt.at(simNb).at(sp).nbStages;

	// has there been a change of simulation number?;
	if (simNb == prev.simNb) { // no
		current.isNewSim = false; 
		current.simLines = prev.simLines + 1;
	}
	else { // yes
		// check for valid simulation number
		current.isNewSim = true; 
		current.simLines = 1;
		if (whichLine > 1 && simNb != prev.simNb + 1) {
			BatchError(whichInputFile, whichLine, 222, " "); 
			current.errors++;
		}
		// check for correct number of lines for previous simulation
		if (mustCheckLines && !(prev.simLines >= prev.reqdSimLines)) {
			BatchError(whichInputFile, whichLine, 0, " "); 
			current.errors++;
			batchLogOfs << "No. of lines for previous Simulation " << prev.simNb
				<< gShouldBeStr << prev.reqdSimLines << endl;
		}
	}
	current.simNb = simNb;

	// validate inStageDep
	if (gUsesStageStruct) {
		if (isStageDep != 0 && isStageDep != 1) {
			BatchError(whichInputFile, whichLine, 1, "StageDep"); 
			current.errors++;
			isStageDep = 1; // to calculate required number of lines
		}
	}
	else if (isStageDep != 0) {
		BatchError(whichInputFile, whichLine, 0, " ");
		current.errors++;
		batchLogOfs << "StageDep must be 0 for non-stage-structured model" << endl;
		isStageDep = 0; // to calculate required number of lines
	}
	// validate inSexDep
	if (nbSexDisp == 2) {
		if (isSexDep != 0 && isSexDep != 1) {
			BatchError(whichInputFile, whichLine, 1, "SexDep"); 
			current.errors++;
			isSexDep = 1; // to calculate required number of lines
		}
	}
	else if (isSexDep != 0) {
		BatchError(whichInputFile, whichLine, 0, " ");
		current.errors++;
		batchLogOfs << "SexDep must be 0 for asexual model" << endl;
		isSexDep = 0; // to calculate required number of lines
	}
	if (current.isNewSim) { // set required number of lines
		if (isStageDep) {
			current.reqdSimLines = nbStg;
			if (isSexDep) current.reqdSimLines *= nbSexDisp;
		}
		else {
			current.reqdSimLines = isSexDep ? nbSexDisp : 1;
		}
	}
	else current.reqdSimLines = prev.reqdSimLines;

	// validate stage
	if (isStageDep) { // there must be 1 or 2 lines for each stage
		if (isSexDep) { // there must be 2 lines for each stage
			expectedStage = current.simLines % 2 ? 
				(current.simLines + 1) / 2: 
				current.simLines / 2;
			if (stage != expectedStage - 1) {
				BatchError(whichInputFile, whichLine, 0, " "); 
				current.errors++;
				batchLogOfs << "Stages must be sequentially numbered from 0" << endl;
			}
		}
		else if (stage != current.simLines - 1) {
			BatchError(whichInputFile, whichLine, 0, " ");
			current.errors++;
			batchLogOfs << "Stages must be sequentially numbered from 0" << endl;
		}
	}
	else if (stage != 0) {
		BatchError(whichInputFile, whichLine, 0, " ");
		current.errors++;
		batchLogOfs << "Stage must be 0 for non-stage-structured model" << endl;
	}
	// validate sex
	if (isSexDep) {
		if (sex != (current.simLines + 1) % 2) {
			BatchError(whichInputFile, whichLine, 0, " "); 
			current.errors++;
			batchLogOfs << "Sex must be alternately 0 and 1 if SexDep is 1" << endl;
		}
	}
	else if (sex != 0) {
		BatchError(whichInputFile, whichLine, 0, " ");
		current.errors++;
		batchLogOfs << "Sex must be 0 if SexDep is 0" << endl;
	}

	// validate inIndVar
	if (isStageDep && !mustCheckStgDepWithIndVar) {
		if (isIndVar != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); 
			current.errors++;
			batchLogOfs << "IndVar must be 0 if stage-dependent" << endl;
		}
	}
	else if (isIndVar < 0 || isIndVar > 1) {
		BatchError(whichInputFile, whichLine, 1, "IndVar");
		current.errors++;
	}
	return current;
}

// Functions to handle and report error conditions
void BatchError(string filename, int line, int option, string fieldname)
{
	if (line == -999) { // message does not cite line number
		batchLogOfs << "*** Error in " << filename << ": ";
	}
	else {
		batchLogOfs << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		batchLogOfs << fieldname << " must be 0 or 1";
		break;
	case 2:
		batchLogOfs << fieldname << " must be 0, 1 or 2";
		break;
	case 3:
		batchLogOfs << fieldname << " must be 0, 1, 2 or 3";
		break;
	case 4:
		batchLogOfs << fieldname << " must be from 0 to 4";
		break;
	case 5:
		batchLogOfs << fieldname << " must be from 0 to 5";
		break;
	case 6:
		batchLogOfs << fieldname << " must be from 0 to 6";
		break;
	case 7:
		batchLogOfs << fieldname << " must be from 0 to 7";
		break;
	case 10:
		batchLogOfs << fieldname << " must be greater than zero";
		break;
	case 11:
		batchLogOfs << fieldname << " must be 1 or more";
		break;
	case 12:
		batchLogOfs << fieldname << " must be 2 or more";
		break;
	case 13:
		batchLogOfs << fieldname << " must be 3 or more";
		break;
	case 18:
		batchLogOfs << fieldname << " must be greater than 1.0";
		break;
	case 19:
		batchLogOfs << fieldname << " must be 0 or more";
		break;
	case 20:
		batchLogOfs << fieldname << " must be between 0 and 1";
		break;
	case 21:
		batchLogOfs << fieldname << " must be greater than 1";
		break;
	case 33:
		batchLogOfs << fieldname << " must be 1, 2 or 3";
		break;
	case 44:
		batchLogOfs << fieldname << " must be from 1 to 4";
		break;
	case 55:
		batchLogOfs << fieldname << " must be from 1 to 5";
		break;
	case 66:
		batchLogOfs << fieldname << " must be from 1 to 6";
		break;
	case 100:
		batchLogOfs << fieldname << " must be between 0 and 100";
		break;
	case 111:
		batchLogOfs << fieldname << " must match the first Simulation in ParameterFile";
		break;
	case 222:
		batchLogOfs << "Simulation numbers must be sequential integers";
		break;
	case 333:
		batchLogOfs << "No. of " << fieldname << " columns must equal max. no. of habitats ("
			<< gMaxNbHab << ") and be sequentially numbered starting from 1";
		break;
	case 444:
		batchLogOfs << "No. of " << fieldname << " columns must be one fewer than no. of stages, and be sequentially numbered starting from 1";
		break;
	case 555:
		batchLogOfs << "No. of " << fieldname << " columns must equal no. of stages, and be sequentially numbered starting from 0";
		break;
	case 666:
		batchLogOfs << fieldname << " must be a unique positive integer";
		break;
	default:
		batchLogOfs << "*** Unspecified error regarding parameter " << fieldname;
	}
	if (option != 0) batchLogOfs << endl;
}

void BatchError(string filename, int line, int option, string fieldname, string fieldname2)
{
	if (line == -999) { // message does not cite line number
		batchLogOfs << "*** Error in " << filename << ": ";
	}
	else {
		batchLogOfs << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		batchLogOfs << fieldname << " must be greater than " << fieldname2;
		break;
	case 2:
		batchLogOfs << fieldname << " must be greater than or equal to " << fieldname2;
		break;
	case 3:
		batchLogOfs << fieldname << " must be less than or equal to " << fieldname2;
		break;
	case 4:
		batchLogOfs << fieldname << " must be less than " << fieldname2;
		break;
	default:
		batchLogOfs << "*** Unspecified error regarding parameters " << fieldname
			<< " and " << fieldname2;
	}
	if (option != 0) batchLogOfs << endl;
}

void printControlFormatError()
{
	cout << "Format error in Control file" << endl;
	batchLogOfs << endl << "***" << endl << "*** Format error in Control file:"
		<< gCaseSensitiveStr << " and file names" << gSpecMustMatchStr
		<< endl
		<< "***" << endl;
}

void FormatError(string filename, int errors)
{
	batchLogOfs << "*** Format error in header line of ";
	if (errors == 0) {
		batchLogOfs << filename << endl;
	}
	else {
		batchLogOfs << filename << ": " << errors << " error";
		if (errors > 1) batchLogOfs << "s";
		batchLogOfs << " detected" << endl;
	}
}

void OpenError(string ftype, string fname)
{
	batchLogOfs << "*** Unable to open " << ftype << " " << fname << endl;
}

void EOFerror(string filename)
{
	batchLogOfs << "*** Failed to read to EOF in " << filename << endl;
}

void FileOK(string ftype, int n, int option)
{
	batchLogOfs << ftype << " OK: total no. of ";
	switch (option) {
	case 0:
		batchLogOfs << "simulations = ";
		break;
	case 1:
		batchLogOfs << "landscapes = ";
		break;
	case 2:
		batchLogOfs << "parameters = ";
		break;
	default:
		batchLogOfs << "PROBLEMS = ";
	}
	batchLogOfs << n << endl;
}

void FileHeadersOK(string filename)
{
	batchLogOfs << filename << " OK" << endl;
}

void SimulnCountError(string filename)
{
	batchLogOfs << "*** No. of simulations in " << filename
		<< " does not match no. in ParameterFile" << endl;
}

//---------------------------------------------------------------------------
int ReadLandFile(Landscape* pLandscape)
{
	landParams ppLand = pLandscape->getLandParams();
	genLandParams ppGenLand = pLandscape->getGenLandParams();

	if (gLandType == 9) { // artificial landscape
		ppLand.rasterType = 9;
		ifsLandFile >> ppLand.landNum >> ppGenLand.isFractal >> ppGenLand.isContinuous
			>> ppLand.dimX >> ppLand.dimY >> ppGenLand.minPct >> ppGenLand.maxPct
			>> ppGenLand.propSuit >> ppGenLand.hurst;
		ppLand.maxX = ppLand.dimX - 1; 
		ppLand.maxY = ppLand.dimY - 1;

		if (ppGenLand.isFractal && ppLand.maxX > ppLand.maxY) {
			return -901;
		}
		if (ppGenLand.isFractal) {
			if ((ppLand.dimX < 3 || ppLand.dimX % 2 != 1)
				|| (ppLand.dimY < 3 || ppLand.dimY % 2 != 1)) {
				return -902;
			}
		}
		// SCFP 26/9/13 - min and max habitat percentages need to be set for all types of
		// fractal landscape (including discrete), as they are passed to the fractal generator
		// NOTE that will not have been checked for a discrete landscape
		if (ppGenLand.isFractal && !ppGenLand.isContinuous) { 
			ppGenLand.minPct = 1; 
			ppGenLand.maxPct = 100;
		}
		if (ppGenLand.isContinuous)
			ppLand.nHab = 2;
		else 
			ppLand.nHab = 1;
	}
	else { // imported raster map
		string inNbHab;
		ifsLandFile >> ppLand.landNum >> inNbHab >> gHabMapName >> gSpLandName;
		ifsLandFile >> gDynLandFileName;
		if (gLandType == 2) 
			ppLand.nHab = 1; // habitat quality landscape has one habitat class
	}

	pLandscape->setLandParams(ppLand, true);
	pLandscape->setGenLandParams(ppGenLand);

	return ppLand.landNum;
}

void ReadSpLandFile(ifstream& ifsSpLand,
	map<species_id, string>& pathsToPatchMaps,
	map<species_id, string>& pathsToCostMaps,
	map<species_id, string>& pathsToSpDistMaps,
	map<species_id, bool>& whichUseSpDist
) {
	int inSp;
	string patchMap, usesCosts, SpDistMap;
	int nbSpecies = whichUseSpDist.size();

	for (int i = 0; i < nbSpecies; i++) {

		ifsSpLand >> inSp >> patchMap >> usesCosts >> SpDistMap;

		patchMap = patchMap == "NULL" ? " " :
			paramsSim->getDir(1) + patchMap;
		pathsToPatchMaps.emplace(inSp, patchMap);

		if (!(usesCosts == "NULL" || usesCosts == "none")) {
			// only populate with species for which costs apply
			usesCosts = paramsSim->getDir(1) + usesCosts;
			pathsToCostMaps.emplace(inSp, usesCosts);
		}

		if (whichUseSpDist.at(inSp))
			pathsToSpDistMaps.emplace(inSp, SpDistMap);
	}
}


//---------------------------------------------------------------------------
int ReadDynLandFile(Landscape* pLandscape) {

	string landChgMap, spLandFile;
	int change, imported;
	int nbChanges = 0;
	bool usesCosts = false;
	landChange chg;
	landParams ppLand = pLandscape->getLandParams();
	string pathToDynLandFile = paramsSim->getDir(1) + gDynLandFileName;

	ifsDynLandFile.open(pathToDynLandFile.c_str());
	if (ifsDynLandFile.is_open()) {
		string header;
		int nheaders = 4;
		for (int i = 0; i < nheaders; i++) 
			ifsDynLandFile >> header;
	}
	else {
		ifsDynLandFile.clear();
		return 72727;
	}

	// read data lines
	change = -98765;
	ifsDynLandFile >> change; // first change number

	while (change != -98765) {
		chg.chgnum = change;
		ifsDynLandFile >> chg.chgyear >> landChgMap >> spLandFile;
		chg.habfile = paramsSim->getDir(1) + landChgMap;
		chg.spLandFile = paramsSim->getDir(1) + spLandFile;
		
		nbChanges++;
		pLandscape->addLandChange(chg);

		// read first field on next line
		change = -98765;
		ifsDynLandFile >> change;
		if (ifsDynLandFile.eof()) change = -98765;
	}

	ifsDynLandFile.close();
	ifsDynLandFile.clear();

	// read landscape change maps
	if (ppLand.usesPatches) {
		pLandscape->createPatchChgMatrix();
	}
	if (usesCosts) {
		pLandscape->createCostsChgMatrix();
	}
	for (int chgIndex = 0; chgIndex < nbChanges; chgIndex++) {
		
		imported = pLandscape->readLandChange(chgIndex, usesCosts);
		if (imported != 0) return imported;

		if (ppLand.usesPatches) {
			pLandscape->recordPatchChanges(chgIndex + 1);
		}
		if (usesCosts) {
			pLandscape->recordCostChanges(chgIndex + 1);
		}
	}
	if (ppLand.usesPatches) {
		// record changes back to original landscape for multiple replicates
		pLandscape->resetPatchChanges();
	}
	if (usesCosts) {
		pLandscape->resetCostChanges();
	}
	return 0;
}

//--------------------------------------------------------------------------

void flushHeaders(ifstream& ifs) {
	string headerLine;
	// Pass the first line (headers) to an empty string...
	std::getline(ifs, headerLine);
	// ... and do nothing with it 
}

int ReadGeneticsFile(speciesMap_t& simSpecies, ifstream& ifs) {

	string indir = paramsSim->getDir(1);
	set<int> patchList;

	if (ifs.is_open()) {
		string line, value;

		// Read 1 line at every call
		std::getline(ifs, line);

		// Convert input parameters to string vector
		stringstream ss(line);
		vector<string> parameters;
		while (std::getline(ss, value, '	'))
			parameters.push_back(value);

		// Assumes all input is correct after errors being handled by CheckGenetics
		species_id sp = stoi(parameters[1]);
		Species* pSpecies = simSpecies.at(sp);
		// not ideal to reset these in here 
		pSpecies->resetGeneticParameters();

		int genomeSize = stoi(parameters[2]);
		set<int> chrEnds = stringToChromosomeEnds(parameters[3], genomeSize);
		float recombinationRate = parameters[3] == "#" ? 0.0 : stof(parameters[4]);

		outputParams out = pSpecies->getOutputParams();
		out.outputGenes = (parameters[5] == "TRUE");
		out.outputWeirCockerham = (parameters[6] == "TRUE");
		out.outputWeirHill = (parameters[7] == "TRUE");
		out.outputStartGenetics = stoi(parameters[8]);
		out.outputGeneticInterval = stoi(parameters[9]);
		pSpecies->setOutputParams(out);

		string inPatches = parameters[10];
		string patchSamplingOption;
		int nPatchesToSample = 0;
		if (inPatches != "all" && inPatches != "random" && inPatches != "random_occupied") {
			// then must be a list of indices
			patchSamplingOption = "list";
			patchList = stringToPatches(inPatches);
			if (patchList.contains(0)) throw logic_error("Patch sampling: ID 0 is reserved for the matrix and should not be sampled.");
		}
		else {
			patchSamplingOption = inPatches;
			if (inPatches == "random" || inPatches == "random_occupied")
				nPatchesToSample = stoi(parameters[11]);
			// patchList remains empty, filled when patches are sampled every gen
		}
		const string strNbInds = parameters[12];
		const int nbStages = pSpecies->getStageParams().nStages;
		set<int> stagesToSampleFrom = stringToStages(parameters[13], nbStages);

		pSpecies->setGeneticParameters(chrEnds, genomeSize, recombinationRate, patchSamplingOption,
			patchList, strNbInds, stagesToSampleFrom, nPatchesToSample);
	}
	else {
		throw runtime_error("GeneticsFile is not open.");
	}
	return 0;
}

int ReadTraitsFile(speciesMap_t& simSpecies, ifstream& ifs, map<species_id, spInputOptions> simOptionsMap) {

	Species* pSpecies;
	int prevsimNb = -998;

	if (ifs.is_open()) {

		//read first header line
		string strLine, entry;

		int nbRowsToRead = 1; // need to read species to get correct number
		for (int i = 0; i < nbRowsToRead; i++) {
			
			// Read input row
			std::getline(ifs, strLine);
			
			// Read input parameters as strings
			stringstream inLine(strLine);
			vector<string> parameters;
			while (std::getline(inLine, entry, '	')) {
				parameters.push_back(entry);
			}
			if (i == 0) {
				species_id sp = stoi(parameters[1]);
				nbRowsToRead = simOptionsMap.at(sp).nbTraitFileRows;
				pSpecies = simSpecies.at(sp);
				pSpecies->clearTraitTable();
			}

			// Create trait from parameters 
			setUpSpeciesTrait(pSpecies, parameters);
		}
	}
	else {
		throw runtime_error("TraitsFile is not open.");
	}
	return 0;
}

// Set up a trait from input parameters and add it Species
void setUpSpeciesTrait(Species* pSpecies, vector<string> parameters) {
	// Assumes all input is correct, errors have been handled by CheckTraits

	const int genomeSize = pSpecies->getGenomeSize();
	TraitType traitType = stringToTraitType(parameters[2]);
	const sex_t sex = stringToSex(parameters[3]);
	if (sex != NA) traitType = addSexDepToTrait(traitType, sex);
	const set<int> positions = stringToLoci(parameters[4], parameters[5], genomeSize);
	const ExpressionType expressionType = stringToExpressionType(parameters[6]);

	// Initial allele distribution parameters
	const DistributionType initDist = stringToDistributionType(parameters[7]);
	const map<GenParamType, float> initParams = stringToParameterMap(parameters[8]);

	// Initial dominance distribution parameters
	const DistributionType initDomDist = stringToDistributionType(parameters[9]);
	const map<GenParamType, float> initDomParams = stringToParameterMap(parameters[10]);

	// Mutation parameters
	bool isInherited = (parameters[11] == "TRUE");
	DistributionType mutationDistribution = isInherited ? 
		stringToDistributionType(parameters[12]) : 
		DistributionType::NONE;
	map<GenParamType, float> mutationParameters;
	if (isInherited) {
		mutationParameters = stringToParameterMap(parameters[13]);
	}

	// Dominance distribution parameters
	const DistributionType dominanceDist = stringToDistributionType(parameters[14]);
	const map<GenParamType, float> dominanceParams = stringToParameterMap(parameters[15]);

	float mutationRate = isInherited ? stof(parameters[16]) : 0.0;
	
	parameters[17].erase(
		// send windows line endings to hell where they belong
		remove(parameters[17].begin(), parameters[17].end(), '\r'),
		parameters[17].end()
	);
	const bool isOutput = parameters[18] == "TRUE";
	int ploidy = pSpecies->getDemogrParams().repType == 0 ? 1 : 2;

	// Create species trait
	unique_ptr<SpeciesTrait> trait(new SpeciesTrait(
		traitType, sex, 
		positions, expressionType, 
		initDist, initParams, 
		initDomDist, initDomParams,
		isInherited, mutationRate, 
		mutationDistribution, mutationParameters,
		dominanceDist, dominanceParams,
		ploidy,
		isOutput
	));
	pSpecies->addTrait(traitType, *trait);
}

// Convert string to corresponding TraitType value, if valid
TraitType stringToTraitType(const std::string& str) {
	// Non-dispersal traits
	if (str == "neutral") return NEUTRAL;
	else if (str == "genetic_load") return GENETIC_LOAD;
	// Sex-invariant dispersal traits
	else if (str == "emigration_d0") return E_D0; // EP uses d0 for trait data
	else if (str == "emigration_alpha") return E_ALPHA;
	else if (str == "emigration_beta") return E_BETA;
	else if (str == "settlement_s0") return S_S0;
	else if (str == "settlement_alpha") return S_ALPHA;
	else if (str == "settlement_beta") return S_BETA;
	else if (str == "kernel_meanDistance1") return KERNEL_MEANDIST_1;
	else if (str == "kernel_meanDistance2") return KERNEL_MEANDIST_2;
	else if (str == "kernel_probability") return KERNEL_PROBABILITY;
	else if (str == "crw_stepLength") return CRW_STEPLENGTH;
	else if (str == "crw_stepCorrelation") return CRW_STEPCORRELATION;
	else if (str == "sms_directionalPersistence") return SMS_DP;
	else if (str == "sms_goalBias") return SMS_GB;
	else if (str == "sms_alphaDB") return SMS_ALPHADB;
	else if (str == "sms_betaDB") return SMS_BETADB;
	else return INVALID_TRAIT;
}

// Convert string to corresponding ExpressionType value, if valid
ExpressionType stringToExpressionType(const std::string& str) {
	if (str == "average") return AVERAGE;
	else if (str == "additive") return ADDITIVE;
	else if (str == "multiplicative") return MULTIPLICATIVE;
	else if (str == "#") return NOTEXPR;
	else throw logic_error(str + " is not a valid gene expression type.");
}

// Convert string to corresponding DistributionType value, if valid
DistributionType stringToDistributionType(const std::string& str) {
	if (str == "#") return NONE;
	else if (str == "uniform") return UNIFORM;
	else if (str == "normal") return NORMAL;
	else if (str == "gamma") return GAMMA;
	else if (str == "scaled") return SCALED;
	else if (str == "negExp") return NEGEXP;
	else if (str == "KAM") return KAM;
	else if (str == "SSM") return SSM;
	else throw logic_error(str + " is not a valid distribution type.");
}

// Convert distribution parameters field into appropriate type
map<GenParamType, float> stringToParameterMap(string parameterString) {

	map<GenParamType, float> paramMap;
	if (parameterString != "#") {
		// drop quotation marks
		parameterString.erase(remove(parameterString.begin(), parameterString.end(), '\"'), parameterString.end());
		stringstream ss(parameterString);

		string singleParamString, valueWithin;
		while (std::getline(ss, singleParamString, ',')) {
			stringstream sss(singleParamString);
			vector<string> paramNameAndVal;
			while (std::getline(sss, valueWithin, '=')) {
				paramNameAndVal.push_back(valueWithin);
			}

			if (paramNameAndVal.size() == 2) {
				GenParamType parameterT = strToGenParamType(paramNameAndVal[0]);
				if (parameterT == INVALID)
					throw logic_error("Invalid genetic parameter name.");
				float value = stof(paramNameAndVal[1]);
				paramMap.emplace(parameterT, value);
			}
			else throw logic_error("Traits file: ERROR - parameter values for a distribution missing, should be e.g. 'mean=0,standard_deviation=0.5' or if not applicable put #");
		}
	}
	return paramMap;
}

// Convert string to corresponding SexType value, if valid
const sex_t stringToSex(const std::string& str) {
	if (str == "female") return FEM;
	else if (str == "male") return MAL;
	else if (str == "#") return NA;
	else return INVALID_SEX;
}

// Convert patches input parameter string into set of patch indices
set<int> stringToPatches(const string& str) {

	set<int> patches;
	stringstream ss(str);
	string strPch;
	int pch;
	// Read comma-separated values
	while (std::getline(ss, strPch, ',')) {
		strPch.erase(remove(strPch.begin(), strPch.end(), '\"'), strPch.end());
		pch = std::stoi(strPch);
		patches.insert(pch);
	}
	return patches;
}

// Convert stages input parameter string into set of stage numbers
set<int> stringToStages(const string& str, const int& nbStages) {
	set<int> stages;
	if (str == "all") {
		for (int stg = 0; stg < nbStages; ++stg) {
			stages.insert(stg);
		}
	}
	else {
		// Parse comma-separated list from input string
		stringstream ss(str);
		string strStg;
		int stg;
		// Read comma-separated values
		while (std::getline(ss, strStg, ',')) {
			strStg.erase(remove(strStg.begin(), strStg.end(), '\"'), strStg.end());
			stg = std::stoi(strStg);
			if (stg > nbStages - 1)
				throw logic_error("Genetics file: ERROR - sampled stage exceeds number of stages.");
			else {
				stages.insert(stg);
			}
		}
	}
	return stages;
}

// Convert ChromosomeEnds input parameter string into set of positions
set<int> stringToChromosomeEnds(string str, const int& genomeSize) {
	set<int> chromosomeEnds;
	if (str == "#")
		chromosomeEnds.insert(genomeSize - 1); // last position in genome
	else {
		// Parse comma-separated list from input string
		// drop quotation marks
		str.erase(remove(str.begin(), str.end(), '\"'), str.end());
		stringstream ss(str);

		string strPos;
		int pos;
		// Read comma-separated positions
		while (std::getline(ss, strPos, ',')) {
			pos = std::stoi(strPos);
			chromosomeEnds.insert(pos);
		}
	}
	return chromosomeEnds;
}

set<int> selectRandomLociPositions(int nbLoci, const int& genomeSize) {
	set<int> positions;
	if (nbLoci > genomeSize) throw logic_error("Number of random loci exceeds genome size.");
	int rndLocus;
	for (int i = 0; i < nbLoci; ++i)
	{
		do {
			rndLocus = pRandom->IRandom(0, genomeSize - 1);
		} while (positions.contains(rndLocus));
		positions.insert(rndLocus);
	}
	return positions;
}

set<int> stringToLoci(string pos, string nLoci, const int& genomeSize) {

	set<int> positions;

	if (pos != "random") {

		// Parse comma-separated list from input string
		stringstream ss(pos);
		string value, valueWithin;
		// Read comma-separated positions
		while (std::getline(ss, value, ',')) {
			stringstream sss(value);
			vector<int> positionRange;
			// Read single positions and dash-separated ranges
			while (std::getline(sss, valueWithin, '-')) {
				valueWithin.erase(remove(valueWithin.begin(), valueWithin.end(), '\"'), valueWithin.end());
				positionRange.push_back(stoi(valueWithin));
			}
			switch (positionRange.size())
			{
			case 1: // single position
				if (positionRange[0] >= genomeSize)
					throw logic_error("Traits file: ERROR - trait positions must not exceed genome size");
				positions.insert(positionRange[0]);
				break;
			case 2: // dash-separated range
				if (positionRange[0] >= genomeSize || positionRange[1] >= genomeSize) {
					throw logic_error("Traits file: ERROR - trait positions must not exceed genome size");
				}
				if (positionRange[0] >= positionRange[1])
					throw logic_error("Position ranges must be in ascending order");
				for (int i = positionRange[0]; i < positionRange[1] + 1; ++i) {
					positions.insert(i);
				}
				break;
			default: // zero or more than 2 values between commas: error
				throw logic_error("Traits file: ERROR - incorrectly formatted position range.");
				break;
			}
		}

		for (auto position : positions) {
			if (position >= genomeSize)
				throw logic_error("Traits file: ERROR - trait positions " + to_string(position) + " must not exceed genome size.\n");
		}
	}
	else { // random
		positions = selectRandomLociPositions(stoi(nLoci), genomeSize);
	}
	return positions;
}

GenParamType strToGenParamType(const string& str) {
	if (str == "mean")
		return MEAN;
	else if (str == "sd")
		return SD;
	else if (str == "min")
		return MIN;
	else if (str == "max")
		return MAX;
	else if (str == "shape")
		return SHAPE;
	else if (str == "scale")
		return SCALE;
	else return INVALID;
}

void ReadSimParameters() {

	string inAbsorbing, inFixRepSeed;
	string inEnvStoch, inEnvStochType;

	simParams sim = paramsSim->getSim();

	ifsSimFile >> sim.simulation >> sim.reps >> sim.years;
	ifsSimFile >> inAbsorbing;
	sim.absorbing = (inAbsorbing == "1");

	ifsSimFile >> inFixRepSeed;
	sim.fixReplicateSeed = inFixRepSeed == "1";

	paramsSim->setSim(sim);

	// Environmental Stochasticity
	envStochParams env;
	ifsSimFile >> inEnvStoch;
	env.usesStoch = inEnvStoch == "1" || inEnvStoch == "2";
	env.stochIsLocal = inEnvStoch == "2";
	ifsSimFile >> inEnvStochType;
	env.inK = (inEnvStochType == "1");
	ifsSimFile >> env.ac >> env.std;
	paramsStoch->setStoch(env);
}

//---------------------------------------------------------------------------
int ReadParameters(const Landscape* pLandscape, speciesMap_t& simSpecies)
{
	int errorCode = 0;
	landParams paramsLand = pLandscape->getLandParams();

	if (!ifsParamFile.is_open()) {
		cout << endl << "ReadParameters(): ERROR - ParameterFile is not open" << endl;
		return 4086534;
	}

	int gradType, shiftBegin, shiftStop;
	float k, grad_inc, opt_y, f, optExt, shift_rate;
	string inAbsorbing, inShifting, inLocalExt, inHeatMaps;
	species_id sp;

	ifsParamFile >> sp;
	Species* pSpecies = simSpecies.at(sp);
	demogrParams dem = pSpecies->getDemogrParams();

	// Environmental gradient
	envGradParams paramsGrad;

	ifsParamFile >> paramsGrad.gradType;
	ifsParamFile >> paramsGrad.gradIncr >> paramsGrad.optY 
		>> paramsGrad.factor >> paramsGrad.extProbOpt >> inShifting;
	ifsParamFile >> shift_rate >> shiftBegin >> shiftStop;
	paramsGrad.doesShift = (inShifting == "1" && gradType != 0);
	paramsGrad.shiftRate = paramsGrad.doesShift ? shift_rate : 0;
	paramsGrad.shiftBegin = paramsGrad.doesShift ? shiftBegin : 0;
	paramsGrad.shiftStop = paramsGrad.doesShift ? shiftStop : 0;
	paramsGrad.optY0 = paramsGrad.optY; // reset between replicates
	pSpecies->setEnvGrad(paramsGrad);

	float minR, maxR, minK, maxK;
	ifsParamFile >> minR >> maxR >> minK >> maxK;
	if (paramsStoch->getStoch().inK) {
		float minKK, maxKK;
		minKK = minK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
		maxKK = maxK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
		pSpecies->setMinMax(minKK, maxKK);
	}
	else pSpecies->setMinMax(minR, maxR);

	// Local extinction
	float locExtProb;
	ifsParamFile >> locExtProb;
	pSpecies->setLocalExtProb(locExtProb);

	// Demographic parameters
	int nbStg;
	ifsParamFile >> nbStg >> dem.repType >> dem.repSeasons;
	ifsParamFile >> dem.propMales >> dem.harem >> dem.bc >> dem.lambda;
	pSpecies->setDemogr(dem);
	if (gUsesStageStruct) pSpecies->setNbStages(nbStg);

	// Artificial landscape
	if (gLandType == 9) {
		// only one value of K is read, but it must be applied as the second habitat if the
		// landscape is discrete (the first is the matrix where K = 0) or as the first 
		// (only) habitat if the landscape is continuous
		genLandParams genland = pLandscape->getGenLandParams();
		int nhab = genland.isContinuous ? 1 : 2;

		pSpecies->createHabK(nhab);
		ifsParamFile >> k;
		k *= (((float)paramsLand.resol * (float)paramsLand.resol)) / 10000.0f;

		if (genland.isContinuous) {
			pSpecies->setHabK(0, k);
		}
		else {
			pSpecies->setHabK(0, 0);
			pSpecies->setHabK(1, k);
		}
	}
	else {
		pSpecies->createHabK(paramsLand.nHabMax);
		for (int i = 0; i < paramsLand.nHabMax; i++) {
			ifsParamFile >> k;
			k *= ((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f;
			pSpecies->setHabK(i, k);
		}
	}

	// Output parameters
	outputParams spParams;
	ifsParamFile >> spParams.outStartPop	>> spParams.outStartInd
				 >> spParams.outStartTraitCell >> spParams.outStartTraitRow 
				 >> spParams.outStartConn >> spParams.outIntRange 
				 >> spParams.outIntOcc >> spParams.outIntPop 
				 >> spParams.outIntInd >> spParams.outIntTraitCell 
				 >> spParams.outIntTraitRow >> spParams.outIntConn;

	spParams.outRange = spParams.outIntRange > 0;
	spParams.outOccup = spParams.outIntOcc > 0;
	spParams.outPop = spParams.outIntPop > 0;
	spParams.outInds = spParams.outIntInd > 0;
	spParams.outTraitsCells = spParams.outIntTraitCell > 0;
	spParams.outTraitsRows = spParams.outIntTraitRow > 0;
	spParams.outConnect = spParams.outIntConn > 0;

	if (paramsLand.usesPatches) {
		if (spParams.outTraitsRows) errorCode = 104;
	}
	else {
		if (spParams.outConnect) errorCode = 105;
	}
	ifsParamFile >> inHeatMaps;
	spParams.saveVisits = inHeatMaps == "1";
	pSpecies->setOutputParams(spParams);

	return errorCode;
}

//---------------------------------------------------------------------------
int ReadStageStructure(speciesMap_t& simSpecies)
{
	int simulation, postDestructn;
	string inputDir = paramsSim->getDir(1);
	species_id sp;

	ifsStageStructFile >> sp;
	Species* pSpecies = simSpecies.at(sp);
	stageParams sstruct = pSpecies->getStageParams();

	ifsStageStructFile >> simulation;
	ifsStageStructFile >> postDestructn >> sstruct.probRep >> sstruct.repInterval >> sstruct.maxAge;
	if (postDestructn == 1) sstruct.disperseOnLoss = true;
	else sstruct.disperseOnLoss = false;

	string StgStructFile;
	ifsStageStructFile >> StgStructFile;
	ifsTransMatrix.open((inputDir + StgStructFile).c_str());
	int nbSexesDem = pSpecies->getDemogrParams().repType == 2 ? 2 : 1;
	ReadTransitionMatrix(pSpecies, sstruct.nStages, nbSexesDem, 0, 0);
	ifsTransMatrix.close(); 
	ifsTransMatrix.clear();
	ifsStageStructFile >> sstruct.survival;

	float devCoeff, survCoeff;
	string fecStgWtsFile;
	ifsStageStructFile >> sstruct.fecDens >> sstruct.fecStageDens >> fecStgWtsFile;
	if (fecStgWtsFile != "NULL") {
		ifsFecDens.open((inputDir + fecStgWtsFile).c_str());
		ReadStageWeights(pSpecies, 1);
		ifsFecDens.close(); 
		ifsFecDens.clear();
	}
	string devStgWtsFile;
	ifsStageStructFile >> sstruct.devDens >> devCoeff >> sstruct.devStageDens >> devStgWtsFile;
	if (devStgWtsFile != "NULL") {
		ifsDevDens.open((inputDir + devStgWtsFile).c_str());
		ReadStageWeights(pSpecies, 2);
		ifsDevDens.close(); 
		ifsDevDens.clear();
	}
	string survStgWtsFile;
	ifsStageStructFile >> sstruct.survDens >> survCoeff >> sstruct.survStageDens >> survStgWtsFile;
	if (survStgWtsFile != "NULL") {
		ifsSurvDens.open((inputDir + survStgWtsFile).c_str());
		ReadStageWeights(pSpecies, 3);
		ifsSurvDens.close(); 
		ifsSurvDens.clear();
	}

	pSpecies->setStage(sstruct);

	if (sstruct.devDens || sstruct.survDens) {
		pSpecies->setDensDep(devCoeff, survCoeff);
	}

	return 0;
}

//---------------------------------------------------------------------------
int ReadTransitionMatrix(Species* pSpecies, short nstages, short nsexesDem, short hab, short season)
{
	int ii;
	int minAge;
	float ss, dd;
	string header;

	// read header line
	for (int i = 0; i < (nstages * nsexesDem) + 2; i++) {
		ifsTransMatrix >> header;
	}

	if (gMatrix != nullptr) {
		for (int j = 0; j < gMatrixSize; j++) delete[] gMatrix[j];
		delete[] gMatrix;
		gMatrix = nullptr; 
		gMatrixSize = 0;
	}

	if (nsexesDem != 2) { // asexual or implicit sexual model
	// create a temporary matrix
		gMatrix = new float* [nstages];
		gMatrixSize = nstages;
		for (int i = 0; i < nstages; i++)
			gMatrix[i] = new float[nstages];

		for (int i = 0; i < nstages; i++) { 
			ifsTransMatrix >> header;
			for (int j = 0; j < nstages; j++) {
				ifsTransMatrix >> gMatrix[j][i];
			}
			ifsTransMatrix >> minAge; 
			pSpecies->setMinAge(i, 0, minAge);
		}

		for (int j = 1; j < nstages; j++)
			pSpecies->setFec(j, 0, gMatrix[j][0]);

		for (int j = 0; j < nstages; j++) {
			ss = 0.0; dd = 0.0;
			for (int i = 0; i < nstages; i++) {
				if (i == j) ss = gMatrix[j][i];
				if (i == (j + 1)) dd = gMatrix[j][i];
			}
			pSpecies->setSurv(j, 0, ss + dd);
			pSpecies->setDev(j, 0, (ss + dd) > 0.0f ? dd / (ss + dd) : 0.0);
		}
	}
	else { // complex sexual model
		gMatrix = new float* [nstages * 2];
		gMatrixSize = nstages * 2;
		for (int j = 0; j < nstages * 2; j++)
			gMatrix[j] = new float[nstages * 2 - 1];

		for (int i = 0; i < nstages * 2 - 1; i++) {
			ifsTransMatrix >> header;
			for (int j = 0; j < nstages * 2; j++) 
				ifsTransMatrix >> gMatrix[j][i];
			if (i == 0) {
				ifsTransMatrix >> minAge; 
				pSpecies->setMinAge(i, 0, minAge); 
				pSpecies->setMinAge(i, 1, minAge);
			}
			else {
				ifsTransMatrix >> minAge;
				if (i % 2) 
					pSpecies->setMinAge((i + 1) / 2, 1, minAge);	// odd lines  - males
				else
					pSpecies->setMinAge(i / 2, 0, minAge);			// even lines - females
			}
		}

		ii = 1;
		for (int j = 2; j < nstages * 2; j++) {
			if (j % 2 == 0)
				pSpecies->setFec(ii, 1, gMatrix[j][0]);
			else {
				pSpecies->setFec(ii, 0, gMatrix[j][0]);
				ii++;
			}
		}
		// survival and development of male juveniles
		pSpecies->setSurv(0, 1, (gMatrix[0][0] + gMatrix[0][1]));
		if ((gMatrix[0][0] + gMatrix[0][1]) > 0.0)
			pSpecies->setDev(0, 1, (gMatrix[0][1] / (gMatrix[0][0] + gMatrix[0][1])));
		else
			pSpecies->setDev(0, 1, 0.0);
		// survival and development of female juveniles
		pSpecies->setSurv(0, 0, (gMatrix[1][0] + gMatrix[1][2]));
		if ((gMatrix[1][0] + gMatrix[1][2]) > 0.0)
			pSpecies->setDev(0, 0, (gMatrix[1][2] / (gMatrix[1][0] + gMatrix[1][2])));
		else
			pSpecies->setDev(0, 0, 0.0);
		// survival and development of stages 1+
		ii = 1;
		for (int j = 2; j < nstages * 2; j++) {
			ss = 0.0; dd = 0.0;
			if (j % 2 == 0) { // males
				for (int i = 0; i < nstages * 2 - 1; i++)
				{
					if (j == i + 1) ss = gMatrix[j][i];
					if (j == i - 1) dd = gMatrix[j][i];
				}
				pSpecies->setSurv(ii, 1, (ss + dd));
				if ((ss + dd) > 0.0)
					pSpecies->setDev(ii, 1, dd / (ss + dd));
				else
					pSpecies->setDev(ii, 1, 0.0);
			}
			else { // females
				for (int i = 0; i < nstages * 2; i++) {
					if (j == i + 1) ss = gMatrix[j][i];
					if (j == i - 1) dd = gMatrix[j][i];
				}
				pSpecies->setSurv(ii, 0, (ss + dd));
				if ((ss + dd) > 0.0)
					pSpecies->setDev(ii, 0, dd / (ss + dd));
				else
					pSpecies->setDev(ii, 0, 0.0);
				ii++;
			}
		}
	}

	if (gMatrix != nullptr) {
		for (int j = 0; j < gMatrixSize; j++)
			delete[] gMatrix[j];
		delete[] gMatrix;
		gMatrix = nullptr;
		gMatrixSize = 0;
	}

	return 0;
}

//---------------------------------------------------------------------------
int ReadStageWeights(Species* pSpecies, int option)
{
	string header;
	int i, j;
	float f;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	int n = sstruct.nStages;
	if (dem.repType == 2) n *= 2;
	
	switch (option) {

	case 1: { // fecundity
		// create stage weights matrix
		pSpecies->createDDwtFec(n);
		for (i = 0; i < n + 1; i++) 
			ifsFecDens >> header;
		// read coefficients
		for (i = 0; i < n; i++) {
			ifsFecDens >> header;
			for (j = 0; j < n; j++) {
				ifsFecDens >> f; 
				pSpecies->setDDwtFec(j, i, f);
			}
		}
		break;
	}

	case 2: { // development
		//create stage weights matrix
		pSpecies->createDDwtDev(n);
		for (i = 0; i < n + 1; i++) 
			ifsDevDens >> header;
		//read coefficients
		for (i = 0; i < n; i++) {
			ifsDevDens >> header;
			for (j = 0; j < n; j++) {
				ifsDevDens >> f; 
				pSpecies->setDDwtDev(j, i, f);
			}
		}
		break;
	}

	case 3: { // sstruct.survival
		//create stage weights matrix
		pSpecies->createDDwtSurv(n);
		for (i = 0; i < n + 1; i++) 
			ifsSurvDens >> header;
		//read coefficients
		for (i = 0; i < n; i++) {
			ifsSurvDens >> header;
			for (j = 0; j < n; j++) {
				ifsSurvDens >> f;
				pSpecies->setDDwtSurv(j, i, f);
			}
		}
		break;
	}

	}

	return 0;
}

//---------------------------------------------------------------------------
int ReadEmigration(speciesMap_t& simSpecies)
{
	int errorCode = 0;
	int inFullKernel, inDensDep, inStgDep, inSexDep, inIndVar;
	int simulationNb, simNbFirstLine = 0, inStage, inSex, inEmigstage;
	float inEp, inD0, inAlpha, inBeta;
	bool isFirstLine = true;
	int nbSexesDisp;
	emigTraits emigrationTraits;
	species_id sp;
	Species* pSpecies;
	demogrParams dem;
	stageParams sstruct;
	emigRules emig;

	int nbLinesToRead = 1; // need to read first line to set correct value
	for (int line = 0; line < nbLinesToRead; line++) {

		ifsEmigrationFile >> simulationNb >> sp >> inDensDep >> inFullKernel
			>> inStgDep >> inSexDep >> inIndVar >> inEmigstage;

		if (isFirstLine) {
			pSpecies = simSpecies.at(sp);
			dem = pSpecies->getDemogrParams();
			sstruct = pSpecies->getStageParams();
			emig = pSpecies->getEmigRules();

			simNbFirstLine = simulationNb;
			emig.densDep = (inDensDep == 1);
			emig.stgDep = (inStgDep == 1);
			emig.indVar = (inIndVar == 1);
			emig.sexDep = (inSexDep == 1);

			if (inEmigstage >= 0 && inEmigstage < sstruct.nStages)
				emig.emigStage = inEmigstage;
			else emig.emigStage = 0;

			// Set nb lines to correct value
			if (emig.stgDep) nbLinesToRead *= sstruct.nStages;
			if (emig.sexDep) nbLinesToRead *= nbSexesDisp;

			pSpecies->setFullKernel(inFullKernel != 0);
			pSpecies->setEmigRules(emig);
			isFirstLine = false;
		}

		if (simulationNb != simNbFirstLine) { // serious problem
			errorCode = 300;
		}
		ifsEmigrationFile >> inStage >> inSex;

		// ERROR MESSAGES SHOULD NEVER BE ACTIVATED ---------------------------------
		if (dem.repType == 0 && emig.sexDep) {
			errorCode = 301;
		} 
		if (!dem.stageStruct && emig.stgDep) {
			errorCode = 303;
		}
		//---------------------------------------------------------------------------

		ifsEmigrationFile >> inEp >> inD0 >> inAlpha >> inBeta;

		emigrationTraits.d0 = emig.densDep ? inD0 : inEp;
		emigrationTraits.alpha = emig.densDep ? inAlpha : 0.0;
		emigrationTraits.beta = emig.densDep ? inBeta : 0.0;
		pSpecies->setSpEmigTraits(
			emig.stgDep ? inStage : 0,
			emig.sexDep ? inSex : 0,
			emigrationTraits
		);
	} // end of Nlines for loop

	return errorCode;
}

//---------------------------------------------------------------------------
int ReadTransferFile(speciesMap_t& simSpecies, landParams paramsLand, int transferType)
{
	int error = 0;
	switch (transferType) {

	case 0: // negative exponential dispersal kernel
		error = ReadTransferKernels(simSpecies, paramsLand);
		break; // end of negative exponential dispersal kernel

	case 1: // SMS
		ReadTransferSMS(simSpecies, paramsLand);
		break; // end of SMS

	case 2: // CRW
		error = ReadTransferCRW(simSpecies, paramsLand);
		break; // end of CRW

	default:
		error = 440;
		break;
	} // end of switch (TransferType)

	if (transferType > 0) {
		int nbHab = paramsLand.isArtificial ?
			paramsLand.nHab : paramsLand.nHabMax;
		for (auto& [sp, pSpecies] : simSpecies)
			pSpecies->createHabCostMort(nbHab);
	}

	return error;
}

int ReadTransferKernels(speciesMap_t& simSpecies, landParams paramsLand) {

	int inKernelType, inDistMort, inIndVar, simNb, inStageDep, inSexDep, inStage, inSex;
	float flushMort;
	int simNbFirstLine = 0;
	stageParams stageStruct;
	demogrParams dem;
	transferRules trfr;
	trfrKernelParams kernParams;
	int sexKernels = 0;
	species_id sp;
	Species* pSpecies;
	bool isFirstLine = true;
	int errorCode = 0;

	int nbLinesToRead = 1;
	for (int line = 0; line < nbLinesToRead; line++) {

		ifsTransferFile >> simNb >> sp >> inStageDep >> inSexDep 
			>> inKernelType >> inDistMort >> inIndVar;
		
		if (isFirstLine) {
			simNbFirstLine = simNb;
			pSpecies = simSpecies.at(sp);
			stageStruct = pSpecies->getStageParams();
			dem = pSpecies->getDemogrParams();
			trfr = pSpecies->getTransferRules();

			trfr.twinKern = (inKernelType == 1);
			trfr.distMort = (inDistMort == 1);
			sexKernels = 2 * inStageDep + inSexDep;
			trfr.indVar = (inIndVar == 1);
			trfr.sexDep = (inSexDep == 1);
			trfr.stgDep = (inStageDep == 1);

			// Set expected nb of lines of input
			if (trfr.stgDep) nbLinesToRead *= stageStruct.nStages;
			if (trfr.sexDep) nbLinesToRead *= dem.repType == 0 ? 1 : 2;
			
			pSpecies->setTrfrRules(trfr);
		}
		if (simNb != simNbFirstLine) { // serious problem
			errorCode = 400;
		}
		ifsTransferFile >> inStage >> inSex;

		if (dem.repType == 0) {
			if (sexKernels == 1 || sexKernels == 3) 
				errorCode = 401;
		}
		else if (sexKernels == 2 || sexKernels == 3) 
				errorCode = 403;

		switch (sexKernels) {

		case 0: // no sex / stage dependence
			ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			pSpecies->setSpKernTraits(0, 0, kernParams, paramsLand.resol);
			break;

		case 1: // sex-dependent
			if (trfr.twinKern) {
				ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				ifsTransferFile >> kernParams.meanDist1; 
				kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(0, inSex, kernParams, paramsLand.resol);
			break;

		case 2: // stage-dependent
			if (trfr.twinKern) {
				ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				ifsTransferFile >> kernParams.meanDist1;
				kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(inStage, 0, kernParams, paramsLand.resol);
			break;

		case 3: // sex- & stage-dependent
			if (trfr.twinKern) {
				ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				ifsTransferFile >> kernParams.meanDist1; kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(inStage, inSex, kernParams, paramsLand.resol);
			break;
		} // end of switch (sexkernels)

		// mortality
		if (inStage == 0 && inSex == 0) {
			trfrMortParams mort;
			ifsTransferFile >> mort.fixedMort >> mort.mortAlpha >> mort.mortBeta;
			pSpecies->setMortParams(mort);
		}
		else for (int i = 0; i < 3; i++) 
			ifsTransferFile >> flushMort;

		if (isFirstLine) pSpecies->setTrfrRules(trfr);
		isFirstLine = false;

	} // end of lines for loop

	return errorCode;
}

void ReadTransferSMS(speciesMap_t& simSpecies, const landParams& paramsLand) {

	int inIndVar, inSMType, inAlphaDB, inBetaDB, inStraightenPath, simNb;
	float inHabMort, flushHabMort, inMortHabitat, inMortMatrix;
	int inCostHab, inCostMatrix;
	trfrMovtParams move;
	species_id sp;

	ifsTransferFile >> simNb >> sp >> inIndVar >> move.pr >> move.prMethod >> move.dp
		>> move.memSize >> move.gb >> move.goalType >> inAlphaDB >> inBetaDB
		>> inStraightenPath >> inSMType >> move.stepMort;
	Species* pSpecies = simSpecies.at(sp);
	transferRules trfr = pSpecies->getTransferRules();

	trfr.indVar = (inIndVar == 1);
	if (move.goalType == 2) { // dispersal bias
		move.alphaDB = inAlphaDB;
		move.betaDB = inBetaDB;
	}
	trfr.habMort = (inSMType == 1);
	move.straightenPath = (inStraightenPath == 1);

	if (!paramsLand.isArtificial) { // imported landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			if (trfr.habMort) { // habitat-dependent step mortality
				for (int i = 0; i < paramsLand.nHabMax; i++) {
					ifsTransferFile >> inHabMort;
					pSpecies->setHabMort(i, inHabMort);
				}
			}
			else { // constant step mortality
				for (int i = 0; i < paramsLand.nHabMax; i++) 
					ifsTransferFile >> flushHabMort;
			}
		}
	}
	else { // artificial landscape
		if (trfr.habMort) { // habitat-dependent step mortality
			// values are for habitat (hab=1) then for matrix (hab=0)
			ifsTransferFile >> inMortHabitat >> inMortMatrix;
			pSpecies->setHabMort(1, inMortHabitat);
			pSpecies->setHabMort(0, inMortMatrix);
		}
		else { // constant step mortality
			ifsTransferFile >> flushHabMort >> flushHabMort;
		}
	}

	trfr.usesCosts = gUseSMSCosts.at(sp);

	if (!paramsLand.isArtificial) { // imported landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			for (int i = 0; i < paramsLand.nHabMax; i++) {
				ifsTransferFile >> inCostHab;
				if (!trfr.usesCosts)
					pSpecies->setHabCost(i, inCostHab);
			}
		}
	}
	else { // artificial landscape
		ifsTransferFile >> inCostHab >> inCostMatrix;
		if (!trfr.usesCosts) {
			// costs are for habitat (hab=1) then for matrix (hab=0)
			pSpecies->setHabCost(1, inCostHab);
			pSpecies->setHabCost(0, inCostMatrix);
		}
	}
	pSpecies->setTrfrRules(trfr);
	pSpecies->setSpMovtTraits(move);
}

int ReadTransferCRW(speciesMap_t& simSpecies, const landParams& paramsLand) {

	int inIndVar, inStraightenPath, inSMconst, simNb;
	float inHabMort, flushHabMort;
	species_id sp;
	int error = 0;
	trfrMovtParams move;
	ifsTransferFile >> simNb >> sp >> inIndVar;
	Species* pSpecies = simSpecies.at(sp);
	transferRules trfr = pSpecies->getTransferRules();
	trfr.indVar = inIndVar != 0;

	ifsTransferFile >> move.stepLength >> move.rho;
	ifsTransferFile >> inStraightenPath >> inSMconst >> move.stepMort;

	trfr.habMort = inSMconst != 0;
	move.straightenPath = inStraightenPath != 0;

	//Habitat-dependent per step mortality
	if (trfr.habMort && paramsLand.rasterType != 0)
		error = 434;

	if (!paramsLand.isArtificial && paramsLand.rasterType == 0) { // imported habitat codes landscape
		if (trfr.habMort) { // habitat-dependent step mortality
			for (int i = 0; i < paramsLand.nHabMax; i++) {
				ifsTransferFile >> inHabMort;
				pSpecies->setHabMort(i, inHabMort);
			}
		}
		else { // constant step mortality
			for (int i = 0; i < paramsLand.nHabMax; i++) 
				ifsTransferFile >> flushHabMort;
		}
	}
	pSpecies->setTrfrRules(trfr);
	pSpecies->setSpMovtTraits(move);
	return error;
}

//---------------------------------------------------------------------------
int ReadSettlement(speciesMap_t& simSpecies)
{
	int simNb, simNbFirstLine = 0, inStageDep, inSexDep, inStage, inSex;
	bool isFirstline = true;
	bool mustFindMate;
	int errorCode = 0;
	demogrParams dem;
	stageParams sstruct;
	transferRules trfr;
	settleType sett;
	settleRules srules;
	settleSteps ssteps;
	settleTraits settleDD;
	int sexSettle = 0, inSettleType = 0, inDensDep, inIndVar, inFindMate;
	species_id sp;
	Species* pSpecies;

	isFirstline = true;

	int nbLinesToRead = 1;
	for (int line = 0; line < nbLinesToRead; line++) {

		ifsSettlementFile >> simNb >> sp >> inStageDep >> inSexDep >> inStage >> inSex;
		
		if (!trfr.usesMovtProc) { // dispersal kernel
			ifsSettlementFile >> inSettleType >> inFindMate;
		}
		else {
			ifsSettlementFile >> inDensDep >> inIndVar >> inFindMate;
		}
		mustFindMate = (inFindMate == 1);

		if (isFirstline) {

			simNbFirstLine = simNb;
			pSpecies = simSpecies.at(sp);
			dem = pSpecies->getDemogrParams();
			sstruct = pSpecies->getStageParams();
			trfr = pSpecies->getTransferRules();
			sett = pSpecies->getSettle();

			sett.stgDep = (inStageDep == 1);
			sett.sexDep = (inSexDep == 1);
			sett.indVar = trfr.usesMovtProc ? inIndVar == 1 : false; // no ind var for kernels
			pSpecies->setSettle(sett);

			// update no.of lines according to known stage- and sex-dependency
			int nbSexesDisp = dem.repType == 0 ? 1 : 2;
			if (sett.sexDep) nbLinesToRead *= nbSexesDisp;
			if (sett.stgDep) nbLinesToRead *= sstruct.nStages;
		}

		if (simNb != simNbFirstLine) { // serious problem
			errorCode = 500;
		}

		if (trfr.usesMovtProc) {
			// Movement process
			bool hasMales = dem.repType > 0;
			if (!hasMales && sett.sexDep) 
				errorCode = 508;
			if (!dem.stageStruct && sett.stgDep) 
				errorCode = 509;

			ifsSettlementFile >> ssteps.minSteps >> ssteps.maxSteps >> ssteps.maxStepsYr;
			ifsSettlementFile >> settleDD.s0 >> settleDD.alpha >> settleDD.beta;

			int stageToSet = sett.stgDep ? inStage : 0;
			int sexToSet = sett.sexDep ? inSex : 0;
			srules = pSpecies->getSettRules(stageToSet, sexToSet);
			srules.densDep = (inDensDep == 1);
			srules.findMate = (inFindMate == 1);

			pSpecies->setSettRules(stageToSet, sexToSet, srules);
			pSpecies->setSteps(stageToSet, sexToSet, ssteps);

			if (srules.densDep) {
				pSpecies->setSpSettTraits(stageToSet, sexToSet, settleDD);
			}

			if (!sett.stgDep) {
				if (!sett.sexDep) {
					if (dem.stageStruct) { // model is structured - also set parameters for all stages
						for (int stg = 1; stg < sstruct.nStages; stg++) {
							pSpecies->setSettRules(stg, 0, srules);
							pSpecies->setSteps(stg, 0, ssteps);
							pSpecies->setSpSettTraits(stg, 0, settleDD);
							if (hasMales) { // model is sexual - also set parameters for males
								pSpecies->setSettRules(stg, 1, srules);
								pSpecies->setSteps(stg, 1, ssteps);
								if (srules.densDep && !sett.indVar) 
									pSpecies->setSpSettTraits(stg, 1, settleDD);
							}
						}
					}
					else {
						if (hasMales) { // model is sexual - also set parameters for males
							pSpecies->setSettRules(0, 1, srules);
							pSpecies->setSteps(0, 1, ssteps);
							if (srules.densDep) {
								pSpecies->setSpSettTraits(0, 1, settleDD);
							}
						}
					}
				}
				else { // stage-dep but not sex-dep
					if (dem.stageStruct) { // model is structured - also set parameters for all stages
						for (int stg = 1; stg < sstruct.nStages; stg++) {
							pSpecies->setSettRules(stg, sexToSet, srules);
							pSpecies->setSteps(stg, sexToSet, ssteps);
							if (srules.densDep && !sett.indVar) 
								pSpecies->setSpSettTraits(stg, sexToSet, settleDD);
						}
					}
				}
			}
			else { // not stage-dep
				if (!sett.sexDep) {
					if (hasMales) { // model is sexual - also set parameters for males
						pSpecies->setSettRules(stageToSet, 1, srules);
						pSpecies->setSteps(stageToSet, 1, ssteps);
						if (srules.densDep) {
							pSpecies->setSpSettTraits(stageToSet, 1, settleDD);
						}
					}
				}
			}
		} // end of movement model
		else { // dispersal kernel

			bool hasMales = dem.repType > 0;
			if (!hasMales && sett.sexDep)
				errorCode = 501;
			if (!dem.stageStruct && sett.stgDep)
				errorCode = 502;
			if (!sett.stgDep && (inSettleType == 1 || inSettleType == 3) && !dem.stageStruct)
				errorCode = 503;
			if (!sett.sexDep && mustFindMate && !hasMales)
				errorCode = 504;

			int stageToSet = sett.stgDep ? inStage : 0;
			int sexToSet = sett.sexDep ? inSex : 0;
			srules = pSpecies->getSettRules(stageToSet, sexToSet);

			srules.wait = inSettleType == 1 || inSettleType == 3;
			srules.goToNeighbourLocn = inSettleType == 2 || inSettleType == 3;
			srules.findMate = mustFindMate;
			pSpecies->setSettRules(stageToSet, sexToSet, srules);

			if (!sett.stgDep && dem.stageStruct) {
				// Must set other stages
				if (!sett.sexDep) {
					for (int stg = 0; stg < sstruct.nStages; stg++) {
						pSpecies->setSettRules(stg, 0, srules);
						if (hasMales) { // model is sexual - also set parameters for males
							pSpecies->setSettRules(stg, 1, srules);
						}
					}
				}
				else {
					for (int stg = 1; stg < sstruct.nStages; stg++) {
						pSpecies->setSettRules(stg, sexToSet, srules);
					}
				}
			}
			if (!sett.sexDep && hasMales) { // Must set males
				pSpecies->setSettRules(sett.stgDep ? stageToSet : 0, 1, srules);
			}
		} // end of dispersal kernel

		isFirstline = false;

	} // end of for line loop

	return errorCode;
}

//---------------------------------------------------------------------------
int ReadInitialisation(const landParams& paramsLand, speciesMap_t& simSpecies)
{
	string inputDir = paramsSim->getDir(1);

	int simNb, maxcells;
	float totalProps;
	int errorCode = 0;

	species_id sp;
	ifsInitFile >> sp;
	Species* pSpecies = simSpecies.at(sp);
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	initParams init;
	ifsInitFile >> simNb >> init.seedType >> init.freeType >> init.spDistType;

	if (init.seedType == 1 && !gUseSpeciesDist.at(sp)) 
		errorCode = 601;

	if (paramsLand.usesPatches) 
		ifsInitFile >> init.initDens >> init.indsHa;
	else 
		ifsInitFile >> init.initDens >> init.indsCell;

	ifsInitFile >> init.minSeedX >> init.maxSeedX 
		>> init.minSeedY >> init.maxSeedY
		>> init.nSeedPatches >> init.nSpDistPatches
		>> init.initFrzYr >> init.restrictRows
		>> init.restrictFreq >> init.finalFrzYr
		>> init.indsFile;

	init.restrictRange = (init.seedType == 0 && init.restrictRows > 0);

	if (dem.stageStruct) {
		float propStage;
		ifsInitFile >> init.initAge;
		totalProps = 0.0;
		for (int stg = 1; stg < sstruct.nStages; stg++) {
			ifsInitFile >> propStage;
			if (init.seedType != 2) {
				totalProps += propStage;
				pSpecies->setProp(stg, propStage);
			}
		}
		if (init.seedType!=2 && totalProps != 1.0) { 
			throw logic_error("The proportion of initial individuals in each stage doesn not sum to 1.");
		}
	}

	pSpecies->setInitParams(init);

	switch (init.seedType) {
	case 0: // free initialisation

		// Set initial distribution parameters
		if (init.minSeedX == gEmptyVal)
			init.minSeedX = 0;
		if (init.minSeedY == gEmptyVal)
			init.minSeedY = 0;
		if (init.maxSeedX == gEmptyVal) 
			init.maxSeedX = paramsLand.maxX;
		if (init.maxSeedY == gEmptyVal) 
			init.maxSeedY = paramsLand.maxY;
		if (init.minSeedY > init.maxSeedY || init.minSeedX > init.maxSeedX) {
			errorCode = 603;
		}
		maxcells = (init.maxSeedY - init.minSeedY) * (init.maxSeedX - init.minSeedX);
		if (init.freeType == 0 && init.nSeedPatches > maxcells) 
			errorCode = 602;

		// Done, initialisation takes place in Community::initialise()
		break;

	case 1: // from species distribution
		
		// Nothing, input is processed in ReadSpLand()
		break;

	case 2: // from initial individuals file
		if (init.indsFile != prevInitialIndsFile) {
			// read and store the list of individuals to be initialised
			ReadInitIndsFile(pSpecies, 0, paramsLand, (inputDir + init.indsFile));
			prevInitialIndsFile = init.indsFile;
		}
		break;
	default:
		throw logic_error("SeedType must be 0, 1, or 2.");
		break;
	}
	return errorCode;
}

//---------------------------------------------------------------------------
int ReadInitIndsFile(Species* pSpecies, int option, const landParams& paramsLand, string indsfile) {
	string header;
	demogrParams dem = pSpecies->getDemogrParams();
	initParams init = pSpecies->getInitParams();

	if (option == 0) { // open file and read header line
		ifsInitIndsFile.open(indsfile.c_str());
		string header;
		int nheaders = 3;
		if (paramsLand.usesPatches) nheaders++;
		else nheaders += 2;
		if (dem.repType > 0) nheaders++;
		if (dem.stageStruct) nheaders += 2;
		for (int i = 0; i < nheaders; i++) ifsInitIndsFile >> header;
		pSpecies->resetInitInds();
		//	return 0;
	}

	if (option == 9) { // close file
		if (ifsInitIndsFile.is_open()) {
			ifsInitIndsFile.close(); 
			ifsInitIndsFile.clear();
		}
		return 0;
	}

	// Read data lines;
	initInd iind;
	int ninds;
	int totinds = 0;

	iind.year = gEmptyVal;
	ifsInitIndsFile >> iind.year;
	bool must_stop = (iind.year == gEmptyVal);

	while (!must_stop) {
		ifsInitIndsFile >> iind.speciesID;

		if (paramsLand.usesPatches) {
			ifsInitIndsFile >> iind.patchID;
			iind.x = iind.y = 0;
		}
		else {
			ifsInitIndsFile >> iind.x >> iind.y; 
			iind.patchID = 0;
		}
		ifsInitIndsFile >> ninds;

		if (dem.repType > 0) 
			ifsInitIndsFile >> iind.sex;
		else 
			iind.sex = 0;

		if (dem.stageStruct) {
			ifsInitIndsFile >> iind.age >> iind.stage;
		}
		else {
			iind.age = iind.stage = 0;
		}
		for (int i = 0; i < ninds; i++) {
			totinds++;
			pSpecies->addInitInd(iind);
		}

		iind.year = gEmptyVal;
		ifsInitIndsFile >> iind.year;
		if (iind.year == gEmptyVal || ifsInitIndsFile.eof())
			must_stop = true;

	} // end of while loop

	if (ifsInitIndsFile.is_open()) ifsInitIndsFile.close();
	ifsInitIndsFile.clear();

	return totinds;
}

//---------------------------------------------------------------------------
void RunBatch()
{
	int land_nr;
	int read_error;
	bool areParamsOk;
	simParams sim = paramsSim->getSim();

	// Create species
	speciesMap_t allSpecies;
	for (species_id sp : gSpeciesNames) {
		allSpecies.emplace(sp, new Species);
	}

	Landscape* pLandscape = nullptr; 

	// Open landscape batch file and read header record
	ifsLandFile.open(landFile);
	if (!ifsLandFile.is_open()) {
		cout << endl << "Error opening landFile - aborting batch run" << endl;
		return;
	}
	flushHeaders(ifsLandFile);

	for (int j = 0; j < gNbLandscapes; j++) {

		// Create new landscape
		if (pLandscape != nullptr) delete pLandscape;
		pLandscape = new Landscape(gSpeciesNames);
		bool landOK = true;

		land_nr = ReadLandFile(pLandscape);
		if (land_nr <= 0) { // error condition
			string msg = "Error code " + to_string(-land_nr)
				+ " returned from reading LandFile - aborting batch run";
			cout << endl << msg << endl;
			ifsLandFile.close();  
			ifsLandFile.clear();
			return;
		}

		landParams paramsLand = pLandscape->getLandParams();
		paramsLand.usesPatches = gUsesPatches;
		paramsLand.resol = gResol;
		paramsLand.rasterType = gLandType;
		if (gLandType == 9) {
			paramsLand.isArtificial = true;
			paramsLand.nHab = 2;
		}
		else {
			paramsLand.isArtificial = false;
			paramsLand.isDynamic = gDynLandFileName != "NULL";
		}
		paramsLand.nHabMax = gMaxNbHab;
		pLandscape->setLandParams(paramsLand, true);

		if (gLandType != 9) { // imported landscape
			string pathToHabMap = paramsSim->getDir(1) + gHabMapName;			
			map<species_id, string> pathsToPatchMaps, pathsToCostMaps, pathsToSpDistMaps;
			ReadSpLandFile(
				ifsSpLandFile,
				pathsToPatchMaps,
				pathsToCostMaps,
				pathsToSpDistMaps,
				gUseSpeciesDist
			);

			if (pLandscape->readLandscape(0, pathToHabMap, pathsToPatchMaps) != 0) {
				cout << "Error reading landscape" << endl;
				landOK = false;
			}
			if (sim.batchMode) {
				if (pLandscape->readCosts(pathsToCostMaps) < 0) {
					cout << "Error reading landscape" << endl;
					landOK = false;
				}
			}

			if (paramsLand.isDynamic) {
				if (ReadDynLandFile(pLandscape) != 0) {
					cout << "Error reading dynamic landscape" << endl;
					landOK = false;
				}
			}
			if (gLandType == 0) pLandscape->updateHabitatIndices();

			// Species Distribution
			for (auto& [sp, pathToMap] : pathsToSpDistMaps) {
				string distname = paramsSim->getDir(1) + pathToMap;
				if (pLandscape->newDistribution(sp, distname) != 0) {
					cout << endl << "Error reading initial distribution for landscape "
						<< land_nr << " - aborting" << endl;
					landOK = false;
				}
			}
		} // end of imported landscape

		if (landOK) {

			// Open all other batch files and read headers
			{
				ifsSimFile.open(gSimFile);
				if (!ifsSimFile.is_open()) {
					cout << endl << "Error opening SimFile - aborting batch run" << endl;
					return;
				}
				flushHeaders(ifsSimFile);

				flushHeaders(ifsParamFile);
				ifsParamFile.open(gParametersFile);
				if (!ifsParamFile.is_open()) {
					cout << endl << "Error opening ParameterFile - aborting batch run" << endl;
					return;
				}
				flushHeaders(ifsParamFile);

				if (gUsesStageStruct) {
					ifsStageStructFile.open(stageStructFile);
					flushHeaders(ifsStageStructFile);
				}

				ifsEmigrationFile.open(emigrationFile);
				flushHeaders(ifsEmigrationFile);

				ifsTransferFile.open(transferFile);
				flushHeaders(ifsTransferFile);

				ifsSettlementFile.open(settleFile);
				flushHeaders(ifsSettlementFile);

				ifsInitFile.open(initialFile);
				flushHeaders(ifsInitFile);

				if (gAnyUsesGenetics) {
					ifsGeneticsFile.open(geneticsFile.c_str());
					flushHeaders(ifsGeneticsFile);
					ifsTraitsFile.open(traitsFile.c_str());
					flushHeaders(ifsTraitsFile);
				}
			}

			for (auto& thisSimulation : gSpInputOpt) {

				int simNb = thisSimulation.first;
				auto& simOptionsMap = thisSimulation.second;
				// Subset species that are used in this simulation
				speciesMap_t simSpecies;
				for (auto& sp : views::keys(simOptionsMap)) {
					simSpecies.emplace(sp, allSpecies.at(sp));
				}

				// Load parameters for this simulation
				areParamsOk = true;
				ReadSimParameters();

				// Read one line of input per simulation and species
				for (int s = 0; s < simSpecies.size(); s++) {
					// species don't have to be read in order
					read_error = ReadParameters(pLandscape, simSpecies);
					if (read_error) areParamsOk = false;
					if (gUsesStageStruct) ReadStageStructure(simSpecies);
					read_error = ReadEmigration(simSpecies);
					if (read_error) areParamsOk = false;
					read_error = ReadTransferFile(simSpecies, paramsLand, gTransferType);
					if (read_error) areParamsOk = false;
					read_error = ReadSettlement(simSpecies);
					if (read_error) areParamsOk = false;
					read_error = ReadInitialisation(paramsLand, simSpecies);
					if (read_error) areParamsOk = false;

					if (gAnyUsesGenetics) {
						read_error = ReadGeneticsFile(simSpecies, ifsGeneticsFile);
						if (read_error) areParamsOk = false;
						read_error = ReadTraitsFile(simSpecies, ifsTraitsFile, simOptionsMap);
						if (read_error) areParamsOk = false;
					}
				}

				if (areParamsOk) {

					cout << endl << "Running simulation nr. "
						<< to_string(paramsSim->getSim().simulation)
						<< " on landscape no. " << to_string(land_nr) << endl;

					// for batch processing, include landscape number in parameter file name
					OutParameters(pLandscape, simSpecies);

					RunModel(pLandscape, simNb, simSpecies);

				}
				else {
					cout << endl << "Error in reading parameter file(s)" << endl;
				}

				// Empty species map for next simulation
				simSpecies.clear();

			} // end of loop through simulations

			// Close input files
			{
				ifsSimFile.close();
				ifsSimFile.clear();
				ifsParamFile.close();
				ifsParamFile.clear();
				if (gUsesStageStruct) {
					ifsStageStructFile.close();
					ifsStageStructFile.clear();
				}
				ifsEmigrationFile.close();
				ifsEmigrationFile.clear();
				ifsTransferFile.close();
				ifsTransferFile.clear();
				ifsSettlementFile.close();
				ifsSettlementFile.clear();
				ifsInitFile.close();
				ifsInitFile.clear();

				if (gAnyUsesGenetics) {
					ifsGeneticsFile.close();
					ifsGeneticsFile.clear();
					ifsTraitsFile.close();
					ifsTraitsFile.clear();
				}
			}
			if (pLandscape != nullptr) {
				delete pLandscape; 
				pLandscape = nullptr;
			}

		} // end of landOK condition

	} // end of nLandscapes loop

	for (auto& [sp, pSpecies] : allSpecies)
		delete pSpecies;

	ifsLandFile.close();  
	ifsLandFile.clear();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


