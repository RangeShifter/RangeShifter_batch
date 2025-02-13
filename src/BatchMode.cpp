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
ifstream ifsParamFile, ifsLandFile, ifsDynLandFile;
ifstream ifsSpDistFile, ifsStageStructFile, ifsTransMatrix;
ifstream ifsStageWeightsFile;
ifstream ifsEmigrationFile, ifsTransferFile, ifsSettlementFile;
ifstream ifsTraitsFile, ifsGeneticsFile;
ifstream ifsInitFile, ifsInitIndsFile;
ifstream ifsFecDens, ifsDevDens, ifsSurvDens;

ofstream batchLogOfs;


// global variables passed between parsing functions...
// should be removed eventually, maybe share variables through members of a class
int gBatchNb;
int gIsPatchModel, gResol, gLandType, gMaxNbHab, gUseSpeciesDist, gDistResol;
int gReproType;
int gNbRepSeasons;
int gStageStruct, gNbStages, gTransferType;
int gNbSexesDem;		// no. of explicit sexes for demographic model
int gNbSexesDisp;	// no. of explicit sexes for dispersal model
int gFirstSimNb = 0; // not great, globals should not be modified.
bool gHasGenetics = true;

set<int> gSimNbs; // record of simulation numbers to check input file use the same numbers

// Track trait-relevant options to check for coherency across input files, 
// e.g. if emig file says emigration is indvar, trait file should have d0 entry
map<int, TraitInputOptions> gTraitOptions;
vector<int> gNbTraitFileRows;

rasterdata landraster;
// ...including names of the input files
string parameterFile;
string landFile;
string gHabMapName, gPatchMapName, gDynLandFileName, gSpDistFileName, gNameCostFile;
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

float** matrix = NULL;	// temporary matrix used in batch mode
int matrixsize = 0; 		// size of temporary matrix

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
batchfiles ParseControlAndCheckInputFiles(string pathToControlFile, string inputDir, string outputDir)
{
	batchfiles b;
	int lines, nSimuls;
	int nbErrors = 0;
	string paramName, filename, fname, batchLogPath, header;
	string whichInputFile = "Control file";
	bool anyFormatError = false;
	
	// Open batch log
	batchLogPath = outputDir + "BatchLog.txt";
	batchLogOfs.open(batchLogPath.c_str());
	if (!batchLogOfs.is_open()) {
		cout << "Error opening batch output log file " << batchLogPath << endl;
		b.ok = false;
		return b;
	}

	// Open control file
	ifstream controlIfs{ pathToControlFile.c_str() };
	if (!controlIfs.is_open()) {
		cout << "Error opening Control file: " << pathToControlFile << endl;
		batchLogOfs << "Error opening Control file: " << pathToControlFile << endl;
		b.ok = false;
		if (batchLogOfs.is_open()) { 
			batchLogOfs.close(); 
			batchLogOfs.clear(); 
		}
		return b;
	}
	else {
		batchLogOfs << "Checking Control file " << pathToControlFile << endl;
	}

	// Check batch parameters

	controlIfs >> paramName >> gBatchNb;
	if (paramName == "BatchNum") {
		if (gBatchNb < 0) {
			BatchError(whichInputFile, -999, 19, "BatchNum"); nbErrors++;
		}
		else b.batchNum = gBatchNb;
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gIsPatchModel;
	if (paramName == "PatchModel") {
		if (gIsPatchModel != 0 && gIsPatchModel != 1) {
			BatchError(whichInputFile, -999, 1, "PatchModel"); nbErrors++;
		}
		else b.isPatchModel = gIsPatchModel;
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gResol;
	if (paramName == "Resolution") {
		if (gResol < 1) {
			BatchError(whichInputFile, -999, 11, "Resolution"); nbErrors++;
		}
		else b.resolution = gResol;
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
			if (gLandType == 9 && gIsPatchModel) {
				BatchError(whichInputFile, -999, 0, "LandType");
				batchLogOfs << "LandType may not be 9 for a patch-based model" << endl;
				nbErrors++;
			}
			else b.landType = gLandType;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gMaxNbHab;
	if (paramName == "MaxHabitats") {
		if (gLandType == 0) { // raster with unique habitat codes
			if (gMaxNbHab < 2) {
				BatchError(whichInputFile, -999, 12, "MaxHabitats"); nbErrors++;
			}
			else b.maxNbHab = gMaxNbHab;
		}
		else { // raster with habitat quality OR artificial landscape
			if (gMaxNbHab != 1) {
				BatchError(whichInputFile, -999, 0, " "); nbErrors++;
				batchLogOfs << "MaxHabitats must be 1 for LandType = " << gLandType << endl;
			}
			else {
				if (gLandType == 9) // artificial landscape
					// although the user enters 1, the actual number of habitats is 2
					b.maxNbHab = 2;
				else
					b.maxNbHab = gMaxNbHab;
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gUseSpeciesDist;
	if (paramName == "SpeciesDist") {
		if (gUseSpeciesDist != 0 && gUseSpeciesDist != 1) {
			BatchError(whichInputFile, -999, 1, "SpeciesDist"); nbErrors++;
		}
		else {
			if (gUseSpeciesDist != 0 && gLandType == 9) {
				BatchError(whichInputFile, -999, 0, "SpeciesDist");
				batchLogOfs << "SpeciesDist must be 0 for an artificial landscape" << endl;
				nbErrors++;

			}
			else b.speciesDist = gUseSpeciesDist;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gDistResol;
	if (paramName == "DistResolution") {
		if (gUseSpeciesDist == 1) { // distribution resolution is required
			if (gDistResol < gResol) {
				BatchError(whichInputFile, -999, 0, "DistResolution");
				batchLogOfs << "DistResolution may not be less than Resolution" << endl;
				nbErrors++;
			}
			else {
				if (gDistResol % gResol) {
					BatchError(whichInputFile, -999, 0, "DistResolution");
					batchLogOfs << "DistResolution must be an integer multiple of Resolution" << endl;
					nbErrors++;
				}
				else b.distResol = gDistResol;
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gReproType;
	gNbSexesDem = gNbSexesDisp = 0;
	if (paramName == "Reproduction") {
		if (gReproType != 0 && gReproType != 1 && gReproType != 2) {
			BatchError(whichInputFile, -999, 2, "Reproduction"); nbErrors++;
		}
		else {
			switch (gReproType) {
			case 0: { gNbSexesDem = 1; gNbSexesDisp = 1; break; }
			case 1: { gNbSexesDem = 1; gNbSexesDisp = 2; break; }
			case 2: { gNbSexesDem = 2; gNbSexesDisp = 2; break; }
			}
			b.reproType = gReproType; 
			b.sexesDem = gNbSexesDem; 
			b.nbSexesDisp = gNbSexesDisp;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gNbRepSeasons;
	if (paramName == "RepSeasons") {
		if (gNbRepSeasons < 1) {
			BatchError(whichInputFile, -999, 11, "RepSeasons"); nbErrors++;
		}
		else b.nbRepSeasons = gNbRepSeasons;
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gStageStruct;
	if (paramName == "StageStruct") {
		if (gStageStruct != 0 && gStageStruct != 1) {
			BatchError(whichInputFile, -999, 1, "StageStruct"); nbErrors++;
		}
		else b.isStageStruct = gStageStruct;
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gNbStages;
	if (paramName == "Stages") {
		if (gStageStruct) {
			if (gNbStages < 2 || gNbStages > 10) {
				BatchError(whichInputFile, -999, 0, " "); nbErrors++;
				batchLogOfs << "Stages must be between 2 and 10" << endl;
			}
			b.nbStages = gNbStages;
		}
		else { // non-stage-structured model must have 2 stages
			b.nbStages = gNbStages = 2;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlIfs >> paramName >> gTransferType;
	if (paramName == "Transfer") {
		if (gTransferType < 0 || gTransferType > 2) {
			BatchError(whichInputFile, -999, 2, "Transfer"); nbErrors++;
		}
		else b.transferType = gTransferType;
	}
	else anyFormatError = true; // wrong control file format

	if (anyFormatError || nbErrors > 0) { // terminate batch error checking
		if (anyFormatError) {
			CtrlFormatError();
		}
		batchLogOfs << endl
			<< "*** Model parameters in Control file must be corrected before further input file checks are conducted"
			<< endl;
		batchLogOfs.close(); 
		batchLogOfs.clear();
		b.ok = false;
		controlIfs.close(); 
		controlIfs.clear();
		return b;
	}

	// Check parameter file
	controlIfs >> paramName >> filename;
	if (paramName == "ParameterFile" && !anyFormatError) {
		fname = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << fname << endl;
		ifsParamFile.open(fname.c_str());
		if (ifsParamFile.is_open()) {
			b.nSimuls = CheckParameterFile();
			if (b.nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramName, b.nSimuls, 0);
				parameterFile = fname;
			}
			ifsParamFile.close();
		}
		else {
			OpenError(paramName, fname); b.ok = false;
			cout << "Unable to open ParameterFile" << endl;
		}
		ifsParamFile.clear();
		if (!b.ok) {
			batchLogOfs << endl
				<< "*** ParameterFile must be corrected before further input file checks are conducted"
				<< endl;
			batchLogOfs.close(); 
			batchLogOfs.clear();
			b.ok = false;
			controlIfs.close(); 
			controlIfs.clear();
			return b;
		}
	}
	else anyFormatError = true; // wrong control file format
	if (ifsParamFile.is_open()) ifsParamFile.close();
	ifsParamFile.clear();

	// Check land file
	controlIfs >> paramName >> filename;
	if (paramName == "LandFile" && !anyFormatError) {
		fname = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << fname << endl;
		ifsLandFile.open(fname.c_str());
		if (ifsLandFile.is_open()) {
			lines = CheckLandFile(gLandType, inputDir);
			if (lines < 0) {
				b.ok = false;
				if (lines < -111)
					batchLogOfs << "*** Format error in " << paramName << endl;
			}
			else {
				FileOK(paramName, lines, 1);
				landFile = fname; 
				b.nLandscapes = lines;
			}
			ifsLandFile.close();
		}
		else {
			OpenError(paramName, fname); b.ok = false;
		}
		ifsLandFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check stage structure file if required file
	controlIfs >> paramName >> filename;
	batchLogOfs << endl;
	if (paramName == "StageStructFile" && !anyFormatError) {
		if (filename == "NULL") {
			if (gStageStruct) {
				batchLogOfs << "*** File name is required for " << paramName << endl;
				b.ok = false;
			}
			else b.stageStructFile = filename;
		}
		else { // filename is not NULL
			if (gStageStruct) { // check file only if it is required
				fname = inputDir + filename;
				batchLogOfs << "Checking " << paramName << " " << fname << endl;
				ifsStageStructFile.open(fname.c_str());
				if (ifsStageStructFile.is_open()) {
					nSimuls = CheckStageFile(inputDir);
					if (nSimuls < 0) {
						b.ok = false;
					}
					else {
						FileOK(paramName, nSimuls, 0);
						if (nSimuls != b.nSimuls) {
							SimulnCountError(filename); b.ok = false;
						}
						else stageStructFile = fname;
					}
					ifsStageStructFile.close();
				}
				else {
					OpenError(paramName, fname); b.ok = false;
				}
				ifsStageStructFile.clear();
			} // end of required
			else { // file is not required, and filename should be NULL
				if (filename != "NULL") {
					batchLogOfs << "*** File name for stageStructFile should be NULL as StageStruct = "
						<< gStageStruct << endl;
					b.ok = false;
				}
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	// Check emigration file
	controlIfs >> paramName >> filename;
	if (paramName == "EmigrationFile" && !anyFormatError) {
		fname = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << fname << endl;
		ifsEmigrationFile.open(fname.c_str());
		if (ifsEmigrationFile.is_open()) {
			nSimuls = CheckEmigFile();
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); 
					b.ok = false;
				}
				else emigrationFile = fname;
			}
			ifsEmigrationFile.close();
		}
		else {
			OpenError(paramName, fname); 
			b.ok = false;
		}
		ifsEmigrationFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check transfer file
	controlIfs >> paramName >> filename;
	if (paramName == "TransferFile" && !anyFormatError) {
		fname = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << fname << endl;
		ifsTransferFile.open(fname.c_str());
		if (ifsTransferFile.is_open()) {
			nSimuls = CheckTransferFile(inputDir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else transferFile = fname;
			}
			ifsTransferFile.close(); ifsTransferFile.clear();
		}
		else {
			OpenError(paramName, fname); b.ok = false;
		}
		ifsTransferFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check settlement file
	controlIfs >> paramName >> filename;
	if (paramName == "SettlementFile" && !anyFormatError) {
		fname = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << fname << endl;
		ifsSettlementFile.open(fname.c_str());
		if (ifsSettlementFile.is_open()) {
			nSimuls = CheckSettleFile();
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); 
					b.ok = false;
				}
				else settleFile = fname;
			}
			ifsSettlementFile.close();
		}
		else {
			OpenError(paramName, fname); 
			b.ok = false;
		}
		ifsSettlementFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check genetics file if required file
	controlIfs >> paramName >> filename;
	batchLogOfs << endl;
	if (paramName == "GeneticsFile" && !anyFormatError) {
		if (filename == "NULL") {
			bool anyIsEmigIndVar = false;
			bool anyIsSettIndVar = false;
			bool anyIsKernTransfIndVar = false;
			bool anyIsSMSTransferIndVar = false;
			for (auto const& [simNb, traitOpt] : gTraitOptions) {
				if (traitOpt.isEmigIndVar) anyIsEmigIndVar = true;
				if (traitOpt.isSettIndVar) anyIsSettIndVar = true;
				if (traitOpt.isKernTransfIndVar) anyIsKernTransfIndVar = true;
				if (traitOpt.isSMSTransfIndVar) anyIsSMSTransferIndVar = true;
			}
			if (anyIsEmigIndVar || anyIsSettIndVar 
				|| anyIsKernTransfIndVar
				|| anyIsSMSTransferIndVar
				)
			{
				batchLogOfs << "Error: GeneticsFile is NULL but one or more dispersal traits has been set to IndVar." << endl;
				b.ok = false;
			}
			else {
				gHasGenetics = false;
				batchLogOfs << "No genetics required " << paramName << endl;
			}
		}
		else {
			gHasGenetics = true;
			fname = inputDir + filename;
			batchLogOfs << "Checking " << paramName << " " << fname << endl;
			ifsGeneticsFile.open(fname.c_str());
			if (ifsGeneticsFile.is_open()) {
				nSimuls = CheckGeneticsFile(inputDir);
				if (nSimuls < 0) {
					b.ok = false;
				}
				else {
					FileOK(paramName, nSimuls, 0);
					geneticsFile = fname;
				}
				ifsGeneticsFile.close();
			}
			else {
				OpenError(paramName, fname); 
				b.ok = false;
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
			if (gHasGenetics)
			{
				batchLogOfs << "Error: Genetics are enabled but no TraitsFile is provided." << endl;
				b.ok = false;
			}
		}
		else {
			fname = inputDir + filename;
			batchLogOfs << "Checking " << paramName << " " << fname << endl;
			ifsTraitsFile.open(fname.c_str());
			if (ifsTraitsFile.is_open()) {
				nSimuls = CheckTraitsFile(inputDir);
				if (nSimuls < 0) {
					b.ok = false;
				}
				else {
					FileOK(paramName, nSimuls, 0);
					traitsFile = fname;
				}
				ifsTraitsFile.close();
			}
			else {
				OpenError(paramName, filename);
				b.ok = false;
			}
			if (ifsTraitsFile.is_open()) ifsTraitsFile.close();
			ifsTraitsFile.clear();
		}
	}

	// Check initialisation file
	controlIfs >> paramName >> filename;
	if (paramName == "InitialisationFile" && !anyFormatError) {
		fname = inputDir + filename;
		batchLogOfs << endl << "Checking " << paramName << " " << fname << endl;
		ifsInitFile.open(fname.c_str());
		if (ifsInitFile.is_open()) {
			nSimuls = CheckInitFile(inputDir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramName, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else initialFile = fname;
			}
			ifsInitFile.close();
		}
		else {
			OpenError(paramName, fname); b.ok = false;
		}
		ifsInitFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	if (anyFormatError) {
		CtrlFormatError();
		b.ok = false;
	}

	if (controlIfs.is_open()) { controlIfs.close(); controlIfs.clear(); }
	if (batchLogOfs.is_open()) { batchLogOfs.close(); batchLogOfs.clear(); }

	return b;
}

//---------------------------------------------------------------------------
int CheckParameterFile()
{
	string header, Kheader, intext;
	int i, simNb, inReplicates, inYears;
	int inAbsorb, inGradient, inShifting, inShiftStart, inShiftEnd, inEnvStoch, inStochType;
	int inOptimum;
	int inLocalExt, inSaveMaps;
	int prevsimul = 0;
	float inMinR, inMaxR, inMinK, inMaxK, sum_K, min_K, max_K;
	float inGradSteep, inGradScalingFactor, inLocalExtOpt, inShiftRate;
	float inStochAC, inStochStD, inLocalExtProb, inPropMales, inHarem;
	float inBc, inRmax, inK;
	int inOutStartPop, inOutStartInd, inOutStartTraitCell, inOutStartTraitRow;
	int inOutStartConn, inOutIntRange, inOutIntOcc, inOutIntPop, inOutIntInd;
	int inOutIntTraitCell, inOutIntTraitRow, inOutIntConn, inMapsInterval;
	int inSMSHeatMap, inDrawLoadedSp, inFixReplicateSeed;
	int nbErrors = 0;
	int nbKerrors = 0;
	string whichFile = "ParameterFile";

	// Parse header line;
	ifsParamFile >> header; if (header != "Simulation") nbErrors++;
	ifsParamFile >> header; if (header != "Replicates") nbErrors++;
	ifsParamFile >> header; if (header != "Years") nbErrors++;
	ifsParamFile >> header; if (header != "Absorbing") nbErrors++;
	ifsParamFile >> header; if (header != "Gradient") nbErrors++;
	ifsParamFile >> header; if (header != "GradSteep") nbErrors++;
	ifsParamFile >> header; if (header != "Optimum") nbErrors++;
	ifsParamFile >> header; if (header != "f") nbErrors++;
	ifsParamFile >> header; if (header != "LocalExtOpt") nbErrors++;
	ifsParamFile >> header; if (header != "Shifting") nbErrors++;
	ifsParamFile >> header; if (header != "ShiftRate") nbErrors++;
	ifsParamFile >> header; if (header != "ShiftStart") nbErrors++;
	ifsParamFile >> header; if (header != "ShiftEnd") nbErrors++;
	ifsParamFile >> header; if (header != "EnvStoch") nbErrors++;
	ifsParamFile >> header; if (header != "EnvStochType") nbErrors++;
	ifsParamFile >> header; if (header != "ac") nbErrors++;
	ifsParamFile >> header; if (header != "std") nbErrors++;
	ifsParamFile >> header; if (header != "minR") nbErrors++;
	ifsParamFile >> header; if (header != "maxR") nbErrors++;
	ifsParamFile >> header; if (header != "minK") nbErrors++;
	ifsParamFile >> header; if (header != "maxK") nbErrors++;
	ifsParamFile >> header; if (header != "LocalExt") nbErrors++;
	ifsParamFile >> header; if (header != "LocalExtProb") nbErrors++;
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
	ifsParamFile >> header; if (header != "FixReplicateSeed") nbErrors++;

	if (nbErrors > 0 || nbKerrors > 0) {
		FormatError(whichFile, nbErrors);
		batchLogOfs << "*** ParameterFile column headers are incorrect." << endl;
		if (nbKerrors > 0) {
			BatchError(whichFile, -999, 333, "K");
		}
		return -111;
	}

	// Parse data lines
	int whichLine = 1;
	int nSimuls = 0;
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
	else {
		prevsimul = gFirstSimNb = simNb; 
		nSimuls++;
	}
	while (simNb != -98765) {

		// Record simulation numbers to cross-check other files
		gSimNbs.insert(simNb);

		// Initialise trait option map with simulation numbers
		gTraitOptions.emplace(simNb, TraitInputOptions());

		ifsParamFile >> inReplicates; 
		if (inReplicates <= 0) { 
			BatchError(whichFile, whichLine, 11, "Replicates"); 
			nbErrors++; 
		}
		ifsParamFile >> inYears; 
		if (inYears <= 0) {
			BatchError(whichFile, whichLine, 11, "Years"); 
			nbErrors++; 
		}
		ifsParamFile >> inAbsorb;
		if (inAbsorb != 0 && inAbsorb != 1) { 
			BatchError(whichFile, whichLine, 1, "Absorbing"); 
			nbErrors++; 
		}
		ifsParamFile >> inGradient;
		if (gIsPatchModel) {
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
		if (inGradient == 4 && (inLocalExtOpt < 0.0 || inLocalExtOpt >= 1.0))
		{
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
		ifsParamFile >> inEnvStoch;
		if (gIsPatchModel == 0) { // cell-based model
			if (inEnvStoch != 0 && inEnvStoch != 1 && inEnvStoch != 2) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "EnvStoch must be 0, 1 or 2 for cell-based model" << endl;
				nbErrors++;
			}
		}
		else { // patch-based model
			if (inEnvStoch != 0 && inEnvStoch != 1) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "EnvStoch must be 0 or 1 for patch-based model" << endl;
				nbErrors++;
			}
		}
		ifsParamFile >> inStochType;
		if (inEnvStoch && (inStochType < 0 || inStochType > 1)) {
			BatchError(whichFile, whichLine, 1, "EnvStochType"); 
			nbErrors++;
		}
		ifsParamFile >> inStochAC;
		if (inEnvStoch && (inStochAC < 0.0 || inStochAC >= 1.0)) {
			BatchError(whichFile, whichLine, 20, "ac"); 
			nbErrors++; 
		}
		ifsParamFile >> inStochStD;
		if (inEnvStoch && (inStochStD <= 0.0 || inStochStD > 1.0)) {
			BatchError(whichFile, whichLine, 20, "std"); 
			nbErrors++; 
		}
		ifsParamFile >> inMinR;
		if (inEnvStoch && inStochType == 0 && inMinR <= 0.0) { 
			BatchError(whichFile, whichLine, 10, "minR"); 
			nbErrors++; 
		}
		ifsParamFile >> inMaxR;
		if (inEnvStoch && inStochType == 0 && inMaxR <= inMinR) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "maxR must be greater than minR" << endl;
			nbErrors++;
		}
		ifsParamFile >> inMinK >> inMaxK;
		if (inEnvStoch && inStochType == 1) {
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
		ifsParamFile >> inLocalExt;
		if (gIsPatchModel == 0) { // cell-based model
			if (inLocalExt < 0 || inLocalExt > 1) {
				BatchError(whichFile, whichLine, 1, "LocalExt");
				nbErrors++;
			}
			else {
				if (inGradient == 4) { // gradient in local extinction probability
					if (inLocalExt != 0) {
						BatchError(whichFile, whichLine, 0, " ");
						batchLogOfs << "LocalExt must be zero if Gradient is 4" << endl;
						nbErrors++;
					}
				}
			}
		}
		else { // patch-based model
			if (inLocalExt != 0) {
				BatchError(whichFile, whichLine, 0, "null");
				batchLogOfs << "LocalExt must be 0 for patch-based model" << endl;
				nbErrors++;
			}
		}
		ifsParamFile >> inLocalExtProb;
		if (gIsPatchModel == 0 && inLocalExt == 1 && (inLocalExtProb <= 0.0 || inLocalExtProb >= 1.0))
		{
			BatchError(whichFile, whichLine, 20, "LocalExtProb"); 
			nbErrors++;
		}
		ifsParamFile >> inPropMales;
		if (gReproType && (inPropMales <= 0.0 || inPropMales >= 1.0)) {
			BatchError(whichFile, whichLine, 0, "");
			batchLogOfs << "PropMales should be above 0 and below 1 for sexual models" << endl;
			nbErrors++;
		}
		ifsParamFile >> inHarem;
		if (gReproType == 2 && inHarem <= 0.0) {
			BatchError(whichFile, whichLine, 10, "Harem"); 
			nbErrors++; 
		}
		ifsParamFile >> inBc;
		if (gStageStruct == 0 && inBc <= 0.0) {
			BatchError(whichFile, whichLine, 10, "bc"); 
			nbErrors++; 
		}
		ifsParamFile >> inRmax;
		if (gStageStruct == 0 && inRmax <= 0.0) {
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
			if (inEnvStoch && inStochType == 1) { // environmental stochasticity in K
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
				if (inReplicates < 2 && inOutIntOcc > 0) {
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
		else {
			if (gIsPatchModel != 1 && inOutIntConn > 0) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLogOfs << "OutIntConn may be >0 only if PatchModel is 1" << endl;
				nbErrors++;
			}
		}
		
		ifsParamFile >> inSMSHeatMap;
		if (inSMSHeatMap != 0 && inSMSHeatMap != 1) {
			BatchError(whichFile, whichLine, 1, "SMSHeatMap");
			nbErrors++;
		}
		ifsParamFile >> inFixReplicateSeed;
		if (inFixReplicateSeed != 0 && inFixReplicateSeed != 1) {
			BatchError(whichFile, whichLine, 1, "FixReplicateSeed");
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
			nSimuls++;
		}
	} // end of while loop
	if (!ifsParamFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}

	if (nbErrors > 0) return -111;
	else return nSimuls;
}

int CheckLandFile(int landtype, string indir)
{
	string fileName, header, inputText, whichInputFile;
	int j, inint, line;
	float infloat;
	rasterdata patchraster, spdistraster, costraster;
	int errors = 0;
	int totlines = 0;
	vector <int> landlist;
	string filetype = "LandFile";

	if (landtype == 0 || landtype == 2) { // real landscape
		// Parse header line;
		ifsLandFile >> header; if (header != "LandNum") errors++;
		ifsLandFile >> header; if (header != "Nhabitats") errors++;
		ifsLandFile >> header; if (header != "LandscapeFile") errors++;
		ifsLandFile >> header; if (header != "PatchFile") errors++;
		ifsLandFile >> header; if (header != "CostMapFile") errors++;
		ifsLandFile >> header; if (header != "DynLandFile") errors++;
		ifsLandFile >> header; if (header != "SpDistFile") errors++;
		if (errors > 0) {
			FormatError(filetype, 0);
			batchLogOfs << "*** Ensure format is correct for real landscape" << endl;
			return -111;
		}
		// Parse data lines
		line = 1;
		inint = -98765;
		ifsLandFile >> inint;
		while (inint != -98765) {
			//		 batchlog << "ParseLandFile(): Landscape no. = " << inint << endl;
			if (inint < 1) {
				BatchError(filetype, line, 11, "LandNum"); errors++;
			}
			else {
				// landscape number must be unique - retain in list to check
				for (j = 0; j < (int)landlist.size(); j++) {
					if (inint == landlist[j]) {
						BatchError(filetype, line, 666, "LandNum"); j = (int)landlist.size() + 1; errors++;
					}
				}
				landlist.push_back(inint);
			}
			ifsLandFile >> inint;
			if (landtype == 0) { // raster map with unique habitat codes
				if (inint < 0) {
					BatchError(filetype, line, 10, "Nhabitats"); errors++;
				}
				if (inint > gMaxNbHab) {
					BatchError(filetype, line, 0, " ");
					batchLogOfs << "Nhabitats may not exceed MaxHabitats in Control file" << endl;
					errors++;
				}
			}

			// check landscape filename
			whichInputFile = "LandscapeFile";
			ifsLandFile >> inputText;
			fileName = indir + inputText;
			landraster = CheckRasterFile(fileName);
			if (landraster.ok) {
				if (landraster.cellsize == gResol)
					batchLogOfs << whichInputFile << " headers OK: " << fileName << endl;
				else {
					errors++;
					batchLogOfs << gResolOfStr << whichInputFile << " " << fileName
						<< gResolNotMatchStr << endl;
				}
			}
			else {
				errors++;
				if (landraster.errors == -111)
					OpenError(whichInputFile, fileName);
				else
					FormatError(fileName, landraster.errors);
			}

			// check patch map filename
			whichInputFile = "PatchFile";
			ifsLandFile >> inputText;
			if (inputText == "NULL") {
				if (gIsPatchModel) {
					BatchError(filetype, line, 0, " "); errors++;
					batchLogOfs << whichInputFile << gPatchReqdStr << endl;
				}
			}
			else {
				if (gIsPatchModel) {
					fileName = indir + inputText;
					patchraster = CheckRasterFile(fileName);
					if (patchraster.ok) {
						if (patchraster.cellsize == gResol) {
							if (patchraster.ncols == landraster.ncols
								&& patchraster.nrows == landraster.nrows
								&& patchraster.cellsize == landraster.cellsize
								&& (int)patchraster.xllcorner == (int)landraster.xllcorner
								&& (int)patchraster.yllcorner == (int)landraster.yllcorner) {
								batchLogOfs << whichInputFile << " headers OK: " << fileName << endl;
							}
							else {
								batchLogOfs << gHeadersOfStr << whichInputFile << " " << fileName
									<< gHeadersNotMatchStr << endl;
								errors++;
							}
						}
						else {
							batchLogOfs << gResolOfStr << whichInputFile << " " << fileName
								<< gResolNotMatchStr << endl;
							errors++;
						}
					}
					else {
						errors++;
						if (patchraster.errors == -111)
							OpenError(whichInputFile, fileName);
						else
							FormatError(fileName, patchraster.errors);
					}
				}
			}

			// check cost map filename
			whichInputFile = "CostMapFile";
			ifsLandFile >> gNameCostFile;
			if (gNameCostFile == "NULL") {
				if (gTransferType == 1) { // SMS
					if (landtype == 2) {
						BatchError(filetype, line, 0, " "); errors++;
						batchLogOfs << whichInputFile << " is required for a habitat quality landscape" << endl;
					}
				}
			}
			else {
				if (gTransferType == 1) { // SMS
					fileName = indir + gNameCostFile;
					costraster = CheckRasterFile(fileName);
					if (costraster.ok) {
						if (costraster.cellsize == gResol) {
							if (costraster.ncols == landraster.ncols
								&& costraster.nrows == landraster.nrows
								&& costraster.cellsize == landraster.cellsize
								&& (int)costraster.xllcorner == (int)landraster.xllcorner
								&& (int)costraster.yllcorner == (int)landraster.yllcorner) {
								batchLogOfs << whichInputFile << " headers OK: " << fileName << endl;
							}
							else {
								batchLogOfs << gHeadersOfStr << whichInputFile << " " << fileName
									<< gHeadersNotMatchStr << endl;
								errors++;
							}
						}
						else {
							batchLogOfs << gResolOfStr << whichInputFile << " " << fileName
								<< gResolNotMatchStr << endl;
							errors++;
						}
					}
					else {
						errors++;
						if (costraster.errors == -111)
							OpenError(whichInputFile, fileName);
						else
							FormatError(fileName, costraster.errors);
					}
				}
				else {
					BatchError(filetype, line, 0, " "); errors++;
					batchLogOfs << whichInputFile << " must be NULL if transfer model is not SMS" << endl;
				}
			}

			// check dynamic landscape filename
			whichInputFile = "DynLandFile";
			ifsLandFile >> inputText;
			if (inputText != "NULL") { // landscape is dynamic
				fileName = indir + inputText;
				batchLogOfs << "Checking " << whichInputFile << " " << fileName << endl;
				ifsDynLandFile.open(fileName.c_str());
				if (ifsDynLandFile.is_open()) {
					int something = CheckDynamicFile(indir, gNameCostFile);
					if (something < 0) {
						errors++;
					}
					ifsDynLandFile.close(); ifsDynLandFile.clear();
				}
				else {
					ifsDynLandFile.clear();
					errors++;
					OpenError(whichInputFile, fileName);
				}
			}

			// check initial distribution map filename
			whichInputFile = "SpDistFile";
			ifsLandFile >> inputText;
			if (inputText == "NULL") {
				if (gUseSpeciesDist) {
					BatchError(filetype, line, 0, " "); 
					errors++;
					batchLogOfs << whichInputFile << " is required as SpeciesDist is 1 in Control file" << endl;
				}
			}
			else {
				if (gUseSpeciesDist) {
					fileName = indir + inputText;
					spdistraster = CheckRasterFile(fileName);
					if (spdistraster.ok) {
						if (spdistraster.cellsize == gDistResol) {
							if (spdistraster.cellsize == landraster.cellsize) {
								// check that extent matches landscape extent
								if (spdistraster.ncols != landraster.ncols
									|| spdistraster.nrows != landraster.nrows) {
									batchLogOfs << "*** Extent of " << whichInputFile
										<< " does not match extent of LandscapeFile" << endl;
									errors++;
								}
								else {
									// check origins match
									if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
										&& (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
										batchLogOfs << whichInputFile << " headers OK: " << fileName << endl;
									}
									else {
										batchLogOfs << "*** Origin co-ordinates of " << whichInputFile
											<< " do not match those of LandscapeFile" << endl;
										errors++;
									}
								}
							}
							else { // not able to check extents match
								// check origins match
								if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
									&& (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
									batchLogOfs << whichInputFile << " headers OK: " << fileName << endl;
								}
								else {
									batchLogOfs << "*** Origin co-ordinates of " << whichInputFile
										<< " do not match those of LandscapeFile" << endl;
									errors++;
								}
							}
						}
						else {
							batchLogOfs << "*** Resolution of " << whichInputFile << " " << fileName
								<< " does not match DistResolution in Control file" << endl;
							errors++;
						}
					}
					else {
						errors++;
						if (spdistraster.errors == -111)
							OpenError(whichInputFile, fileName);
						else
							FormatError(fileName, spdistraster.errors);
					}
				}
			}
			totlines++; line++;
			// read first field on next line
			inint = -98765;
			ifsLandFile >> inint;
		} // end of while loop
		landlist.clear();
	} // end of real landscape
	else {
		if (landtype == 9) { // artificial landscape
			int fractal, type, Xdim, Ydim;
			float minhab, maxhab;
			// Parse header line;
			ifsLandFile >> header; if (header != "LandNum") errors++;
			ifsLandFile >> header; if (header != "Fractal") errors++;
			ifsLandFile >> header; if (header != "Type") errors++;
			ifsLandFile >> header; if (header != "Xdim") errors++;
			ifsLandFile >> header; if (header != "Ydim") errors++;
			ifsLandFile >> header; if (header != "MinHab") errors++;
			ifsLandFile >> header; if (header != "MaxHab") errors++;
			ifsLandFile >> header; if (header != "Psuit") errors++;
			ifsLandFile >> header; if (header != "H") errors++;
			if (errors > 0) {
				FormatError(filetype, 0);
				batchLogOfs << "*** Ensure format is correct for artificial landscape" << endl;
				return -111;
			}
			// Parse data lines
			line = 1;
			inint = -98765;
			ifsLandFile >> inint;
			while (inint != -98765) {
				for (j = 0; j < (int)landlist.size(); j++) {
					if (inint < 1 || inint == landlist[j]) {
						BatchError(filetype, line, 666, "LandNum"); j = (int)landlist.size() + 1; errors++;
					}
				}
				landlist.push_back(inint);
				ifsLandFile >> fractal;
				if (fractal < 0 || fractal > 1) {
					BatchError(filetype, line, 1, "Fractal"); errors++;
				}
				ifsLandFile >> type;
				if (type < 0 || type > 1) {
					BatchError(filetype, line, 1, "Type"); errors++;
				}
				ifsLandFile >> Xdim >> Ydim;
				if (fractal == 1) {
					if (Xdim < 3) {
						BatchError(filetype, line, 13, "Xdim"); errors++;
					}
					if (Ydim < 3) {
						BatchError(filetype, line, 13, "Ydim"); errors++;
					}
				}
				else {
					if (Xdim < 1) {
						BatchError(filetype, line, 11, "Xdim"); errors++;
					}
					if (Ydim < 1) {
						BatchError(filetype, line, 11, "Ydim"); errors++;
					}
				}
				if (fractal == 1) {
					if (Ydim < Xdim) {
						BatchError(filetype, line, 0, " ");
						batchLogOfs << "Y dimension may not be less than X dimension" << endl; errors++;
					}
					if ((Xdim > 2 && !isValidFractalDim(Xdim - 1))
						|| (Ydim > 2 && !isValidFractalDim(Ydim - 1))) {
						BatchError(filetype, line, 0, " ");
						batchLogOfs << "X and Y dimensions must be a power of 2 plus 1" << endl; errors++;
					}
				}
				ifsLandFile >> minhab >> maxhab;
				if (type == 1) { // continuous landscape
					if (minhab <= 0.0 || minhab >= 100.0) {
						BatchError(filetype, line, 100, "MinHab"); errors++;
					}
					if (maxhab <= 0.0 || maxhab > 100.0) {
						BatchError(filetype, line, 100, "MaxHab"); errors++;
					}
					if (maxhab <= minhab) {
						BatchError(filetype, line, 0, " ");
						batchLogOfs << "MaxHab must exceed MinHab" << endl; errors++;
					}
				}
				ifsLandFile >> infloat;
				if (infloat < 0.0 || infloat > 1.0) {
					BatchError(filetype, line, 20, "Psuit"); errors++;
				}
				ifsLandFile >> infloat;
				if (fractal == 1) {
					if (infloat <= 0.0 || infloat >= 1.0) {
						BatchError(filetype, line, 20, "H"); errors++;
					}
				}
				totlines++; line++;
				// read first field on next line
				inint = -98765;
				ifsLandFile >> inint;
			} // end of while loop
		} // end of artificial landscape
		else { // ERROR condition which should not occur
			batchLogOfs << "*** Critical error in land file. "
				<< "Invalid value of landscape type passed to function ParseLandFile()" << endl;
			errors++;
		}
	}
	if (!ifsLandFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return totlines;

}

int CheckDynamicFile(string indir, string costfile) {

	string header, filename, fname, ftype, intext;
	int change, prevchange, year, prevyear = 0;
	rasterdata landchgraster, patchchgraster, costchgraster;
	int errors = 0;
	string filetype = "DynLandFile";
	//int totlines = 0;

	ifsDynLandFile >> header; if (header != "Change") errors++;
	ifsDynLandFile >> header; if (header != "Year") errors++;
	ifsDynLandFile >> header; if (header != "LandChangeFile") errors++;
	ifsDynLandFile >> header; if (header != "PatchChangeFile") errors++;
	ifsDynLandFile >> header; if (header != "CostChangeFile") errors++;

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	change = -98765;
	ifsDynLandFile >> change; // first change number
	if (change != 1) {
		batchLogOfs << "*** Error in DynLandFile - first change number must be 1" << endl;
		errors++;
	}
	else {
		prevchange = change;
	}
	while (change != -98765) {

		ifsDynLandFile >> year; if (year <= 0) { BatchError(filetype, line, 10, "Year"); errors++; }
		if (line > 1) {
			if (year <= prevyear) {
				BatchError(filetype, line, 1, "Year", "previous Year"); errors++;
			}
		}
		prevyear = year;

		// check landscape filename
		ftype = "LandChangeFile";
		ifsDynLandFile >> intext;
		//batchlog << "***** indir=" << indir << " intext=" << intext << endl;
		fname = indir + intext;
		landchgraster = CheckRasterFile(fname);
		if (landchgraster.ok) {
			if (landchgraster.cellsize == gResol)
				if (landchgraster.ncols == landraster.ncols
					&& landchgraster.nrows == landraster.nrows
					&& landchgraster.cellsize == landraster.cellsize
					&& (int)landchgraster.xllcorner == (int)landraster.xllcorner
					&& (int)landchgraster.yllcorner == (int)landraster.yllcorner) {
					batchLogOfs << ftype << " headers OK: " << fname << endl;
				}
				else {
					batchLogOfs << gHeadersOfStr << ftype << " " << fname
						<< gHeadersNotMatchStr << endl;
					errors++;
				}
			else {
				errors++;
				batchLogOfs << gResolOfStr << ftype << " " << fname << gResolNotMatchStr << endl;
			}
		}
		else {
			errors++;
			if (landchgraster.errors == -111)
				OpenError(ftype, fname);
			else
				FormatError(fname, landchgraster.errors);
		}

		// check patch filename
		ftype = "PatchChangeFile";
		ifsDynLandFile >> intext;
		if (intext == "NULL") {
			if (gIsPatchModel) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << ftype << gPatchReqdStr << endl;
			}
		}
		else {
			if (gIsPatchModel) {
				fname = indir + intext;
				patchchgraster = CheckRasterFile(fname);
				if (patchchgraster.ok) {
					if (patchchgraster.cellsize == gResol) {
						if (patchchgraster.ncols == landraster.ncols
							&& patchchgraster.nrows == landraster.nrows
							&& patchchgraster.cellsize == landraster.cellsize
							&& (int)patchchgraster.xllcorner == (int)landraster.xllcorner
							&& (int)patchchgraster.yllcorner == (int)landraster.yllcorner) {
							batchLogOfs << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchLogOfs << gHeadersOfStr << ftype << " " << fname
								<< gHeadersNotMatchStr << endl;
							errors++;
						}
					}
					else {
						batchLogOfs << gResolOfStr << ftype << " " << fname
							<< gResolNotMatchStr << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (patchchgraster.errors == -111)
						OpenError(ftype, fname);
					else
						FormatError(fname, patchchgraster.errors);
				}
			}
		}

		// check costs change filename
		ftype = "CostChangeFile";
		ifsDynLandFile >> intext;
		if (intext == "NULL") {
			if (costfile != "NULL") {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << ftype << " must be supplied " << endl;
			}
		}
		else {
			if (costfile == "NULL") {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << ftype << " must be NULL to match LandFile " << endl;
			}
			else {
				fname = indir + intext;
				costchgraster = CheckRasterFile(fname);
				if (costchgraster.ok) {
					if (costchgraster.cellsize == gResol) {
						if (costchgraster.ncols == landraster.ncols
							&& costchgraster.nrows == landraster.nrows
							&& costchgraster.cellsize == landraster.cellsize
							&& (int)costchgraster.xllcorner == (int)landraster.xllcorner
							&& (int)costchgraster.yllcorner == (int)landraster.yllcorner) {
							batchLogOfs << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchLogOfs << gHeadersOfStr << ftype << " " << fname
								<< gHeadersNotMatchStr << endl;
							errors++;
						}
					}
					else {
						batchLogOfs << gResolOfStr << ftype << " " << fname
							<< gResolNotMatchStr << endl;
						errors++;
					}
				}
				else {
					errors++;
					if (costchgraster.errors == -111)
						OpenError(ftype, fname);
					else
						FormatError(fname, costchgraster.errors);
				}
			}
		}

		line++;
		// read first field on next line
		change = -98765;
		ifsDynLandFile >> change;
		if (ifsDynLandFile.eof()) {
			change = -98765;
		}
		else { // check for valid change number
			if (change != prevchange + 1) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "Change numbers must be sequential integers" << endl;
				errors++;
			}
			prevchange = change;
		}
	}

	if (errors > 0) return -111;
	else return 0;

}

//---------------------------------------------------------------------------
int CheckStageFile(string indir)
{
	string header, filename, fname, ftype2;
	int inint, i, err, fecdensdep, fecstagewts, devdensdep, devstagewts, survdensdep, survstagewts;
	float infloat;
	int errors = 0;
	int simuls = 0;
	int prevsimul;
	bool checkfile;
	vector <string> transfiles, wtsfiles;
	string filetype = "StageStructFile";

	// Parse header line;
	ifsStageStructFile >> header; if (header != "Simulation") errors++;
	ifsStageStructFile >> header; if (header != "PostDestructn") errors++;
	ifsStageStructFile >> header; if (header != "PRep") errors++;
	ifsStageStructFile >> header; if (header != "RepInterval") errors++;
	ifsStageStructFile >> header; if (header != "MaxAge") errors++;
	ifsStageStructFile >> header; if (header != "TransMatrixFile") errors++;
	ifsStageStructFile >> header; if (header != "SurvSched") errors++;
	ifsStageStructFile >> header; if (header != "FecDensDep") errors++;
	ifsStageStructFile >> header; if (header != "FecStageWts") errors++;
	ifsStageStructFile >> header; if (header != "FecStageWtsFile") errors++;
	ifsStageStructFile >> header; if (header != "DevDensDep") errors++;
	ifsStageStructFile >> header; if (header != "DevDensCoeff") errors++;
	ifsStageStructFile >> header; if (header != "DevStageWts") errors++;
	ifsStageStructFile >> header; if (header != "DevStageWtsFile") errors++;
	ifsStageStructFile >> header; if (header != "SurvDensDep") errors++;
	ifsStageStructFile >> header; if (header != "SurvDensCoeff") errors++;
	ifsStageStructFile >> header; if (header != "SurvStageWts") errors++;
	ifsStageStructFile >> header; if (header != "SurvStageWtsFile") errors++;
	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	inint = -98765;
	ifsStageStructFile >> inint;
	// first simulation number must match first one in parameterFile
	if (inint != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	prevsimul = inint;
	while (inint != -98765) {

		if (!gSimNbs.contains(inint)) {
			BatchError(filetype, line, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			errors++;
		}

		simuls++;
		ifsStageStructFile >> inint;
		if (inint < 0 || inint > 1) { BatchError(filetype, line, 1, "PostDestructn"); errors++; }
		ifsStageStructFile >> infloat;
		if (infloat <= 0 || infloat > 1.0) { BatchError(filetype, line, 20, "PRep"); errors++; }
		ifsStageStructFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "RepInterval"); errors++; }
		ifsStageStructFile >> inint;
		if (inint < 2) { BatchError(filetype, line, 12, "MaxAge"); errors++; }

		ifsStageStructFile >> filename;
		// transition matrix file - compulsory
		ftype2 = "TransMatrixFile";
		checkfile = true;
		for (i = 0; i < (int)transfiles.size(); i++) {
			if (filename == transfiles[i]) { // file has already been checked
				checkfile = false;
			}
		}
		if (checkfile) {
			if (filename == "NULL") {
				batchLogOfs << "*** " << ftype2 << " is compulsory for stage-structured model" << endl;
				errors++;
			}
			else {
				fname = indir + filename;
				batchLogOfs << "Checking " << ftype2 << " " << fname << endl;
				ifsTransMatrix.open(fname.c_str());
				if (ifsTransMatrix.is_open()) {
					err = CheckTransitionFile(gNbStages, gNbSexesDem);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					ifsTransMatrix.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (ifsTransMatrix.is_open()) ifsTransMatrix.close();
				ifsTransMatrix.clear();
			}
		}
		transfiles.push_back(filename);

		ifsStageStructFile >> inint;
		if (inint < 0 || inint > 2) { BatchError(filetype, line, 2, "SurvSched"); errors++; }
		ifsStageStructFile >> fecdensdep;
		if (fecdensdep < 0 || fecdensdep > 1)
		{
			BatchError(filetype, line, 1, "FecDensDep"); errors++; fecdensdep = 1;
		}
		ifsStageStructFile >> fecstagewts;
		if (fecdensdep) {
			if (fecstagewts < 0 || fecstagewts > 1)
			{
				BatchError(filetype, line, 1, "FecStageWts"); errors++; fecstagewts = 1;
			}
		}
		else {
			if (fecstagewts != 0) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "FecStageWts must be 0 if FecDensDep is 0" << endl; errors++;
				errors++; fecstagewts = 1;
			}
		}

		// fecundity stage weights file - optional
		ftype2 = "FecStageWtsFile";
		ifsStageStructFile >> filename;
		if (filename == "NULL") {
			if (fecstagewts) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << ftype2 << " is compulsory unless FecStageWts is 0" << endl;
				errors++;
			}
		}
		else {
			checkfile = true;
			for (i = 0; i < (int)wtsfiles.size(); i++) {
				if (filename == wtsfiles[i]) checkfile = false; // file has already been checked
			}
			if (checkfile) {
				fname = indir + filename;
				batchLogOfs << "Checking " << ftype2 << " " << fname << endl;
				ifsStageWeightsFile.open(fname.c_str());
				if (ifsStageWeightsFile.is_open()) {
					err = CheckWeightsFile(ftype2);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					ifsStageWeightsFile.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (ifsStageWeightsFile.is_open()) ifsStageWeightsFile.close();
				ifsStageWeightsFile.clear();
			}
			wtsfiles.push_back(filename);
		}

		ifsStageStructFile >> devdensdep;
		if (devdensdep < 0 || devdensdep > 1)
		{
			BatchError(filetype, line, 1, "DevDensDep"); errors++; devdensdep = 1;
		}
		ifsStageStructFile >> infloat >> devstagewts;
		if (devdensdep) {
			if (infloat <= 0.0) {
				BatchError(filetype, line, 10, "DevDensCoeff"); errors++;
			}
			if (devstagewts < 0 || devstagewts > 1) {
				BatchError(filetype, line, 1, "DevStageWts"); errors++; devstagewts = 1;
			}
		}
		else {
			if (devstagewts != 0) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "DevStageWts must be 0 if DevDensDep is 0" << endl; errors++;
				errors++; devstagewts = 1;
			}
		}

		// development stage weights file - optional
		ftype2 = "DevStageWtsFile";
		ifsStageStructFile >> filename;
		if (filename == "NULL") {
			if (devstagewts) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << ftype2 << " is compulsory unless DevStageWts is 0" << endl;
				errors++;
			}
		}
		else {
			checkfile = true;
			for (i = 0; i < (int)wtsfiles.size(); i++) {
				if (filename == wtsfiles[i]) checkfile = false; // file has already been checked
			}
			if (checkfile) {
				fname = indir + filename;
				batchLogOfs << "Checking " << ftype2 << " " << fname << endl;
				ifsStageWeightsFile.open(fname.c_str());
				if (ifsStageWeightsFile.is_open()) {
					err = CheckWeightsFile(ftype2);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					ifsStageWeightsFile.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (ifsStageWeightsFile.is_open()) ifsStageWeightsFile.close();
				ifsStageWeightsFile.clear();
			}
			wtsfiles.push_back(filename);
		}

		ifsStageStructFile >> survdensdep;
		if (survdensdep < 0 || survdensdep > 1)
		{
			BatchError(filetype, line, 1, "SurvDensDep"); errors++; survdensdep = 1;
		}
		ifsStageStructFile >> infloat >> survstagewts;
		if (survdensdep) {
			if (infloat <= 0.0) {
				BatchError(filetype, line, 10, "SurvDensCoeff"); errors++;
			}
			if (survstagewts < 0 || survstagewts > 1) {
				BatchError(filetype, line, 1, "SurvStageWts"); errors++; survstagewts = 1;
			}
		}
		else {
			if (survstagewts != 0) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "SurvStageWts must be 0 if SurvDensDep is 0" << endl; errors++;
				errors++; survstagewts = 1;
			}
		}

		// survival stage weights file - optional
		ftype2 = "SurvStageWtsFile";
		ifsStageStructFile >> filename;
		if (filename == "NULL") {
			if (survstagewts) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << ftype2 << " is compulsory unless SurvStageWts is 0" << endl;
				errors++;
			}
		}
		else {
			checkfile = true;
			for (i = 0; i < (int)wtsfiles.size(); i++) {
				if (filename == wtsfiles[i]) checkfile = false; // file has already been checked
			}
			if (checkfile) {
				fname = indir + filename;
				batchLogOfs << "Checking " << ftype2 << " " << fname << endl;
				ifsStageWeightsFile.open(fname.c_str());
				if (ifsStageWeightsFile.is_open()) {
					err = CheckWeightsFile(ftype2);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					ifsStageWeightsFile.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (ifsStageWeightsFile.is_open()) ifsStageWeightsFile.close();
				ifsStageWeightsFile.clear();
			}
			wtsfiles.push_back(filename);
		}

		// read next simulation
		line++;
		inint = -98765;
		ifsStageStructFile >> inint;
		if (ifsStageStructFile.eof()) {
			inint = -98765;
		}
		else { // check for valid simulation number
			if (inint != prevsimul + 1) {
				BatchError(filetype, line, 222, " ");
				errors++;
			}
			prevsimul = inint;
		}
	}
	if (!ifsStageStructFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	transfiles.clear();
	wtsfiles.clear();

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
// Check transition matrix file
int CheckTransitionFile(short nstages, short nsexesDem)
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
	//prevminage = minage;
	//				batchlog << "MINAGE = " << minage << endl;
	if (minage != 0) {
		BatchError(filetype, line, 0, " ");
		batchLogOfs << "MinAge must be zero for juvenile stage" << endl; errors++;
	}

	// one row for each stage/sex combination
	//				batchlog << "HEADER = " << header << endl;
	for (stage = 1; stage < nstages; stage++) {
		for (sex = 0; sex < nsexesDem; sex++) {
			line++;
			// row header
			ifsTransMatrix >> header;
			if (nsexesDem == 1) expectedHeader = to_string(stage);
			else {
				if (sex == 0) expectedHeader = to_string(stage) + "m"; else expectedHeader = to_string(stage) + "f";
			}
			if (header != expectedHeader) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "Invalid row header" << endl; errors++;
			}
			for (iStage = 0; iStage < nstages; iStage++) {
				for (iSex = 0; iSex < nsexesDem; iSex++) {
					ifsTransMatrix >> infloat;
					//				batchlog << "TRANS PROB = " << infloat << endl;
					if (infloat < 0.0 || infloat > 1) {
						BatchError(filetype, line, 20, "Transition probability"); errors++;
					}
				}
			}
			//		 prevminage = minage;
			ifsTransMatrix >> minage;
			//				batchlog << "MINAGE = " << minage << endl;
			if (stage == 1 && minage != 0) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "MinAge must be zero for stage 1" << endl; errors++;
			}
			if (stage > 1) {
				if (minage < 0) {
					BatchError(filetype, line, 19, "MinAge"); errors++;
				}
				// SCFP 30/9/13 - IDEALLY OUGHT TO TEST THAT MINAGE IS NO LESS THAN PREVIOUS MINAGE
				// BUT WOULD NEED TO BE PREV MINAGE FOR THE SAME SEX
				// HOWEVER, IT IS NOT CRITICAL, AS A MINAGE OF LESS THAN PREVIOUS CANNOT CAUSE ANY
				// PROBLEM, AS PREVIOUS MINAGE GETS APPLIED EARLIER
	//			if (minage < prevminage) {
	//				BatchError(filetype,line,0," ");
	//				batchlog << "MinAge may not be less than MinAge of previous stage" << endl; errors++;
	//			}
			}
		}
	}
	// final read should hit EOF
	ifsTransMatrix >> header;

	if (!ifsTransMatrix.eof()) {
		EOFerror(filetype);
		errors++;
	}

	return errors;

}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Check stage weights matrix file
int CheckWeightsFile(string filetype)
{
	string header, hhh;
	int i, j, stage, sex, line;
	float infloat;
	int errors = 0;

	// check header records
	ifsStageWeightsFile >> header; if (header != "StageWts") errors++;
	for (i = 0; i < gNbStages; i++) {
		for (j = 0; j < gNbSexesDem; j++) {
			ifsStageWeightsFile >> header;
			if (gNbSexesDem == 1) hhh = to_string(i);
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
	for (stage = 0; stage < gNbStages; stage++) {
		for (sex = 0; sex < gNbSexesDem; sex++) {
			line++;
			// row header
			ifsStageWeightsFile >> header;
			if (gNbSexesDem == 1) hhh = to_string(stage);
			else {
				if (sex == 0) hhh = to_string(stage) + "m"; else hhh = to_string(stage) + "f";
			}
			if (header != hhh) {
				BatchError(filetype, line, 0, " ");
				batchLogOfs << "Invalid row header" << endl; errors++;
			}
			for (i = 0; i < gNbStages; i++) {
				for (j = 0; j < gNbSexesDem; j++) {
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

	return errors;

}

//---------------------------------------------------------------------------
int CheckEmigFile()
{
	string header;
	int simNb;
	int inDensDep, inUseFullKern, inStgDep, inSexDep, inIndVar, inEmigStg, inStage, inSex;
	bool isDensDep, isIndVar;
	float inEP, inD0, inAlpha, inBeta;
	int nbErrors = 0;
	int nbSims = 0;
	string whichInputFile = "EmigrationFile";

	isDensDep = false;
	isIndVar = false;
	inEP = 0.0;

	// Parse header line;
	ifsEmigrationFile >> header; if (header != "Simulation") nbErrors++;
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
	simNb = gFirstSimNb + 1; // that is, NOT first sim number
	prevLine.simNb = -999;
	prevLine.simLines = prevLine.reqdSimLines = 0;
	ifsEmigrationFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichInputFile, lineNb, 111, "Simulation"); 
		nbErrors++;
		readNextLine = false;
	}

	while (readNextLine) {

		if (!gSimNbs.contains(simNb)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}

		// read and validate columns relating to stage and sex-dependency and to IIV
		ifsEmigrationFile >> inDensDep >> inUseFullKern >> inStgDep >> inSexDep;
		ifsEmigrationFile >> inIndVar >> inEmigStg >> inStage >> inSex;
		currentLine = CheckStageSex(whichInputFile, lineNb, simNb, prevLine, inStgDep, inSexDep, inStage, inSex, inIndVar, true, false);
		if (currentLine.isNewSim) nbSims++;
		nbErrors += currentLine.errors;
		prevLine = currentLine;

		// validate density dependency
		if (inDensDep != 0 && inDensDep != 1) {
			BatchError(whichInputFile, lineNb, 1, "DensDep"); 
			nbErrors++;
		}
		else {
			gTraitOptions.at(simNb).isEmigDensDep = (inDensDep == 1);
		}

		// validate individual variation
		if (inIndVar != 0 && inIndVar != 1) {
			BatchError(whichInputFile, lineNb, 1, "IndVar");
			nbErrors++;
		}
		else {
			gTraitOptions.at(simNb).isEmigIndVar = (inIndVar == 1);
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
		if (gStageStruct && !inStgDep && inIndVar == 1
			&& inStage == 0 && inSex == 0
			&& (inEmigStg < 0 || inEmigStg >= gNbStages)) {
			BatchError(whichInputFile, lineNb, 0, "EmigStage");
			nbErrors++;
			batchLogOfs << "EmigStage must be from 0 to " << to_string(gNbStages - 1) << endl;
		}
		if (inSexDep != 0 && inSexDep != 1) {
			BatchError(whichInputFile, lineNb, 1, "SexDep");
			nbErrors++;
		} 
		else {
			gTraitOptions.at(simNb).isEmigSexDep = (inSexDep == 1);
		}

		if (inStage == 0 && inSex == 0) { // first line of a simulation
			// record whether density dependence and individual variability are applied
			isDensDep = (inDensDep == 1);
			isIndVar = (inIndVar == 1);
		}

		// read remaining columns of the current record
		ifsEmigrationFile >> inEP >> inD0 >> inAlpha >> inBeta;

		if (gTraitOptions.at(simNb).isEmigIndVar) {
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
		if (simNb == errSimNb || ifsEmigrationFile.eof()) readNextLine = false;
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
	int i, simNb, inStageDep, inSexDep, inKernelType, inDistMort, inIndVar, inStage, inSex;
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
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichFile, whichLine, 111, "Simulation"); 
		errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		if (!gSimNbs.contains(simNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			errors++;
		}

		switch (gTransferType) {

		case 0: { // negative exponential dispersal kernel
			// read and validate columns relating to stage and sex-dependency and to IIV
			ifsTransferFile >> inStageDep >> inSexDep >> inKernelType >> inDistMort;
			ifsTransferFile >> inIndVar >> inStage >> inSex;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, inStageDep, inSexDep, inStage, inSex, inIndVar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			// validate kernel type
			if (inKernelType != 0 && inKernelType != 1) {
				BatchError(whichFile, whichLine, 1, "KernelType"); errors++;
			}
			else {
				gTraitOptions.at(simNb).usesTwoKernels = (inKernelType == 1);
			}
			// validate mortality
			if (inDistMort != 0 && inDistMort != 1) {
				BatchError(whichFile, whichLine, 1, "DistMort"); errors++;
			}
			// read remaining columns of the current record
			ifsTransferFile >> meanDistI >> meanDistII >> ProbKernelI;
			ifsTransferFile >> mortProb >> slope >> inflPoint;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				gTraitOptions.at(simNb).isKernTransfIndVar = (inIndVar == 1);
			}

			if (inSexDep != 0 && inSexDep != 1) {
				BatchError(whichFile, whichLine, 1, "SexDep"); errors++;
			}
			else {
				gTraitOptions.at(simNb).isKernTransfSexDep = (inSexDep == 1);
			}

			// validate mortality
			if (inDistMort != 0 && inDistMort != 1) {
				BatchError(whichFile, whichLine, 1, "DistMort"); errors++;
			}

			if (gTraitOptions.at(simNb).isKernTransfIndVar) {
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
			current = CheckStageSex(whichFile, whichLine, simNb, prev, 0, 0, 0, 0, 0, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				gTraitOptions.at(simNb).isSMSTransfIndVar = (inIndVar == 1);
			}

			// validate SMS movement parameters
			if (inPerceptualRange < 1) {
				BatchError(whichFile, whichLine, 11, "PR"); errors++;
			}
			if (inPercRangeMethod < 1 || inPercRangeMethod > 3) {
				BatchError(whichFile, whichLine, 33, "PRmethod"); errors++;
			}
			if (gTraitOptions.at(simNb).isSMSTransfIndVar) {
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
				gTraitOptions.at(simNb).usesSMSGoalBias = (inGoalType == 2);
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
			else {
				if (inSMType != 0 && inSMType != 1) {
					BatchError(whichFile, whichLine, 1, "SMtype"); errors++;
				}
			}
			if (inSMType == 0)
			{
				if (inSMConst < 0.0 || inSMConst >= 1.0) {
					BatchError(whichFile, whichLine, 20, "SMconst"); errors++;
				}
			}
			switch (gLandType) {

			case 0: { // raster map with unique habitat codes
				for (i = 0; i < gMaxNbHab; i++) {
					ifsTransferFile >> morthab;
					if (inSMType == 1)
					{
						if (morthab < 0.0 || morthab >= 1.0) {
							colheader = "MortHab" + to_string(i + 1);
							BatchError(whichFile, whichLine, 20, colheader); errors++;
						}
					}
				}
				for (i = 0; i < gMaxNbHab; i++) {
					ifsTransferFile >> costhab;
					if (gNameCostFile == "NULL") {
						if (costhab < 1) {
							colheader = "CostHab" + to_string(i + 1);
							BatchError(whichFile, whichLine, 11, colheader); errors++;
						}
					}
				}
				break;
			} // end of raster map with unique habitat codes

			case 2: { // raster map with habitat quality

				break;
			} // end of raster map with habitat quality

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
			current = CheckStageSex(whichFile, whichLine, simNb, prev, 0, 0, 0, 0, inIndVar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				gTraitOptions.at(simNb).isCRWTransfIndVar = (inIndVar == 1);
			}

			if (gTraitOptions.at(simNb).isCRWTransfIndVar) {
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
			if (gLandType == 0) { // real landscape with habitat types
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
			else { // real landscape with quality OR artificial landscape
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
		BatchError(whichFile, whichLine, 0, " "); errors++;
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
	int simNb, inStageDep, inSexDep, inStage, inSex, inSettleType;
	int inDensDep, inIndVar, inFindMate, inMinSteps, inMaxSteps, inMaxStepsYear;
	float inS0, inAlphaS, inBetaS;
	int nbErrors = 0;
	int nbSims = 0;
	string whichFile = "SettlementFile";

	// Parse header line;
	ifsSettlementFile >> header; if (header != "Simulation") nbErrors++;
	ifsSettlementFile >> header; if (header != "StageDep") nbErrors++;
	ifsSettlementFile >> header; if (header != "SexDep") nbErrors++;
	ifsSettlementFile >> header; if (header != "Stage") nbErrors++;
	ifsSettlementFile >> header; if (header != "Sex") nbErrors++;
	if (gTransferType == 0)
	{ // dispersal kernel
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
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichFile, whichLine, 111, "Simulation"); nbErrors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		if (!gSimNbs.contains(simNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}

		if (gTransferType == 0)
		{ // dispersal kernel
			// read and validate columns relating to stage and sex-dependency (NB no IIV here)
			ifsSettlementFile >> inStageDep >> inSexDep >> inStage >> inSex >> inSettleType >> inFindMate;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, inStageDep, inSexDep, inStage, inSex, 0, true, false);
			if (current.isNewSim) nbSims++;
			nbErrors += current.errors;
			prev = current;
			if (inSettleType < 0 || inSettleType > 3) {
				BatchError(whichFile, whichLine, 3, "SettleType"); nbErrors++;
			}
			if (!gStageStruct && (inSettleType == 1 || inSettleType == 3)) {
				BatchError(whichFile, whichLine, 0, " "); nbErrors++;
				batchLogOfs << "Invalid SettleType for a non-stage-structured population" << endl;
			}
			if (gNbSexesDisp > 1) {
				if (inFindMate < 0 || inFindMate > 1) {
					BatchError(whichFile, whichLine, 1, "FindMate"); nbErrors++;
				}
			}
		}
		else { // movement method
			// read and validate columns relating to stage and sex-dependency (IIV psossible)
			ifsSettlementFile >> inStageDep >> inSexDep >> inStage >> inSex >> inDensDep >> inIndVar >> inFindMate;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, inStageDep, inSexDep, inStage, inSex, inIndVar, true, false);
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
				gTraitOptions.at(simNb).isSettIndVar = inIndVar == 1;
			}

			if (inSexDep != 0 && inSexDep != 1) {
				BatchError(whichFile, whichLine, 1, "SexDep");
				nbErrors++;
			}
			else {
				gTraitOptions.at(simNb).isSettSexDep = inSexDep == 1;
			}

			if (gReproType != 0 && gNbSexesDisp > 1) {
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

			if (gTraitOptions.at(simNb).isSettIndVar) {
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
		if (ifsSettlementFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation

	if (current.simLines != current.reqdSimLines) {
		BatchError(whichFile, whichLine, 0, " "); nbErrors++;
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
	int simNb, nextLineSimNb;
	string filename, inTraitType, inSex, inInitDist, inInitParams, inInitDomDist, inInitDomParams,
		inDominanceDist, inDominanceParams, inIsInherited, inMutationDist, 
		inMutationParams, inPositions, inNbPositions, inExpressionType, inMutationRate, inIsOutput;
	int nbErrors = 0;
	int nbSims = 0;
	int nbGenLoadTraits = 0;
	vector <string> archfiles;
	const string whichInputFile = "TraitsFile";
	vector <TraitType> allReadTraits;

	// Parse header line
	ifsTraitsFile >> header; if (header != "Simulation") nbErrors++;
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

	ifsTraitsFile >> simNb;

	bool stopReading = (simNb == simNbNotRead);
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichInputFile, lineNb, 111, "Simulation"); 
		nbErrors++;
	}
	int nbRowsToRead = 0;

	while (!stopReading) {

		if (!gSimNbs.contains(simNb)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}

		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		ifsTraitsFile >> inTraitType >> inSex >> inPositions >> inNbPositions 
			>> inExpressionType >> inInitDist >> inInitParams >> inInitDomDist >> inInitDomParams
			>> inIsInherited >> inMutationDist >> inMutationParams >> inDominanceDist >> inDominanceParams
			>> inMutationRate >> inIsOutput;

		current = CheckStageSex(whichInputFile, lineNb, simNb, prev, 0, 0, 0, 0, 0, true, false);
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
		ifsTraitsFile >> nextLineSimNb;
		if (nextLineSimNb == simNbNotRead
			// Exit loop
			|| ifsTraitsFile.eof()) {
			stopReading = true;
			nbErrors += checkTraitSetCoherency(allReadTraits, simNb);
			gNbTraitFileRows.push_back(nbRowsToRead);
		}
		else if (nextLineSimNb != simNb) {
			// About to change sim, conduct checks of all read traits
			nbErrors += checkTraitSetCoherency(allReadTraits, simNb);
			// Store nb of rows to help reading file later on
			gNbTraitFileRows.push_back(nbRowsToRead);
			nbRowsToRead = 0; // reset for next sim
			nbGenLoadTraits = 0;
			allReadTraits.clear();
			simNb = nextLineSimNb;
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

int checkTraitSetCoherency(const vector <TraitType>& allReadTraits, const int& simNb) {
	int nbErrors = 0;
	const string whichInputFile = "TraitsFile";

	if (gTraitOptions.at(simNb).anyNeutral && !traitExists(NEUTRAL, allReadTraits)) {
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

	if (gTraitOptions.at(simNb).isEmigIndVar) {
		if (!hasD0) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "EP or d0 is missing." << endl;
			nbErrors++;
		}
		if (gTraitOptions.at(simNb).isEmigSexDep) {
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

		if (gTraitOptions.at(simNb).isEmigDensDep) {
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
			if (gTraitOptions.at(simNb).isEmigSexDep) {
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

	if (gTraitOptions.at(simNb).isKernTransfIndVar) {
		if (!hasKern1) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "(First) kernel mean is missing." << endl;
			nbErrors++;
		}
		if (gTraitOptions.at(simNb).isKernTransfSexDep) {
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
		if (gTraitOptions.at(simNb).usesTwoKernels) {
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
			if (gTraitOptions.at(simNb).isKernTransfSexDep) {
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
	if (gTraitOptions.at(simNb).isSMSTransfIndVar) {
		if (!hasDP) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "SMS directional persistence trait is missing." << endl;
			nbErrors++;
		}
		if (gTraitOptions.at(simNb).usesSMSGoalBias) {
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
	if (gTraitOptions.at(simNb).isCRWTransfIndVar) {
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

	if (gTraitOptions.at(simNb).isSettIndVar) {
		if (!hasS0) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLogOfs << "Settlement probability trait is missing." << endl;
			nbErrors++;
		}
		if (gTraitOptions.at(simNb).isSettSexDep) {
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
		if (gTraitOptions.at(simNb).isSettSexDep) {
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
	int simNb, prevSimNb, errCode;
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
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichFile, whichLine, 111, "Simulation"); 
		nbErrors++;
	}
	while (simNb != -98765) {

		if (!gSimNbs.contains(simNb)) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			nbErrors++;
		}

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
		if (gNbSexesDisp == 1 && inRecombinationRate != "#") {
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
		gTraitOptions.at(simNb).anyNeutral = inOutWeirCockerham == "TRUE"
			|| inOutWeirHill == "TRUE";
		bool anyGeneticsOutput = inOutGeneValues == "TRUE" 
			|| gTraitOptions.at(simNb).anyNeutral;

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
	int i, simNb;
	int seedtype, freetype, sptype, initdens, indscell = 0, minX, maxX, minY, maxY;
	int nCells, nSpCells, initAge;
	int initFreezeYear, restrictRows, restrictFreq, finalFreezeYear;
	float inds_per_ha;

	int errors = 0; int propnerrors = 0;
	int simuls = 0;
	string filetype = "InitialisationFile";

	// Parse header line;
	ifsInitFile >> header; if (header != "Simulation") errors++;
	ifsInitFile >> header; if (header != "SeedType") errors++;
	ifsInitFile >> header; if (header != "FreeType") errors++;
	ifsInitFile >> header; if (header != "SpType") errors++;
	ifsInitFile >> header; if (header != "InitDens") errors++;
	ifsInitFile >> header;
	if (gIsPatchModel) { if (header != "IndsHa") errors++; }
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
	if (gStageStruct) {
		ifsInitFile >> header; if (header != "InitAge") errors++;
		for (i = 1; i < gNbStages; i++) {
			colheader = "PropStage" + to_string(i);
			ifsInitFile >> header; if (header != colheader) propnerrors++;
		}
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
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		if (!gSimNbs.contains(simNb)) {
			BatchError(filetype, line, 0, " ");
			batchLogOfs << "Simulation number doesn't match those in ParametersFile" << endl;
			errors++;
		}

		current = CheckStageSex(filetype, line, simNb, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) simuls++;
		errors += current.errors;
		prev = current;

		ifsInitFile >> seedtype >> freetype >> sptype >> initdens >> inds_per_ha;
		if (!gIsPatchModel) indscell = (int)inds_per_ha;
		if (seedtype < 0 || seedtype > 2) {
			BatchError(filetype, line, 2, "SeedType"); errors++;
		}
		if (gLandType == 9 && seedtype != 0) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLogOfs << "SeedType must be 0 for an artificial landscape"
				<< endl;
		}
		if (!gUseSpeciesDist && seedtype == 1) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLogOfs << "SeedType may not be 1 if there is no initial species distribution map"
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
				if (gIsPatchModel) {
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
						err = CheckInitIndsFile();
						if (err == 0) FileHeadersOK(ftype2); else errors++;
						ifsInitIndsFile.close();
					}
					else {
						OpenError(ftype2, fname); errors++;
					}
					if (ifsInitIndsFile.is_open()) ifsInitIndsFile.close();
					ifsInitIndsFile.clear();
					indsfiles.push_back(filename);
				}
			}
			else {
				BatchError(filetype, line, 0, " "); errors++;
				batchLogOfs << ftype2 << " must be NULL for SeedType "
					<< seedtype << endl;
			}
		}

		if (gStageStruct) {
			ifsInitFile >> initAge;
			if (seedtype != 2 && (initAge < 0 || initAge > 2)) {
				BatchError(filetype, line, 2, "initAge"); errors++;
			}
			float propstage;
			float cumprop = 0.0;
			for (i = 1; i < gNbStages; i++) {
				ifsInitFile >> propstage;
				cumprop += propstage;
				if (seedtype != 2 && (propstage < 0.0 || propstage > 1.0)) {
					colheader = "PropStage" + to_string(i);
					BatchError(filetype, line, 20, colheader); errors++;
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
int CheckInitIndsFile() {
	string header;
	int year, species, patchID, x, y, ninds, sex, age, stage, prevyear;

	int errors = 0;
	string filetype = "InitIndsFile";

	// Parse header line
	ifsInitIndsFile >> header; 
	if (header != "Year") errors++;
	ifsInitIndsFile >> header; 
	if (header != "Species") errors++;
	if (gIsPatchModel) {
		ifsInitIndsFile >> header; 
		if (header != "PatchID") errors++;
	}
	else {
		ifsInitIndsFile >> header; 
		if (header != "X") errors++;
		ifsInitIndsFile >> header;
		if (header != "Y") errors++;
	}
	ifsInitIndsFile >> header; 
	if (header != "Ninds") errors++;
	if (gReproType > 0) {
		ifsInitIndsFile >> header; 
		if (header != "Sex") errors++;
	}
	if (gStageStruct) {
		ifsInitIndsFile >> header; 
		if (header != "Age") errors++;
		ifsInitIndsFile >> header; 
		if (header != "Stage") errors++;
	}

	// Report any errors in headers, and if so, terminate validation
	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	string filename, ftype2, fname;
	year = prevyear = -98765;
	ifsInitIndsFile >> year;
	while (year != -98765) {
		if (year < 0) {
			BatchError(filetype, line, 19, "Year"); errors++;
		}
		else {
			if (year < prevyear) {
				BatchError(filetype, line, 2, "Year", "previous Year"); errors++;
			}
		}
		prevyear = year;
		ifsInitIndsFile >> species;
		if (species != 0) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLogOfs << "Species must be 0" << endl;
		}
		if (gIsPatchModel) {
			ifsInitIndsFile >> patchID;
			if (patchID < 1) {
				BatchError(filetype, line, 11, "PatchID"); errors++;
			}
		}
		else {
			ifsInitIndsFile >> x >> y;
			if (x < 0 || y < 0) {
				BatchError(filetype, line, 19, "X and Y"); errors++;
			}
		}
		ifsInitIndsFile >> ninds;
		if (ninds < 1) {
			BatchError(filetype, line, 11, "Ninds"); errors++;
		}
		if (gReproType > 0) {
			ifsInitIndsFile >> sex;
			if (sex < 0 || sex > 1) {
				BatchError(filetype, line, 1, "Sex"); errors++;
			}
		}
		if (gStageStruct) {
			ifsInitIndsFile >> age >> stage;
			if (age < 1) {
				BatchError(filetype, line, 11, "Age"); errors++;
			}
			if (stage < 1) {
				BatchError(filetype, line, 11, "Stage"); errors++;
			}
			if (stage >= gNbStages) {
				BatchError(filetype, line, 4, "Stage", "no. of stages"); errors++;
			}
		}
		line++;
		year = -98765;
		ifsInitIndsFile >> year;
		if (ifsInitIndsFile.eof()) year = -98765;
	} // end of while loop
	if (!ifsInitIndsFile.eof()) {
		EOFerror(filetype);
		errors++;
	}
	return errors;
}

//---------------------------------------------------------------------------
/*
Check stage- and sex-dependency fields for any of the dispersal files.
Check that the number of records for a simulation matches the stage-
and sex-dependency settings (unless checklines is false).
Validate the IIV field (if present).
*/
simCheck CheckStageSex(string whichInputFile, int whichLine, int simNb, simCheck prev,
	int isStageDep, int isSexDep, int stage, int sex, int isIndVar,
	bool mustCheckLines, bool mustCheckStgDepWithIndVar)
{
	simCheck current;
	current.errors = 0;
	int expectedStage;

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
			BatchError(whichInputFile, whichLine, 222, " "); current.errors++;
		}
		// check for correct number of lines for previous simulation
		if (mustCheckLines && !(prev.simLines >= prev.reqdSimLines)) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "No. of lines for previous Simulation " << prev.simNb
				<< gShouldBeStr << prev.reqdSimLines << endl;
		}
	}
	current.simNb = simNb;

	// validate inStageDep
	if (gStageStruct) {
		if (isStageDep != 0 && isStageDep != 1) {
			BatchError(whichInputFile, whichLine, 1, "StageDep"); current.errors++;
			isStageDep = 1; // to calculate required number of lines
		}
	}
	else {
		if (isStageDep != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "StageDep must be 0 for non-stage-structured model" << endl;
			isStageDep = 0; // to calculate required number of lines
		}
	}
	// validate inSexDep
	if (gNbSexesDisp == 2) {
		if (isSexDep != 0 && isSexDep != 1) {
			BatchError(whichInputFile, whichLine, 1, "SexDep"); current.errors++;
			isSexDep = 1; // to calculate required number of lines
		}
	}
	else {
		if (isSexDep != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "SexDep must be 0 for asexual model" << endl;
			isSexDep = 0; // to calculate required number of lines
		}
	}
	if (current.isNewSim) { // set required number of lines
		if (isStageDep) {
			current.reqdSimLines = gNbStages;
			if (isSexDep) current.reqdSimLines *= gNbSexesDisp;
		}
		else {
			current.reqdSimLines = isSexDep ? gNbSexesDisp : 1;
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
		else { // there must be 1 line for each stage
			if (stage != current.simLines - 1) {
				BatchError(whichInputFile, whichLine, 0, " "); 
				current.errors++;
				batchLogOfs << "Stages must be sequentially numbered from 0" << endl;
			}
		}
	}
	else { // no stage-dependent emigration
		if (stage != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "Stage must be 0 for non-stage-structured model" << endl;
		}
	}
	// validate sex
	if (isSexDep) {
		if (sex != (current.simLines + 1) % 2) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "Sex must be alternately 0 and 1 if SexDep is 1" << endl;
		}
	}
	else {
		if (sex != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "Sex must be 0 if SexDep is 0" << endl;
		}
	}

	// validate inIndVar
	if (isStageDep && !mustCheckStgDepWithIndVar) {
		if (isIndVar != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLogOfs << "IndVar must be 0 if stage-dependent" << endl;
		}
	}
	else {
		if (isIndVar < 0 || isIndVar > 1) {
			BatchError(whichInputFile, whichLine, 1, "IndVar"); current.errors++;
		}
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
		batchLogOfs << "No. of " << fieldname << " columns must be one fewer than no. of stages, i.e. "
			<< gNbStages - 1 << ", and be sequentially numbered starting from 1";
		break;
	case 555:
		batchLogOfs << "No. of " << fieldname << " columns must equal no. of stages, i.e. "
			<< gNbStages << ", and be sequentially numbered starting from 0";
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

void CtrlFormatError(void)
{
	cout << "Format error in Control file" << endl;
	batchLogOfs << endl << "***" << endl << "*** Format error in Control file:"
		<< gCaseSensitiveStr << " and file names" << gSpecMustMatchStr
		<< endl
		<< "***" << endl;
}

void ArchFormatError(void)
{
	batchLogOfs << "*** Format error in ArchFile:" << gCaseSensitiveStr << gSpecMustMatchStr << endl;
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
		ifsLandFile >> ppLand.landNum >> ppGenLand.fractal >> ppGenLand.continuous
			>> ppLand.dimX >> ppLand.dimY >> ppGenLand.minPct >> ppGenLand.maxPct
			>> ppGenLand.propSuit >> ppGenLand.hurst;
		ppLand.maxX = ppLand.dimX - 1; 
		ppLand.maxY = ppLand.dimY - 1;

		if (ppGenLand.fractal && ppLand.maxX > ppLand.maxY) {
			return -901;
		}
		if (ppGenLand.fractal) {
			if ((ppLand.dimX < 3 || ppLand.dimX % 2 != 1)
				|| (ppLand.dimY < 3 || ppLand.dimY % 2 != 1)) {
				return -902;
			}
		}
		// SCFP 26/9/13 - min and max habitat percentages need to be set for all types of
		// fractal landscape (including discrete), as they are passed to the fractal generator
		// NOTE that will not have been checked for a discrete landscape
		if (ppGenLand.fractal && !ppGenLand.continuous) { 
			ppGenLand.minPct = 1; 
			ppGenLand.maxPct = 100;
		}
		if (ppGenLand.continuous)
			ppLand.nHab = 2;
		else 
			ppLand.nHab = 1;
	}
	else { // imported raster map
		string inNbHab;
		ifsLandFile >> ppLand.landNum >> inNbHab >> gHabMapName >> gPatchMapName;
		ifsLandFile >> gNameCostFile >> gDynLandFileName >> gSpDistFileName;
		if (gLandType == 2) 
			ppLand.nHab = 1; // habitat quality landscape has one habitat class
	}

	pLandscape->setLandParams(ppLand, true);
	pLandscape->setGenLandParams(ppGenLand);

	return ppLand.landNum;
}

//---------------------------------------------------------------------------
int ReadDynLandFile(Landscape* pLandscape) {

	string landchangefile, patchchangefile, costchangefile;
	int change, imported;
	int nbChanges = 0;
	bool usesCosts = false;
	landChange chg;
	landParams ppLand = pLandscape->getLandParams();
	string fname = paramsSim->getDir(1) + gDynLandFileName;

	ifsDynLandFile.open(fname.c_str());
	if (ifsDynLandFile.is_open()) {
		string header;
		int nheaders = 5;
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
		ifsDynLandFile >> chg.chgyear >> landchangefile >> patchchangefile >> costchangefile;
		chg.habfile = paramsSim->getDir(1) + landchangefile;
		chg.pchfile = paramsSim->getDir(1) + patchchangefile;
		usesCosts = costchangefile != "NULL";
		chg.costfile = usesCosts ? paramsSim->getDir(1) + costchangefile : "none";

		nbChanges++;
		pLandscape->addLandChange(chg);

		// read first field on next line
		change = -98765;
		ifsDynLandFile >> change;
		if (ifsDynLandFile.eof()) {
			change = -98765;
		}
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
		if (imported != 0) {
			return imported;
		}

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

int ReadGeneticsFile(ifstream& ifs, Landscape* pLandscape) {

	string indir = paramsSim->getDir(1);
	bool outputGeneValues, outputWeirCockerham, outputWeirHill;
	int outputStartGenetics, outputGeneticInterval;
	set<int> patchList;

	//not ideal to reset these in here 
	pSpecies->resetGeneticParameters();

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
		int genomeSize = stoi(parameters[1]);
		set<int> chrEnds = stringToChromosomeEnds(parameters[2], genomeSize);
		float recombinationRate = parameters[3] == "#" ? 0.0 : stof(parameters[3]);
		outputGeneValues = (parameters[4] == "TRUE");
		outputWeirCockerham = (parameters[5] == "TRUE");
		outputWeirHill = (parameters[6] == "TRUE");
		outputStartGenetics = stoi(parameters[7]);
		outputGeneticInterval = stoi(parameters[8]);

		string inPatches = parameters[9];
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
				nPatchesToSample = stoi(parameters[10]);
			// patchList remains empty, filled when patches are sampled every gen
		}
		const string strNbInds = parameters[11];
		const int nbStages = pSpecies->getStageParams().nStages;
		set<int> stagesToSampleFrom = stringToStages(parameters[12], nbStages);

		pSpecies->setGeneticParameters(chrEnds, genomeSize, recombinationRate,
			patchList, strNbInds, stagesToSampleFrom, nPatchesToSample);
		paramsSim->setGeneticSim(patchSamplingOption, outputGeneValues, outputWeirCockerham, outputWeirHill, outputStartGenetics, outputGeneticInterval);
	}
	else {
		throw runtime_error("GeneticsFile is not open.");
	}
	return 0;
}

int ReadTraitsFile(ifstream& ifs, const int& nbRowsToRead) {

	pSpecies->clearTraitTable();
	int prevsimNb = -998;

	if (ifs.is_open()) {
		//read first header line
		string strLine, entry;

		for (int i = 0; i < nbRowsToRead; i++) {
			
			// Read input row
			std::getline(ifs, strLine);

			// Read input parameters as strings
			stringstream inLine(strLine);
			vector<string> parameters;
			while (std::getline(inLine, entry, '	'))
			{
				parameters.push_back(entry);
			}

			// Create trait from parameters 
			setUpSpeciesTrait(parameters);
		}
	}
	else {
		throw runtime_error("TraitsFile is not open.");
	}
	return 0;
}

// Set up a trait from input parameters and add it Species
void setUpSpeciesTrait(vector<string> parameters) {
	// Assumes all input is correct, errors have been handled by CheckTraits

	const int genomeSize = pSpecies->getGenomeSize();
	TraitType traitType = stringToTraitType(parameters[1]);
	const sex_t sex = stringToSex(parameters[2]);
	if (sex != NA) traitType = addSexDepToTrait(traitType, sex);
	const set<int> positions = stringToLoci(parameters[3], parameters[4], genomeSize);
	const ExpressionType expressionType = stringToExpressionType(parameters[5]);

	// Initial allele distribution parameters
	const DistributionType initDist = stringToDistributionType(parameters[6]);
	const map<GenParamType, float> initParams = stringToParameterMap(parameters[7]);

	// Initial dominance distribution parameters
	const DistributionType initDomDist = stringToDistributionType(parameters[8]);
	const map<GenParamType, float> initDomParams = stringToParameterMap(parameters[9]);

	// Mutation parameters
	bool isInherited = (parameters[10] == "TRUE");
	DistributionType mutationDistribution = isInherited ? 
		stringToDistributionType(parameters[11]) : 
		DistributionType::NONE;
	map<GenParamType, float> mutationParameters;
	if (isInherited) {
		mutationParameters = stringToParameterMap(parameters[12]);
	}

	// Dominance distribution parameters
	const DistributionType dominanceDist = stringToDistributionType(parameters[13]);
	const map<GenParamType, float> dominanceParams = stringToParameterMap(parameters[14]);

	float mutationRate = isInherited ? stof(parameters[15]) : 0.0;
	
	parameters[16].erase(
		// send windows line endings to hell where they belong
		remove(parameters[16].begin(), parameters[16].end(), '\r'),
		parameters[16].end()
	);
	const bool isOutput = parameters[16] == "TRUE";

	// Create species trait
	int ploidy = gNbSexesDisp;
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

//---------------------------------------------------------------------------
int ReadParameters(Landscape* pLandscape)
{
	int errorCode = 0;
	landParams paramsLand = pLandscape->getLandParams();

	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	simParams sim = paramsSim->getSim();

	if (!ifsParamFile.is_open()) {
		cout << endl << "ReadParameters(): ERROR - ParameterFile is not open" << endl;
		return 4086534;
	}

	int gradType, shiftBegin, shiftStop;
	float k, grad_inc, opt_y, f, optExt, shift_rate;
	string inAbsorbing, inShifting, inEnvStoch, inEnvStochType, inLocalExt,
		inSaveMaps, inHeatMaps, inDrawLoaded, inFixRepSeed;

	ifsParamFile >> sim.simulation >> sim.reps >> sim.years;
	ifsParamFile >> inAbsorbing;
	sim.absorbing = (inAbsorbing == "1");

	// Environmental gradient
	ifsParamFile >> gradType;
	ifsParamFile >> grad_inc >> opt_y >> f >> optExt >> inShifting >> shift_rate;

	bool isShifting = (inShifting == "1" && gradType != 0);
	ifsParamFile >> shiftBegin >> shiftStop;
	paramGrad paramsGrad;
	// also set grad.optY0 for species
	paramsGrad.setGradient(gradType, grad_inc, opt_y, f, optExt);
	if (isShifting) paramsGrad.setShifting(shift_rate, shiftBegin, shiftStop);
	else paramsGrad.noShifting();
	pSpecies->setGrad(paramsGrad);

	// Environmental Stochasticity
	envStochParams env;
	ifsParamFile >> inEnvStoch;
	env.usesStoch = inEnvStoch == "1" || inEnvStoch == "2";
	env.stochIsLocal = inEnvStoch == "2";
	if (paramsLand.usesPatches && env.stochIsLocal) errorCode = 101;

	ifsParamFile >> inEnvStochType;
	env.inK = (inEnvStochType == "1");
	float minR, maxR, minK, maxK;
	ifsParamFile >> env.ac >> env.std >> minR >> maxR >> minK >> maxK;
	if (env.inK) {
		float minKK, maxKK;
		minKK = minK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
		maxKK = maxK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
		pSpecies->setMinMax(minKK, maxKK);
	}
	else pSpecies->setMinMax(minR, maxR);

	// Local extinction
	ifsParamFile >> inLocalExt;
	env.usesLocalExt = (inLocalExt == "1");
	if (paramsLand.usesPatches && env.usesLocalExt) errorCode = 102;

	ifsParamFile >> env.locExtProb;
	pSpecies->setStoch(env);

	// Demographic parameters
	ifsParamFile >> dem.propMales >> dem.harem >> dem.bc >> dem.lambda;
	pSpecies->setDemogr(dem);

	// Artificial landscape
	if (gLandType == 9) {
		// only one value of K is read, but it must be applied as the second habitat if the
		// landscape is discrete (the first is the matrix where K = 0) or as the first 
		// (only) habitat if the landscape is continuous
		genLandParams genland = pLandscape->getGenLandParams();
		int nhab = genland.continuous ? 1 : 2;

		pSpecies->createHabK(nhab);
		ifsParamFile >> k;
		k *= (((float)paramsLand.resol * (float)paramsLand.resol)) / 10000.0f;

		if (genland.continuous) {
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
	ifsParamFile	>> spParams.outStartPop	>> spParams.outStartInd
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

	if (spParams.outOccup && sim.reps < 2) errorCode = 103;
	if (paramsLand.usesPatches) {
		if (spParams.outTraitsRows) errorCode = 104;
	}
	else {
		if (spParams.outConnect) errorCode = 105;
	}
	ifsParamFile >> inHeatMaps;
	spParams.saveVisits = inHeatMaps == "1";
	pSpecies->setSpeciesParams(spParams);

	ifsParamFile >> inFixRepSeed;
	sim.fixReplicateSeed = inFixRepSeed == "1";
	paramsSim->setSim(sim);
	return errorCode;
}

//---------------------------------------------------------------------------
int ReadStageStructure()
{
	string name;
	int simulation, postDestructn;
	stageParams sstruct = pSpecies->getStageParams();
	string Inputs = paramsSim->getDir(1);

	ifsStageStructFile >> simulation;
	ifsStageStructFile >> postDestructn >> sstruct.probRep >> sstruct.repInterval >> sstruct.maxAge;
	if (postDestructn == 1) sstruct.disperseOnLoss = true;
	else sstruct.disperseOnLoss = false;

	ifsStageStructFile >> name;
	// 'name' is TransMatrixFile
	ifsTransMatrix.open((Inputs + name).c_str());
	ReadTransitionMatrix(sstruct.nStages, gNbSexesDem, 0, 0);
	ifsTransMatrix.close(); ifsTransMatrix.clear();
	ifsStageStructFile >> sstruct.survival;

	float devCoeff, survCoeff;
	ifsStageStructFile >> sstruct.fecDens >> sstruct.fecStageDens >> name; // 'name' is FecStageWtsFile
	if (name != "NULL") {
		ifsFecDens.open((Inputs + name).c_str());
		ReadStageWeights(1);
		ifsFecDens.close(); 
		ifsFecDens.clear();
	}
	ifsStageStructFile >> sstruct.devDens >> devCoeff >> sstruct.devStageDens >> name; // 'name' is DevStageWtsFile
	if (name != "NULL") {
		ifsDevDens.open((Inputs + name).c_str());
		ReadStageWeights(2);
		ifsDevDens.close(); 
		ifsDevDens.clear();
	}
	ifsStageStructFile >> sstruct.survDens >> survCoeff >> sstruct.survStageDens >> name; // 'name' is SurvStageWtsFile
	if (name != "NULL") {
		ifsSurvDens.open((Inputs + name).c_str());
		ReadStageWeights(3);
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
int ReadTransitionMatrix(short nstages, short nsexesDem, short hab, short season)
{
	int ii;
	int minAge;
	float ss, dd;
	string header;
	demogrParams dem = pSpecies->getDemogrParams();

	// read header line
	for (int i = 0; i < (nstages * nsexesDem) + 2; i++)
	{
		ifsTransMatrix >> header;
	}

	if (matrix != NULL) {
		for (int j = 0; j < matrixsize; j++) delete[] matrix[j];
		delete[] matrix;
		matrix = NULL; matrixsize = 0;
	}

	if (dem.repType != 2) { // asexual or implicit sexual model
	// create a temporary matrix
		matrix = new float* [nstages];
		matrixsize = nstages;
		for (int i = 0; i < nstages; i++)
			matrix[i] = new float[nstages];

		for (int i = 0; i < nstages; i++)
		{ // i = row; j = coloumn
			ifsTransMatrix >> header;
			for (int j = 0; j < nstages; j++)
			{
				ifsTransMatrix >> matrix[j][i];
			}
			ifsTransMatrix >> minAge; pSpecies->setMinAge(i, 0, minAge);
		}

		for (int j = 1; j < nstages; j++)
			pSpecies->setFec(j, 0, matrix[j][0]);
		for (int j = 0; j < nstages; j++)
		{
			ss = 0.0; dd = 0.0;
			for (int i = 0; i < nstages; i++)
			{
				if (i == j) ss = matrix[j][i];
				if (i == (j + 1)) dd = matrix[j][i];
			}
			pSpecies->setSurv(j, 0, ss + dd);
			if ((ss + dd) > 0.0f)
				pSpecies->setDev(j, 0, dd / (ss + dd));
			else
				pSpecies->setDev(j, 0, 0.0);
		}
	}
	else { // complex sexual model
		matrix = new float* [nstages * 2];
		matrixsize = nstages * 2;
		for (int j = 0; j < nstages * 2; j++)
			matrix[j] = new float[nstages * 2 - 1];

		for (int i = 0; i < nstages * 2 - 1; i++)
		{ // i = row; j = coloumn
			ifsTransMatrix >> header;
			for (int j = 0; j < nstages * 2; j++) 
				ifsTransMatrix >> matrix[j][i];
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
		for (int j = 2; j < nstages * 2; j++)
		{
			if (j % 2 == 0)
				pSpecies->setFec(ii, 1, matrix[j][0]);
			else {
				pSpecies->setFec(ii, 0, matrix[j][0]);
				ii++;
			}
		}
		// survival and development of male juveniles
		pSpecies->setSurv(0, 1, (matrix[0][0] + matrix[0][1]));
		if ((matrix[0][0] + matrix[0][1]) > 0.0)
			pSpecies->setDev(0, 1, (matrix[0][1] / (matrix[0][0] + matrix[0][1])));
		else
			pSpecies->setDev(0, 1, 0.0);
		// survival and development of female juveniles
		pSpecies->setSurv(0, 0, (matrix[1][0] + matrix[1][2]));
		if ((matrix[1][0] + matrix[1][2]) > 0.0)
			pSpecies->setDev(0, 0, (matrix[1][2] / (matrix[1][0] + matrix[1][2])));
		else
			pSpecies->setDev(0, 0, 0.0);
		// survival and development of stages 1+
		ii = 1;
		//	for (int j = 2; j < sstruct.nStages*2; j++) 
		for (int j = 2; j < nstages * 2; j++)
		{
			ss = 0.0; dd = 0.0;
			if (j % 2 == 0) { // males
				for (int i = 0; i < nstages * 2 - 1; i++)
				{
					if (j == i + 1) ss = matrix[j][i];
					if (j == i - 1) dd = matrix[j][i];
				}
				pSpecies->setSurv(ii, 1, (ss + dd));
				if ((ss + dd) > 0.0)
					pSpecies->setDev(ii, 1, dd / (ss + dd));
				else
					pSpecies->setDev(ii, 1, 0.0);
			}
			else { // females
				for (int i = 0; i < nstages * 2; i++)
				{
					if (j == i + 1) ss = matrix[j][i];
					if (j == i - 1) dd = matrix[j][i];
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

	if (matrix != NULL) {
		for (int j = 0; j < matrixsize; j++)
			delete[] matrix[j];
		delete[] matrix;
		matrix = NULL; 
		matrixsize = 0;
	}

	return 0;
}

//---------------------------------------------------------------------------
int ReadStageWeights(int option)
{
	string header;
	int i, j, n;
	float f;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	if (dem.repType != 2) 
		n = sstruct.nStages;
	else 
		n = sstruct.nStages * gMaxNbSexes;

	switch (option) {

	case 1: { // fecundity
		// create stage weights matrix
		pSpecies->createDDwtFec(n);
		for (i = 0; i < n + 1; i++) ifsFecDens >> header;
		// read coefficients
		for (i = 0; i < n; i++) {
			ifsFecDens >> header;
			for (j = 0; j < n; j++) {
				ifsFecDens >> f; pSpecies->setDDwtFec(j, i, f);
			}
		}
		break;
	}

	case 2: { // development
		//create stage weights matrix
		pSpecies->createDDwtDev(n);
		for (i = 0; i < n + 1; i++) ifsDevDens >> header;
		//read coefficients
		for (i = 0; i < n; i++) {
			ifsDevDens >> header;
			for (j = 0; j < n; j++) {
				ifsDevDens >> f; pSpecies->setDDwtDev(j, i, f);
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
				ifsSurvDens >> f; pSpecies->setDDwtSurv(j, i, f);
			}
		}
		break;
	}

	}

	return 0;
}

//---------------------------------------------------------------------------
int ReadEmigration()
{
	int errorCode = 0;
	int inFullKernel, inDensDep, inStgDep, inSexDep, inIndVar;
	int Nlines, simulationNb, simNbFirstLine = 0, inStage, inSex, inEmigstage;
	float inEp, inD0, inAlpha, inBeta;
	bool isFirstLine = true;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	emigTraits emigrationTraits;

	// set no.of lines assuming maximum stage- and sex-dependency
	if (sstruct.nStages == 0) Nlines = gNbSexesDisp;
	else Nlines = sstruct.nStages * gNbSexesDisp;

	for (int line = 0; line < Nlines; line++) {

		ifsEmigrationFile >> simulationNb >> inDensDep >> inFullKernel 
				 >> inStgDep >> inSexDep >> inIndVar >> inEmigstage;

		if (isFirstLine) {
			simNbFirstLine = simulationNb;
			emig.densDep = (inDensDep == 1);
			emig.stgDep = (inStgDep == 1);
			emig.indVar = (inIndVar == 1);
			emig.sexDep = (inSexDep == 1);
			if (inEmigstage >= 0 && inEmigstage < sstruct.nStages)
				emig.emigStage = inEmigstage;
			else emig.emigStage = 0;
			// update no.of lines according to known stage- and sex-dependency
			if (emig.stgDep) {
				if (emig.sexDep) Nlines = sstruct.nStages * gNbSexesDisp;
				else Nlines = sstruct.nStages;
			}
			else {
				if (emig.sexDep) Nlines = gNbSexesDisp;
				else Nlines = 1;
			}

			if (inFullKernel == 0) pSpecies->setFullKernel(false); 
			else pSpecies->setFullKernel(true);
			pSpecies->setEmigRules(emig);
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

		if (emig.sexDep) {
			if (emig.stgDep) {
				if (emig.densDep) {
					emigrationTraits.d0 = inD0; 
					emigrationTraits.alpha = inAlpha;
					emigrationTraits.beta = inBeta;
				}
				else {
					emigrationTraits.d0 = inEp; 
					emigrationTraits.alpha = emigrationTraits.beta = 0.0;
				}
				pSpecies->setSpEmigTraits(inStage, inSex, emigrationTraits);
			}
			else { // !emig.stgDep

				if (emig.densDep) {
					emigrationTraits.d0 = inD0;
					emigrationTraits.alpha = inAlpha; 
					emigrationTraits.beta = inBeta;
				}
				else {
					emigrationTraits.d0 = inEp; 
					emigrationTraits.alpha = emigrationTraits.beta = 0.0;
				}
				pSpecies->setSpEmigTraits(0, inSex, emigrationTraits);

			}
		}
		else { // !emig.sexDep
			if (emig.stgDep) {
				if (emig.densDep) {
					emigrationTraits.d0 = inD0; 
					emigrationTraits.alpha = inAlpha; 
					emigrationTraits.beta = inBeta;
					pSpecies->setSpEmigTraits(inStage, 0, emigrationTraits);
				}
				else {
					emigrationTraits.d0 = inEp; 
					emigrationTraits.alpha = emigrationTraits.beta = 0.0;
					pSpecies->setSpEmigTraits(inStage, 0, emigrationTraits);
				}
			}
			else { // !emig.stgDep
				if (emig.densDep) {
					emigrationTraits.d0 = inD0; 
					emigrationTraits.alpha = inAlpha; 
					emigrationTraits.beta = inBeta;
				}
				else {
					emigrationTraits.d0 = inEp; 
					emigrationTraits.alpha = emigrationTraits.beta = 0.0;
				}
				pSpecies->setSpEmigTraits(0, 0, emigrationTraits);
			}
		}

		isFirstLine = false;

	} // end of Nlines for loop

	return errorCode;
}

//---------------------------------------------------------------------------
int ReadTransferFile(Landscape* pLandscape)
{
	int error = 0;
	landParams paramsLand = pLandscape->getLandParams();
	transferRules trfr = pSpecies->getTransferRules();

	// new local variable to replace former global variable
	int TransferType = trfr.usesMovtProc ? trfr.moveType : 0; 

	switch (TransferType) {

	case 0: // negative exponential dispersal kernel
		error = ReadTransferKernels(trfr, paramsLand);
		break; // end of negative exponential dispersal kernel

	case 1: // SMS
		ReadTransferSMS(trfr, paramsLand);
		break; // end of SMS

	case 2: // CRW
		error = ReadTransferCRW(trfr, paramsLand);
		break; // end of CRW

	default:
		error = 440;
		break;
	} // end of switch (TransferType)

	return error;
}

int ReadTransferKernels(transferRules trfr, const landParams& paramsLand) {

	int inKernelType, inDistMort, inIndVar, simNb, inStageDep, inSexDep, inStage, inSex;
	float flushMort;
	int simNbFirstLine = 0;
	stageParams stageStruct = pSpecies->getStageParams();
	demogrParams dem = pSpecies->getDemogrParams();
	trfrKernelParams kernParams;
	int sexKernels = 0;
	bool isFirstLine = true;
	int errorCode = 0;

	// set no.of lines assuming maximum stage- and sex-dependency
	int Nlines = stageStruct.nStages == 0 ? gNbSexesDisp : gNbSexesDisp * stageStruct.nStages;

	for (int line = 0; line < Nlines; line++) {

		ifsTransferFile >> simNb >> inStageDep >> inSexDep >> inKernelType >> inDistMort >> inIndVar;
		if (isFirstLine) {
			simNbFirstLine = simNb;
			trfr.twinKern = (inKernelType == 1);
			trfr.distMort = (inDistMort == 1);
			sexKernels = 2 * inStageDep + inSexDep;
			trfr.indVar = (inIndVar == 1);
			trfr.sexDep = (inSexDep == 1);
			// update no.of lines according to known stage- and sex-dependency
			trfr.stgDep = (inStageDep == 1);
			if (trfr.stgDep) {
				Nlines = inSexDep ? stageStruct.nStages * gNbSexesDisp : stageStruct.nStages;
			}
			else {
				Nlines = inSexDep ? gNbSexesDisp : 1;
			}
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
		else {
			if (sexKernels == 2 || sexKernels == 3) 
				errorCode = 403;
		}

		switch (sexKernels) {

		case 0: // no sex / stage dependence
			ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			pSpecies->setSpKernTraits(0, 0, kernParams, paramsLand.resol);
			break;

		case 1: // sex-dependent
			if (trfr.twinKern)
			{
				ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				ifsTransferFile >> kernParams.meanDist1; kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(0, inSex, kernParams, paramsLand.resol);

			break;

		case 2: // stage-dependent
			if (trfr.twinKern)
			{
				ifsTransferFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				ifsTransferFile >> kernParams.meanDist1; kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(inStage, 0, kernParams, paramsLand.resol);
			break;

		case 3: // sex- & stage-dependent
			if (trfr.twinKern)
			{
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

void ReadTransferSMS(transferRules trfr, const landParams& paramsLand) {

	int inIndVar, inSMType, inAlphaDB, inBetaDB, inStraightenPath, simNb;
	float inHabMort, flushHabMort, inMortHabitat, inMortMatrix;
	int inCostHab, flushCostHab, inCostMatrix;
	trfrMovtParams move;

	ifsTransferFile >> simNb >> inIndVar >> move.pr >> move.prMethod >> move.dp
		>> move.memSize >> move.gb >> move.goalType >> inAlphaDB >> inBetaDB
		>> inStraightenPath >> inSMType >> move.stepMort;

	trfr.indVar = (inIndVar == 1);
	if (move.goalType == 2) { // dispersal bias
		move.alphaDB = inAlphaDB;
		move.betaDB = inBetaDB;
	}
	trfr.habMort = (inSMType == 1);
	move.straightenPath = (inStraightenPath == 1);

	if (!paramsLand.generated) { // imported landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			if (trfr.habMort)
			{ // habitat-dependent step mortality
				for (int i = 0; i < paramsLand.nHabMax; i++)
				{
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
		if (trfr.habMort)
		{ // habitat-dependent step mortality
			// values are for habitat (hab=1) then for matrix (hab=0)
			ifsTransferFile >> inMortHabitat >> inMortMatrix;
			pSpecies->setHabMort(1, inMortHabitat);
			pSpecies->setHabMort(0, inMortMatrix);
		}
		else { // constant step mortality
			ifsTransferFile >> flushHabMort >> flushHabMort;
		}
	}
	trfr.costMap = (gNameCostFile != "NULL") ? true : false;

	if (!paramsLand.generated) { // imported landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			if (trfr.costMap)
			{
				for (int i = 0; i < paramsLand.nHabMax; i++) 
					ifsTransferFile >> flushCostHab;
			}
			else { // not costMap
				for (int i = 0; i < paramsLand.nHabMax; i++) {
					ifsTransferFile >> inCostHab; 
					pSpecies->setHabCost(i, inCostHab);
				}
			}
		}
	}
	else { // artificial landscape
		if (trfr.costMap) // should not occur 
		{
			ifsTransferFile >> flushCostHab >> flushCostHab;
		}
		else { // not costMap
			// costs are for habitat (hab=1) then for matrix (hab=0)
			ifsTransferFile >> inCostHab >> inCostMatrix;
			pSpecies->setHabCost(1, inCostHab);
			pSpecies->setHabCost(0, inCostMatrix);
		}
	}
	pSpecies->setTrfrRules(trfr);
	pSpecies->setSpMovtTraits(move);
}

int ReadTransferCRW(transferRules trfr, const landParams& paramsLand) {

	int inIndVar, inStraightenPath, inSMconst, simNb;
	float inHabMort, flushHabMort;

	int error = 0;
	trfrMovtParams move;
	ifsTransferFile >> simNb >> inIndVar;
	if (inIndVar == 0) trfr.indVar = false;
	else trfr.indVar = true;

	ifsTransferFile >> move.stepLength >> move.rho;
	ifsTransferFile >> inStraightenPath >> inSMconst >> move.stepMort;

	if (inSMconst == 0) trfr.habMort = false;
	else trfr.habMort = true;

	if (inStraightenPath == 0) move.straightenPath = false;
	else move.straightenPath = true;

	//Habitat-dependent per step mortality
	if (trfr.habMort && paramsLand.rasterType != 0)
		error = 434;

	if (!paramsLand.generated && paramsLand.rasterType == 0) { // imported habitat codes landscape
		if (trfr.habMort)
		{ // habitat-dependent step mortality
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
int ReadSettlement()
{
	int Nlines, simNb, simNbFirstLine = 0, inStageDep, inSexDep, inStage, inSex;
	bool isFirstline = true;
	bool mustFindMate;
	int errorCode = 0;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	settleRules srules;
	settleSteps ssteps;
	settleTraits settleDD;
	int sexSettle = 0, inSettleType = 0, inDensDep, inIndVar, inFindMate;

	isFirstline = true;

	// set no.of lines assuming maximum stage- and sex-dependency
	if (sstruct.nStages == 0) Nlines = gNbSexesDisp;
	else Nlines = sstruct.nStages * gNbSexesDisp;

	for (int line = 0; line < Nlines; line++) {

		ifsSettlementFile >> simNb >> inStageDep >> inSexDep >> inStage >> inSex;
		if (!trfr.usesMovtProc)
		{ // dispersal kernel
			ifsSettlementFile >> inSettleType >> inFindMate;
		}
		else {
			ifsSettlementFile >> inDensDep >> inIndVar >> inFindMate;
		}
		mustFindMate = (inFindMate == 1);

		if (isFirstline) {
			simNbFirstLine = simNb;
			sett.stgDep = (inStageDep == 1);
			sett.sexDep = (inSexDep == 1);
			if (trfr.usesMovtProc) {// no ind var for kernels
				sett.indVar = (inIndVar == 1);
			}
			else {
				sett.indVar = false;
			}
			pSpecies->setSettle(sett);

			// update no.of lines according to known stage- and sex-dependency
			Nlines = sett.sexDep ? gNbSexesDisp : 1;
			if (sett.stgDep) Nlines *= sstruct.nStages;
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

			switch (inSettleType) {
			case 0:
				srules.wait = false;
				srules.goToNeighbourLocn = false;
				break;
			case 1:
				srules.wait = true;
				srules.goToNeighbourLocn = false;
				break;
			case 2:
				srules.wait = false;
				srules.goToNeighbourLocn = true;
				break;
			case 3:
				srules.wait = true;
				srules.goToNeighbourLocn = true;
				break;
			}
			srules.findMate = mustFindMate;
			pSpecies->setSettRules(stageToSet, sexToSet, srules);

			// Surely this can be simplified further
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
			if (!sett.sexDep && hasMales) {
				// Must set males
				if (!sett.stgDep) {
					pSpecies->setSettRules(0, 1, srules);
					// males of other stages already set above
					// actually stage 0 males too?
				}
				else {
					pSpecies->setSettRules(stageToSet, 1, srules);
				}
			}
		} // end of dispersal kernel

		isFirstline = false;

	} // end of for line loop

	return errorCode;
}

//---------------------------------------------------------------------------
int ReadInitialisation(Landscape* pLandscape)
{
	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	initParams init = paramsInit->getInit();
	string inputDir = paramsSim->getDir(1);

	int simNb, maxcells;
	float totalProps;
	int errorCode = 0;

	ifsInitFile >> simNb >> init.seedType >> init.freeType >> init.spDistType;

	if (init.seedType == 1 && !paramsLand.useSpDist) 
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
			if(init.seedType!=2){
				totalProps += propStage;
				paramsInit->setProp(stg, propStage);
			}
		}
		if (init.seedType!=2 && totalProps != 1.0) { 
			throw logic_error("The proportion of initial individuals in each stage doesn not sum to 1.");
		}
	}

	paramsInit->setInit(init);

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
		// Nothing, initialisation takes place in Community::initialise()
		break;
	case 2: // from initial individuals file
		if (init.indsFile != prevInitialIndsFile) {
			// read and store the list of individuals to be initialised
			ReadInitIndsFile(0, pLandscape, (inputDir + init.indsFile));
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
int ReadInitIndsFile(int option, Landscape* pLandscape, string indsfile) {
	string header;
	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	initParams init = paramsInit->getInit();

	if (option == 0) { // open file and read header line
		ifsInitIndsFile.open(indsfile.c_str());
		string header;
		int nheaders = 3;
		if (paramsLand.usesPatches) nheaders++;
		else nheaders += 2;
		if (dem.repType > 0) nheaders++;
		if (dem.stageStruct) nheaders += 2;
		for (int i = 0; i < nheaders; i++) ifsInitIndsFile >> header;
		paramsInit->resetInitInds();
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
			paramsInit->addInitInd(iind);
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
void RunBatch(int nSimuls, int nLandscapes, Species* pSpecies)
{
	int land_nr;
	int read_error;
	bool params_ok;
	simParams sim = paramsSim->getSim();

	Landscape* pLandscape = nullptr;  		// pointer to landscape

	speciesMap_t allSpecies{ {0, pSpecies} }; // only one for now

	// Open landscape batch file and read header record
	ifsLandFile.open(landFile);
	if (!ifsLandFile.is_open()) {
		cout << endl << "Error opening landFile - aborting batch run" << endl;
		return;
	}
	flushHeaders(ifsLandFile);

	for (int j = 0; j < nLandscapes; j++) {

		// Create new landscape
		if (pLandscape != nullptr) delete pLandscape;
		pLandscape = new Landscape(allSpecies);
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
		paramsLand.usesPatches = gIsPatchModel;
		paramsLand.resol = gResol;
		paramsLand.rasterType = gLandType;
		if (gLandType == 9) {
			paramsLand.generated = true;
			paramsLand.nHab = 2;
		}
		else {
			paramsLand.generated = false;
			if (gDynLandFileName == "NULL") paramsLand.dynamic = false;
			else paramsLand.dynamic = true;
		}
		paramsLand.nHabMax = gMaxNbHab;
		paramsLand.useSpDist = gUseSpeciesDist;
		paramsLand.spResol = gDistResol;
		pLandscape->setLandParams(paramsLand, true);

		if (gLandType != 9) { // imported landscape
			string pathToHabMap = paramsSim->getDir(1) + gHabMapName;
			int landcode;
			string pathToCostMap;
			if (gNameCostFile == "NULL" || gNameCostFile == "none") pathToCostMap = "NULL";
			else pathToCostMap = paramsSim->getDir(1) + gNameCostFile;

			string pathToPatchMap = paramsLand.usesPatches ? 
				paramsSim->getDir(1) + gPatchMapName : " ";
			landcode = pLandscape->readLandscape(0, pathToHabMap, pathToPatchMap, pathToCostMap);
			if (landcode != 0) {
				cout << "Error reading landscape" << endl;
				landOK = false;
			}

			if (paramsLand.dynamic) {
				landcode = ReadDynLandFile(pLandscape);
				if (landcode != 0) {
					cout << "Error reading dynamic landscape" << endl;
					landOK = false;
				}
			}
			if (gLandType == 0) pLandscape->updateHabitatIndices();

			// Species Distribution
			if (paramsLand.useSpDist) { // read initial species distribution
				string distname = paramsSim->getDir(1) + gSpDistFileName;
				landcode = pLandscape->newDistribution(pSpecies, distname);
				if (landcode != 0) {
					cout << endl << "Error reading initial distribution for landscape "
						<< land_nr << " - aborting" << endl;
					landOK = false;
				}
			}
		} // end of imported landscape

		if (landOK) {

			// Open all other batch files and read header records
			ifsParamFile.open(parameterFile);
			if (!ifsParamFile.is_open()) {
				cout << endl << "Error opening ParameterFile - aborting batch run" << endl;
				return;
			}
			flushHeaders(ifsParamFile);

			if (gStageStruct) {
				ifsStageStructFile.open(stageStructFile);
				flushHeaders(ifsStageStructFile);
			}

			ifsEmigrationFile.open(emigrationFile);
			flushHeaders(ifsEmigrationFile);

			ifsTransferFile.open(transferFile);
			flushHeaders(ifsTransferFile);
			if (pSpecies->getTransferRules().usesMovtProc) {
				if (paramsLand.generated)
					pSpecies->createHabCostMort(paramsLand.nHab);
				else pSpecies->createHabCostMort(paramsLand.nHabMax);
			}
			ifsSettlementFile.open(settleFile);
			flushHeaders(ifsSettlementFile);

			ifsInitFile.open(initialFile);
			flushHeaders(ifsInitFile);

			if (gHasGenetics) {
				ifsGeneticsFile.open(geneticsFile.c_str());
				flushHeaders(ifsGeneticsFile);
				ifsTraitsFile.open(traitsFile.c_str());
				flushHeaders(ifsTraitsFile);
			}

			// nSimuls is the total number of lines (simulations) in
			// the batch and is set in the control function
			string msgsim = "Simulation,";
			string msgerr = ",ERROR CODE,";
			string msgabt = ",simulation aborted";

			for (int i = 0; i < nSimuls; i++) {

				params_ok = true;
				read_error = ReadParameters(pLandscape);
				if (read_error) {
					params_ok = false;
				}
				if (gStageStruct) {
					ReadStageStructure();
				}
				read_error = ReadEmigration();
				if (read_error) {
					params_ok = false;
				}
				read_error = ReadTransferFile(pLandscape);
				if (read_error) {
					params_ok = false;
				}
				read_error = ReadSettlement();
				if (read_error) {
					params_ok = false;
				}
				read_error = ReadInitialisation(pLandscape);
				if (read_error) {
					params_ok = false;
				}

				if (gHasGenetics) {
					read_error = ReadGeneticsFile(ifsGeneticsFile, pLandscape);
					if (read_error) {
						params_ok = false;
					}
					read_error = ReadTraitsFile(ifsTraitsFile, gNbTraitFileRows[i]);
					if (read_error) {
						params_ok = false;
					}
				}
				
				if (params_ok) {

					cout << endl << "Running simulation nr. " 
						<< to_string(paramsSim->getSim().simulation)
						<< " on landscape no. " << to_string(land_nr) << endl;

					// for batch processing, include landscape number in parameter file name
					OutParameters(pLandscape, allSpecies);

					RunModel(pLandscape, i, allSpecies);

				} // end of if (params_ok)
				else {
					cout << endl << "Error in reading parameter file(s)" << endl;
				}
			} // end of nSimuls for loop

			// close input files
			ifsParamFile.close();
			ifsParamFile.clear();
			if (gStageStruct) {
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

			if (gHasGenetics) {
				ifsGeneticsFile.close();
				ifsGeneticsFile.clear();
				ifsTraitsFile.close();
				ifsTraitsFile.clear();
			}

			if (pLandscape != nullptr) {
				delete pLandscape; 
				pLandscape = nullptr;
			}

		} // end of landOK condition

	} // end of nLandscapes loop

	ifsLandFile.close();  
	ifsLandFile.clear();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


