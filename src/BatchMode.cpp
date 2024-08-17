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

ifstream controlFile;
// Note - all batch files are prefixed 'b' here for reasons concerned with RS v1.0
ifstream bParamFile, bLandFile, bDynLandFile;
ifstream bSpDistFile, bStageStructFile, bTransMatrix;
ifstream bStageWeightsFile;
ifstream bEmigrationFile, bTransferFile, bSettlementFile;
ifstream bTraitsFile, bGeneticsFile;
ifstream bInitFile, bInitIndsFile;

ofstream batchLog;

ofstream rsLog; // performance log for recording simulation times, etc.

// NOTE: THE STREAMS USED TO READ THE DATA AT RUN TIME COULD TAKE THE SAME NAMES AS
// USED DURING PARSING (ABOVE)
ifstream parameters;
ifstream ssfile, tmfile, fdfile, ddfile, sdfile;
ifstream emigFile, transFile, settFile, initFile, initIndsFile;
ifstream landfile, dynlandfile;
ifstream ifsGenetics, ifsTraits;

// global variables passed between parsing functions...
// should be removed eventually, maybe share variables through members of a class
int batchnum;
int patchmodel, resolution, landtype, maxNhab, speciesdist, distresolution;
int reproductn;
int repseasons;
int stagestruct, stages, gTransferType;
int sexesDem;		// no. of explicit sexes for demographic model
int gNbSexesDisp;	// no. of explicit sexes for dispersal model
int gFirstSimNb = 0; // not great, globals should not be modified.
int fileNtraits; // no. of traits defined in genetic architecture file
bool gHasGenetics = true;

DispersalTraitInputOptions gDispTraitOpt;
vector<int> gNbTraitFileRows;

rasterdata landraster;
// ...including names of the input files
string parameterFile;
string landFile;
string name_landscape, name_patch, name_dynland, name_sp_dist, gNameCostFile;
string stageStructFile, transMatrix;
string emigrationFile, transferFile, settleFile, geneticsFile, gPathToTraitsFile, initialFile;
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
// Returns input value less next highest power of 2 (for x > 2)
int power2check(int x) {
	if (x < 2) return 0;
	int r = x % 2;
	while (r == 0) {
		x /= 2; r = x % 2;
	}
	return x;
}

//---------------------------------------------------------------------------
batchfiles ParseControlAndCheckInputFiles(string pathToControlFile, string indir, string outdir)
{
	batchfiles b;
	int lines, nSimuls;
	int nbErrors = 0;
	string paramname, filename, fname, batchLogPath, header;
	string whichInputFile = "Control file";
	bool anyFormatError = false;
	b.ok = true; 
	b.nSimuls = 0; 
	b.nLandscapes = 0;

	// open batch log file
	batchLogPath = outdir + "BatchLog.txt";
	batchLog.open(batchLogPath.c_str());
	if (!batchLog.is_open()) {
		cout << "Error opening batch output log file " << batchLogPath << endl;
		b.ok = false;
		return b;
	}

	controlFile.open(pathToControlFile.c_str());

	if (!controlFile.is_open()) {
		cout << "Error opening Control file: " << pathToControlFile << endl;
		batchLog << "Error opening Control file: " << pathToControlFile << endl;
		b.ok = false;
		if (batchLog.is_open()) { 
			batchLog.close(); 
			batchLog.clear(); 
		}
		return b;
	}
	else {
		batchLog << "Checking Control file " << pathToControlFile << endl;
	}

	// Check fixed model parameters

	controlFile >> paramname >> batchnum;
	if (paramname == "BatchNum") {
		if (batchnum < 0) {
			BatchError(whichInputFile, -999, 19, "BatchNum"); nbErrors++;
		}
		else b.batchNum = batchnum;
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> patchmodel;
	if (paramname == "PatchModel") {
		if (patchmodel != 0 && patchmodel != 1) {
			BatchError(whichInputFile, -999, 1, "PatchModel"); nbErrors++;
		}
		else b.patchmodel = patchmodel;
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> resolution;
	if (paramname == "Resolution") {
		if (resolution < 1) {
			BatchError(whichInputFile, -999, 11, "Resolution"); nbErrors++;
		}
		else b.resolution = resolution;
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> landtype;
	if (paramname == "LandType") {
		if (landtype != 0 && landtype != 2 && landtype != 9) {
			BatchError(whichInputFile, -999, 0, "LandType");
			batchLog << "LandType must be 0, 2 or 9" << endl;
			nbErrors++;
		}
		else {
			if (landtype == 9 && patchmodel) {
				BatchError(whichInputFile, -999, 0, "LandType");
				batchLog << "LandType may not be 9 for a patch-based model" << endl;
				nbErrors++;
			}
			else b.landtype = landtype;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> maxNhab;
	if (paramname == "MaxHabitats") {
		if (landtype == 0) { // raster with unique habitat codes
			if (maxNhab < 2) {
				BatchError(whichInputFile, -999, 12, "MaxHabitats"); nbErrors++;
			}
			else b.maxNhab = maxNhab;
		}
		else { // raster with habitat quality OR artificial landscape
			if (maxNhab != 1) {
				BatchError(whichInputFile, -999, 0, " "); nbErrors++;
				batchLog << "MaxHabitats must be 1 for LandType = " << landtype << endl;
			}
			else {
				if (landtype == 9) // artificial landscape
					// although the user enters 1, the actual number of habitats is 2
					b.maxNhab = 2;
				else
					b.maxNhab = maxNhab;
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> speciesdist;
	if (paramname == "SpeciesDist") {
		if (speciesdist != 0 && speciesdist != 1) {
			BatchError(whichInputFile, -999, 1, "SpeciesDist"); nbErrors++;
		}
		else {
			if (speciesdist != 0 && landtype == 9) {
				BatchError(whichInputFile, -999, 0, "SpeciesDist");
				batchLog << "SpeciesDist must be 0 for an artificial landscape" << endl;
				nbErrors++;

			}
			else b.speciesdist = speciesdist;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> distresolution;
	if (paramname == "DistResolution") {
		if (speciesdist == 1) { // distribution resolution is required
			if (distresolution < resolution) {
				BatchError(whichInputFile, -999, 0, "DistResolution");
				batchLog << "DistResolution may not be less than Resolution" << endl;
				nbErrors++;
			}
			else {
				if (distresolution % resolution) {
					BatchError(whichInputFile, -999, 0, "DistResolution");
					batchLog << "DistResolution must be an integer multiple of Resolution" << endl;
					nbErrors++;
				}
				else b.distresolution = distresolution;
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> reproductn;
	sexesDem = gNbSexesDisp = 0;
	if (paramname == "Reproduction") {
		if (reproductn != 0 && reproductn != 1 && reproductn != 2) {
			BatchError(whichInputFile, -999, 2, "Reproduction"); nbErrors++;
		}
		else {
			switch (reproductn) {
			case 0: { sexesDem = 1; gNbSexesDisp = 1; break; }
			case 1: { sexesDem = 1; gNbSexesDisp = 2; break; }
			case 2: { sexesDem = 2; gNbSexesDisp = 2; break; }
			}
			b.reproductn = reproductn; 
			b.sexesDem = sexesDem; 
			b.nbSexesDisp = gNbSexesDisp;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> repseasons;
	if (paramname == "RepSeasons") {
		if (repseasons < 1) {
			BatchError(whichInputFile, -999, 11, "RepSeasons"); nbErrors++;
		}
		else b.repseasons = repseasons;
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> stagestruct;
	if (paramname == "StageStruct") {
		if (stagestruct != 0 && stagestruct != 1) {
			BatchError(whichInputFile, -999, 1, "StageStruct"); nbErrors++;
		}
		else b.stagestruct = stagestruct;
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> stages;
	if (paramname == "Stages") {
		if (stagestruct) {
			if (stages < 2 || stages > 10) {
				BatchError(whichInputFile, -999, 0, " "); nbErrors++;
				batchLog << "Stages must be between 2 and 10" << endl;
			}
			b.stages = stages;
		}
		else { // non-stage-structured model must have 2 stages
			b.stages = stages = 2;
		}
	}
	else anyFormatError = true; // wrong control file format

	controlFile >> paramname >> gTransferType;
	if (paramname == "Transfer") {
		if (gTransferType < 0 || gTransferType > 2) {
			BatchError(whichInputFile, -999, 2, "Transfer"); nbErrors++;
		}
		else b.transfer = gTransferType;
	}
	else anyFormatError = true; // wrong control file format

	if (anyFormatError || nbErrors > 0) { // terminate batch error checking
		if (anyFormatError) {
			CtrlFormatError();
		}
		batchLog << endl
			<< "*** Model parameters in Control file must be corrected before further input file checks are conducted"
			<< endl;
		batchLog.close(); 
		batchLog.clear();
		b.ok = false;
		controlFile.close(); 
		controlFile.clear();
		return b;
	}

	// Check parameter file
	controlFile >> paramname >> filename;
	if (paramname == "ParameterFile" && !anyFormatError) {
		fname = indir + filename;
		batchLog << endl << "Checking " << paramname << " " << fname << endl;
		bParamFile.open(fname.c_str());
		if (bParamFile.is_open()) {
			b.nSimuls = CheckParameterFile();
			if (b.nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname, b.nSimuls, 0);
				b.parameterFile = fname;
			}
			bParamFile.close();
		}
		else {
			OpenError(paramname, fname); b.ok = false;
			cout << "Unable to open ParameterFile" << endl;
		}
		bParamFile.clear();
		if (!b.ok) {
			batchLog << endl
				<< "*** ParameterFile must be corrected before further input file checks are conducted"
				<< endl;
			batchLog.close(); batchLog.clear();
			b.ok = false;
			controlFile.close(); controlFile.clear();
			return b;
		}
	}
	else anyFormatError = true; // wrong control file format
	if (bParamFile.is_open()) bParamFile.close();
	bParamFile.clear();

	// Check land file
	controlFile >> paramname >> filename;
	if (paramname == "LandFile" && !anyFormatError) {
		fname = indir + filename;
		batchLog << endl << "Checking " << paramname << " " << fname << endl;
		bLandFile.open(fname.c_str());
		if (bLandFile.is_open()) {
			lines = CheckLandFile(landtype, indir);
			if (lines < 0) {
				b.ok = false;
				if (lines < -111)
					batchLog << "*** Format error in " << paramname << endl;
			}
			else {
				FileOK(paramname, lines, 1);
				b.landFile = fname; b.nLandscapes = lines;
			}
			bLandFile.close();
		}
		else {
			OpenError(paramname, fname); b.ok = false;
		}
		bLandFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check stage structure file if required file
	controlFile >> paramname >> filename;
	batchLog << endl;
	if (paramname == "StageStructFile" && !anyFormatError) {
		if (filename == "NULL") {
			if (stagestruct) {
				batchLog << "*** File name is required for " << paramname << endl;
				b.ok = false;
			}
			else b.stageStructFile = filename;
		}
		else { // filename is not NULL
			if (stagestruct) { // check file only if it is required
				fname = indir + filename;
				batchLog << "Checking " << paramname << " " << fname << endl;
				bStageStructFile.open(fname.c_str());
				if (bStageStructFile.is_open()) {
					nSimuls = CheckStageFile(indir);
					if (nSimuls < 0) {
						b.ok = false;
					}
					else {
						FileOK(paramname, nSimuls, 0);
						if (nSimuls != b.nSimuls) {
							SimulnCountError(filename); b.ok = false;
						}
						else b.stageStructFile = fname;
					}
					bStageStructFile.close();
				}
				else {
					OpenError(paramname, fname); b.ok = false;
				}
				bStageStructFile.clear();
			} // end of required
			else { // file is not required, and filename should be NULL
				if (filename != "NULL") {
					batchLog << "*** File name for stageStructFile should be NULL as StageStruct = "
						<< stagestruct << endl;
					b.ok = false;
				}
			}
		}
	}
	else anyFormatError = true; // wrong control file format

	// Check emigration file
	controlFile >> paramname >> filename;
	if (paramname == "EmigrationFile" && !anyFormatError) {
		fname = indir + filename;
		batchLog << endl << "Checking " << paramname << " " << fname << endl;
		bEmigrationFile.open(fname.c_str());
		if (bEmigrationFile.is_open()) {
			nSimuls = CheckEmigFile();
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); 
					b.ok = false;
				}
				else b.emigrationFile = fname;
			}
			bEmigrationFile.close();
		}
		else {
			OpenError(paramname, fname); 
			b.ok = false;
		}
		bEmigrationFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check transfer file
	controlFile >> paramname >> filename;
	if (paramname == "TransferFile" && !anyFormatError) {
		fname = indir + filename;
		batchLog << endl << "Checking " << paramname << " " << fname << endl;
		bTransferFile.open(fname.c_str());
		if (bTransferFile.is_open()) {
			nSimuls = CheckTransferFile(indir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else b.transferFile = fname;
			}
			bTransferFile.close(); bTransferFile.clear();
		}
		else {
			OpenError(paramname, fname); b.ok = false;
		}
		bTransferFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check settlement file
	controlFile >> paramname >> filename;
	if (paramname == "SettlementFile" && !anyFormatError) {
		fname = indir + filename;
		batchLog << endl << "Checking " << paramname << " " << fname << endl;
		bSettlementFile.open(fname.c_str());
		if (bSettlementFile.is_open()) {
			nSimuls = CheckSettleFile();
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); 
					b.ok = false;
				}
				else b.settleFile = fname;
			}
			bSettlementFile.close();
		}
		else {
			OpenError(paramname, fname); 
			b.ok = false;
		}
		bSettlementFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	// Check genetics file if required file
	controlFile >> paramname >> filename;
	batchLog << endl;
	if (paramname == "GeneticsFile" && !anyFormatError) {
		if (filename == "NULL") {
			if (gDispTraitOpt.isEmigIndVar
				|| gDispTraitOpt.isSettIndVar
				|| gDispTraitOpt.isKernTransfIndVar
				|| gDispTraitOpt.isSMSTransfIndVar
				)
			{
				batchLog << "Error: one or more dispersal traits has been set to IndVar." << endl;
				b.ok = false;
			}
			else {
				gHasGenetics = false;
				batchLog << "No genetics required " << paramname << endl;
			}
		}
		else {
			gHasGenetics = true;
			fname = indir + filename;
			batchLog << "Checking " << paramname << " " << fname << endl;
			bGeneticsFile.open(fname.c_str());
			if (bGeneticsFile.is_open()) {
				nSimuls = CheckGeneticsFile(indir);
				if (nSimuls < 0) {
					b.ok = false;
				}
				else {
					FileOK(paramname, nSimuls, 0);
					b.geneticsFile = fname;
				}
				bGeneticsFile.close();
			}
			else {
				OpenError(paramname, fname); b.ok = false;
			}
			bGeneticsFile.clear();
		}
	}
	else anyFormatError = true; // wrong control file format

	// Check initialisation file
	controlFile >> paramname >> filename;
	if (paramname == "InitialisationFile" && !anyFormatError) {
		fname = indir + filename;
		batchLog << endl << "Checking " << paramname << " " << fname << endl;
		bInitFile.open(fname.c_str());
		if (bInitFile.is_open()) {
			nSimuls = CheckInitFile(indir);
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else b.initFile = fname;
			}
			bInitFile.close();
		}
		else {
			OpenError(paramname, fname); b.ok = false;
		}
		bInitFile.clear();
	}
	else anyFormatError = true; // wrong control file format

	if (anyFormatError) {
		CtrlFormatError();
		b.ok = false;
	}

	if (controlFile.is_open()) { controlFile.close(); controlFile.clear(); }
	if (batchLog.is_open()) { batchLog.close(); batchLog.clear(); }

	// NOTE: THE FOLLOWING ELEMENTS COULD BE REMOVED FROM b ...
	parameterFile = b.parameterFile;
	landFile = b.landFile;
	stageStructFile = b.stageStructFile;
	emigrationFile = b.emigrationFile;
	transferFile = b.transferFile;
	settleFile = b.settleFile;
	geneticsFile = b.geneticsFile;
	initialFile = b.initFile;

	return b;
}

//---------------------------------------------------------------------------
int CheckParameterFile(void)
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
	bParamFile >> header; if (header != "Simulation") nbErrors++;
	bParamFile >> header; if (header != "Replicates") nbErrors++;
	bParamFile >> header; if (header != "Years") nbErrors++;
	bParamFile >> header; if (header != "Absorbing") nbErrors++;
	bParamFile >> header; if (header != "Gradient") nbErrors++;
	bParamFile >> header; if (header != "GradSteep") nbErrors++;
	bParamFile >> header; if (header != "Optimum") nbErrors++;
	bParamFile >> header; if (header != "f") nbErrors++;
	bParamFile >> header; if (header != "LocalExtOpt") nbErrors++;
	bParamFile >> header; if (header != "Shifting") nbErrors++;
	bParamFile >> header; if (header != "ShiftRate") nbErrors++;
	bParamFile >> header; if (header != "ShiftStart") nbErrors++;
	bParamFile >> header; if (header != "ShiftEnd") nbErrors++;
	bParamFile >> header; if (header != "EnvStoch") nbErrors++;
	bParamFile >> header; if (header != "EnvStochType") nbErrors++;
	bParamFile >> header; if (header != "ac") nbErrors++;
	bParamFile >> header; if (header != "std") nbErrors++;
	bParamFile >> header; if (header != "minR") nbErrors++;
	bParamFile >> header; if (header != "maxR") nbErrors++;
	bParamFile >> header; if (header != "minK") nbErrors++;
	bParamFile >> header; if (header != "maxK") nbErrors++;
	bParamFile >> header; if (header != "LocalExt") nbErrors++;
	bParamFile >> header; if (header != "LocalExtProb") nbErrors++;
	bParamFile >> header; if (header != "PropMales") nbErrors++;
	bParamFile >> header; if (header != "Harem") nbErrors++;
	bParamFile >> header; if (header != "bc") nbErrors++;
	bParamFile >> header; if (header != "Rmax") nbErrors++;
	for (i = 0; i < maxNhab; i++) {
		Kheader = "K" + to_string(i + 1);
		bParamFile >> header; 
		if (header != Kheader) nbKerrors++;
	}
	bParamFile >> header; if (header != "OutStartPop") nbErrors++;
	bParamFile >> header; if (header != "OutStartInd") nbErrors++;
	bParamFile >> header; if (header != "OutStartTraitCell") nbErrors++;
	bParamFile >> header; if (header != "OutStartTraitRow") nbErrors++;
	bParamFile >> header; if (header != "OutStartConn") nbErrors++;
	bParamFile >> header; if (header != "OutIntRange") nbErrors++;
	bParamFile >> header; if (header != "OutIntOcc") nbErrors++;
	bParamFile >> header; if (header != "OutIntPop") nbErrors++;
	bParamFile >> header; if (header != "OutIntInd") nbErrors++;
	bParamFile >> header; if (header != "OutIntTraitCell") nbErrors++;
	bParamFile >> header; if (header != "OutIntTraitRow") nbErrors++;
	bParamFile >> header; if (header != "OutIntConn") nbErrors++;
	bParamFile >> header; if (header != "SaveMaps") nbErrors++;
	bParamFile >> header; if (header != "MapsInterval") nbErrors++;
	bParamFile >> header; if (header != "SMSHeatMap") nbErrors++;
	bParamFile >> header; if (header != "DrawLoadedSp") nbErrors++;
	bParamFile >> header; if (header != "FixReplicateSeed") nbErrors++;

	if (nbErrors > 0 || nbKerrors > 0) {
		FormatError(whichFile, nbErrors);
		batchLog << "*** ParameterFile column headers are incorrect." << endl;
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
	bParamFile >> simNb; // first simulation number
	if (simNb == errSimNb) {
		batchLog << "*** Error in ParameterFile - first simulation number could not be read." << endl;
		nbErrors++;
	}
	else if (simNb < 0) {
		batchLog << "*** Error in ParameterFile - first simulation number must be >= 0" << endl;
		nbErrors++;
	}
	else {
		prevsimul = gFirstSimNb = simNb; 
		nSimuls++;
	}
	while (simNb != -98765) {
		bParamFile >> inReplicates; 
		if (inReplicates <= 0) { 
			BatchError(whichFile, whichLine, 11, "Replicates"); 
			nbErrors++; 
		}
		bParamFile >> inYears; 
		if (inYears <= 0) {
			BatchError(whichFile, whichLine, 11, "Years"); 
			nbErrors++; 
		}
		bParamFile >> inAbsorb;
		if (inAbsorb != 0 && inAbsorb != 1) { 
			BatchError(whichFile, whichLine, 1, "Absorbing"); 
			nbErrors++; 
		}
		bParamFile >> inGradient;
		if (patchmodel) {
			if (inGradient != 0) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "Gradient must be 0 for patch-based model" << endl;
				nbErrors++;
				inGradient = 0; // to prevent checking of subsequent fields
			}
			inGradient = 0; // to prevent unnecessary checking of subsequent fields
		}
		else { // cell-based model
			if (inGradient < 0 || inGradient > 3) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "Gradient must be between 0 and 3 for cell-based model" << endl;
				nbErrors++;
			}
		}
		bParamFile >> inGradSteep;
		if (inGradient && inGradSteep < 0.0) {
			BatchError(whichFile, whichLine, 19, "GradSteep"); 
			nbErrors++; 
		}
		bParamFile >> inOptimum;
		if (inGradient && inOptimum < 0) {
			BatchError(whichFile, whichLine, 19, "Optimum"); 
			nbErrors++; 
		}
		bParamFile >> inGradScalingFactor;
		if (inGradient && inGradScalingFactor < 0.0) {
			BatchError(whichFile, whichLine, 19, "f"); 
			nbErrors++; 
		}
		bParamFile >> inLocalExtOpt;
		if (inGradient == 4 && (inLocalExtOpt < 0.0 || inLocalExtOpt >= 1.0))
		{
			BatchError(whichFile, whichLine, 20, "LocalExtOpt"); 
			nbErrors++;
		}
		bParamFile >> inShifting;
		if (inGradient && (inShifting != 0 && inShifting != 1)) { 
			BatchError(whichFile, whichLine, 1, "Shifting");
			nbErrors++; 
		}
		bParamFile >> inShiftRate;
		if (inGradient && inShifting && inShiftRate <= 0.0) {
			BatchError(whichFile, whichLine, 10, "ShiftRate"); 
			nbErrors++; 
		}
		bParamFile >> inShiftStart;
		if (inGradient && inShifting && inShiftStart <= 0) {
			BatchError(whichFile, whichLine, 10, "ShiftStart");
			nbErrors++;
		}
		bParamFile >> inShiftEnd;
		if (inGradient && inShifting && inShiftEnd <= inShiftStart) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "ShiftEnd must be greater than ShiftStart" << endl;
			nbErrors++;
		}
		bParamFile >> inEnvStoch;
		if (patchmodel == 0) { // cell-based model
			if (inEnvStoch != 0 && inEnvStoch != 1 && inEnvStoch != 2) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "EnvStoch must be 0, 1 or 2 for cell-based model" << endl;
				nbErrors++;
			}
		}
		else { // patch-based model
			if (inEnvStoch != 0 && inEnvStoch != 1) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "EnvStoch must be 0 or 1 for patch-based model" << endl;
				nbErrors++;
			}
		}
		bParamFile >> inStochType;
		if (inEnvStoch && (inStochType < 0 || inStochType > 1)) {
			BatchError(whichFile, whichLine, 1, "EnvStochType"); 
			nbErrors++;
		}
		bParamFile >> inStochAC;
		if (inEnvStoch && (inStochAC < 0.0 || inStochAC >= 1.0)) {
			BatchError(whichFile, whichLine, 20, "ac"); 
			nbErrors++; 
		}
		bParamFile >> inStochStD;
		if (inEnvStoch && (inStochStD <= 0.0 || inStochStD > 1.0)) {
			BatchError(whichFile, whichLine, 20, "std"); 
			nbErrors++; 
		}
		bParamFile >> inMinR;
		if (inEnvStoch && inStochType == 0 && inMinR <= 0.0) { 
			BatchError(whichFile, whichLine, 10, "minR"); 
			nbErrors++; 
		}
		bParamFile >> inMaxR;
		if (inEnvStoch && inStochType == 0 && inMaxR <= inMinR) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "maxR must be greater than minR" << endl;
			nbErrors++;
		}
		bParamFile >> inMinK >> inMaxK;
		if (inEnvStoch && inStochType == 1) {
			if (inMinK <= 0.0) { 
				BatchError(whichFile, whichLine, 10, "minK"); 
				nbErrors++; 
			}
			if (inMaxK <= inMinK) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "maxK must be greater than minK" << endl;
				nbErrors++;
			}
		}
		bParamFile >> inLocalExt;
		if (patchmodel == 0) { // cell-based model
			if (inLocalExt < 0 || inLocalExt > 1) {
				BatchError(whichFile, whichLine, 1, "LocalExt");
				nbErrors++;
			}
			else {
				if (inGradient == 4) { // gradient in local extinction probability
					if (inLocalExt != 0) {
						BatchError(whichFile, whichLine, 0, " ");
						batchLog << "LocalExt must be zero if Gradient is 4" << endl;
						nbErrors++;
					}
				}
			}
		}
		else { // patch-based model
			if (inLocalExt != 0) {
				BatchError(whichFile, whichLine, 0, "null");
				batchLog << "LocalExt must be 0 for patch-based model" << endl;
				nbErrors++;
			}
		}
		bParamFile >> inLocalExtProb;
		if (patchmodel == 0 && inLocalExt == 1 && (inLocalExtProb <= 0.0 || inLocalExtProb >= 1.0))
		{
			BatchError(whichFile, whichLine, 20, "LocalExtProb"); 
			nbErrors++;
		}
		bParamFile >> inPropMales;
		if (reproductn && (inPropMales <= 0.0 || inPropMales >= 1.0)) {
			BatchError(whichFile, whichLine, 20, "PropMales");
			nbErrors++;
		}
		bParamFile >> inHarem;
		if (reproductn == 2 && inHarem <= 0.0) {
			BatchError(whichFile, whichLine, 10, "Harem"); 
			nbErrors++; 
		}
		bParamFile >> inBc;
		if (stagestruct == 0 && inBc <= 0.0) {
			BatchError(whichFile, whichLine, 10, "bc"); 
			nbErrors++; 
		}
		bParamFile >> inRmax;
		if (stagestruct == 0 && inRmax <= 0.0) {
			BatchError(whichFile, whichLine, 10, "Rmax");
			nbErrors++; 
		}
		sum_K = 0.0; 
		min_K = 9999999.0; 
		max_K = 0.0;
		for (i = 0; i < maxNhab; i++) {
			bParamFile >> inK;
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
			batchLog << "At least one K column must be non-zero" << endl;
		}
		else {
			if (inEnvStoch && inStochType == 1) { // environmental stochasticity in K
				if (min_K < inMinK || max_K > inMaxK) {
					BatchError(whichFile, whichLine, 0, " "); 
					nbErrors++;
					batchLog << "Non-zero K values must lie between minK and maxK" << endl;
				}
			}
		}

		bParamFile >> inOutStartPop;
		if (inOutStartPop < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartPop"); 
			nbErrors++; 
		}
		bParamFile >> inOutStartInd;
		if (inOutStartInd < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartInd"); 
			nbErrors++; 
		}
		bParamFile >> inOutStartTraitCell;
		if (inOutStartTraitCell < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartTraitCell"); 
			nbErrors++; 
		}
		bParamFile >> inOutStartTraitRow;
		if (inOutStartTraitRow < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartTraitRow"); 
			nbErrors++; 
		}
		bParamFile >> inOutStartConn;
		if (inOutStartConn < 0) {
			BatchError(whichFile, whichLine, 19, "OutStartConn"); 
			nbErrors++; 
		}
		bParamFile >> inOutIntRange;
		if (inOutIntRange < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntRange"); 
			nbErrors++;
		}
		bParamFile >> inOutIntOcc;
		if (inOutIntOcc < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntOcc"); 
			nbErrors++; 
		}
		else {
			if (landtype == 9) {
				if (inOutIntOcc > 0) {
					BatchError(whichFile, whichLine, 0, " "); 
					nbErrors++;
					batchLog << "OutIntOcc must be zero for a generated landscape" << endl;
				}
			}
			else {
				if (inReplicates < 2 && inOutIntOcc > 0) {
					BatchError(whichFile, whichLine, 0, " "); 
					nbErrors++;
					batchLog << "OutIntOcc may be non-zero only if Replicates >= 2" << endl;
				}
			}
		}
		bParamFile >> inOutIntPop;
		if (inOutIntPop < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntPop");
			nbErrors++; 
		}
		bParamFile >> inOutIntInd;
		if (inOutIntInd < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntInd");
			nbErrors++;
		}
		bParamFile >> inOutIntTraitCell;
		if (inOutIntTraitCell < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntTraitCell");
			nbErrors++; 
		}
		bParamFile >> inOutIntTraitRow;
		if (inOutIntTraitRow < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntTraitRow"); 
			nbErrors++;
		}
		bParamFile >> inOutIntConn;
		if (inOutIntConn < 0) {
			BatchError(whichFile, whichLine, 19, "OutIntConn"); 
			nbErrors++; 
		}
		else {
			if (patchmodel != 1 && inOutIntConn > 0) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "OutIntConn may be >0 only if PatchModel is 1" << endl;
				nbErrors++;
			}
		}
		bParamFile >> inSaveMaps; 
		if (inSaveMaps != 0 && inSaveMaps != 1)
		{
			BatchError(whichFile, whichLine, 1, "SaveMaps"); 
			nbErrors++;
		}
		bParamFile >> inMapsInterval; 
		if (inSaveMaps == 1 && inMapsInterval < 1) {
			BatchError(whichFile, whichLine, 11, "MapsInterval");
			nbErrors++;
		}
		bParamFile >> inSMSHeatMap; 
		if (inSMSHeatMap != 0 && inSMSHeatMap != 1) {
			BatchError(whichFile, whichLine, 1, "SMSHeatMap");
			nbErrors++;
		}
		bParamFile >> inDrawLoadedSp; 
		if (inSaveMaps == 1 && (inDrawLoadedSp != 0 && inDrawLoadedSp != 1)) {
			BatchError(whichFile, whichLine, 1, "DrawLoadedSp");
			nbErrors++;
		}
		bParamFile >> inFixReplicateSeed; 
		if (inFixReplicateSeed != 0 && inFixReplicateSeed != 1) {
			BatchError(whichFile, whichLine, 1, "FixReplicateSeed");
			nbErrors++;
		}

		whichLine++;
		// read next simulation number
		simNb = -98765;
		bParamFile >> simNb;
		if (bParamFile.eof()) {
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
	if (!bParamFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}

	if (nbErrors > 0) return -111;
	else return nSimuls;
}

int CheckLandFile(int landtype, string indir)
{
	string fname, header, intext, ftype;
	int j, inint, line;
	float infloat;
	rasterdata patchraster, spdistraster, costraster;
	int errors = 0;
	int totlines = 0;
	vector <int> landlist;
	string filetype = "LandFile";

	if (landtype == 0 || landtype == 2) { // real landscape
		// Parse header line;
		bLandFile >> header; if (header != "LandNum") errors++;
		bLandFile >> header; if (header != "Nhabitats") errors++;
		bLandFile >> header; if (header != "LandscapeFile") errors++;
		bLandFile >> header; if (header != "PatchFile") errors++;
		bLandFile >> header; if (header != "CostMapFile") errors++;
		bLandFile >> header; if (header != "DynLandFile") errors++;
		bLandFile >> header; if (header != "SpDistFile") errors++;
		if (errors > 0) {
			FormatError(filetype, 0);
			batchLog << "*** Ensure format is correct for real landscape" << endl;
			return -111;
		}
		// Parse data lines
		line = 1;
		inint = -98765;
		bLandFile >> inint;
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
			bLandFile >> inint;
			if (landtype == 0) { // raster map with unique habitat codes
				if (inint < 0) {
					BatchError(filetype, line, 10, "Nhabitats"); errors++;
				}
				if (inint > maxNhab) {
					BatchError(filetype, line, 0, " ");
					batchLog << "Nhabitats may not exceed MaxHabitats in Control file" << endl;
					errors++;
				}
			}

			// check landscape filename
			ftype = "LandscapeFile";
			bLandFile >> intext;
			fname = indir + intext;
			landraster = CheckRasterFile(fname);
			if (landraster.ok) {
				if (landraster.cellsize == resolution)
					batchLog << ftype << " headers OK: " << fname << endl;
				else {
					errors++;
					batchLog << gResolOfStr << ftype << " " << fname
						<< gResolNotMatchStr << endl;
				}
			}
			else {
				errors++;
				if (landraster.errors == -111)
					OpenError(ftype, fname);
				else
					FormatError(fname, landraster.errors);
			}

			// check patch map filename
			ftype = "PatchFile";
			bLandFile >> intext;
			if (intext == "NULL") {
				if (patchmodel) {
					BatchError(filetype, line, 0, " "); errors++;
					batchLog << ftype << gPatchReqdStr << endl;
				}
			}
			else {
				if (patchmodel) {
					fname = indir + intext;
					patchraster = CheckRasterFile(fname);
					if (patchraster.ok) {
						if (patchraster.cellsize == resolution) {
							if (patchraster.ncols == landraster.ncols
								&& patchraster.nrows == landraster.nrows
								&& patchraster.cellsize == landraster.cellsize
								&& (int)patchraster.xllcorner == (int)landraster.xllcorner
								&& (int)patchraster.yllcorner == (int)landraster.yllcorner) {
								batchLog << ftype << " headers OK: " << fname << endl;
							}
							else {
								batchLog << gHeadersOfStr << ftype << " " << fname
									<< gHeadersNotMatchStr << endl;
								errors++;
							}
						}
						else {
							batchLog << gResolOfStr << ftype << " " << fname
								<< gResolNotMatchStr << endl;
							errors++;
						}
					}
					else {
						errors++;
						if (patchraster.errors == -111)
							OpenError(ftype, fname);
						else
							FormatError(fname, patchraster.errors);
					}
				}
			}

			// check cost map filename
			ftype = "CostMapFile";
			bLandFile >> gNameCostFile;
			if (gNameCostFile == "NULL") {
				if (gTransferType == 1) { // SMS
					if (landtype == 2) {
						BatchError(filetype, line, 0, " "); errors++;
						batchLog << ftype << " is required for a habitat quality landscape" << endl;
					}
				}
			}
			else {
				if (gTransferType == 1) { // SMS
					fname = indir + gNameCostFile;
					costraster = CheckRasterFile(fname);
					if (costraster.ok) {
						if (costraster.cellsize == resolution) {
							if (costraster.ncols == landraster.ncols
								&& costraster.nrows == landraster.nrows
								&& costraster.cellsize == landraster.cellsize
								&& (int)costraster.xllcorner == (int)landraster.xllcorner
								&& (int)costraster.yllcorner == (int)landraster.yllcorner) {
								batchLog << ftype << " headers OK: " << fname << endl;
							}
							else {
								batchLog << gHeadersOfStr << ftype << " " << fname
									<< gHeadersNotMatchStr << endl;
								errors++;
							}
						}
						else {
							batchLog << gResolOfStr << ftype << " " << fname
								<< gResolNotMatchStr << endl;
							errors++;
						}
					}
					else {
						errors++;
						if (costraster.errors == -111)
							OpenError(ftype, fname);
						else
							FormatError(fname, costraster.errors);
					}
				}
				else {
					BatchError(filetype, line, 0, " "); errors++;
					batchLog << ftype << " must be NULL if transfer model is not SMS" << endl;
				}
			}

			// check dynamic landscape filename
			ftype = "DynLandFile";
			bLandFile >> intext;
			if (intext != "NULL") { // landscape is dynamic
				fname = indir + intext;
				batchLog << "Checking " << ftype << " " << fname << endl;
				bDynLandFile.open(fname.c_str());
				if (bDynLandFile.is_open()) {
					int something = CheckDynamicFile(indir, gNameCostFile);
					if (something < 0) {
						errors++;
					}
					bDynLandFile.close(); bDynLandFile.clear();
				}
				else {
					bDynLandFile.clear();
					errors++;
					OpenError(ftype, fname);
				}
			}

			// check initial distribution map filename
			ftype = "SpDistFile";
			bLandFile >> intext;
			if (intext == "NULL") {
				if (speciesdist) {
					BatchError(filetype, line, 0, " "); errors++;
					batchLog << ftype << " is required as SpeciesDist is 1 in Control file" << endl;
				}
			}
			else {
				if (speciesdist) {
					fname = indir + intext;
					spdistraster = CheckRasterFile(fname);
					if (spdistraster.ok) {
						if (spdistraster.cellsize == distresolution) {
							if (spdistraster.cellsize == landraster.cellsize) {
								// check that extent matches landscape extent
								if (spdistraster.ncols != landraster.ncols
									|| spdistraster.nrows != landraster.nrows) {
									batchLog << "*** Extent of " << ftype
										<< " does not match extent of LandscapeFile" << endl;
									errors++;
								}
								else {
									// check origins match
									if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
										&& (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
										batchLog << ftype << " headers OK: " << fname << endl;
									}
									else {
										batchLog << "*** Origin co-ordinates of " << ftype
											<< " do not match those of LandscapeFile" << endl;
										errors++;
									}
								}
							}
							else { // not able to check extents match
								// check origins match
								if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
									&& (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
									batchLog << ftype << " headers OK: " << fname << endl;
								}
								else {
									batchLog << "*** Origin co-ordinates of " << ftype
										<< " do not match those of LandscapeFile" << endl;
									errors++;
								}
							}
						}
						else {
							batchLog << "*** Resolution of " << ftype << " " << fname
								<< " does not match DistResolution in Control file" << endl;
							errors++;
						}
					}
					else {
						errors++;
						if (spdistraster.errors == -111)
							OpenError(ftype, fname);
						else
							FormatError(fname, spdistraster.errors);
					}
				}
			}

			totlines++; line++;
			// read first field on next line
			inint = -98765;
			bLandFile >> inint;
			//		batchlog << "ParseLandFile(): first item of next line = " << inint << endl;
		} // end of while loop
		landlist.clear();
	} // end of real landscape
	else {
		if (landtype == 9) { // artificial landscape
			int fractal, type, Xdim, Ydim;
			float minhab, maxhab;
			// Parse header line;
			bLandFile >> header; if (header != "LandNum") errors++;
			bLandFile >> header; if (header != "Fractal") errors++;
			bLandFile >> header; if (header != "Type") errors++;
			bLandFile >> header; if (header != "Xdim") errors++;
			bLandFile >> header; if (header != "Ydim") errors++;
			bLandFile >> header; if (header != "MinHab") errors++;
			bLandFile >> header; if (header != "MaxHab") errors++;
			bLandFile >> header; if (header != "Psuit") errors++;
			bLandFile >> header; if (header != "H") errors++;
			if (errors > 0) {
				FormatError(filetype, 0);
				batchLog << "*** Ensure format is correct for artificial landscape" << endl;
				return -111;
			}
			// Parse data lines
			line = 1;
			inint = -98765;
			bLandFile >> inint;
			while (inint != -98765) {
				for (j = 0; j < (int)landlist.size(); j++) {
					if (inint < 1 || inint == landlist[j]) {
						BatchError(filetype, line, 666, "LandNum"); j = (int)landlist.size() + 1; errors++;
					}
				}
				//		batchlog << "ParseLandFile(): Adding landscape no. " << inint
				//			<< " to landscape list" << endl;
				landlist.push_back(inint);
				bLandFile >> fractal;
				if (fractal < 0 || fractal > 1) {
					BatchError(filetype, line, 1, "Fractal"); errors++;
				}
				bLandFile >> type;
				if (type < 0 || type > 1) {
					BatchError(filetype, line, 1, "Type"); errors++;
				}
				bLandFile >> Xdim >> Ydim;
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
						batchLog << "Y dimension may not be less than X dimension" << endl; errors++;
					}
					if ((Xdim > 2 && power2check(Xdim - 1) != 1)
						|| (Ydim > 2 && power2check(Ydim - 1) != 1)) {
						BatchError(filetype, line, 0, " ");
						batchLog << "X and Y dimensions must be a power of 2 plus 1" << endl; errors++;
					}
				}
				bLandFile >> minhab >> maxhab;
				if (type == 1) { // continuous landscape
					if (minhab <= 0.0 || minhab >= 100.0) {
						BatchError(filetype, line, 100, "MinHab"); errors++;
					}
					if (maxhab <= 0.0 || maxhab > 100.0) {
						BatchError(filetype, line, 100, "MaxHab"); errors++;
					}
					if (maxhab <= minhab) {
						BatchError(filetype, line, 0, " ");
						batchLog << "MaxHab must exceed MinHab" << endl; errors++;
					}
				}
				bLandFile >> infloat;
				if (infloat < 0.0 || infloat > 1.0) {
					BatchError(filetype, line, 20, "Psuit"); errors++;
				}
				bLandFile >> infloat;
				if (fractal == 1) {
					if (infloat <= 0.0 || infloat >= 1.0) {
						BatchError(filetype, line, 20, "H"); errors++;
					}
				}
				totlines++; line++;
				// read first field on next line
				inint = -98765;
				bLandFile >> inint;
			} // end of while loop
		} // end of artificial landscape
		else { // ERROR condition which should not occur
			batchLog << "*** Critical error in land file. "
				<< "Invalid value of landscape type passed to function ParseLandFile()" << endl;
			errors++;
		}
	}
	if (!bLandFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return totlines;

}

int CheckDynamicFile(string indir, string costfile) {
#if RSDEBUG
	DEBUGLOG << "ParseDynamicFile(): costfile=" << costfile << endl;
#endif
	string header, filename, fname, ftype, intext;
	int change, prevchange, year, prevyear = 0;
	rasterdata landchgraster, patchchgraster, costchgraster;
	int errors = 0;
	string filetype = "DynLandFile";
	//int totlines = 0;

	bDynLandFile >> header; if (header != "Change") errors++;
	bDynLandFile >> header; if (header != "Year") errors++;
	bDynLandFile >> header; if (header != "LandChangeFile") errors++;
	bDynLandFile >> header; if (header != "PatchChangeFile") errors++;
	bDynLandFile >> header; if (header != "CostChangeFile") errors++;

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	change = -98765;
	bDynLandFile >> change; // first change number
	if (change != 1) {
		batchLog << "*** Error in DynLandFile - first change number must be 1" << endl;
		errors++;
	}
	else {
		prevchange = change;
	}
	while (change != -98765) {

		bDynLandFile >> year; if (year <= 0) { BatchError(filetype, line, 10, "Year"); errors++; }
		if (line > 1) {
			if (year <= prevyear) {
				BatchError(filetype, line, 1, "Year", "previous Year"); errors++;
			}
		}
		prevyear = year;

		// check landscape filename
		ftype = "LandChangeFile";
		bDynLandFile >> intext;
		//batchlog << "***** indir=" << indir << " intext=" << intext << endl;
		fname = indir + intext;
		landchgraster = CheckRasterFile(fname);
		if (landchgraster.ok) {
			if (landchgraster.cellsize == resolution)
				if (landchgraster.ncols == landraster.ncols
					&& landchgraster.nrows == landraster.nrows
					&& landchgraster.cellsize == landraster.cellsize
					&& (int)landchgraster.xllcorner == (int)landraster.xllcorner
					&& (int)landchgraster.yllcorner == (int)landraster.yllcorner) {
					batchLog << ftype << " headers OK: " << fname << endl;
				}
				else {
					batchLog << gHeadersOfStr << ftype << " " << fname
						<< gHeadersNotMatchStr << endl;
					errors++;
				}
			else {
				errors++;
				batchLog << gResolOfStr << ftype << " " << fname << gResolNotMatchStr << endl;
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
		bDynLandFile >> intext;
		if (intext == "NULL") {
			if (patchmodel) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLog << ftype << gPatchReqdStr << endl;
			}
		}
		else {
			if (patchmodel) {
				fname = indir + intext;
				patchchgraster = CheckRasterFile(fname);
				if (patchchgraster.ok) {
					if (patchchgraster.cellsize == resolution) {
						if (patchchgraster.ncols == landraster.ncols
							&& patchchgraster.nrows == landraster.nrows
							&& patchchgraster.cellsize == landraster.cellsize
							&& (int)patchchgraster.xllcorner == (int)landraster.xllcorner
							&& (int)patchchgraster.yllcorner == (int)landraster.yllcorner) {
							batchLog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchLog << gHeadersOfStr << ftype << " " << fname
								<< gHeadersNotMatchStr << endl;
							errors++;
						}
					}
					else {
						batchLog << gResolOfStr << ftype << " " << fname
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
		bDynLandFile >> intext;
		if (intext == "NULL") {
			if (costfile != "NULL") {
				BatchError(filetype, line, 0, " "); errors++;
				batchLog << ftype << " must be supplied " << endl;
			}
		}
		else {
			if (costfile == "NULL") {
				BatchError(filetype, line, 0, " "); errors++;
				batchLog << ftype << " must be NULL to match LandFile " << endl;
			}
			else {
				fname = indir + intext;
				costchgraster = CheckRasterFile(fname);
				if (costchgraster.ok) {
					if (costchgraster.cellsize == resolution) {
						if (costchgraster.ncols == landraster.ncols
							&& costchgraster.nrows == landraster.nrows
							&& costchgraster.cellsize == landraster.cellsize
							&& (int)costchgraster.xllcorner == (int)landraster.xllcorner
							&& (int)costchgraster.yllcorner == (int)landraster.yllcorner) {
							batchLog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchLog << gHeadersOfStr << ftype << " " << fname
								<< gHeadersNotMatchStr << endl;
							errors++;
						}
					}
					else {
						batchLog << gResolOfStr << ftype << " " << fname
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
		bDynLandFile >> change;
		if (bDynLandFile.eof()) {
			change = -98765;
		}
		else { // check for valid change number
			if (change != prevchange + 1) {
				BatchError(filetype, line, 0, " ");
				batchLog << "Change numbers must be sequential integers" << endl;
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
	bStageStructFile >> header; if (header != "Simulation") errors++;
	bStageStructFile >> header; if (header != "PostDestructn") errors++;
	bStageStructFile >> header; if (header != "PRep") errors++;
	bStageStructFile >> header; if (header != "RepInterval") errors++;
	bStageStructFile >> header; if (header != "MaxAge") errors++;
	bStageStructFile >> header; if (header != "TransMatrixFile") errors++;
	bStageStructFile >> header; if (header != "SurvSched") errors++;
	bStageStructFile >> header; if (header != "FecDensDep") errors++;
	bStageStructFile >> header; if (header != "FecStageWts") errors++;
	bStageStructFile >> header; if (header != "FecStageWtsFile") errors++;
	bStageStructFile >> header; if (header != "DevDensDep") errors++;
	bStageStructFile >> header; if (header != "DevDensCoeff") errors++;
	bStageStructFile >> header; if (header != "DevStageWts") errors++;
	bStageStructFile >> header; if (header != "DevStageWtsFile") errors++;
	bStageStructFile >> header; if (header != "SurvDensDep") errors++;
	bStageStructFile >> header; if (header != "SurvDensCoeff") errors++;
	bStageStructFile >> header; if (header != "SurvStageWts") errors++;
	bStageStructFile >> header; if (header != "SurvStageWtsFile") errors++;
	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	inint = -98765;
	bStageStructFile >> inint;
	// first simulation number must match first one in parameterFile
	if (inint != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	prevsimul = inint;
	while (inint != -98765) {
		simuls++;
		bStageStructFile >> inint;
		if (inint < 0 || inint > 1) { BatchError(filetype, line, 1, "PostDestructn"); errors++; }
		bStageStructFile >> infloat;
		if (infloat <= 0 || infloat > 1.0) { BatchError(filetype, line, 20, "PRep"); errors++; }
		bStageStructFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "RepInterval"); errors++; }
		bStageStructFile >> inint;
		if (inint < 2) { BatchError(filetype, line, 12, "MaxAge"); errors++; }

		bStageStructFile >> filename;
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
				batchLog << "*** " << ftype2 << " is compulsory for stage-structured model" << endl;
				errors++;
			}
			else {
				fname = indir + filename;
				batchLog << "Checking " << ftype2 << " " << fname << endl;
				bTransMatrix.open(fname.c_str());
				if (bTransMatrix.is_open()) {
					err = CheckTransitionFile(stages, sexesDem);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bTransMatrix.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (bTransMatrix.is_open()) bTransMatrix.close();
				bTransMatrix.clear();
			}
		}
		transfiles.push_back(filename);

		bStageStructFile >> inint;
		if (inint < 0 || inint > 2) { BatchError(filetype, line, 2, "SurvSched"); errors++; }
		bStageStructFile >> fecdensdep;
		if (fecdensdep < 0 || fecdensdep > 1)
		{
			BatchError(filetype, line, 1, "FecDensDep"); errors++; fecdensdep = 1;
		}
		bStageStructFile >> fecstagewts;
		if (fecdensdep) {
			if (fecstagewts < 0 || fecstagewts > 1)
			{
				BatchError(filetype, line, 1, "FecStageWts"); errors++; fecstagewts = 1;
			}
		}
		else {
			if (fecstagewts != 0) {
				BatchError(filetype, line, 0, " ");
				batchLog << "FecStageWts must be 0 if FecDensDep is 0" << endl; errors++;
				errors++; fecstagewts = 1;
			}
		}

		// fecundity stage weights file - optional
		ftype2 = "FecStageWtsFile";
		bStageStructFile >> filename;
		if (filename == "NULL") {
			if (fecstagewts) {
				BatchError(filetype, line, 0, " ");
				batchLog << ftype2 << " is compulsory unless FecStageWts is 0" << endl;
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
				batchLog << "Checking " << ftype2 << " " << fname << endl;
				bStageWeightsFile.open(fname.c_str());
				if (bStageWeightsFile.is_open()) {
					err = CheckWeightsFile(ftype2);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bStageWeightsFile.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (bStageWeightsFile.is_open()) bStageWeightsFile.close();
				bStageWeightsFile.clear();
			}
			wtsfiles.push_back(filename);
		}

		bStageStructFile >> devdensdep;
		if (devdensdep < 0 || devdensdep > 1)
		{
			BatchError(filetype, line, 1, "DevDensDep"); errors++; devdensdep = 1;
		}
		bStageStructFile >> infloat >> devstagewts;
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
				batchLog << "DevStageWts must be 0 if DevDensDep is 0" << endl; errors++;
				errors++; devstagewts = 1;
			}
		}

		// development stage weights file - optional
		ftype2 = "DevStageWtsFile";
		bStageStructFile >> filename;
		if (filename == "NULL") {
			if (devstagewts) {
				BatchError(filetype, line, 0, " ");
				batchLog << ftype2 << " is compulsory unless DevStageWts is 0" << endl;
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
				batchLog << "Checking " << ftype2 << " " << fname << endl;
				bStageWeightsFile.open(fname.c_str());
				if (bStageWeightsFile.is_open()) {
					err = CheckWeightsFile(ftype2);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bStageWeightsFile.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (bStageWeightsFile.is_open()) bStageWeightsFile.close();
				bStageWeightsFile.clear();
			}
			wtsfiles.push_back(filename);
		}

		bStageStructFile >> survdensdep;
		if (survdensdep < 0 || survdensdep > 1)
		{
			BatchError(filetype, line, 1, "SurvDensDep"); errors++; survdensdep = 1;
		}
		bStageStructFile >> infloat >> survstagewts;
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
				batchLog << "SurvStageWts must be 0 if SurvDensDep is 0" << endl; errors++;
				errors++; survstagewts = 1;
			}
		}

		// survival stage weights file - optional
		ftype2 = "SurvStageWtsFile";
		bStageStructFile >> filename;
		if (filename == "NULL") {
			if (survstagewts) {
				BatchError(filetype, line, 0, " ");
				batchLog << ftype2 << " is compulsory unless SurvStageWts is 0" << endl;
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
				batchLog << "Checking " << ftype2 << " " << fname << endl;
				bStageWeightsFile.open(fname.c_str());
				if (bStageWeightsFile.is_open()) {
					err = CheckWeightsFile(ftype2);
					if (err == 0) FileHeadersOK(ftype2); else errors++;
					bStageWeightsFile.close();
				}
				else {
					OpenError(ftype2, fname); errors++;
				}
				if (bStageWeightsFile.is_open()) bStageWeightsFile.close();
				bStageWeightsFile.clear();
			}
			wtsfiles.push_back(filename);
		}

		// read next simulation
		line++;
		inint = -98765;
		bStageStructFile >> inint;
		if (bStageStructFile.eof()) {
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
	if (!bStageStructFile.eof()) {
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
	string header, hhh;
	int i, j, stage, sex, line, minage;
	//int prevminage;
	float infloat;
	int errors = 0;
	string filetype = "TransMatrixFile";

	// check header records
	bTransMatrix >> header; if (header != "Transition") errors++;
	for (i = 0; i < nstages; i++) {
		for (j = 0; j < nsexesDem; j++) {
			bTransMatrix >> header;
			if (nsexesDem == 1) hhh = to_string(i);
			else {
				if (j == 0) hhh = to_string(i) + "m"; else hhh = to_string(i) + "f";
			}
			if (header != hhh) errors++;
			//		batchlog << "i = " << i << " j = " << j << " hhh = " << hhh << " header = " << header
			//			<< " errors = " << errors << endl;
		}
	}
	bTransMatrix >> header; if (header != "MinAge") errors++;

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// check matrix, including row headers

	// single row for juveniles
	line = 1;
	bTransMatrix >> header;
	if (header != "0") {
		BatchError(filetype, line, 0, " ");
		batchLog << "Invalid row header" << endl; errors++;
	}
	float totfecundity = 0.0;
	for (i = 0; i < nstages; i++) {
		for (j = 0; j < nsexesDem; j++) {
			bTransMatrix >> infloat;
			if (i > 0) {
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
	bTransMatrix >> minage;
	//prevminage = minage;
	//				batchlog << "MINAGE = " << minage << endl;
	if (minage != 0) {
		BatchError(filetype, line, 0, " ");
		batchLog << "MinAge must be zero for juvenile stage" << endl; errors++;
	}

	// one row for each stage/sex combination
	//				batchlog << "HEADER = " << header << endl;
	for (stage = 1; stage < nstages; stage++) {
		for (sex = 0; sex < nsexesDem; sex++) {
			line++;
			// row header
			bTransMatrix >> header;
			if (nsexesDem == 1) hhh = to_string(stage);
			else {
				if (sex == 0) hhh = to_string(stage) + "m"; else hhh = to_string(stage) + "f";
			}
			if (header != hhh) {
				BatchError(filetype, line, 0, " ");
				batchLog << "Invalid row header" << endl; errors++;
			}
			for (i = 0; i < nstages; i++) {
				for (j = 0; j < nsexesDem; j++) {
					bTransMatrix >> infloat;
					//				batchlog << "TRANS PROB = " << infloat << endl;
					if (infloat < 0.0 || infloat > 1) {
						BatchError(filetype, line, 20, "Transition probability"); errors++;
					}
				}
			}
			//		 prevminage = minage;
			bTransMatrix >> minage;
			//				batchlog << "MINAGE = " << minage << endl;
			if (stage == 1 && minage != 0) {
				BatchError(filetype, line, 0, " ");
				batchLog << "MinAge must be zero for stage 1" << endl; errors++;
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
	bTransMatrix >> header;

	if (!bTransMatrix.eof()) {
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
	bStageWeightsFile >> header; if (header != "StageWts") errors++;
	for (i = 0; i < stages; i++) {
		for (j = 0; j < sexesDem; j++) {
			bStageWeightsFile >> header;
			if (sexesDem == 1) hhh = to_string(i);
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
	for (stage = 0; stage < stages; stage++) {
		for (sex = 0; sex < sexesDem; sex++) {
			line++;
			// row header
			bStageWeightsFile >> header;
			if (sexesDem == 1) hhh = to_string(stage);
			else {
				if (sex == 0) hhh = to_string(stage) + "m"; else hhh = to_string(stage) + "f";
			}
			if (header != hhh) {
				BatchError(filetype, line, 0, " ");
				batchLog << "Invalid row header" << endl; errors++;
			}
			for (i = 0; i < stages; i++) {
				for (j = 0; j < sexesDem; j++) {
					bStageWeightsFile >> infloat;
					// NOTE - any real number is acceptable - no check required
				}
			}
		}
	}
	// final read should hit EOF
	bStageWeightsFile >> header;

	if (!bStageWeightsFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	return errors;

}

//---------------------------------------------------------------------------
int CheckEmigFile(void)
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
	bEmigrationFile >> header; if (header != "Simulation") nbErrors++;
	bEmigrationFile >> header; if (header != "DensDep") nbErrors++;
	bEmigrationFile >> header; if (header != "UseFullKern") nbErrors++;
	bEmigrationFile >> header; if (header != "StageDep") nbErrors++;
	bEmigrationFile >> header; if (header != "SexDep") nbErrors++;
	bEmigrationFile >> header; if (header != "IndVar") nbErrors++;
	bEmigrationFile >> header; if (header != "EmigStage") nbErrors++;
	bEmigrationFile >> header; if (header != "Stage") nbErrors++;
	bEmigrationFile >> header; if (header != "Sex") nbErrors++;
	bEmigrationFile >> header; if (header != "EP") nbErrors++;
	bEmigrationFile >> header; if (header != "D0") nbErrors++;
	bEmigrationFile >> header; if (header != "alpha") nbErrors++;
	bEmigrationFile >> header; if (header != "beta") nbErrors++;

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
	bEmigrationFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichInputFile, lineNb, 111, "Simulation"); 
		nbErrors++;
		readNextLine = false;
	}

	while (readNextLine) {
		// read and validate columns relating to stage and sex-dependency and to IIV
		bEmigrationFile >> inDensDep >> inUseFullKern >> inStgDep >> inSexDep;
		bEmigrationFile >> inIndVar >> inEmigStg >> inStage >> inSex;
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
			gDispTraitOpt.isEmigDensDep = (inDensDep == 1);
		}

		// validate individual variation
		if (inIndVar != 0 && inIndVar != 1) {
			BatchError(whichInputFile, lineNb, 1, "IndVar");
			nbErrors++;
		}
		else {
			gDispTraitOpt.isEmigIndVar = (inIndVar == 1);
		}


		// validate use full kernel
		if (inUseFullKern != 0 && inUseFullKern != 1) {
			BatchError(whichInputFile, lineNb, 1, "UseFullKern"); 
			nbErrors++;
		}
		if (inDensDep == 1 && inUseFullKern != 0) {
				BatchError(whichInputFile, lineNb, 0, "UseFullKern"); 
				nbErrors++;
				batchLog << "UseFullKern must be 0 if there is density-dependent emigration" << endl;
		}
		// validate emigration stage
		if (stagestruct && !inStgDep && inIndVar == 1
			&& inStage == 0 && inSex == 0
			&& (inEmigStg < 0 || inEmigStg >= stages)) {
			BatchError(whichInputFile, lineNb, 0, "EmigStage");
			nbErrors++;
			batchLog << "EmigStage must be from 0 to " << to_string(stages - 1) << endl;
		}
		if (inSexDep != 0 && inSexDep != 1) {
			BatchError(whichInputFile, lineNb, 1, "SexDep");
			nbErrors++;
		} 
		else {
			gDispTraitOpt.isEmigSexDep = (inSexDep == 1);
		}

		if (inStage == 0 && inSex == 0) { // first line of a simulation
			// record whether density dependence and individual variability are applied
			isDensDep = (inDensDep == 1);
			isIndVar = (inIndVar == 1);
		}

		// read remaining columns of the current record
		bEmigrationFile >> inEP >> inD0 >> inAlpha >> inBeta;

		if (isDensDep) {
			if (inEP != gEmptyVal) {
				batchLog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is enabled EP must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inD0 < 0.0 || inD0 > 1.0) {
				BatchError(whichInputFile, lineNb, 20, "D0"); 
				nbErrors++;
			}
			// if (inAlpha < 0.0) {
			//	BatchError(whichInputFile, lineNb, 10, "alpha");
			//	nbErrors++;
			// }
			if (inBeta < 0.0) {
				BatchError(whichInputFile, lineNb, 10, "beta");
				nbErrors++;
			}
		}
		else { // !densdepset
			if (inEP < 0.0 || inEP > 1.0) {
				BatchError(whichInputFile, lineNb, 20, "EP"); 
				nbErrors++;
			}
			if (inD0 != gEmptyVal) {
				batchLog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled D0 must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inAlpha != gEmptyVal) {
				batchLog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled alpha must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inBeta != gEmptyVal) {
				batchLog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled beta must be " << gEmptyVal << endl;
				nbErrors++;
			}
		}

		// read next simulation
		lineNb++;
		const int errSimNb = -98765;
		simNb = errSimNb;
		bEmigrationFile >> simNb;
		if (simNb == errSimNb || bEmigrationFile.eof()) readNextLine = false;
	} // end of while loop

	// check for correct number of lines for previous simulation
	if (currentLine.simLines != currentLine.reqdSimLines) {
		BatchError(whichInputFile, lineNb, 0, " "); 
		nbErrors++;
		batchLog << gNbLinesStr << currentLine.simNb
			<< gShouldBeStr << currentLine.reqdSimLines << endl;
	}
	if (!bEmigrationFile.eof()) {
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
	bTransferFile >> header; 
	if (header != "Simulation") errors++;

	switch (gTransferType) {

	case 0: { // negative exponential dispersal kernel
		batchLog << "Checking dispersal kernel format file" << endl;
		bTransferFile >> header; if (header != "StageDep") errors++;
		bTransferFile >> header; if (header != "SexDep") errors++;
		bTransferFile >> header; if (header != "KernelType") errors++;
		bTransferFile >> header; if (header != "DistMort") errors++;
		bTransferFile >> header; if (header != "IndVar") errors++;
		bTransferFile >> header; if (header != "Stage") errors++;
		bTransferFile >> header; if (header != "Sex") errors++;
		bTransferFile >> header; if (header != "meanDistI") errors++;
		bTransferFile >> header; if (header != "meanDistII") errors++;
		bTransferFile >> header; if (header != "ProbKernelI") errors++;
		bTransferFile >> header; if (header != "MortProb") errors++;
		bTransferFile >> header; if (header != "Slope") errors++;
		bTransferFile >> header; if (header != "InflPoint") errors++;
		break;
	} // end of negative exponential dispersal kernel

	case 1: { // SMS
		batchLog << "Checking SMS format file ";
		bTransferFile >> header; if (header != "IndVar") errors++;
		bTransferFile >> header; if (header != "PR") errors++;
		bTransferFile >> header; if (header != "PRMethod") errors++;
		bTransferFile >> header; if (header != "DP") errors++;
		bTransferFile >> header; if (header != "MemSize") errors++;
		bTransferFile >> header; if (header != "GB") errors++;
		bTransferFile >> header; if (header != "GoalType") errors++;
		bTransferFile >> header; if (header != "AlphaDB") errors++;
		bTransferFile >> header; if (header != "BetaDB") errors++;
		bTransferFile >> header; if (header != "StraightenPath") errors++;
		bTransferFile >> header; if (header != "SMtype") errors++;
		bTransferFile >> header; if (header != "SMconst") errors++;
		switch (landtype) {
		case 0: { // raster map with unique habitat codes
			batchLog << "for LandType = 0" << endl;
			for (i = 0; i < maxNhab; i++) {
				colheader = "MortHab" + to_string(i + 1);
				bTransferFile >> header; if (header != colheader) morthaberrors++;
			}
			for (i = 0; i < maxNhab; i++) {
				colheader = "CostHab" + to_string(i + 1);
				bTransferFile >> header; if (header != colheader) costerrors++;
			}
			break;
		} // end of raster map with unique habitat codes
		case 2: { // raster map with habitat quality
			batchLog << "for LandType = 2" << endl;
			break;
		} // end of raster map with habitat quality
		case 9: { // artificial landscape
			batchLog << "for LandType = 9" << endl;
			bTransferFile >> header; if (header != "MortHabitat") errors++;
			bTransferFile >> header; if (header != "MortMatrix") errors++;
			bTransferFile >> header; if (header != "CostHabitat") errors++;
			bTransferFile >> header; if (header != "CostMatrix") errors++;
			break;
		} // end of artificial landscape
		} // end of switch (landtype)
		break;
	} // end of SMS

	case 2: { // CRW
		batchLog << "Checking CRW format file" << endl;
		bTransferFile >> header; if (header != "IndVar") errors++;
		bTransferFile >> header; if (header != "SL") errors++;
		bTransferFile >> header; if (header != "Rho") errors++;
		bTransferFile >> header; if (header != "StraightenPath") errors++;
		bTransferFile >> header; if (header != "SMtype") errors++;
		bTransferFile >> header; if (header != "SMconst") errors++;
		if (landtype == 0) {
			for (i = 0; i < maxNhab; i++) {
				colheader = "MortHab" + to_string(i + 1);
				bTransferFile >> header; if (header != colheader) morthaberrors++;
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
	bTransferFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichFile, whichLine, 111, "Simulation"); errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		switch (gTransferType) {

		case 0: { // negative exponential dispersal kernel
			// read and validate columns relating to stage and sex-dependency and to IIV
			bTransferFile >> inStageDep >> inSexDep >> inKernelType >> inDistMort;
			bTransferFile >> inIndVar >> inStage >> inSex;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, inStageDep, inSexDep, inStage, inSex, inIndVar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			// validate kernel type
			if (inKernelType != 0 && inKernelType != 1) {
				BatchError(whichFile, whichLine, 1, "KernelType"); errors++;
			}
			else {
				gDispTraitOpt.usesTwoKernels = (inKernelType == 1);
			}
			// validate mortality
			if (inDistMort != 0 && inDistMort != 1) {
				BatchError(whichFile, whichLine, 1, "DistMort"); errors++;
			}
			// read remaining columns of the current record
			bTransferFile >> meanDistI >> meanDistII >> ProbKernelI;
			bTransferFile >> mortProb >> slope >> inflPoint;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				gDispTraitOpt.isKernTransfIndVar = (inIndVar == 1);
			}

			if (inSexDep != 0 && inSexDep != 1) {
				BatchError(whichFile, whichLine, 1, "SexDep"); errors++;
			}
			else {
				gDispTraitOpt.isKernTransfSexDep = (inSexDep == 1);
			}

			// validate mortality
			if (inDistMort != 0 && inDistMort != 1) {
				BatchError(whichFile, whichLine, 1, "DistMort"); errors++;
			}

			if (!inIndVar) {
				if (meanDistI < resolution) {
					// NOTE - should also check whether emigration prob is constant and equal to 1
					//but checks across diffferent input files are not yet implemented
					BatchError(whichFile, whichLine, 2, "meanDistI", "Resolution"); errors++;
				}
				if (inKernelType != 0) {
					if (meanDistII < resolution) {
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
			bTransferFile >> inIndVar;
			bTransferFile >> inPerceptualRange >> inPercRangeMethod >> inDirPersistence;
			bTransferFile >> inMemSize >> inGoalBias >> inGoalType >> inAlphaDispBias >> inBetaDispBias;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, 0, 0, 0, 0, 0, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				gDispTraitOpt.isSMSTransfIndVar = (inIndVar == 1);
			}

			// validate SMS movement parameters
			if (inPerceptualRange < 1) {
				BatchError(whichFile, whichLine, 11, "PR"); errors++;
			}
			if (inPercRangeMethod < 1 || inPercRangeMethod > 3) {
				BatchError(whichFile, whichLine, 33, "PRmethod"); errors++;
			}
			if (!inIndVar && inDirPersistence < 1.0) {
				BatchError(whichFile, whichLine, 11, "DP"); errors++;
			}
			if (inMemSize < 1 || inMemSize > 14) {
				BatchError(whichFile, whichLine, 0, "MemSize"); errors++;
				batchLog << "MemSize must be from 1 to 14" << endl;
			}
			if (!inIndVar && inGoalBias < 1.0) {
				BatchError(whichFile, whichLine, 11, "GB"); errors++;
			}
			if (inGoalType != 0 && inGoalType != 2) {
				BatchError(whichFile, whichLine, 2, "GoalType"); errors++;
			}
			else {
				gDispTraitOpt.usesSMSGoalBias = (inGoalType == 2);
			}

			if (!inIndVar && inGoalType == 2) { // dispersal bias
				if (inAlphaDispBias <= 0.0) {
					BatchError(whichFile, whichLine, 10, "AlphaDB"); errors++;
				}
				if (inBetaDispBias <= 0.0) {
					BatchError(whichFile, whichLine, 10, "BetaDB"); errors++;
				}
			}
			bTransferFile >> inStraightenPath >> inSMType >> inSMConst;
			if (inStraightenPath != 0 && inStraightenPath != 1) {
				BatchError(whichFile, whichLine, 1, "StraightenPath"); errors++;
			}
			if (landtype == 2) // habitat quality landscape 
			{ // must have constant mortality
				if (inSMType != 0) {
					BatchError(whichFile, whichLine, 0, " "); errors++;
					batchLog << "SMtype must be 0 for LandType 2" << endl;
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
			switch (landtype) {

			case 0: { // raster map with unique habitat codes
				//			batchlog << "for LandType = 0" << endl;
				for (i = 0; i < maxNhab; i++) {
					bTransferFile >> morthab;
					if (inSMType == 1)
					{
						if (morthab < 0.0 || morthab >= 1.0) {
							colheader = "MortHab" + to_string(i + 1);
							BatchError(whichFile, whichLine, 20, colheader); errors++;
						}
					}
				}
				for (i = 0; i < maxNhab; i++) {
					bTransferFile >> costhab;
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
				bTransferFile >> morthab >> mortmatrix;
				bTransferFile >> costhab >> costmatrix;
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
			bTransferFile >> inIndVar >> inStepLength >> inStepCorr >> inStraightenPath >> inSMType >> inSMConst;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, 0, 0, 0, 0, inIndVar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (inIndVar != 0 && inIndVar != 1) {
				BatchError(whichFile, whichLine, 1, "IndVar"); errors++;
			}
			else {
				gDispTraitOpt.isCRWTransfIndVar = (inIndVar == 1);
			}

			if (inStepLength <= 0.0) {
				BatchError(whichFile, whichLine, 10, "SL"); errors++;
			}
			if (inStepCorr <= 0.0 || inStepCorr >= 1.0) {
				BatchError(whichFile, whichLine, 20, "Rho"); errors++;
			}
			if (inStraightenPath != 0 && inStraightenPath != 1) {
				BatchError(whichFile, whichLine, 1, "StraightenPath"); errors++;
			}
			if (landtype == 0) { // real landscape with habitat types
				if (inSMType != 0 && inSMType != 1) {
					BatchError(whichFile, whichLine, 1, "SMtype"); errors++;
				}
				if (!inSMType) {
					if (inSMConst < 0.0 || inSMConst >= 1.0) {
						BatchError(whichFile, whichLine, 20, "SMconst"); errors++;
					}
				}
				for (int i = 0; i < maxNhab; i++) {
					bTransferFile >> morthab;
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
					batchLog << "SMtype must be 0 for LandType 2 or 9" << endl;
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
		bTransferFile >> simNb;
		if (bTransferFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (gTransferType == 0 // no. of lines checked for dispersal kernel transfer method only
		&& current.simLines != current.reqdSimLines) {
		BatchError(whichFile, whichLine, 0, " "); errors++;
		batchLog << gNbLinesStr << current.simNb
			<< gShouldBeStr << current.reqdSimLines << endl;
	}
	if (!bTransferFile.eof()) {
		EOFerror(whichFile);
		errors++;
	}
	costsfiles.clear();

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int CheckSettleFile(void)
{
	string header;
	int simNb, inStageDep, inSexDep, inStage, inSex, inSettleType;
	int inDensDep, inIndVar, inFindMate, inMinSteps, inMaxSteps, inMaxStepsYear;
	float inS0, inAlphaS, inBetaS;
	int nbErrors = 0;
	int nbSims = 0;
	string whichFile = "SettlementFile";

	// Parse header line;
	bSettlementFile >> header; if (header != "Simulation") nbErrors++;
	bSettlementFile >> header; if (header != "StageDep") nbErrors++;
	bSettlementFile >> header; if (header != "SexDep") nbErrors++;
	bSettlementFile >> header; if (header != "Stage") nbErrors++;
	bSettlementFile >> header; if (header != "Sex") nbErrors++;
	if (gTransferType == 0)
	{ // dispersal kernel
		bSettlementFile >> header; if (header != "SettleType") nbErrors++;
		bSettlementFile >> header; if (header != "FindMate") nbErrors++;
	}
	else { // movement method
		bSettlementFile >> header; if (header != "DensDep") nbErrors++;
		bSettlementFile >> header; if (header != "IndVar") nbErrors++;
		bSettlementFile >> header; if (header != "FindMate") nbErrors++;
		bSettlementFile >> header; if (header != "MinSteps") nbErrors++;
		bSettlementFile >> header; if (header != "MaxSteps") nbErrors++;
		bSettlementFile >> header; if (header != "MaxStepsYear") nbErrors++;
		bSettlementFile >> header; if (header != "S0") nbErrors++;
		bSettlementFile >> header; if (header != "AlphaS") nbErrors++;
		bSettlementFile >> header; if (header != "BetaS") nbErrors++;
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
	bSettlementFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichFile, whichLine, 111, "Simulation"); nbErrors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {
		if (gTransferType == 0)
		{ // dispersal kernel
			// read and validate columns relating to stage and sex-dependency (NB no IIV here)
			bSettlementFile >> inStageDep >> inSexDep >> inStage >> inSex >> inSettleType >> inFindMate;
			current = CheckStageSex(whichFile, whichLine, simNb, prev, inStageDep, inSexDep, inStage, inSex, 0, true, false);
			if (current.isNewSim) nbSims++;
			nbErrors += current.errors;
			prev = current;
			if (inSettleType < 0 || inSettleType > 3) {
				BatchError(whichFile, whichLine, 3, "SettleType"); nbErrors++;
			}
			if (!stagestruct && (inSettleType == 1 || inSettleType == 3)) {
				BatchError(whichFile, whichLine, 0, " "); nbErrors++;
				batchLog << "Invalid SettleType for a non-stage-structured population" << endl;
			}
			if (gNbSexesDisp > 1) {
				if (inFindMate < 0 || inFindMate > 1) {
					BatchError(whichFile, whichLine, 1, "FindMate"); nbErrors++;
				}
			}
		}
		else { // movement method
			// read and validate columns relating to stage and sex-dependency (IIV psossible)
			bSettlementFile >> inStageDep >> inSexDep >> inStage >> inSex >> inDensDep >> inIndVar >> inFindMate;
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
				batchLog << "IndVar must be 0 if DensDep is 0" << endl;
			}
			else {
				gDispTraitOpt.isSettIndVar = inIndVar == 1;
			}

			if (inSexDep != 0 && inSexDep != 1) {
				BatchError(whichFile, whichLine, 1, "SexDep");
				nbErrors++;
			}
			else {
				gDispTraitOpt.isSettSexDep = inSexDep == 1;
			}

			if (reproductn != 0 && gNbSexesDisp > 1) {
				if (inFindMate != 0 && inFindMate != 1) {
					BatchError(whichFile, whichLine, 1, "FindMate"); 
					nbErrors++;
				}
			}
			bSettlementFile >> inMinSteps >> inMaxSteps >> inMaxStepsYear;
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
			bSettlementFile >> inS0 >> inAlphaS >> inBetaS;

			if (inDensDep == 1) {

				if (inS0 <= 0.0 || inS0 > 1.0) {
					BatchError(whichFile, whichLine, 20, "S0"); 
					nbErrors++;
				}
				// NOTE: alphaS and betaS can take any value
			}
		}
		// read next simulation
		whichLine++;
		simNb = -98765;
		bSettlementFile >> simNb;
		if (bSettlementFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation

	if (current.simLines != current.reqdSimLines) {
		BatchError(whichFile, whichLine, 0, " "); nbErrors++;
		batchLog << gNbLinesStr << current.simNb
			<< gShouldBeStr << current.reqdSimLines << endl;
	}
	if (!bSettlementFile.eof()) {
		EOFerror(whichFile);
		nbErrors++;
	}

	if (nbErrors > 0) return -111;
	else return nbSims;

}

//---------------------------------------------------------------------------
int CheckTraitsFile(string indir, const bool& anyNeutralGenetics)
{
	string header, colheader;
	int simNb, nextLineSimNb;
	string filename, inTraitType, inSex, inInitDist, inInitParams,
		inDominanceDist, inDominanceParams, inIsInherited, inMutationDist, 
		inMutationParams, inPositions, inNbPositions, inExpressionType, inMutationRate;
	int nbErrors = 0;
	int nbSims = 0;
	int nbGenLoadTraits = 0;
	vector <string> archfiles;
	const string whichInputFile = "TraitsFile";
	vector <TraitType> allReadTraits;

	// Parse header line
	bTraitsFile >> header; if (header != "Simulation") nbErrors++;
	bTraitsFile >> header; if (header != "TraitType") nbErrors++;
	bTraitsFile >> header; if (header != "ExprSex") nbErrors++;
	bTraitsFile >> header; if (header != "Positions") nbErrors++;
	bTraitsFile >> header; if (header != "NbrOfPositions") nbErrors++;
	bTraitsFile >> header; if (header != "ExpressionType") nbErrors++;
	bTraitsFile >> header; if (header != "InitialDistribution") nbErrors++;
	bTraitsFile >> header; if (header != "InitialParameters") nbErrors++;
	bTraitsFile >> header; if (header != "DominanceDistribution") nbErrors++;
	bTraitsFile >> header; if (header != "DominanceParameters") nbErrors++;
	bTraitsFile >> header; if (header != "IsInherited") nbErrors++;
	bTraitsFile >> header; if (header != "MutationDistribution") nbErrors++;
	bTraitsFile >> header; if (header != "MutationParameters") nbErrors++;
	bTraitsFile >> header; if (header != "MutationRate") nbErrors++;

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

	bTraitsFile >> simNb;

	bool stopReading = (simNb == simNbNotRead);
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichInputFile, lineNb, 111, "Simulation"); 
		nbErrors++;
	}
	while (!stopReading) {
		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		bTraitsFile >> inTraitType >> inSex >> inPositions >> inNbPositions >> inExpressionType >> inInitDist >> inInitParams
			>> inDominanceDist >> inDominanceParams >> inIsInherited >> inMutationDist >> inMutationParams
			>> inMutationRate;

		current = CheckStageSex(whichInputFile, lineNb, simNb, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) nbSims++;
		nbErrors += current.errors;
		prev = current;

		////  Validate parameters

		// Check sex is valid
		sex_t sex = stringToSex(inSex);
		if (sex == sex_t::INVALID_SEX) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << inSex << " is invalid: ExprSex must be either female, male, or # (if not applicable)." << endl;
			nbErrors++;
		}

		// Check trait type is legal
		TraitType tr = stringToTraitType(inTraitType);
		if (tr == TraitType::INVALID_TRAIT) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << inTraitType << " is not a valid TraitType." << endl;
			nbErrors++;
		}
		// Can trait be sex-dependent?
		const bool canBeSexDep = tr == E_D0 || tr == E_ALPHA || tr == E_BETA
			|| tr == S_S0 || tr == S_ALPHA || tr == S_BETA
			|| tr == KERNEL_MEANDIST_1 || tr == KERNEL_MEANDIST_2
			|| tr == KERNEL_PROBABILITY;
		if (!canBeSexDep && (sex == FEM || sex == MAL)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << inTraitType << " cannot be sex-dependent so ExprSex must be left blank (#)." << endl;
			nbErrors++;
		}
		if (sex != NA) // add sex to trait if present
			tr = addSexDepToTrait(tr, sex);

		// There can be up to 5 genetic load traits
		if (tr == GENETIC_LOAD) {
			nbGenLoadTraits++;
			if (nbGenLoadTraits > 5) {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "There cannot be more than 5 genetic load traits." << endl;
				nbErrors++;
			}
		}
		else if (traitExists(tr, allReadTraits)) {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "Trait " << to_string(tr) << " is supplied multiple times." << endl;
			nbErrors++;
		}
		allReadTraits.push_back(tr);

		// Check Positions and NbrOfPositions
		const regex patternPositions{ "^\"?(([0-9]+-)?[0-9]+,)*([0-9]+-)?[0-9]+\"?$" };
		bool isMatch = regex_search(inPositions, patternPositions);
		if (!isMatch && inPositions != "random") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "Positions must be either a comma-separated list of integer ranges, or random." << endl;
			nbErrors++;
		}
		if (inPositions == "random") {
			if (stoi(inNbPositions) <= 0) {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "NbrOfPositions must be a strictly positive integrer." << endl;
				nbErrors++;
			}
		}
		else if (inNbPositions != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "If Positions is not random NbrOfPositions must be blank (#)." << endl;
			nbErrors++;
		}

		// Check ExpressionType
		if (tr == NEUTRAL && inExpressionType != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "ExpressionType must be left blank (#) for the neutral trait." << endl;
			nbErrors++;
		}
		if (tr == GENETIC_LOAD && inExpressionType != "multiplicative") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "ExpressionType must be \"multiplicative\" for genetic load traits." << endl;
			nbErrors++;
		}
		const bool isDisp = tr != NEUTRAL && tr != GENETIC_LOAD && tr != INVALID_TRAIT;
		if (isDisp && inExpressionType != "additive" && inExpressionType != "average") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "ExpressionType must be \"additive\" or \"average\" for dispersal traits." << endl;
			nbErrors++;
		}

		// Check InitialDistribution
		if (tr == NEUTRAL && inInitDist != "#" && inInitDist != "uniform") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "InitialDistribution must be either uniform or left blank (#) for the neutral trait." << endl;
			nbErrors++;
		}
		if (tr == GENETIC_LOAD && inInitDist != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "InitialDistribution must be blank (#) for genetic load traits." << endl;
			nbErrors++;
		}
		if (isDisp && inInitDist != "normal" && inInitDist != "uniform") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "InitialDistribution must be either normal or uniform for dispersal traits." << endl;
			nbErrors++;
		}

		// Check InitialParameters
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
					batchLog << "For neutral trait with uniform initialisation, InitialParameters must have form max=int" << endl;
					nbErrors++;
				}
				else {
					const int maxVal = stoi(inInitParams.substr(4));
					if (maxVal > 255) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLog << "For neutral trait with uniform initialisation, max parameter must be between 0 and 255." << endl;
						nbErrors++;
					}
				}
			}
			// if not uniform then initDist must be blank, no params
			else if (inInitParams != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "For neutral trait with no initialisation, InitialParameters must be blank (#)" << endl;
				nbErrors++;
			}
		}
		if (tr == GENETIC_LOAD && inInitParams != "#") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "For genetic load traits, InitialParameters must be blank (#)" << endl;
			nbErrors++;
		}
		if (isDisp) {
			if (inInitDist == "uniform") {
				isMatch = regex_search(inInitParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For dispersal trait uniform initialisation, InitialParameters must have form min=float,max=float" << endl;
					nbErrors++;
				}
			}
			else if (inInitDist == "normal") {
				isMatch = regex_search(inInitParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For normal initialisation, InitialParameters must have form mean=float,sd=float" << endl;
					nbErrors++;
				}
			}
		}

		// Check DominanceDistribution and DominanceParameters
		if (tr == NEUTRAL) {
			if (inDominanceDist != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "DominanceDistribution must be left blank (#) for the neutral trait." << endl;
				nbErrors++;
			}
			if (inDominanceParams != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "DominanceParameters must be left blank (#) for the neutral trait." << endl;
				nbErrors++;
			}
		}
		if (isDisp) {
			if (inDominanceDist != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "DominanceDistribution must be left blank (#) for dispersal traits." << endl;
				nbErrors++;
			}
			if (inDominanceParams != "#") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "DominanceParameters must be left blank (#) for dispersal traits." << endl;
				nbErrors++;
			}
		}
		if (tr == GENETIC_LOAD) {
			if (inDominanceDist == "normal") {
				isMatch = regex_search(inDominanceParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a normal dominance distribution, DominanceParams must have form mean=float,sd=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "gamma") {
				isMatch = regex_search(inDominanceParams, patternParamsGamma);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a Gamma dominance distribution, DominanceParams must have form shape=float,scale=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "uniform") {
				isMatch = regex_search(inDominanceParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a uniform dominance distribution, DominanceParams must have form min=float,max=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "negExp") {
				isMatch = regex_search(inDominanceParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a negative exponential dominance distribution, DominanceParams must have form mean=float" << endl;
					nbErrors++;
				}
			}
			else if (inDominanceDist == "scaled") {
				isMatch = regex_search(inDominanceParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a scaled dominance distribution, DominanceParams must have form mean=float" << endl;
					nbErrors++;
				}
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "DominanceDistribution must be either normal, gamma, uniform, negExp or scaled for genetic load traits." << endl;
				nbErrors++;
			}
		}

		// Check isInherited and MutationRate
		if ((tr == NEUTRAL || tr == GENETIC_LOAD) && inIsInherited != "TRUE") {
			BatchError(whichInputFile, lineNb, 0, " ");
			batchLog << "isInherited must always be TRUE for neutral and genetic load traits." << endl;
			nbErrors++;
		}
		else if (isDisp) {
			if (inIsInherited != "TRUE" && inIsInherited != "FALSE") {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "IsInherited must be either TRUE or FALSEfor dispersal traits." << endl;
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
			batchLog << "If isInherited if off, mutationRate must be blank (#)." << endl;
			nbErrors++;
		}

		// Check MutationDistribution and MutationParameters
		if (tr == NEUTRAL) {
			if (inMutationDist == "KAM" || inMutationDist == "SSM") {
				isMatch = regex_search(inMutationParams, patternParamsNeutral);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a neutral trait, mutationParams must have form max=int." << endl;
					nbErrors++;
				}
				else {
					const int maxVal = stoi(inInitParams.substr(4));
					if (maxVal > 255) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLog << "For the neutral trait mutation max parameter must be between 0 and 255." << endl;
						nbErrors++;
					}
				}
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "For a neutral trait, mutationDistribution must be either KAM or SSM." << endl;
				nbErrors++;
			}
		}
		if (isDisp) {
			if (inIsInherited == "TRUE") {
				if (inMutationDist == "uniform") {
					isMatch = regex_search(inMutationParams, patternParamsUnif);
					if (!isMatch) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLog << "For a uniform distribution, mutationParams must have form min=float,max=float." << endl;
						nbErrors++;
					}
				}
				else if (inMutationDist == "normal") {
					isMatch = regex_search(inMutationParams, patternParamsNormal);
					if (!isMatch) {
						BatchError(whichInputFile, lineNb, 0, " ");
						batchLog << "For a normal distribution, mutationParams must have form mean=float,sd=float." << endl;
						nbErrors++;
					}
				}
				else {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For dispersal traits, mutationDistribution must be either uniform or normal" << endl;
					nbErrors++;
				}
			}
			else { // not inherited
				if (inMutationDist != "#" || inMutationParams != "#") {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "If isInherited is turned off, mutationDistribution and mutationParameters must be left blank (#)." << endl;
					nbErrors++;
				}
			}
		}
		if (tr == GENETIC_LOAD) {
			if (inMutationDist == "uniform") {
				isMatch = regex_search(inMutationParams, patternParamsUnif);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a uniform distribution, mutationParams must have form min=float,max=float." << endl;
					nbErrors++;
				}
			}
			else if (inMutationDist == "normal") {
				isMatch = regex_search(inMutationParams, patternParamsNormal);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a normal distribution, mutationParams must have form mean=float,sd=float." << endl;
					nbErrors++;
				}
			}
			else if (inMutationDist == "gamma") {
				isMatch = regex_search(inMutationParams, patternParamsGamma);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a Gamma distribution, mutationParams must have form shape=float,scale=float." << endl;
					nbErrors++;
				}
			}
			else if (inMutationDist == "negExp") {
				isMatch = regex_search(inMutationParams, patternParamsMean);
				if (!isMatch) {
					BatchError(whichInputFile, lineNb, 0, " ");
					batchLog << "For a negative exponential distribution, mutationParams must have form mean=float." << endl;
					nbErrors++;
				}
			}
			else {
				BatchError(whichInputFile, lineNb, 0, " ");
				batchLog << "For genetic load traits, mutationDistribution must be either uniform, gamma, negExp or normal" << endl;
				nbErrors++;
			}
		}

		// Preview next line
		nextLineSimNb = simNbNotRead;
		bTraitsFile >> nextLineSimNb;
		if (nextLineSimNb == simNbNotRead
			|| bTraitsFile.eof()) {
			stopReading = true;
			gNbTraitFileRows.push_back(lineNb);
		}
		else if (nextLineSimNb != simNb) {
			// about to change sim, conduct checks of all read traits
			nbErrors += checkTraitSetCoherency(allReadTraits, anyNeutralGenetics);
			// Store nb of rows to help reading file later on
			gNbTraitFileRows.push_back(lineNb);
			allReadTraits.clear();
			simNb = nextLineSimNb;
		} // else continue reading traits for same sim
		lineNb++; 
	} // end of while loop

	if (!bTraitsFile.eof()) {
		EOFerror(whichInputFile);
		nbErrors++;
	}

	if (nbErrors > 0) 
		return -111;
	else return 0;
}

int checkTraitSetCoherency(const vector <TraitType>& allReadTraits, const bool& anyNeutralGenetics) {
	int nbErrors = 0;
	const string whichInputFile = "TraitsFile";

	// If genetic output is enabled, a neutral trait must exists
	bool hasNeutral = traitExists(NEUTRAL, allReadTraits);
	if (anyNeutralGenetics && !hasNeutral) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLog << "A neutral stats output option is turned on in genetics file but no neutral trait is specified in traits file." << endl;
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

	if (gDispTraitOpt.isEmigIndVar) {
		if (!hasD0) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "EP or d0 is missing." << endl;
			nbErrors++;
		}
		if (gDispTraitOpt.isEmigSexDep) {
			if (anyEmigNeitherSex) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Emigration SexDep is on but a trait has been supplied without a sex." << endl;
				nbErrors++;
			}
			if (!bothSexesD0) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Either sex is missing for D0 trait." << endl;
				nbErrors++;
			}
		}
		else if (anyEmigSexDep) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "Emigration SexDep is off but a trait has been supplied with a sex." << endl;
			nbErrors++;
		}

		if (gDispTraitOpt.isEmigDensDep) {
			if (!hasEmigAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Emigration alpha is missing." << endl;
				nbErrors++;
			}
			if (!hasEmigBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Emigration beta is missing." << endl;
				nbErrors++;
			}
			if (gDispTraitOpt.isEmigSexDep) {
				if (!bothSexesEmigAlpha) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLog << "Either sex is missing for emigration alpha trait." << endl;
					nbErrors++;
				}
				if (!bothSexesEmigBeta) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLog << "Either sex is missing for emigration beta trait." << endl;
					nbErrors++;
				}
			}
		}
		else {
			if (hasEmigAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Specified emigration alpha, but emigration is not density-dependent." << endl;
				nbErrors++;
			}
			if (hasEmigBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Specified emigration beta, but emigration is not density-dependent." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasD0 || hasEmigAlpha || hasEmigBeta) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLog << "Specified emigration trait, but emigration is not variable." << endl;
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

	if (gDispTraitOpt.isKernTransfIndVar) {
		if (!hasKern1) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "(First) kernel mean is missing." << endl;
			nbErrors++;
		}
		if (gDispTraitOpt.isKernTransfSexDep) {
			if (anyKernelNeitherSex) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Kernel SexDep is on but a trait has been supplied without a sex." << endl;
				nbErrors++;
			}
			if (!bothSexesMeanDist1) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Either sex is missing for first kernel mean trait." << endl;
				nbErrors++;
			}
		}
		else if (anyKernelSexDep) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "Kernel SexDep is off but a trait has been supplied with a sex." << endl;
			nbErrors++;
		}
		if (gDispTraitOpt.usesTwoKernels) {
			if (!hasKern2) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Second kernel mean is missing." << endl;
				nbErrors++;
			}
			if (!hasKernProb) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Kernel probability is missing." << endl;
				nbErrors++;
			}
			if (gDispTraitOpt.isKernTransfSexDep) {
				if (!bothSexesMeanDist2) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLog << "Either sex is missing for second kernel mean trait." << endl;
					nbErrors++;
				}
				if (!bothSexesKernProb) {
					BatchError(whichInputFile, -999, 0, " ");
					batchLog << "Either sex is missing for kernel probability trait." << endl;
					nbErrors++;
				}
			}
		}
		else {
			if (hasKern2) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Specified second kernel, but only one kernel is used." << endl;
				nbErrors++;
			}
			if (hasKernProb) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Specified kernel probability, but only one kernel is used." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasKern1 || hasKern2 || hasKernProb) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLog << "Specified kernel transfer trait, but kernel transfer is not variable." << endl;
		nbErrors++;
	}

	/// SMS
	bool hasDP = traitExists(SMS_DP, allReadTraits);
	bool hasGB = traitExists(SMS_GB, allReadTraits);
	bool hasSMSAlpha = traitExists(SMS_ALPHADB, allReadTraits);
	bool hasSMSBeta = traitExists(SMS_BETADB, allReadTraits);
	if (gDispTraitOpt.isSMSTransfIndVar) {
		if (!hasDP) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "SMS directional persistence trait is missing." << endl;
			nbErrors++;
		}
		if (gDispTraitOpt.usesSMSGoalBias) {
			if (!hasGB) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "SMS goal bias trait is missing." << endl;
				nbErrors++;
			}
			if (!hasSMSAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "SMS alpha direction bias trait is missing." << endl;
				nbErrors++;
			}
			if (!hasSMSBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "SMS beta direction bias trait is missing." << endl;
				nbErrors++;
			}
		}
		else {
			if (hasGB) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "SMS goal bias trait supplied, but SMS GoalType not set to option 2." << endl;
				nbErrors++;
			}
			if (hasSMSAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "SMS alpha direction bias trait supplied, but SMS GoalType not set to option 2." << endl;
				nbErrors++;
			}
			if (hasSMSBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "SMS beta direction bias trait supplied, but SMS GoalType not set to option 2." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasDP || hasGB || hasSMSAlpha || hasSMSBeta) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLog << "Specified SMS trait, but SMS not set to be variable." << endl;
		nbErrors++;
	}

	/// CRW
	bool hasStepLen = traitExists(CRW_STEPLENGTH, allReadTraits);
	bool hasRho = traitExists(CRW_STEPCORRELATION, allReadTraits);
	if (gDispTraitOpt.isCRWTransfIndVar) {
		if (!hasStepLen) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "CRW step length trait is missing." << endl;
			nbErrors++;
		}
		if (!hasRho) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "CRW step correlation trait is missing." << endl;
			nbErrors++;
		}
	}
	else if (hasStepLen || hasRho) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLog << "Specified CRW trait, but CRW not set to be variable." << endl;
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

	if (gDispTraitOpt.isSettIndVar) {
		if (!hasS0) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "Settlement probability trait is missing." << endl;
			nbErrors++;
		}
		if (gDispTraitOpt.isSettSexDep) {
			if (anySettNeitherSex) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Settlement SexDep is on but a trait has been supplied without a sex." << endl;
				nbErrors++;
			}
			if (!bothSexesS0) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Either sex is missing for settlement probabibility trait." << endl;
				nbErrors++;
			}
		}
		else if (anySettSexDep) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "Settlement SexDep is off but a trait has been supplied with a sex." << endl;
			nbErrors++;
		}
		// if settlement is IndVar, it is always density-dependent
		if (!hasSettAlpha) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "Settlement alpha trait is missing." << endl;
			nbErrors++;
		}
		if (!hasSMSBeta) {
			BatchError(whichInputFile, -999, 0, " ");
			batchLog << "Settlement beta trait is missing." << endl;
			nbErrors++;
		}
		if (gDispTraitOpt.isSettSexDep) {
			if (!bothSexesSettAlpha) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Either sex is missing for settlement alpha trait." << endl;
				nbErrors++;
			}
			if (!bothSexesSettBeta) {
				BatchError(whichInputFile, -999, 0, " ");
				batchLog << "Either sex is missing for settlement beta trait." << endl;
				nbErrors++;
			}
		}
	}
	else if (hasS0 || hasSettAlpha || hasSettBeta) {
		BatchError(whichInputFile, -999, 0, " ");
		batchLog << "Specified settlement trait, but settlement not set to be variable." << endl;
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

	string header, traitFileName, traitFileStr;
	int simNb, prevSimNb, errCode;
	string inChromosomeEnds, inRecombinationRate, inTraitsFile, inPatchList, inStages,
		inOutGeneValues, inOutputNeutralStatistics, inOutputPerLocusWCFstat, inOutputPairwiseFst,
		inOutStartGenetics, inOutputInterval, inNbrPatchesToSample, inNIndsToSample;
	int inGenomeSize;
	int nbErrors = 0;
	int nbSims = 0;
	string whichFile = "GeneticsFile";

	const regex patternIntList{ "^\"?([0-9]+,)*[0-9]+\"?$" }; // comma-separated integer list
	bool isMatch = false;

	// Parse header line;
	bGeneticsFile >> header; if (header != "Simulation") nbErrors++;
	bGeneticsFile >> header; if (header != "GenomeSize") nbErrors++;
	bGeneticsFile >> header; if (header != "ChromosomeEnds") nbErrors++;
	bGeneticsFile >> header; if (header != "RecombinationRate") nbErrors++;
	bGeneticsFile >> header; if (header != "OutputGeneValues") nbErrors++;
	bGeneticsFile >> header; if (header != "OutputNeutralStatistics") nbErrors++;
	bGeneticsFile >> header; if (header != "OutputPerLocusWCFstat") nbErrors++;
	bGeneticsFile >> header; if (header != "OutputPairwiseFst") nbErrors++;
	bGeneticsFile >> header; if (header != "OutputStartGenetics") nbErrors++;
	bGeneticsFile >> header; if (header != "OutputInterval") nbErrors++;
	bGeneticsFile >> header; if (header != "PatchList") nbErrors++;
	bGeneticsFile >> header; if (header != "NbrPatchesToSample") nbErrors++;
	bGeneticsFile >> header; if (header != "nIndividualsToSample") nbErrors++;
	bGeneticsFile >> header; if (header != "Stages") nbErrors++;
	bGeneticsFile >> header; if (header != "TraitsFile") nbErrors++;

	if (nbErrors > 0) {
		FormatError(whichFile, nbErrors);
		return -111;
	}

	// Parse data lines
	int whichLine = 1;
	simNb = -98765;
	bGeneticsFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(whichFile, whichLine, 111, "Simulation"); 
		nbErrors++;
	}
	while (simNb != -98765) {
		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		bGeneticsFile >> inGenomeSize >> inChromosomeEnds >> inRecombinationRate >> inOutGeneValues >> inOutputNeutralStatistics >>
			inOutputPerLocusWCFstat >> inOutputPairwiseFst >> inOutStartGenetics >> inOutputInterval >> inPatchList >> inNbrPatchesToSample
			>> inNIndsToSample >> inStages >> inTraitsFile;

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
			batchLog << "ChromosomeEnds must be either a comma-separated list of integers, or blank (#)." << endl;
			nbErrors++;
		}
		set<int> chrEnds = stringToChromosomeEnds(inChromosomeEnds, inGenomeSize);
		const int maxVal = *chrEnds.rbegin();
		if (maxVal >= inGenomeSize) {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "Positions for ChromosomeEnds cannot exceed GenomeSize." << endl;
			nbErrors++;
		}

		// Check RecombinationRate
		if (gNbSexesDisp && inRecombinationRate != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "Do not specify a recombination rate for haploid/asexual systems." << endl;
			nbErrors++;
		}
		else if (inRecombinationRate != "#") {
			float recombinationRate = stof(inRecombinationRate);
			if (recombinationRate < 0.0 || recombinationRate > 0.5) {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "RecombinationRate must be positive cannot exceed 0.5." << endl;
			}
		}

		// Check genetic output fields
		if (inOutGeneValues != "TRUE" && inOutGeneValues != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "OutGeneValues must be either TRUE or FALSE" << endl;
			nbErrors++;
		}
		if (inOutputNeutralStatistics != "TRUE" && inOutputNeutralStatistics != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "OutputNeutralStatistics must be either TRUE or FALSE" << endl;
			nbErrors++;
		}
		if (inOutputPerLocusWCFstat != "TRUE" && inOutputPerLocusWCFstat != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "OutputPerLocusWCFstat must be either TRUE or FALSE" << endl;
			nbErrors++;
		}
		if (inOutputPairwiseFst != "TRUE" && inOutputPairwiseFst != "FALSE") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "OutputPairwiseFst must be either TRUE or FALSE" << endl;
			nbErrors++;
		}

		bool anyNeutralGenetics = inOutputNeutralStatistics == "TRUE"
			|| inOutputPerLocusWCFstat == "TRUE"
			|| inOutputPairwiseFst == "TRUE";
		bool anyGeneticsOutput = inOutGeneValues == "TRUE" || anyNeutralGenetics;

		if (anyGeneticsOutput) {
			if (inOutStartGenetics == "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "OutStartGenetics cannot be left blank (#) if any genetic output option is TRUE." << endl;
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
				batchLog << "OutputInterval cannot be left blank (#) or 0 if any genetic output option is TRUE." << endl;
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
				batchLog << "OutStartGenetics should be blank (#) if all genetic output options are FALSE." << endl;
				nbErrors++;
			}
			if (inOutputInterval != "#" && inOutputInterval != "0") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "OutputInterval should be blank (#) or 0 if all genetic output options are FALSE." << endl;
				nbErrors++;
			}
		}

		// Check PatchList
		if (anyGeneticsOutput) {
			if (inPatchList == "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "PatchList cannot be left blank (#) if any genetic output option is TRUE." << endl;
				nbErrors++;
			}
			else {
				isMatch = regex_search(inPatchList, patternIntList);
				if (!isMatch && inPatchList != "random" && inPatchList != "all" && inPatchList != "random_occupied") {
					BatchError(whichFile, whichLine, 0, " ");
					batchLog << "PatchList must be either a comma-separated list of integers, random, random_occupied or all." << endl;
					nbErrors++;
				}
			}
		}
		else if (inPatchList != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "PatchList should be blank (#) if all genetic output options are FALSE." << endl;
			nbErrors++;
		}

		// Check NbrPatchesToSample
		if (inPatchList == "random" || inPatchList == "random_occupied") {
			if (inNbrPatchesToSample == "#" || inNbrPatchesToSample == "0") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "NbrPatchesToSample cannot be blank (#) or 0 if PatchList is random or random_occupied." << endl;
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
			batchLog << "NbrPatchesToSample must be blank (#) or zero if PatchList is not random or random_occupied." << endl;
			nbErrors++;
		}

		// Check IndividualsToSample
		if (anyGeneticsOutput) {
			if (inNIndsToSample == "#" || inNIndsToSample == "0") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "NIndsToSample cannot be blank (#) or zero if any genetics output option is TRUE." << endl;
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
			batchLog << "NIndsToSample must be blank (#) or zero if all genetics output options are FALSE." << endl;
			nbErrors++;
		}

		// Check Stages
		if (anyGeneticsOutput) {
			if (inStages == "#") {
				BatchError(whichFile, whichLine, 0, " ");
				batchLog << "Stages cannot be blank (#) if any genetic output option is TRUE." << endl;
				nbErrors++;
			}
			else {
				isMatch = regex_search(inStages, patternIntList);
				if (!isMatch && inStages != "all") {
					BatchError(whichFile, whichLine, 0, " ");
					batchLog << "Stages must be either a comma-separated list of integers, or \"all\"." << endl;
					nbErrors++;
				}
			}
		}
		else if (inStages != "#") {
			BatchError(whichFile, whichLine, 0, " ");
			batchLog << "Stages must be blank (#) if all genetic output options are FALSE." << endl;
			nbErrors++;
		}

		// Check TraitsFile
		if (inTraitsFile == "NULL") {
			batchLog << "*** " << inTraitsFile << " is compulsory for genetic models" << endl;
			nbErrors++;
		}
		else {
			traitFileName = inputDirectory + inTraitsFile;
			traitFileStr = "Traits file";
			batchLog << "Checking " << traitFileStr << " " << traitFileName << endl;
			bTraitsFile.open(traitFileName.c_str());
			if (bTraitsFile.is_open()) {
				errCode = CheckTraitsFile(inputDirectory, anyNeutralGenetics);
				if (errCode >= 0) 
					FileHeadersOK(traitFileStr); 
				else 
					nbErrors++;
				bTraitsFile.close();
			}
			else {
				OpenError(traitFileStr, traitFileName); 
				nbErrors++;
			}
			if (bTraitsFile.is_open()) bTraitsFile.close();
			bTraitsFile.clear();
		}

		// read next simulation
		whichLine++;
		simNb = -98765;
		bGeneticsFile >> simNb;
		if (bGeneticsFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation

	if (!bGeneticsFile.eof()) {
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
	bInitFile >> header; if (header != "Simulation") errors++;
	bInitFile >> header; if (header != "SeedType") errors++;
	bInitFile >> header; if (header != "FreeType") errors++;
	bInitFile >> header; if (header != "SpType") errors++;
	bInitFile >> header; if (header != "InitDens") errors++;
	bInitFile >> header;
	if (patchmodel) { if (header != "IndsHa") errors++; }
	else { if (header != "IndsCell") errors++; }
	bInitFile >> header; if (header != "minX") errors++;
	bInitFile >> header; if (header != "maxX") errors++;
	bInitFile >> header; if (header != "minY") errors++;
	bInitFile >> header; if (header != "maxY") errors++;
	bInitFile >> header; if (header != "NCells") errors++;
	bInitFile >> header; if (header != "NSpCells") errors++;
	bInitFile >> header; if (header != "InitFreezeYear") errors++;
	bInitFile >> header; if (header != "RestrictRows") errors++;
	bInitFile >> header; if (header != "RestrictFreq") errors++;
	bInitFile >> header; if (header != "FinalFreezeYear") errors++;
	bInitFile >> header; if (header != "InitIndsFile") errors++;
	if (stagestruct) {
		bInitFile >> header; if (header != "InitAge") errors++;
		for (i = 1; i < stages; i++) {
			colheader = "PropStage" + to_string(i);
			bInitFile >> header; if (header != colheader) propnerrors++;
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
	bInitFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {
		current = CheckStageSex(filetype, line, simNb, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) simuls++;
		errors += current.errors;
		prev = current;

		bInitFile >> seedtype >> freetype >> sptype >> initdens >> inds_per_ha;
		if (!patchmodel) indscell = (int)inds_per_ha;
		if (seedtype < 0 || seedtype > 2) {
			BatchError(filetype, line, 2, "SeedType"); errors++;
		}
		if (landtype == 9 && seedtype != 0) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLog << "SeedType must be 0 for an artificial landscape"
				<< endl;
		}
		if (!speciesdist && seedtype == 1) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLog << "SeedType may not be 1 if there is no initial species distribution map"
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
				if (patchmodel) {
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

		bInitFile >> minX >> maxX >> minY >> maxY >> nCells >> nSpCells;
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
				batchLog << "NCells may not be greater than the area specified (i.e. "
					<< range_cells << " cells)" << endl;
			}
		}
		if (seedtype == 1 && sptype == 1 && nSpCells < 1) {
			BatchError(filetype, line, 11, "NSpCells"); errors++;
		}

		bInitFile >> initFreezeYear >> restrictRows >> restrictFreq >> finalFreezeYear;
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

		bInitFile >> filename;
		if (filename == "NULL") {
			if (seedtype == 2) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLog << ftype2 << " is compulsory for SeedType 2" << endl;
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
					batchLog << "Checking " << ftype2 << " " << fname << endl;
					bInitIndsFile.open(fname.c_str());
					if (bInitIndsFile.is_open()) {
						err = CheckInitIndsFile();
						if (err == 0) FileHeadersOK(ftype2); else errors++;
						bInitIndsFile.close();
					}
					else {
						OpenError(ftype2, fname); errors++;
					}
					if (bInitIndsFile.is_open()) bInitIndsFile.close();
					bInitIndsFile.clear();
					indsfiles.push_back(filename);
				}
			}
			else {
				BatchError(filetype, line, 0, " "); errors++;
				batchLog << ftype2 << " must be NULL for SeedType "
					<< seedtype << endl;
			}
		}

		if (stagestruct) {
			bInitFile >> initAge;
			if (seedtype != 2 && (initAge < 0 || initAge > 2)) {
				BatchError(filetype, line, 2, "initAge"); errors++;
			}
			float propstage;
			float cumprop = 0.0;
			for (i = 1; i < stages; i++) {
				bInitFile >> propstage;
				cumprop += propstage;
				if (seedtype != 2 && (propstage < 0.0 || propstage > 1.0)) {
					colheader = "PropStage" + to_string(i);
					BatchError(filetype, line, 20, colheader); errors++;
				}
			}
			if (seedtype != 2 && (cumprop < 0.99999 || cumprop > 1.00001)) {
				BatchError(filetype, line, 0, " "); errors++;
				batchLog << "Initial proportions must sum to 1.0" << endl;
			}
		}

		// read next simulation
		line++;
		simNb = -98765;
		bInitFile >> simNb;
		if (bInitFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (current.simLines != current.reqdSimLines) {
		BatchError(filetype, line, 0, " "); errors++;
		batchLog << gNbLinesStr << current.simNb
			<< gShouldBeStr << current.reqdSimLines << endl;
	}
	if (!bInitFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int CheckInitIndsFile(void) {
	string header;
	int year, species, patchID, x, y, ninds, sex, age, stage, prevyear;

	int errors = 0;
	string filetype = "InitIndsFile";

	// Parse header line
	bInitIndsFile >> header; if (header != "Year") errors++;
	bInitIndsFile >> header; if (header != "Species") errors++;
	if (patchmodel) {
		bInitIndsFile >> header; if (header != "PatchID") errors++;
	}
	else {
		bInitIndsFile >> header; if (header != "X") errors++;
		bInitIndsFile >> header; if (header != "Y") errors++;
	}
	bInitIndsFile >> header; if (header != "Ninds") errors++;
	if (reproductn > 0) {
		bInitIndsFile >> header; if (header != "Sex") errors++;
	}
	if (stagestruct) {
		bInitIndsFile >> header; if (header != "Age") errors++;
		bInitIndsFile >> header; if (header != "Stage") errors++;
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
	bInitIndsFile >> year;
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
		bInitIndsFile >> species;
		if (species != 0) {
			BatchError(filetype, line, 0, " "); errors++;
			batchLog << "Species must be 0" << endl;
		}
		if (patchmodel) {
			bInitIndsFile >> patchID;
			if (patchID < 1) {
				BatchError(filetype, line, 11, "PatchID"); errors++;
			}
		}
		else {
			bInitIndsFile >> x >> y;
			if (x < 0 || y < 0) {
				BatchError(filetype, line, 19, "X and Y"); errors++;
			}
		}
		bInitIndsFile >> ninds;
		if (ninds < 1) {
			BatchError(filetype, line, 11, "Ninds"); errors++;
		}
		if (reproductn > 0) {
			bInitIndsFile >> sex;
			if (sex < 0 || sex > 1) {
				BatchError(filetype, line, 1, "Sex"); errors++;
			}
		}
		if (stagestruct) {
			bInitIndsFile >> age >> stage;
			if (age < 1) {
				BatchError(filetype, line, 11, "Age"); errors++;
			}
			if (stage < 1) {
				BatchError(filetype, line, 11, "Stage"); errors++;
			}
			if (stage >= stages) {
				BatchError(filetype, line, 4, "Stage", "no. of stages"); errors++;
			}
		}
		line++;
		year = -98765;
		bInitIndsFile >> year;
		if (bInitIndsFile.eof()) year = -98765;
	} // end of while loop
	if (!bInitIndsFile.eof()) {
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
			batchLog << "No. of lines for previous Simulation " << prev.simNb
				<< gShouldBeStr << prev.reqdSimLines << endl;
		}
	}
	current.simNb = simNb;

	// validate inStageDep
	if (stagestruct) {
		if (isStageDep != 0 && isStageDep != 1) {
			BatchError(whichInputFile, whichLine, 1, "StageDep"); current.errors++;
			isStageDep = 1; // to calculate required number of lines
		}
	}
	else {
		if (isStageDep != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLog << "StageDep must be 0 for non-stage-structured model" << endl;
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
			batchLog << "SexDep must be 0 for asexual model" << endl;
			isSexDep = 0; // to calculate required number of lines
		}
	}
	if (current.isNewSim) { // set required number of lines
		if (isStageDep) {
			current.reqdSimLines = stages;
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
				batchLog << "Stages must be sequentially numbered from 0" << endl;
			}
		}
		else { // there must be 1 line for each stage
			if (stage != current.simLines - 1) {
				BatchError(whichInputFile, whichLine, 0, " "); 
				current.errors++;
				batchLog << "Stages must be sequentially numbered from 0" << endl;
			}
		}
	}
	else { // no stage-dependent emigration
		if (stage != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLog << "Stage must be 0 for non-stage-structured model" << endl;
		}
	}
	// validate sex
	if (isSexDep) {
		if (sex != (current.simLines + 1) % 2) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLog << "Sex must be alternately 0 and 1 if SexDep is 1" << endl;
		}
	}
	else {
		if (sex != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLog << "Sex must be 0 if SexDep is 0" << endl;
		}
	}

	// validate inIndVar
	if (isStageDep && !mustCheckStgDepWithIndVar) {
		if (isIndVar != 0) {
			BatchError(whichInputFile, whichLine, 0, " "); current.errors++;
			batchLog << "IndVar must be 0 if stage-dependent" << endl;
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
		batchLog << "*** Error in " << filename << ": ";
	}
	else {
		batchLog << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		batchLog << fieldname << " must be 0 or 1";
		break;
	case 2:
		batchLog << fieldname << " must be 0, 1 or 2";
		break;
	case 3:
		batchLog << fieldname << " must be 0, 1, 2 or 3";
		break;
	case 4:
		batchLog << fieldname << " must be from 0 to 4";
		break;
	case 5:
		batchLog << fieldname << " must be from 0 to 5";
		break;
	case 6:
		batchLog << fieldname << " must be from 0 to 6";
		break;
	case 7:
		batchLog << fieldname << " must be from 0 to 7";
		break;
	case 10:
		batchLog << fieldname << " must be greater than zero";
		break;
	case 11:
		batchLog << fieldname << " must be 1 or more";
		break;
	case 12:
		batchLog << fieldname << " must be 2 or more";
		break;
	case 13:
		batchLog << fieldname << " must be 3 or more";
		break;
	case 18:
		batchLog << fieldname << " must be greater than 1.0";
		break;
	case 19:
		batchLog << fieldname << " must be 0 or more";
		break;
	case 20:
		batchLog << fieldname << " must be between 0 and 1";
		break;
	case 21:
		batchLog << fieldname << " must be greater than 1";
		break;
	case 33:
		batchLog << fieldname << " must be 1, 2 or 3";
		break;
	case 44:
		batchLog << fieldname << " must be from 1 to 4";
		break;
	case 55:
		batchLog << fieldname << " must be from 1 to 5";
		break;
	case 66:
		batchLog << fieldname << " must be from 1 to 6";
		break;
	case 100:
		batchLog << fieldname << " must be between 0 and 100";
		break;
	case 111:
		batchLog << fieldname << " must match the first Simulation in ParameterFile";
		break;
	case 222:
		batchLog << "Simulation numbers must be sequential integers";
		break;
	case 333:
		batchLog << "No. of " << fieldname << " columns must equal max. no. of habitats ("
			<< maxNhab << ") and be sequentially numbered starting from 1";
		break;
	case 444:
		batchLog << "No. of " << fieldname << " columns must be one fewer than no. of stages, i.e. "
			<< stages - 1 << ", and be sequentially numbered starting from 1";
		break;
	case 555:
		batchLog << "No. of " << fieldname << " columns must equal no. of stages, i.e. "
			<< stages << ", and be sequentially numbered starting from 0";
		break;
	case 666:
		batchLog << fieldname << " must be a unique positive integer";
		break;
	default:
		batchLog << "*** Unspecified error regarding parameter " << fieldname;
	}
	if (option != 0) batchLog << endl;
}

void BatchError(string filename, int line, int option, string fieldname, string fieldname2)
{
	if (line == -999) { // message does not cite line number
		batchLog << "*** Error in " << filename << ": ";
	}
	else {
		batchLog << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		batchLog << fieldname << " must be greater than " << fieldname2;
		break;
	case 2:
		batchLog << fieldname << " must be greater than or equal to " << fieldname2;
		break;
	case 3:
		batchLog << fieldname << " must be less than or equal to " << fieldname2;
		break;
	case 4:
		batchLog << fieldname << " must be less than " << fieldname2;
		break;
	default:
		batchLog << "*** Unspecified error regarding parameters " << fieldname
			<< " and " << fieldname2;
	}
	if (option != 0) batchLog << endl;
}

void CtrlFormatError(void)
{
	cout << "Format error in Control file" << endl;
	batchLog << endl << "***" << endl << "*** Format error in Control file:"
		<< gCaseSensitiveStr << " and file names" << gSpecMustMatchStr
		<< endl
		<< "***" << endl;
}

void ArchFormatError(void)
{
	batchLog << "*** Format error in ArchFile:" << gCaseSensitiveStr << gSpecMustMatchStr << endl;
}

void FormatError(string filename, int errors)
{
	batchLog << "*** Format error in header line of ";
	if (errors == 0) {
		batchLog << filename << endl;
	}
	else {
		batchLog << filename << ": " << errors << " error";
		if (errors > 1) batchLog << "s";
		batchLog << " detected" << endl;
	}
}

void OpenError(string ftype, string fname)
{
	batchLog << "*** Unable to open " << ftype << " " << fname << endl;
}

void EOFerror(string filename)
{
	batchLog << "*** Failed to read to EOF in " << filename << endl;
}

void FileOK(string ftype, int n, int option)
{
	batchLog << ftype << " OK: total no. of ";
	switch (option) {
	case 0:
		batchLog << "simulations = ";
		break;
	case 1:
		batchLog << "landscapes = ";
		break;
	case 2:
		batchLog << "parameters = ";
		break;
	default:
		batchLog << "PROBLEMS = ";
	}
	batchLog << n << endl;
}

void FileHeadersOK(string filename)
{
	batchLog << filename << " OK" << endl;
}

void SimulnCountError(string filename)
{
	batchLog << "*** No. of simulations in " << filename
		<< " does not match no. in ParameterFile" << endl;
}

//---------------------------------------------------------------------------
int ReadLandFile(int option)
{
	if (option == 0) { // open file and read header line
		landfile.open(landFile.c_str());
		if (landfile.is_open()) {
			string header;
			int nheaders;
			if (landtype == 9) nheaders = 9; // artificial landscape
			else { // imported raster map
				nheaders = 7;
			}
			for (int i = 0; i < nheaders; i++) landfile >> header;
		}
		else return 1;
	}

	if (option == 9) { // close file
		if (landfile.is_open()) {
			landfile.close();  landfile.clear();
		}
	}
	return 0;
}

//---------------------------------------------------------------------------
int ReadLandFile(int option, Landscape* pLandscape)
{
	landParams ppLand = pLandscape->getLandParams();
	genLandParams ppGenLand = pLandscape->getGenLandParams();
	simParams sim = paramsSim->getSim();

	if (landtype == 9) { //artificial landscape
		ppLand.rasterType = 9;
		landfile >> ppLand.landNum >> ppGenLand.fractal >> ppGenLand.continuous
			>> ppLand.dimX >> ppLand.dimY >> ppGenLand.minPct >> ppGenLand.maxPct
			>> ppGenLand.propSuit >> ppGenLand.hurst;
		ppLand.maxX = ppLand.dimX - 1; ppLand.maxY = ppLand.dimY - 1;

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
		string inNHabPlaceholder; // no longer necessary to read no. of habitats from landFile
		landfile >> ppLand.landNum >> inNHabPlaceholder >> name_landscape >> name_patch;
		landfile >> gNameCostFile >> name_dynland >> name_sp_dist;
		if (landtype == 2) 
			ppLand.nHab = 1; // habitat quality landscape has one habitat class
	}

	pLandscape->setLandParams(ppLand, sim.batchMode);
	pLandscape->setGenLandParams(ppGenLand);

	return ppLand.landNum;
}

//---------------------------------------------------------------------------
int ReadDynLandFile(Landscape* pLandscape) {

	string landchangefile, patchchangefile, costchangefile;
	int change, imported;
	int nchanges = 0;
	landChange chg;
	landParams ppLand = pLandscape->getLandParams();
	string fname = paramsSim->getDir(1) + name_dynland;

	dynlandfile.open(fname.c_str());
	if (dynlandfile.is_open()) {
		string header;
		int nheaders = 5;
		for (int i = 0; i < nheaders; i++) 
			dynlandfile >> header;
	}
	else {
		dynlandfile.clear();
		return 72727;
	}

	// read data lines
	change = -98765;
	dynlandfile >> change; // first change number
	while (change != -98765) {
		chg.chgnum = change;
		dynlandfile >> chg.chgyear >> landchangefile >> patchchangefile >> costchangefile;
		chg.habfile = paramsSim->getDir(1) + landchangefile;
		chg.pchfile = paramsSim->getDir(1) + patchchangefile;
		if (costchangefile == "NULL") 
			chg.costfile = "none";
		else 
			chg.costfile = paramsSim->getDir(1) + costchangefile;
		nchanges++;
		pLandscape->addLandChange(chg);
		// read first field on next line
		change = -98765;
		dynlandfile >> change;
		if (dynlandfile.eof()) {
			change = -98765;
		}
	}

	dynlandfile.close();
	dynlandfile.clear();

	// read landscape change maps
	if (ppLand.patchModel) {
		pLandscape->createPatchChgMatrix();
	}
	if (costchangefile != "NULL") {
		pLandscape->createCostsChgMatrix();
	}
	for (int i = 0; i < nchanges; i++) {
		if (costchangefile == "NULL") 
			imported = pLandscape->readLandChange(i, false);
		else 
			imported = pLandscape->readLandChange(i, true);
		if (imported != 0) {
			return imported;
		}
		if (ppLand.patchModel) {
			pLandscape->recordPatchChanges(i + 1);
		}
		if (costchangefile != "NULL") {
			pLandscape->recordCostChanges(i + 1);
		}
	}
	if (ppLand.patchModel) {
		// record changes back to original landscape for multiple replicates
		pLandscape->recordPatchChanges(0);
		pLandscape->deletePatchChgMatrix();
	}
	if (costchangefile != "NULL") {
		pLandscape->recordCostChanges(0);
		pLandscape->deleteCostsChgMatrix();
	}
	return 0;
}

//--------------------------------------------------------------------------

void flushHeader(ifstream& ifs) {
	string headerLine;
	// Pass the first line (headers) to an e;pty string...
	std::getline(ifs, headerLine);
	// ... and do nothing with it 
}

int ReadGeneticsFile(ifstream& ifs, Landscape* pLandscape) {

	string indir = paramsSim->getDir(1);
	bool outputGeneValues, outputWCFstat, outputPerLocusWCFstat, outputPairwiseFst;
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
		outputWCFstat = (parameters[5] == "TRUE");
		outputPerLocusWCFstat = (parameters[6] == "TRUE");
		outputPairwiseFst = (parameters[7] == "TRUE");
		outputStartGenetics = stoi(parameters[8]);
		outputGeneticInterval = stoi(parameters[9]);

		string inPatches = parameters[10];
		string patchSamplingOption;
		int nPatchesToSample = stoi(parameters[11]);
		if (inPatches != "all" && inPatches != "random" && inPatches != "random_occupied") {
			// then must be a list of indices
			patchSamplingOption = "list";
			patchList = stringToPatches(inPatches);
			if (patchList.contains(0)) throw logic_error("Patch sampling: ID 0 is reserved for the matrix and should not be sampled.");
		}
		else {
			patchSamplingOption = inPatches;
			// patchList remains empty, filled when patches are sampled every gen
		}
		const string strNbInds = parameters[12];
		const int nbStages = pSpecies->getStageParams().nStages;
		set<int> stagesToSampleFrom = stringToStages(parameters[13], nbStages);

		pSpecies->setGeneticParameters(chrEnds, genomeSize, recombinationRate,
			patchList, strNbInds, stagesToSampleFrom, nPatchesToSample);
		paramsSim->setGeneticSim(patchSamplingOption, outputGeneValues, outputWCFstat, outputPerLocusWCFstat, outputPairwiseFst, outputStartGenetics, outputGeneticInterval);

		gPathToTraitsFile = indir + parameters[14];
	}
	return 0;
}

int ReadTraitsFile(ifstream& ifs, const int& whichSim) {

	pSpecies->clearTraitTable();
	int prevsimNb = -998, simNb;
	int nbRowsToRead = gNbTraitFileRows[whichSim];

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

	// Initial distribution parameters
	const DistributionType initDist = stringToDistributionType(parameters[6]);
	const map<GenParamType, float> initParams = stringToParameterMap(parameters[7]);

	// Dominance distribution parameters
	const DistributionType dominanceDist = stringToDistributionType(parameters[8]);
	const map<GenParamType, float> dominanceParams = stringToParameterMap(parameters[9]);

	// Mutation parameters
	bool isInherited = (parameters[10] == "TRUE");
	DistributionType mutationDistribution = isInherited ? 
		stringToDistributionType(parameters[11]) : 
		DistributionType::NONE;
	map<GenParamType, float> mutationParameters;
	float mutationRate = isInherited ? stof(parameters[13]) : 0.0;
	if (isInherited) {
		mutationParameters = stringToParameterMap(parameters[12]);
	}

	int ploidy = gNbSexesDisp;

	// Create species trait
	SpeciesTrait* trait = new SpeciesTrait(
		traitType, sex, 
		positions, expressionType, 
		initDist, initParams, 
		dominanceDist, dominanceParams, 
		isInherited, mutationRate, 
		mutationDistribution, mutationParameters,
		ploidy
	);
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
				positionRange.push_back(stoi(valueWithin));
			}
			switch (positionRange.size())
			{
			case 1: // single position
				if (positionRange[0] > genomeSize)
					throw logic_error("Traits file: ERROR - trait positions must not exceed genome size");
				positions.insert(positionRange[0]);
				break;
			case 2: // dash-separated range
				if (positionRange[0] > genomeSize || positionRange[1] > genomeSize) {
					throw logic_error("Traits file: ERROR - trait positions must not exceed genome size");
				}
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
			if (position > genomeSize)
				cout << endl << "Traits file: ERROR - trait positions " << position << " must not exceed genome size" << endl;
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
int ReadParameters(int option, Landscape* pLandscape)
{
	int iiii;
	int error = 0;
	landParams paramsLand = pLandscape->getLandParams();

	if (option == 0) { // open file and read header line
		parameters.open(parameterFile.c_str());
		if (parameters.is_open()) {
			string header;
			int nheaders = 44 + paramsLand.nHabMax;
			for (int i = 0; i < nheaders; i++) parameters >> header;
			return 0;
		}
		else return 1;
	}

	if (option == 9) { // close file
		if (parameters.is_open()) {
			parameters.close();  parameters.clear();
		}
		return 0;
	}

	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	simParams sim = paramsSim->getSim();
	simView v = paramsSim->getViews();

	if (!parameters.is_open()) {
		cout << endl << "ReadParameters(): ERROR - ParameterFile is not open" << endl;
		return 4086534;
	}

	int gradType, shift_begin, shift_stop;
	float grad_inc, opt_y, f, optEXT, shift_rate;
	bool shifting;

	parameters >> sim.simulation >> sim.reps >> sim.years;
	parameters >> iiii;
	if (iiii == 1) sim.absorbing = true; else sim.absorbing = false;
	parameters >> gradType;
	parameters >> grad_inc >> opt_y >> f >> optEXT >> iiii >> shift_rate;
	if (iiii == 1 && gradType != 0) shifting = true; else shifting = false;
	parameters >> shift_begin >> shift_stop;
	paramsGrad->setGradient(gradType, grad_inc, opt_y, f, optEXT);
	if (shifting) paramsGrad->setShifting(shift_rate, shift_begin, shift_stop);
	else paramsGrad->noShifting();

	parameters >> iiii;
	if (iiii == 0) env.stoch = false;
	else {
		env.stoch = true;
		if (iiii == 2) env.local = true; else env.local = false;
	}
	if (paramsLand.patchModel && env.local) error = 101;
	parameters >> iiii;
	if (iiii == 1) env.inK = true; else env.inK = false;
	// as from v1.1, there is just one pair of min & max values,
	// which are attributes of the species
	// ULTIMATELY, THE PARAMETER FILE SHOULD HAVE ONLY TWO COLUMNS ...
	float minR, maxR, minK, maxK;
	parameters >> env.ac >> env.std >> minR >> maxR >> minK >> maxK;
	if (env.inK) {
		float minKK, maxKK;
		minKK = minK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
		maxKK = maxK * (((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f);
		pSpecies->setMinMax(minKK, maxKK);
	}
	else pSpecies->setMinMax(minR, maxR);
	parameters >> iiii;
	if (iiii == 1) env.localExt = true; else env.localExt = false;
	if (paramsLand.patchModel && env.localExt) error = 102;
	parameters >> env.locExtProb;
	paramsStoch->setStoch(env);

	parameters >> dem.propMales >> dem.harem >> dem.bc >> dem.lambda;
	pSpecies->setDemogr(dem);

	float k;

	if (landtype == 9) { // artificial landscape
		// only one value of K is read, but it must be applied as the second habitat if the
		// landscape is discrete (the first is the matrix where K = 0) or as the first 
		// (only) habitat if the landscape is continuous
		genLandParams genland = pLandscape->getGenLandParams();
		int nhab;
		if (genland.continuous) nhab = 1;
		else nhab = 2;
		pSpecies->createHabK(nhab);
		parameters >> k;
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
			parameters >> k;
			k *= ((float)paramsLand.resol * (float)paramsLand.resol) / 10000.0f;
			pSpecies->setHabK(i, k);
		}
	}

	parameters >> sim.outStartPop >> sim.outStartInd
		>> sim.outStartTraitCell >> sim.outStartTraitRow >> sim.outStartConn;
	parameters >> sim.outIntRange >> sim.outIntOcc >> sim.outIntPop >> sim.outIntInd
		>> sim.outIntTraitCell >> sim.outIntTraitRow >> sim.outIntConn;
	if (sim.outIntRange > 0)     sim.outRange = true; else sim.outRange = false;
	if (sim.outIntOcc > 0)       sim.outOccup = true; else sim.outOccup = false;
	if (sim.outIntPop > 0)       sim.outPop = true; else sim.outPop = false;
	if (sim.outIntInd > 0)       sim.outInds = true; else sim.outInds = false;
	if (sim.outIntRange > 0)     sim.outRange = true; else sim.outRange = false;
	if (sim.outIntTraitCell > 0) sim.outTraitsCells = true; else sim.outTraitsCells = false;
	if (sim.outIntTraitRow > 0)  sim.outTraitsRows = true; else sim.outTraitsRows = false;
	if (sim.outIntConn > 0)      sim.outConnect = true; else sim.outConnect = false;
	if (sim.outOccup && sim.reps < 2) error = 103;
	if (paramsLand.patchModel) {
		if (sim.outTraitsRows) error = 104;
	}
	else {
		if (sim.outConnect) error = 105;
	}
#if RSDEBUG
	DEBUGLOG << "ReadParameters(): outRange=" << sim.outRange << " outInt=" << sim.outIntRange
		<< endl;
#endif
	parameters >> iiii >> sim.mapInt;
	if (iiii == 0) sim.saveMaps = false;
	else sim.saveMaps = true;
	parameters >> iiii;
	if (iiii == 0) sim.saveVisits = false;
	else sim.saveVisits = true;
	parameters >> iiii;
	if (iiii == 0) sim.drawLoaded = false; else sim.drawLoaded = true;
	parameters >> iiii;
	if (iiii == 1) sim.fixReplicateSeed = true; else sim.fixReplicateSeed = false;

	paramsSim->setSim(sim);
	paramsSim->setViews(v);

	return error;
}

//---------------------------------------------------------------------------
int ReadStageStructure(int option)
{
	string name;
	int simulation, postDestructn;
	stageParams sstruct = pSpecies->getStageParams();
	string Inputs = paramsSim->getDir(1);

	if (option == 0) { // open file and read header line
		ssfile.open(stageStructFile.c_str());
		string header;
		int nheaders = 18;
		for (int i = 0; i < nheaders; i++) ssfile >> header;
		return 0;
	}

	if (option == 9) { // close file
		if (ssfile.is_open()) {
			ssfile.close(); ssfile.clear();
		}
		return 0;
	}

	ssfile >> simulation;
	ssfile >> postDestructn >> sstruct.probRep >> sstruct.repInterval >> sstruct.maxAge;
	if (postDestructn == 1) sstruct.disperseOnLoss = true;
	else sstruct.disperseOnLoss = false;

	ssfile >> name;
	// 'name' is TransMatrixFile
	tmfile.open((Inputs + name).c_str());
	ReadTransitionMatrix(sstruct.nStages, sexesDem, 0, 0);
	tmfile.close(); tmfile.clear();
	ssfile >> sstruct.survival;

	float devCoeff, survCoeff;
	ssfile >> sstruct.fecDens >> sstruct.fecStageDens >> name; // 'name' is FecStageWtsFile
	if (name != "NULL") {
		fdfile.open((Inputs + name).c_str());
		ReadStageWeights(1);
		fdfile.close(); fdfile.clear();
	}
	ssfile >> sstruct.devDens >> devCoeff >> sstruct.devStageDens >> name; // 'name' is DevStageWtsFile
	if (name != "NULL") {
		ddfile.open((Inputs + name).c_str());
		ReadStageWeights(2);
		ddfile.close(); ddfile.clear();
	}
	ssfile >> sstruct.survDens >> survCoeff >> sstruct.survStageDens >> name; // 'name' is SurvStageWtsFile
	if (name != "NULL") {
		sdfile.open((Inputs + name).c_str());
		ReadStageWeights(3);
		sdfile.close(); sdfile.clear();
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
		tmfile >> header;
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
			tmfile >> header;
			for (int j = 0; j < nstages; j++)
			{
				tmfile >> matrix[j][i];
			}
			tmfile >> minAge; pSpecies->setMinAge(i, 0, minAge);
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
			tmfile >> header;
			for (int j = 0; j < nstages * 2; j++) 
				tmfile >> matrix[j][i];
			if (i == 0) {
				tmfile >> minAge; 
				pSpecies->setMinAge(i, 0, minAge); 
				pSpecies->setMinAge(i, 1, minAge);
			}
			else {
				tmfile >> minAge;
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
		for (i = 0; i < n + 1; i++) fdfile >> header;
		// read coefficients
		for (i = 0; i < n; i++) {
			fdfile >> header;
			for (j = 0; j < n; j++) {
				fdfile >> f; pSpecies->setDDwtFec(j, i, f);
			}
		}
		break;
	}

	case 2: { // development
		//create stage weights matrix
		pSpecies->createDDwtDev(n);
		for (i = 0; i < n + 1; i++) ddfile >> header;
		//read coefficients
		for (i = 0; i < n; i++) {
			ddfile >> header;
			for (j = 0; j < n; j++) {
				ddfile >> f; pSpecies->setDDwtDev(j, i, f);
			}
		}
		break;
	}

	case 3: { // sstruct.survival
		//create stage weights matrix
		pSpecies->createDDwtSurv(n);
		for (i = 0; i < n + 1; i++) 
			sdfile >> header;
		//read coefficients
		for (i = 0; i < n; i++) {
			sdfile >> header;
			for (j = 0; j < n; j++) {
				sdfile >> f; pSpecies->setDDwtSurv(j, i, f);
			}
		}
		break;
	}

	}

	return 0;
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
int ReadEmigration(int option)
{
	int errorCode = 0;

	if (option == 0) { // open file and read header line
		emigFile.open(emigrationFile.c_str());
		string header;
		for (int i = 0; i < nHeadersEmig; i++) emigFile >> header;
		return 0;
	}
	if (option == 9) { // close file
		if (emigFile.is_open()) {
			emigFile.close(); emigFile.clear();
		}
		return 0;
	}

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

		emigFile >> simulationNb >> inDensDep >> inFullKernel 
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
		emigFile >> inStage >> inSex;

		// ERROR MESSAGES SHOULD NEVER BE ACTIVATED ---------------------------------
		if (dem.repType == 0 && emig.sexDep) {
			errorCode = 301;
		}
		if (!dem.stageStruct && emig.stgDep) {
			errorCode = 303;
		}
		//---------------------------------------------------------------------------

		emigFile >> inEp >> inD0 >> inAlpha >> inBeta;

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
int ReadTransferFile(int option, Landscape* pLandscape)
{
	int error = 0;
	landParams paramsLand = pLandscape->getLandParams();
	transferRules trfr = pSpecies->getTransferRules();

	if (option == 0) { // open file and read header line
		transFile.open(transferFile.c_str());
		string header;
		int nheaders = 0;
		if (trfr.usesMovtProc) {

			if (paramsLand.generated)
				pSpecies->createHabCostMort(paramsLand.nHab);
			else
				pSpecies->createHabCostMort(paramsLand.nHabMax);

			if (trfr.moveType == 1) { // SMS
				int standardcols = 10;
				if (paramsLand.generated) {
					nheaders = standardcols + 6; // artificial landscape
				}
				else { // real landscape
					if (paramsLand.rasterType == 0)
						nheaders = standardcols + 3 + 2 * paramsLand.nHabMax; // habitat codes
					else nheaders = standardcols + 3; // habitat quality
				}
			}
			else { // CRW
				if (paramsLand.generated) {
					nheaders = 7;
				}
				else {
					if (paramsLand.rasterType == 0) nheaders = 7 + paramsLand.nHabMax;
					else nheaders = 7;
				}
			}
		} else { // dispersal kernel
			nheaders = 14;
		}
		for (int i = 0; i < nheaders; i++) 
			transFile >> header;
		return 0;
	}

	if (option == 9) { // close file
		if (transFile.is_open()) {
			transFile.close(); transFile.clear();
		}
		return 0;
	}

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

		transFile >> simNb >> inStageDep >> inSexDep >> inKernelType >> inDistMort >> inIndVar;
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
		transFile >> inStage >> inSex;

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
			transFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			pSpecies->setSpKernTraits(0, 0, kernParams, paramsLand.resol);
			break;

		case 1: // sex-dependent
			if (trfr.twinKern)
			{
				transFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				transFile >> kernParams.meanDist1; kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(0, inSex, kernParams, paramsLand.resol);

			break;

		case 2: // stage-dependent
			if (trfr.twinKern)
			{
				transFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				transFile >> kernParams.meanDist1; kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(inStage, 0, kernParams, paramsLand.resol);
			break;

		case 3: // sex- & stage-dependent
			if (trfr.twinKern)
			{
				transFile >> kernParams.meanDist1 >> kernParams.meanDist2 >> kernParams.probKern1;
			}
			else {
				transFile >> kernParams.meanDist1; kernParams.meanDist2 = kernParams.meanDist1; 
				kernParams.probKern1 = 1.0;
			}
			pSpecies->setSpKernTraits(inStage, inSex, kernParams, paramsLand.resol);
			break;
		} // end of switch (sexkernels)

		// mortality
		if (inStage == 0 && inSex == 0) {
			trfrMortParams mort;
			transFile >> mort.fixedMort >> mort.mortAlpha >> mort.mortBeta;
			pSpecies->setMortParams(mort);
		}
		else for (int i = 0; i < 3; i++) 
			transFile >> flushMort;

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

	transFile >> simNb >> inIndVar >> move.pr >> move.prMethod >> move.dp
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
					transFile >> inHabMort;
					pSpecies->setHabMort(i, inHabMort);
				}
			}
			else { // constant step mortality
				for (int i = 0; i < paramsLand.nHabMax; i++) 
					transFile >> flushHabMort;
			}
		}
	}
	else { // artificial landscape
		if (trfr.habMort)
		{ // habitat-dependent step mortality
			// values are for habitat (hab=1) then for matrix (hab=0)
			transFile >> inMortHabitat >> inMortMatrix;
			pSpecies->setHabMort(1, inMortHabitat);
			pSpecies->setHabMort(0, inMortMatrix);
		}
		else { // constant step mortality
			transFile >> flushHabMort >> flushHabMort;
		}
	}
	trfr.costMap = (gNameCostFile != "NULL") ? true : false;

	if (!paramsLand.generated) { // imported landscape
		if (paramsLand.rasterType == 0) { // habitat codes
			if (trfr.costMap)
			{
				for (int i = 0; i < paramsLand.nHabMax; i++) 
					transFile >> flushCostHab;
			}
			else { // not costMap
				for (int i = 0; i < paramsLand.nHabMax; i++) {
					transFile >> inCostHab; 
					pSpecies->setHabCost(i, inCostHab);
				}
			}
		}
	}
	else { // artificial landscape
		if (trfr.costMap) // should not occur 
		{
			transFile >> flushCostHab >> flushCostHab;
		}
		else { // not costMap
			// costs are for habitat (hab=1) then for matrix (hab=0)
			transFile >> inCostHab >> inCostMatrix;
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
	transFile >> simNb >> inIndVar;
	if (inIndVar == 0) trfr.indVar = false;
	else trfr.indVar = true;

	transFile >> move.stepLength >> move.rho;
	transFile >> inStraightenPath >> inSMconst >> move.stepMort;

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
				transFile >> inHabMort;
				pSpecies->setHabMort(i, inHabMort);
			}
		}
		else { // constant step mortality
			for (int i = 0; i < paramsLand.nHabMax; i++) 
				transFile >> flushHabMort;
		}
	}
	pSpecies->setTrfrRules(trfr);
	pSpecies->setSpMovtTraits(move);
	return error;
}

//---------------------------------------------------------------------------
int ReadSettlement(int option)
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

	if (option == 0) { // open file and read header line
		settFile.open(settleFile.c_str());
		string header;
		int nheaders = 0;
		if (trfr.usesMovtProc) nheaders = 14;
		else nheaders = 7;
		for (int i = 0; i < nheaders; i++) {
			settFile >> header;
		}
		return 0;
	}
	if (option == 9) { // close file
		if (settFile.is_open()) {
			settFile.close(); settFile.clear();
		}
		return 0;
	}

	isFirstline = true;

	// set no.of lines assuming maximum stage- and sex-dependency
	if (sstruct.nStages == 0) Nlines = gNbSexesDisp;
	else Nlines = sstruct.nStages * gNbSexesDisp;

	for (int line = 0; line < Nlines; line++) {

		settFile >> simNb >> inStageDep >> inSexDep >> inStage >> inSex;
		if (!trfr.usesMovtProc)
		{ // dispersal kernel
			settFile >> inSettleType >> inFindMate;
		}
		else {
			settFile >> inDensDep >> inIndVar >> inFindMate;
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

			settFile >> ssteps.minSteps >> ssteps.maxSteps >> ssteps.maxStepsYr;
			settFile >> settleDD.s0 >> settleDD.alpha >> settleDD.beta;

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
				srules.go2nbrLocn = false;
				break;
			case 1:
				srules.wait = true;
				srules.go2nbrLocn = false;
				break;
			case 2:
				srules.wait = false;
				srules.go2nbrLocn = true;
				break;
			case 3:
				srules.wait = true;
				srules.go2nbrLocn = true;
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
int ReadInitialisation(int option, Landscape* pLandscape)
{
	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	initParams init = paramsInit->getInit();
	string inputDir = paramsSim->getDir(1);

	int simNb, maxcells;
	float totalProps;
	int error = 0;

	if (option == 0) { // open file and read header line
		initFile.open(initialFile.c_str());
		string header;
		int nheaders = 17;
		if (dem.stageStruct) nheaders += 1 + sstruct.nStages - 1;
		for (int i = 0; i < nheaders; i++) initFile >> header;
		return 0;
	}
	if (option == 9) { // close file
		if (initFile.is_open()) {
			initFile.close(); initFile.clear();
		}
		return 0;
	}

	initFile >> simNb >> init.seedType >> init.freeType >> init.spDistType;

	if (init.seedType == 1 && !paramsLand.spDist) 
		error = 601;

	if (paramsLand.patchModel) 
		initFile >> init.initDens >> init.indsHa;
	else 
		initFile >> init.initDens >> init.indsCell;

	initFile >> init.minSeedX >> init.maxSeedX 
		>> init.minSeedY >> init.maxSeedY
		>> init.nSeedPatches >> init.nSpDistPatches
		>> init.initFrzYr >> init.restrictRows
		>> init.restrictFreq >> init.finalFrzYr
		>> init.indsFile;

	init.restrictRange = (init.seedType == 0 && init.restrictRows > 0);

	if (dem.stageStruct) {
		float propStage;
		initFile >> init.initAge;
		totalProps = 0.0;
		for (int stg = 1; stg < sstruct.nStages; stg++) {
			initFile >> propStage;
			if(init.seedType!=2){
				totalProps += propStage;
				paramsInit->setProp(stg, propStage);
			}
		}
		if (init.seedType!=2 && totalProps != 1.0)
		{ 
			throw logic_error("The proportion of initial individuals in each stage doesn not sum to 1.");
		}
	}

	paramsInit->setInit(init);

	switch (init.seedType) {
	case 0: // free initialisation
		if (init.minSeedX == gEmptyVal)
			init.minSeedX = 0;
		if (init.minSeedY == gEmptyVal)
			init.minSeedY = 0;
		if (init.maxSeedX == gEmptyVal) 
			init.maxSeedX = paramsLand.maxX;
		if (init.maxSeedY == gEmptyVal) 
			init.maxSeedY = paramsLand.maxY;
		if (init.minSeedY > init.maxSeedY || init.minSeedX > init.maxSeedX) {
			error = 603;
		}
		maxcells = (init.maxSeedY - init.minSeedY) * (init.maxSeedX - init.minSeedX);
		if (init.freeType == 0 && init.nSeedPatches > maxcells) 
			error = 602;
		break;
	case 1: // from species distribution
		// nothing to do here
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
	return error;
}

//---------------------------------------------------------------------------
int ReadInitIndsFile(int option, Landscape* pLandscape, string indsfile) {
	string header;
	landParams paramsLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	initParams init = paramsInit->getInit();

	if (option == 0) { // open file and read header line
		initIndsFile.open(indsfile.c_str());
		string header;
		int nheaders = 3;
		if (paramsLand.patchModel) nheaders++;
		else nheaders += 2;
		if (dem.repType > 0) nheaders++;
		if (dem.stageStruct) nheaders += 2;
		for (int i = 0; i < nheaders; i++) initIndsFile >> header;
		paramsInit->resetInitInds();
		//	return 0;
	}

	if (option == 9) { // close file
		if (initIndsFile.is_open()) {
			initIndsFile.close(); initIndsFile.clear();
		}
		return 0;
	}

	// Read data lines;
	initInd iind;
	int ninds;
	int totinds = 0;

	iind.year = gEmptyVal;
	initIndsFile >> iind.year;
	bool must_stop = (iind.year == gEmptyVal);

	while (!must_stop) {
		initIndsFile >> iind.species;

		if (paramsLand.patchModel) {
			initIndsFile >> iind.patchID;
			iind.x = iind.y = 0;
		}
		else {
			initIndsFile >> iind.x >> iind.y; 
			iind.patchID = 0;
		}
		initIndsFile >> ninds;

		if (dem.repType > 0) 
			initIndsFile >> iind.sex;
		else 
			iind.sex = 0;

		if (dem.stageStruct) {
			initIndsFile >> iind.age >> iind.stage;
		}
		else {
			iind.age = iind.stage = 0;
		}
		for (int i = 0; i < ninds; i++) {
			totinds++;
			paramsInit->addInitInd(iind);
		}

		iind.year = gEmptyVal;
		initIndsFile >> iind.year;
		if (iind.year == gEmptyVal || initIndsFile.eof())
			must_stop = true;
	} // end of while loop

	if (initIndsFile.is_open()) initIndsFile.close();
	initIndsFile.clear();

	return totinds;
}

//---------------------------------------------------------------------------
void RunBatch(int nSimuls, int nLandscapes)
{
	int land_nr;
	int t0, t1, t00, t01;
	int read_error;
	bool params_ok;
	simParams sim = paramsSim->getSim();

	Landscape* pLandscape = NULL;  		// pointer to landscape

	t0 = (int)time(0);

	string name = paramsSim->getDir(2) + "Batch" + to_string(sim.batchNum) + "_RS_log.csv";
	if (rsLog.is_open()) {
		rsLog.close(); rsLog.clear();
	}
	rsLog.open(name.c_str());
	if (!rsLog.is_open()) {
		cout << endl << "Error - unable to open Batch" << sim.batchNum
			<< "_RS_log.csv file - aborting batch run" << endl;
		return;
	}
	rsLog << "Event,Number,Reps,Years,Time" << endl;
	rsLog << "RANDOM SEED," << RS_random_seed << ",,," << endl;

	// Open landscape batch file and read header record
	if (ReadLandFile(0)) {
		cout << endl << "Error opening landFile - aborting batch run" << endl;
		return;
	}

	for (int j = 0; j < nLandscapes; j++) {
		// create new landscape
		if (pLandscape != NULL) delete pLandscape;
		pLandscape = new Landscape;
		bool landOK = true;

		t00 = (int)time(0);
		land_nr = ReadLandFile(1, pLandscape);
		if (land_nr <= 0) { // error condition
			string msg = "Error code " + to_string(-land_nr)
				+ " returned from reading LandFile - aborting batch run";
			cout << endl << msg << endl;
			ReadLandFile(9); // close the landscape file
			return;
		}
		landParams paramsLand = pLandscape->getLandParams();
		paramsLand.patchModel = patchmodel;
		paramsLand.resol = resolution;
		paramsLand.rasterType = landtype;
		if (landtype == 9) {
			paramsLand.generated = true;
			paramsLand.nHab = 2;
		}
		else {
			paramsLand.generated = false;
			if (name_dynland == "NULL") paramsLand.dynamic = false;
			else paramsLand.dynamic = true;
		}
		paramsLand.nHabMax = maxNhab;
		paramsLand.spDist = speciesdist;
		paramsLand.spResol = distresolution;
		pLandscape->setLandParams(paramsLand, sim.batchMode);

		if (landtype != 9) { // imported landscape
			string hname = paramsSim->getDir(1) + name_landscape;
			int landcode;
			string cname;
			if (gNameCostFile == "NULL" || gNameCostFile == "none") cname = "NULL";
			else cname = paramsSim->getDir(1) + gNameCostFile;
			if (paramsLand.patchModel) {
				string pname = paramsSim->getDir(1) + name_patch;
#if RSDEBUG
				time_t t02a = time(0);
#endif
				landcode = pLandscape->readLandscape(0, hname, pname, cname);
#if RSDEBUG
				time_t t02b = time(0);
				DEBUGLOG << "RunBatch(): TIME for readLandscape() " << t02b - t02a << endl;
#endif
			}
			else {
				landcode = pLandscape->readLandscape(0, hname, " ", cname);
			}
			if (landcode != 0) {
				rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
				cout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
				landOK = false;
			}
			if (paramsLand.dynamic) {
#if RSDEBUG
				time_t t03a = time(0);
#endif
				landcode = ReadDynLandFile(pLandscape);
#if RSDEBUG
				time_t t03b = time(0);
				DEBUGLOG << "RunBatch(): TIME for ReadDynLandFile() " << t03b - t03a << endl;
#endif
				if (landcode != 0) {
					rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
					cout << endl << "Error reading landscape " << land_nr << " - aborting" << endl;
					landOK = false;
				}
			}
			if (landtype == 0) {
				pLandscape->updateHabitatIndices();
			}

			// species distribution

			if (paramsLand.spDist) { // read initial species distribution
				string distname = paramsSim->getDir(1) + name_sp_dist;
				landcode = pLandscape->newDistribution(pSpecies, distname);
				if (landcode == 0) {
				}
				else {
					rsLog << "Landscape," << land_nr << ",ERROR,CODE," << landcode << endl;
					cout << endl << "Error reading initial distribution for landscape "
						<< land_nr << " - aborting" << endl;
					landOK = false;
				}
			}
			paramsSim->setSim(sim);

			if (landOK) {
				t01 = static_cast<int>(time(0));
				rsLog << "Landscape," << land_nr << ",,," << t01 - t00 << endl;

			} // end of landOK condition

		} // end of imported landscape

		if (landOK) {

			// Open all other batch files and read header records
			if (ReadParameters(0, pLandscape)) {
				cout << endl << "Error opening ParameterFile - aborting batch run" << endl;
				return;
			}
			if (stagestruct) {
				ReadStageStructure(0);
			}
			ReadEmigration(0);
			ReadTransferFile(0, pLandscape);
			ReadSettlement(0);
			ReadInitialisation(0, pLandscape);

			if (gHasGenetics) {
				ifsGenetics.open(geneticsFile.c_str());
				flushHeader(ifsGenetics);
			}

			// nSimuls is the total number of lines (simulations) in
			// the batch and is set in the control function
			string msgsim = "Simulation,";
			string msgerr = ",ERROR CODE,";
			string msgabt = ",simulation aborted";

			for (int i = 0; i < nSimuls; i++) {

				t00 = (int)time(0);
				params_ok = true;
				read_error = ReadParameters(1, pLandscape);
				simParams sim = paramsSim->getSim();
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				if (stagestruct) {
					ReadStageStructure(1);
				}
				read_error = ReadEmigration(1);
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				read_error = ReadTransferFile(1, pLandscape);
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				read_error = ReadSettlement(1);
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}
				read_error = ReadInitialisation(1, pLandscape);
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}

				if (gHasGenetics) {
					read_error = ReadGeneticsFile(ifsGenetics, pLandscape);
					if (read_error) {
						rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
						params_ok = false;
					}

					if (!ifsTraits.is_open()) {
						// First simulation
						ifsTraits.open(gPathToTraitsFile.c_str());
						flushHeader(ifsTraits);
					}
					read_error = ReadTraitsFile(ifsTraits, i);
					if (read_error) {
						rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
						params_ok = false;
					}
				}
				
				if (params_ok) {
					simParams sim = paramsSim->getSim();

					cout << endl << "Running simulation nr. " << to_string(sim.simulation)
						<< " on landscape no. " << to_string(land_nr) << endl;

					// for batch processing, include landscape number in parameter file name
					OutParameters(pLandscape);
					RunModel(pLandscape, i);

					t01 = (int)time(0);
					rsLog << msgsim << sim.simulation << "," << sim.reps
						<< "," << sim.years << "," << t01 - t00 << endl;
				} // end of if (params_ok)
				else {
					cout << endl << "Error in reading parameter file(s)" << endl;
				}
			} // end of nSimuls for loop

			// close input files
			ReadParameters(9, pLandscape);
			if (stagestruct) 
				ReadStageStructure(9);
			ReadEmigration(9);
			ReadTransferFile(9, pLandscape);
			ReadSettlement(9);
			ReadInitialisation(9, pLandscape);
			if (gHasGenetics) {
				ifsGenetics.close();
				ifsGenetics.clear();
				ifsTraits.close();
				ifsTraits.clear();
			}

			if (pLandscape != NULL)
			{
				delete pLandscape; 
				pLandscape = NULL;
			}

		} // end of landOK condition

	} // end of nLandscapes loop

	ReadLandFile(9); // close the landFile

	// Write performance data to log file
	t1 = (int)time(0);
	rsLog << endl << "Batch,,,," << t1 - t0 << endl;

	if (rsLog.is_open()) {
		rsLog.close(); rsLog.clear();
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


