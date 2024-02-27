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

ifstream controlfile;
// Note - all batch files are prefixed 'b' here for reasons concerned with RS v1.0
ifstream bParamFile, bLandFile, bDynLandFile;
ifstream bSpDistFile, bStageStructFile, bTransMatrix;
ifstream bStageWeightsFile;
ifstream bEmigrationFile, bTransferFile, bSettlementFile;
ifstream bTraitsFile, bGeneticsFile;
ifstream bInitFile, bInitIndsFile;

ofstream batchlog;

ofstream rsLog; // performance log for recording simulation times, etc.

// NOTE: THE STREAMS USED TO READ THE DATA AT RUN TIME COULD TAKE THE SAME NAMES AS
// USED DURING PARSING (ABOVE)
ifstream parameters;
ifstream ssfile, tmfile, fdfile, ddfile, sdfile;
ifstream emigFile, transFile, settFile, initFile, initIndsFile;
ifstream landfile, dynlandfile;

// global variables passed between parsing functions...
int batchnum;
int patchmodel, resolution, landtype, maxNhab, speciesdist, distresolution;
int reproductn;
int repseasons;
int stagestruct, stages, gTransferType;
int sexesDem;		// no. of explicit sexes for demographic model
int gNbSexesDisp;	// no. of explicit sexes for dispersal model
int gFirstSimNb = 0; // BAD, globals should not be modified.
int fileNtraits; // no. of traits defined in genetic architecture file
//rasterdata landraster,patchraster,spdistraster;
rasterdata landraster;
// ...including names of the input files
string parameterFile;
string landFile;
string name_landscape, name_patch, name_dynland, name_sp_dist, gNameCostFile;
string stageStructFile, transMatrix;
string emigrationFile, transferFile, settleFile, geneticsFile, traitsFile, initialFile;
string prevInitialIndsFile = " ";

string msgnlines = "No. of lines for final Simulation ";
string msgshldbe = " should be ";
string msgresol0 = "*** Resolution of ";
string msgresol1 = " does not match Resolution in Control file ";
string msghdrs0 = "*** Headers of ";
string msghdrs1 = " do not match headers of LandscapeFile";
string msgpatch = " is required for patch-based model";
string msgmatch = " must match the specification exactly";
string msgcase = " case-sensitive parameter names";

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
batchfiles ParseControlFile(string ctrlfile, string indir, string outdir)
{
	batchfiles b;
	int lines, nSimuls;
	int errors = 0;
	string paramname, filename, fname, logname, header;
	string filetype = "Control file";
	bool controlFormatError = false;
	b.ok = true; b.nSimuls = 0; b.nLandscapes = 0;

	// open batch log file
	logname = outdir + "BatchLog.txt";
	batchlog.open(logname.c_str());
	if (!batchlog.is_open()) {
		//	MessageDlg("Error opening batch output log file",mtError, TMsgDlgButtons() << mbOK,0);
		cout << "Error opening batch output log file " << logname << endl;
		b.ok = false;
		return b;
	}

	controlfile.open(ctrlfile.c_str());

	if (!controlfile.is_open()) {
		//	MessageDlg("Error opening Control file",mtError, TMsgDlgButtons() << mbOK,0);
		cout << "Error opening Control file: " << ctrlfile << endl;
		batchlog << "Error opening Control file: " << ctrlfile << endl;
		b.ok = false;
		if (batchlog.is_open()) { batchlog.close(); batchlog.clear(); }
		return b;
	}
	else {
		batchlog << "Checking Control file " << ctrlfile << endl;
	}

	// Check fixed model parameters

	controlfile >> paramname >> batchnum;
	if (paramname == "BatchNum") {
		if (batchnum < 0) {
			BatchError(filetype, -999, 19, "BatchNum"); errors++;
		}
		else b.batchNum = batchnum;
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> patchmodel;
	if (paramname == "PatchModel") {
		if (patchmodel < 0 || patchmodel > 1) {
			BatchError(filetype, -999, 1, "PatchModel"); errors++;
		}
		else b.patchmodel = patchmodel;
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> resolution;
	if (paramname == "Resolution") {
		if (resolution < 1) {
			BatchError(filetype, -999, 11, "Resolution"); errors++;
		}
		else b.resolution = resolution;
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> landtype;
	if (paramname == "LandType") {
		if (landtype != 0 && landtype != 2 && landtype != 9) {
			BatchError(filetype, -999, 0, "LandType");
			batchlog << "LandType must be 0, 2 or 9" << endl;
			errors++;
		}
		else {
			if (landtype == 9 && patchmodel) {
				BatchError(filetype, -999, 0, "LandType");
				batchlog << "LandType may not be 9 for a patch-based model" << endl;
				errors++;
			}
			else b.landtype = landtype;
		}
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> maxNhab;
	if (paramname == "MaxHabitats") {
		if (landtype == 0) { // raster with unique habitat codes
			if (maxNhab < 2) {
				BatchError(filetype, -999, 12, "MaxHabitats"); errors++;
			}
			else b.maxNhab = maxNhab;
		}
		else { // raster with habitat quality OR artificial landscape
			if (maxNhab != 1) {
				BatchError(filetype, -999, 0, " "); errors++;
				batchlog << "MaxHabitats must be 1 for LandType = " << landtype << endl;
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
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> speciesdist;
	if (paramname == "SpeciesDist") {
		if (speciesdist < 0 || speciesdist > 1) {
			BatchError(filetype, -999, 1, "SpeciesDist"); errors++;
		}
		else {
			if (speciesdist != 0 && landtype == 9) {
				BatchError(filetype, -999, 0, "SpeciesDist");
				batchlog << "SpeciesDist must be 0 for an artificial landscape" << endl;
				errors++;

			}
			else b.speciesdist = speciesdist;
		}
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> distresolution;
	if (paramname == "DistResolution") {
		if (speciesdist == 1) { // distribution resolution is required
			if (distresolution < resolution) {
				BatchError(filetype, -999, 0, "DistResolution");
				batchlog << "DistResolution may not be less than Resolution" << endl;
				errors++;
			}
			else {
				if (distresolution % resolution) {
					BatchError(filetype, -999, 0, "DistResolution");
					batchlog << "DistResolution must be an integer multiple of Resolution" << endl;
					errors++;
				}
				else b.distresolution = distresolution;
			}
		}
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> reproductn;
	sexesDem = gNbSexesDisp = 0;
	if (paramname == "Reproduction") {
		if (reproductn < 0 || reproductn > 2) {
			BatchError(filetype, -999, 2, "Reproduction"); errors++;
		}
		else {
			switch (reproductn) {
			case 0: { sexesDem = 1; gNbSexesDisp = 1; break; }
			case 1: { sexesDem = 1; gNbSexesDisp = 2; break; }
			case 2: { sexesDem = 2; gNbSexesDisp = 2; break; }
			}
			b.reproductn = reproductn; b.sexesDem = sexesDem; b.gNbSexesDisp = gNbSexesDisp;
		}
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> repseasons;
	if (paramname == "RepSeasons") {
		if (repseasons < 1) {
			BatchError(filetype, -999, 11, "RepSeasons"); errors++;
		}
		else b.repseasons = repseasons;
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> stagestruct;
	if (paramname == "StageStruct") {
		if (stagestruct < 0 || stagestruct > 1) {
			BatchError(filetype, -999, 1, "StageStruct"); errors++;
		}
		else b.stagestruct = stagestruct;
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> stages;
	if (paramname == "Stages") {
		if (stagestruct) {
			if (stages < 2 || stages > 10) {
				BatchError(filetype, -999, 0, " "); errors++;
				batchlog << "Stages must be between 2 and 10" << endl;
			}
			b.stages = stages;
		}
		else { // non-stage-structured model must have 2 stages
			b.stages = stages = 2;
		}
	}
	else controlFormatError = true; // wrong control file format

	controlfile >> paramname >> gTransferType;
	if (paramname == "Transfer") {
		if (gTransferType < 0 || gTransferType > 2) {
			BatchError(filetype, -999, 2, "Transfer"); errors++;
		}
		else b.transfer = gTransferType;
	}
	else controlFormatError = true; // wrong control file format

	if (controlFormatError || errors > 0) { // terminate batch error checking
		if (controlFormatError) {
			CtrlFormatError();
		}
		batchlog << endl
			<< "*** Model parameters in Control file must be corrected before further input file checks are conducted"
			<< endl;
		batchlog.close(); batchlog.clear();
		b.ok = false;
		controlfile.close(); controlfile.clear();
		return b;
	}

	// Check parameter file
	controlfile >> paramname >> filename;
	if (paramname == "ParameterFile" && !controlFormatError) {
		fname = indir + filename;
		batchlog << endl << "Checking " << paramname << " " << fname << endl;
		bParamFile.open(fname.c_str());
		if (bParamFile.is_open()) {
			b.nSimuls = ParseParameterFile();
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
			batchlog << endl
				<< "*** ParameterFile must be corrected before further input file checks are conducted"
				<< endl;
			batchlog.close(); batchlog.clear();
			b.ok = false;
			controlfile.close(); controlfile.clear();
			return b;
		}
	}
	else controlFormatError = true; // wrong control file format
	if (bParamFile.is_open()) bParamFile.close();
	bParamFile.clear();

	// Check land file
	controlfile >> paramname >> filename;
	if (paramname == "LandFile" && !controlFormatError) {
		fname = indir + filename;
		batchlog << endl << "Checking " << paramname << " " << fname << endl;
		bLandFile.open(fname.c_str());
		if (bLandFile.is_open()) {
			lines = ParseLandFile(landtype, indir);
			if (lines < 0) {
				b.ok = false;
				if (lines < -111)
					batchlog << "*** Format error in " << paramname << endl;
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
	else controlFormatError = true; // wrong control file format

	// Check stage structure file if required file
	controlfile >> paramname >> filename;
	batchlog << endl;
	if (paramname == "StageStructFile" && !controlFormatError) {
		if (filename == "NULL") {
			if (stagestruct) {
				batchlog << "*** File name is required for " << paramname << endl;
				b.ok = false;
			}
			else b.stageStructFile = filename;
		}
		else { // filename is not NULL
			if (stagestruct) { // check file only if it is required
				fname = indir + filename;
				batchlog << "Checking " << paramname << " " << fname << endl;
				bStageStructFile.open(fname.c_str());
				if (bStageStructFile.is_open()) {
					nSimuls = ParseStageFile(indir);
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
					batchlog << "*** File name for stageStructFile should be NULL as StageStruct = "
						<< stagestruct << endl;
					b.ok = false;
				}
			}
		}
	}
	else controlFormatError = true; // wrong control file format

	// Check emigration file
	controlfile >> paramname >> filename;
	if (paramname == "EmigrationFile" && !controlFormatError) {
		fname = indir + filename;
		batchlog << endl << "Checking " << paramname << " " << fname << endl;
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
	else controlFormatError = true; // wrong control file format

	// Check transfer file
	controlfile >> paramname >> filename;
	if (paramname == "TransferFile" && !controlFormatError) {
		fname = indir + filename;
		batchlog << endl << "Checking " << paramname << " " << fname << endl;
		bTransferFile.open(fname.c_str());
		if (bTransferFile.is_open()) {
			nSimuls = ParseTransferFile(indir);
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
	else controlFormatError = true; // wrong control file format

	// Check settlement file
	controlfile >> paramname >> filename;
	if (paramname == "SettlementFile" && !controlFormatError) {
		fname = indir + filename;
		batchlog << endl << "Checking " << paramname << " " << fname << endl;
		bSettlementFile.open(fname.c_str());
		if (bSettlementFile.is_open()) {
			nSimuls = ParseSettleFile();
			if (nSimuls < 0) {
				b.ok = false;
			}
			else {
				FileOK(paramname, nSimuls, 0);
				if (nSimuls != b.nSimuls) {
					SimulnCountError(filename); b.ok = false;
				}
				else b.settleFile = fname;
			}
			bSettlementFile.close();
		}
		else {
			OpenError(paramname, fname); b.ok = false;
		}
		bSettlementFile.clear();
	}
	else controlFormatError = true; // wrong control file format

	// Check genetics file if required file
	controlfile >> paramname >> filename;
	batchlog << endl;
	if (paramname == "GeneticsFile" && !controlFormatError) {
		if (filename == "NULL") {

			batchlog << "No genetics required " << paramname << endl;

		}
		else { // filename is not NULL
			fname = indir + filename;
			batchlog << "Checking " << paramname << " " << fname << endl;
			bGeneticsFile.open(fname.c_str());
			if (bGeneticsFile.is_open()) {
				nSimuls = ParseGeneticsFile(indir);
				if (nSimuls < 0) {
					b.ok = false;
				}
				else {
					FileOK(paramname, nSimuls, 0);
					if (nSimuls != b.nSimuls) {
						SimulnCountError(filename); b.ok = false;
					}
					else b.geneticsFile = fname;
				}
				bGeneticsFile.close();
			}
			else {
				OpenError(paramname, fname); b.ok = false;
			}
			bGeneticsFile.clear();
		}
	}
	else controlFormatError = true; // wrong control file format

	// Check initialisation file
	controlfile >> paramname >> filename;
	if (paramname == "InitialisationFile" && !controlFormatError) {
		fname = indir + filename;
		batchlog << endl << "Checking " << paramname << " " << fname << endl;
		bInitFile.open(fname.c_str());
		if (bInitFile.is_open()) {
			nSimuls = ParseInitFile(indir);
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
	else controlFormatError = true; // wrong control file format

	if (controlFormatError) {
		CtrlFormatError();
		b.ok = false;
	}

	if (controlfile.is_open()) { controlfile.close(); controlfile.clear(); }
	if (batchlog.is_open()) { batchlog.close(); batchlog.clear(); }

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
int ParseParameterFile(void)
{
	string header, Kheader, intext;
	int i, inint, replicates, years;
	int absorb, gradient, shifting, shiftstart, shiftend, envstoch, stochtype;
	int localext, savemaps;
	int prevsimul = 0;
	float infloat, minR, maxR, minK, maxK, sum_K, min_K, max_K;
	int errors = 0;
	int Kerrors = 0;
	string filetype = "ParameterFile";

	//batchlog << "ParseParametersFile(): starting " << endl;
	// Parse header line;
	bParamFile >> header; if (header != "Simulation") errors++;
	bParamFile >> header; if (header != "Replicates") errors++;
	bParamFile >> header; if (header != "Years") errors++;
	bParamFile >> header; if (header != "Absorbing") errors++;
	bParamFile >> header; if (header != "Gradient") errors++;
	bParamFile >> header; if (header != "GradSteep") errors++;
	bParamFile >> header; if (header != "Optimum") errors++;
	bParamFile >> header; if (header != "f") errors++;
	bParamFile >> header; if (header != "LocalExtOpt") errors++;
	bParamFile >> header; if (header != "Shifting") errors++;
	bParamFile >> header; if (header != "ShiftRate") errors++;
	bParamFile >> header; if (header != "ShiftStart") errors++;
	bParamFile >> header; if (header != "ShiftEnd") errors++;
	bParamFile >> header; if (header != "EnvStoch") errors++;
	bParamFile >> header; if (header != "EnvStochType") errors++;
	bParamFile >> header; if (header != "ac") errors++;
	bParamFile >> header; if (header != "std") errors++;
	bParamFile >> header; if (header != "minR") errors++;
	bParamFile >> header; if (header != "maxR") errors++;
	bParamFile >> header; if (header != "minK") errors++;
	bParamFile >> header; if (header != "maxK") errors++;
	bParamFile >> header; if (header != "LocalExt") errors++;
	bParamFile >> header; if (header != "LocalExtProb") errors++;
	bParamFile >> header; if (header != "PropMales") errors++;
	bParamFile >> header; if (header != "Harem") errors++;
	bParamFile >> header; if (header != "bc") errors++;
	bParamFile >> header; if (header != "Rmax") errors++;
	for (i = 0; i < maxNhab; i++) {
		Kheader = "K" + Int2Str(i + 1);
		bParamFile >> header; if (header != Kheader) Kerrors++;
	}
	bParamFile >> header; if (header != "OutStartPop") errors++;
	bParamFile >> header; if (header != "OutStartInd") errors++;
	bParamFile >> header; if (header != "OutStartTraitCell") errors++;
	bParamFile >> header; if (header != "OutStartTraitRow") errors++;
	bParamFile >> header; if (header != "OutStartConn") errors++;
	bParamFile >> header; if (header != "OutIntRange") errors++;
	bParamFile >> header; if (header != "OutIntOcc") errors++;
	bParamFile >> header; if (header != "OutIntPop") errors++;
	bParamFile >> header; if (header != "OutIntInd") errors++;
	bParamFile >> header; if (header != "OutIntTraitCell") errors++;
	bParamFile >> header; if (header != "OutIntTraitRow") errors++;
	bParamFile >> header; if (header != "OutIntConn") errors++;
	bParamFile >> header; if (header != "SaveMaps") errors++;
	bParamFile >> header; if (header != "MapsInterval") errors++;
	bParamFile >> header; if (header != "SMSHeatMap") errors++;
	bParamFile >> header; if (header != "DrawLoadedSp") errors++;
	bParamFile >> header; if (header != "FixReplicateSeed") errors++;
	if (errors > 0 || Kerrors > 0) {
		FormatError(filetype, errors);
		batchlog << "*** Ensure column headers are correct to continue checking data" << endl;
		if (Kerrors > 0) {
			BatchError(filetype, -999, 333, "K");
		}
		return -111;
	}

	// Parse data lines
	int line = 1;
	int nSimuls = 0;
	inint = -98765;
	bParamFile >> inint; // first simulation number
	if (inint < 0) {
		batchlog << "*** Error in ParameterFile - first simulation number must be >= 0" << endl;
		errors++;
	}
	else {
		prevsimul = gFirstSimNb = inint; 
		nSimuls++;
	}
	while (inint != -98765) {
		bParamFile >> replicates; if (replicates <= 0) { BatchError(filetype, line, 11, "Replicates"); errors++; }
		bParamFile >> years; if (years <= 0) { BatchError(filetype, line, 11, "Years"); errors++; }
		bParamFile >> absorb;
		if (absorb < 0 || absorb > 1) { BatchError(filetype, line, 1, "Absorbing"); errors++; }
		bParamFile >> gradient;
		if (patchmodel) {
			if (gradient != 0) {
				BatchError(filetype, line, 0, " ");
				batchlog << "Gradient must be 0 for patch-based model" << endl;
				errors++;
				gradient = 0; // to prevent checking of subsequent fields
			}
			gradient = 0; // to prevent unnecessary checking of subsequent fields
		}
		else { // cell-based model
			if (gradient < 0 || gradient > 3) {
				BatchError(filetype, line, 0, " ");
				batchlog << "Gradient must be between 0 and 3 for cell-based model" << endl;
				errors++;
			}
		}
		bParamFile >> infloat;
		if (gradient && infloat < 0.0) { BatchError(filetype, line, 19, "GradSteep"); errors++; }
		bParamFile >> inint;
		if (gradient && inint < 0) { BatchError(filetype, line, 19, "Optimum"); errors++; }
		bParamFile >> infloat;
		if (gradient && infloat < 0.0) { BatchError(filetype, line, 19, "f"); errors++; }
		bParamFile >> infloat;
		if (gradient == 4 && (infloat < 0.0 || infloat >= 1.0))
		{
			BatchError(filetype, line, 20, "LocalExtOpt"); errors++;
		}
		bParamFile >> shifting;
		if (gradient && (shifting < 0 || shifting > 1)) { BatchError(filetype, line, 1, "Shifting"); errors++; }
		bParamFile >> infloat;
		if (gradient && shifting && infloat <= 0.0) { BatchError(filetype, line, 10, "ShiftRate"); errors++; }
		bParamFile >> shiftstart;
		if (gradient && shifting && shiftstart <= 0) { BatchError(filetype, line, 10, "ShiftStart"); errors++; }
		bParamFile >> shiftend;
		if (gradient && shifting && shiftend <= shiftstart) {
			BatchError(filetype, line, 0, " ");
			batchlog << "ShiftEnd must be greater than ShiftStart" << endl;
			errors++;
		}
		bParamFile >> envstoch;
		if (patchmodel == 0) { // cell-based model
			if (envstoch < 0 || envstoch > 2) {
				BatchError(filetype, line, 0, " ");
				batchlog << "EnvStoch must be 0, 1 or 2 for cell-based model" << endl;
				errors++;
				//			envstoch = 0; // to prevent checking of subsequent fields
			}
		}
		else { // patch-based model
			if (envstoch < 0 || envstoch > 1) {
				BatchError(filetype, line, 0, " ");
				batchlog << "EnvStoch must be 0 or 1 for patch-based model" << endl;
				errors++;
				//			envstoch = 0; // to prevent checking of subsequent fields
			}
		}
		bParamFile >> stochtype;
		if (envstoch && (stochtype < 0 || stochtype > 1)) {
			BatchError(filetype, line, 1, "EnvStochType"); errors++;
		}
		bParamFile >> infloat;
		if (envstoch && (infloat < 0.0 || infloat >= 1.0)) { BatchError(filetype, line, 20, "ac"); errors++; }
		bParamFile >> infloat;
		if (envstoch && (infloat <= 0.0 || infloat > 1.0)) { BatchError(filetype, line, 20, "std"); errors++; }
		bParamFile >> minR;
		if (envstoch && stochtype == 0 && minR <= 0.0) { BatchError(filetype, line, 10, "minR"); errors++; }
		bParamFile >> maxR;
		if (envstoch && stochtype == 0 && maxR <= minR) {
			BatchError(filetype, line, 0, " ");
			batchlog << "maxR must be greater than minR" << endl;
			errors++;
		}
		bParamFile >> minK >> maxK;
		if (envstoch && stochtype == 1) {
			if (minK <= 0.0) { BatchError(filetype, line, 10, "minK"); errors++; }
			if (maxK <= minK) {
				BatchError(filetype, line, 0, " ");
				batchlog << "maxK must be greater than minK" << endl;
				errors++;
			}
		}
		bParamFile >> localext;
		if (patchmodel == 0) { // cell-based model
			if (localext < 0 || localext > 1) {
				BatchError(filetype, line, 1, "LocalExt");
				errors++;
			}
			else {
				if (gradient == 4) { // gradient in local extinction probability
					if (localext != 0) {
						BatchError(filetype, line, 0, " ");
						batchlog << "LocalExt must be zero if Gradient is 4" << endl;
						errors++;
					}
				}
			}
		}
		else { // patch-based model
			if (localext != 0) {
				BatchError(filetype, line, 0, "null");
				batchlog << "LocalExt must be 0 for patch-based model" << endl;
				errors++;
			}
		}
		bParamFile >> infloat;
		if (patchmodel == 0 && localext == 1 && (infloat <= 0.0 || infloat >= 1.0))
		{
			BatchError(filetype, line, 20, "LocalExtProb"); errors++;
		}
		bParamFile >> infloat;
		if (reproductn && (infloat <= 0.0 || infloat >= 1.0)) {
			BatchError(filetype, line, 20, "PropMales"); errors++;
		}
		bParamFile >> infloat;
		if (reproductn == 2 && infloat <= 0.0) { BatchError(filetype, line, 10, "Harem"); errors++; }
		bParamFile >> infloat;
		if (stagestruct == 0 && infloat <= 0.0) { BatchError(filetype, line, 10, "bc"); errors++; }
		bParamFile >> infloat;
		if (stagestruct == 0 && infloat <= 0.0) { BatchError(filetype, line, 10, "Rmax"); errors++; }
		sum_K = 0.0; min_K = 9999999.0; max_K = 0.0;
		for (i = 0; i < maxNhab; i++) {
			bParamFile >> infloat;
			if (infloat < 0.0) {
				Kheader = "K" + Int2Str(i + 1);
				BatchError(filetype, line, 19, Kheader); errors++;
			}
			else {
				sum_K += infloat;
				if (infloat > 0.0) {
					if (infloat < min_K) min_K = infloat;
					if (infloat > max_K) max_K = infloat;
				}
			}
		}
		if (sum_K <= 0.0) {
			BatchError(filetype, line, 0, " "); errors++;
			batchlog << "At least one K column must be non-zero" << endl;
		}
		else {
			if (envstoch && stochtype == 1) { // environmental stochasticity in K
				if (min_K < minK || max_K > maxK) {
					BatchError(filetype, line, 0, " "); errors++;
					batchlog << "Non-zero K values must lie between minK and maxK" << endl;
				}
			}
		}

		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutStartPop"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutStartInd"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutStartTraitCell"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutStartTraitRow"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutStartConn"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntRange"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntOcc"); errors++; }
		else {
			if (landtype == 9) {
				if (inint > 0) {
					BatchError(filetype, line, 0, " "); errors++;
					batchlog << "OutIntOcc must be zero for a generated landscape" << endl;
				}
			}
			else {
				if (replicates < 2 && inint > 0) {
					BatchError(filetype, line, 0, " "); errors++;
					batchlog << "OutIntOcc may be non-zero only if Replicates >= 2" << endl;
				}
			}
		}
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntPop"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntInd"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntTraitCell"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntTraitRow"); errors++; }
		bParamFile >> inint;
		if (inint < 0) { BatchError(filetype, line, 19, "OutIntConn"); errors++; }
		else {
			if (patchmodel != 1 && inint > 0) {
				BatchError(filetype, line, 0, " ");
				batchlog << "OutIntConn may be >0 only if PatchModel is 1" << endl;
				errors++;
			}
		}
		bParamFile >> savemaps; if (savemaps < 0 || savemaps > 1)
		{
			BatchError(filetype, line, 1, "SaveMaps"); errors++;
		}
		bParamFile >> inint; if (savemaps == 1 && inint < 1) {
			BatchError(filetype, line, 11, "MapsInterval");
			errors++;
		}
		bParamFile >> inint; if (inint < 0 || inint > 1) {
			BatchError(filetype, line, 1, "SMSHeatMap");
			errors++;
		}
		bParamFile >> inint; if (savemaps == 1 && (inint < 0 || inint > 1)) {
			BatchError(filetype, line, 1, "DrawLoadedSp");
			errors++;
		}
		bParamFile >> inint; if (inint < 0 || inint > 1) {
			BatchError(filetype, line, 1, "FixReplicateSeed");
			errors++;
		}

		line++;
		// read next simulation number
		inint = -98765;
		bParamFile >> inint;
		if (bParamFile.eof()) {
			inint = -98765;
		}
		else { // check for valid simulation number
			if (inint != prevsimul + 1) {
				BatchError(filetype, line, 222, " ");
				errors++;
			}
			prevsimul = inint; nSimuls++;
		}
		//	batchlog << "ParseParametersFile(): First item of next line = " << inint << endl;
	} // end of while loop
	if (!bParamFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return nSimuls;
}

int ParseLandFile(int landtype, string indir)
{
	//string fname,header,intext,ftype,costfile;
	string fname, header, intext, ftype;
	int j, inint, line;
	float infloat;
	rasterdata patchraster, spdistraster, costraster;
	int errors = 0;
	//int Kerrors = 0;
	int totlines = 0;
	//bool errorshown = false;
	vector <int> landlist;
	string filetype = "LandFile";

	//batchlog << "ParseLandFile(): starting " << endl;
	if (landtype == 0 || landtype == 2) { // real landscape
		// Parse header line;
		bLandFile >> header; if (header != "LandNum") errors++;
		//	batchlog << "ParseLandFile(): header = " << header << endl;
		bLandFile >> header; if (header != "Nhabitats") errors++;
		bLandFile >> header; if (header != "LandscapeFile") errors++;
		bLandFile >> header; if (header != "PatchFile") errors++;
		bLandFile >> header; if (header != "CostMapFile") errors++;
		bLandFile >> header; if (header != "DynLandFile") errors++;
		bLandFile >> header; if (header != "SpDistFile") errors++;
		if (errors > 0) {
			FormatError(filetype, 0);
			batchlog << "*** Ensure format is correct for real landscape" << endl;
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
				//		batchlog << "ParseLandFile(): Adding landscape no. " << inint
				//			<< " to landscape list" << endl;
				landlist.push_back(inint);
			}
			bLandFile >> inint;
			if (landtype == 0) { // raster map with unique habitat codes
				if (inint < 0) {
					BatchError(filetype, line, 10, "Nhabitats"); errors++;
				}
				if (inint > maxNhab) {
					BatchError(filetype, line, 0, " ");
					batchlog << "Nhabitats may not exceed MaxHabitats in Control file" << endl;
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
					batchlog << ftype << " headers OK: " << fname << endl;
				else {
					errors++;
					batchlog << msgresol0 << ftype << " " << fname
						<< msgresol1 << endl;
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
					batchlog << ftype << msgpatch << endl;
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
								batchlog << ftype << " headers OK: " << fname << endl;
							}
							else {
								batchlog << msghdrs0 << ftype << " " << fname
									<< msghdrs1 << endl;
								errors++;
							}
						}
						else {
							batchlog << msgresol0 << ftype << " " << fname
								<< msgresol1 << endl;
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
						batchlog << ftype << " is required for a habitat quality landscape" << endl;
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
								batchlog << ftype << " headers OK: " << fname << endl;
							}
							else {
								batchlog << msghdrs0 << ftype << " " << fname
									<< msghdrs1 << endl;
								errors++;
							}
						}
						else {
							batchlog << msgresol0 << ftype << " " << fname
								<< msgresol1 << endl;
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
					batchlog << ftype << " must be NULL if transfer model is not SMS" << endl;
				}
			}

			// check dynamic landscape filename
			ftype = "DynLandFile";
			bLandFile >> intext;
			if (intext != "NULL") { // landscape is dynamic
				fname = indir + intext;
				batchlog << "Checking " << ftype << " " << fname << endl;
				bDynLandFile.open(fname.c_str());
				if (bDynLandFile.is_open()) {
					int something = ParseDynamicFile(indir, gNameCostFile);
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
					batchlog << ftype << " is required as SpeciesDist is 1 in Control file" << endl;
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
									batchlog << "*** Extent of " << ftype
										<< " does not match extent of LandscapeFile" << endl;
									errors++;
								}
								else {
									// check origins match
									if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
										&& (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
										batchlog << ftype << " headers OK: " << fname << endl;
									}
									else {
										batchlog << "*** Origin co-ordinates of " << ftype
											<< " do not match those of LandscapeFile" << endl;
										errors++;
									}
								}
							}
							else { // not able to check extents match
								// check origins match
								if ((int)spdistraster.xllcorner == (int)landraster.xllcorner
									&& (int)spdistraster.yllcorner == (int)landraster.yllcorner) {
									batchlog << ftype << " headers OK: " << fname << endl;
								}
								else {
									batchlog << "*** Origin co-ordinates of " << ftype
										<< " do not match those of LandscapeFile" << endl;
									errors++;
								}
							}
						}
						else {
							batchlog << "*** Resolution of " << ftype << " " << fname
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
				batchlog << "*** Ensure format is correct for artificial landscape" << endl;
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
						batchlog << "Y dimension may not be less than X dimension" << endl; errors++;
					}
					if ((Xdim > 2 && power2check(Xdim - 1) != 1)
						|| (Ydim > 2 && power2check(Ydim - 1) != 1)) {
						BatchError(filetype, line, 0, " ");
						batchlog << "X and Y dimensions must be a power of 2 plus 1" << endl; errors++;
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
						batchlog << "MaxHab must exceed MinHab" << endl; errors++;
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
			batchlog << "*** Critical error in land file. "
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

int ParseDynamicFile(string indir, string costfile) {
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
		batchlog << "*** Error in DynLandFile - first change number must be 1" << endl;
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
					batchlog << ftype << " headers OK: " << fname << endl;
				}
				else {
					batchlog << msghdrs0 << ftype << " " << fname
						<< msghdrs1 << endl;
					errors++;
				}
			else {
				errors++;
				batchlog << msgresol0 << ftype << " " << fname << msgresol1 << endl;
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
				batchlog << ftype << msgpatch << endl;
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
							batchlog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchlog << msghdrs0 << ftype << " " << fname
								<< msghdrs1 << endl;
							errors++;
						}
					}
					else {
						batchlog << msgresol0 << ftype << " " << fname
							<< msgresol1 << endl;
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
				batchlog << ftype << " must be supplied " << endl;
			}
		}
		else {
			if (costfile == "NULL") {
				BatchError(filetype, line, 0, " "); errors++;
				batchlog << ftype << " must be NULL to match LandFile " << endl;
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
							batchlog << ftype << " headers OK: " << fname << endl;
						}
						else {
							batchlog << msghdrs0 << ftype << " " << fname
								<< msghdrs1 << endl;
							errors++;
						}
					}
					else {
						batchlog << msgresol0 << ftype << " " << fname
							<< msgresol1 << endl;
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
				batchlog << "Change numbers must be sequential integers" << endl;
				errors++;
			}
			prevchange = change;
		}
	}

	if (errors > 0) return -111;
	else return 0;

}

//---------------------------------------------------------------------------
int ParseStageFile(string indir)
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
				//			batchlog << "*** line = " << line << " i = " << i << " filename = " << filename
				//				<< " transfiles[i] = " << transfiles[i] << endl;
				checkfile = false;
			}
		}
		if (checkfile) {
			if (filename == "NULL") {
				batchlog << "*** " << ftype2 << " is compulsory for stage-structured model" << endl;
				errors++;
			}
			else {
				fname = indir + filename;
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bTransMatrix.open(fname.c_str());
				if (bTransMatrix.is_open()) {
					err = ParseTransitionFile(stages, sexesDem);
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
				batchlog << "FecStageWts must be 0 if FecDensDep is 0" << endl; errors++;
				errors++; fecstagewts = 1;
			}
		}

		// fecundity stage weights file - optional
		ftype2 = "FecStageWtsFile";
		bStageStructFile >> filename;
		if (filename == "NULL") {
			if (fecstagewts) {
				BatchError(filetype, line, 0, " ");
				batchlog << ftype2 << " is compulsory unless FecStageWts is 0" << endl;
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
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bStageWeightsFile.open(fname.c_str());
				if (bStageWeightsFile.is_open()) {
					err = ParseWeightsFile(ftype2);
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
				batchlog << "DevStageWts must be 0 if DevDensDep is 0" << endl; errors++;
				errors++; devstagewts = 1;
			}
		}

		// development stage weights file - optional
		ftype2 = "DevStageWtsFile";
		bStageStructFile >> filename;
		if (filename == "NULL") {
			if (devstagewts) {
				BatchError(filetype, line, 0, " ");
				batchlog << ftype2 << " is compulsory unless DevStageWts is 0" << endl;
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
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bStageWeightsFile.open(fname.c_str());
				if (bStageWeightsFile.is_open()) {
					err = ParseWeightsFile(ftype2);
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
				batchlog << "SurvStageWts must be 0 if SurvDensDep is 0" << endl; errors++;
				errors++; survstagewts = 1;
			}
		}

		// survival stage weights file - optional
		ftype2 = "SurvStageWtsFile";
		bStageStructFile >> filename;
		if (filename == "NULL") {
			if (survstagewts) {
				BatchError(filetype, line, 0, " ");
				batchlog << ftype2 << " is compulsory unless SurvStageWts is 0" << endl;
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
				batchlog << "Checking " << ftype2 << " " << fname << endl;
				bStageWeightsFile.open(fname.c_str());
				if (bStageWeightsFile.is_open()) {
					err = ParseWeightsFile(ftype2);
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
int ParseTransitionFile(short nstages, short nsexesDem)
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
			if (nsexesDem == 1) hhh = Int2Str(i);
			else {
				if (j == 0) hhh = Int2Str(i) + "m"; else hhh = Int2Str(i) + "f";
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
		batchlog << "Invalid row header" << endl; errors++;
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
		batchlog << "MinAge must be zero for juvenile stage" << endl; errors++;
	}

	// one row for each stage/sex combination
	//				batchlog << "HEADER = " << header << endl;
	for (stage = 1; stage < nstages; stage++) {
		for (sex = 0; sex < nsexesDem; sex++) {
			line++;
			// row header
			bTransMatrix >> header;
			if (nsexesDem == 1) hhh = Int2Str(stage);
			else {
				if (sex == 0) hhh = Int2Str(stage) + "m"; else hhh = Int2Str(stage) + "f";
			}
			if (header != hhh) {
				BatchError(filetype, line, 0, " ");
				batchlog << "Invalid row header" << endl; errors++;
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
				batchlog << "MinAge must be zero for stage 1" << endl; errors++;
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
int ParseWeightsFile(string filetype)
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
			if (sexesDem == 1) hhh = Int2Str(i);
			else {
				if (j == 0) hhh = Int2Str(i) + "m"; else hhh = Int2Str(i) + "f";
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
			if (sexesDem == 1) hhh = Int2Str(stage);
			else {
				if (sex == 0) hhh = Int2Str(stage) + "m"; else hhh = Int2Str(stage) + "f";
			}
			if (header != hhh) {
				BatchError(filetype, line, 0, " ");
				batchlog << "Invalid row header" << endl; errors++;
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
		// validate use full kernel
		if (inUseFullKern != 0 && inUseFullKern != 1) {
			BatchError(whichInputFile, lineNb, 1, "UseFullKern"); 
			nbErrors++;
		}
		if (inDensDep == 1 && inUseFullKern != 0) {
				BatchError(whichInputFile, lineNb, 0, "UseFullKern"); 
				nbErrors++;
				batchlog << "UseFullKern must be 0 if there is density-dependent emigration" << endl;
		}
		// validate emigration stage
		if (stagestruct && !inStgDep && inIndVar == 1
			&& inStage == 0 && inSex == 0
			&& (inEmigStg < 0 || inEmigStg >= stages)) {
			BatchError(whichInputFile, lineNb, 0, "EmigStage");
			nbErrors++;
			batchlog << "EmigStage must be from 0 to " << Int2Str(stages - 1) << endl;
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
				batchlog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is enabled EP must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inD0 < 0.0 || inD0 > 1.0) {
				BatchError(whichInputFile, lineNb, 20, "D0"); 
				nbErrors++;
			}
			if (inAlpha < 0.0) {
				BatchError(whichInputFile, lineNb, 10, "alpha");
				nbErrors++;
			}
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
				batchlog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled D0 must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inAlpha != gEmptyVal) {
				batchlog << "*** Error in " << whichInputFile << ": "
					<< "if density-dependence is disabled alpha must be " << gEmptyVal << endl;
				nbErrors++;
			}
			if (inAlpha != gEmptyVal) {
				batchlog << "*** Error in " << whichInputFile << ": "
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
		batchlog << msgnlines << currentLine.simNb
			<< msgshldbe << currentLine.reqdSimLines << endl;
	}
	if (!bEmigrationFile.eof()) {
		EOFerror(whichInputFile);
		nbErrors++;
	}
	if (nbErrors > 0) return -111;
	else return nbSims;
}

//---------------------------------------------------------------------------
int ParseTransferFile(string indir)
{
	string header, colheader, intext, fname, ftype;
	int i, simNb, stagedep, sexdep, kerneltype, distmort, indvar, stage, sex;
	int	prMethod, smtype, inStraightenPath;
	float pr, dp, smconst;
	int goaltype, memsize, betaDB; float gb, alphaDB;
	float meanDistI, meanDistII, ProbKernelI;
	float mortProb, slope, inflPoint;
	float morthab, mortmatrix;
	int costhab, costmatrix;
	float SL, rho;

	vector <string> costsfiles;

	int errors = 0; int morthaberrors = 0; int costerrors = 0; int hrerrors = 0;
	int simuls = 0;
	string filetype = "TransferFile";

	// Parse header line;
	bTransferFile >> header; if (header != "Simulation") errors++;
	switch (gTransferType) {

	case 0: { // negative exponential dispersal kernel
		batchlog << "Checking dispersal kernel format file" << endl;
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
		batchlog << "Checking SMS format file ";
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
			batchlog << "for LandType = 0" << endl;
			for (i = 0; i < maxNhab; i++) {
				colheader = "MortHab" + Int2Str(i + 1);
				bTransferFile >> header; if (header != colheader) morthaberrors++;
			}
			for (i = 0; i < maxNhab; i++) {
				colheader = "CostHab" + Int2Str(i + 1);
				bTransferFile >> header; if (header != colheader) costerrors++;
			}
			break;
		} // end of raster map with unique habitat codes
		case 2: { // raster map with habitat quality
			batchlog << "for LandType = 2" << endl;
			break;
		} // end of raster map with habitat quality
		case 9: { // artificial landscape
			batchlog << "for LandType = 9" << endl;
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
		batchlog << "Checking CRW format file" << endl;
		bTransferFile >> header; if (header != "IndVar") errors++;
		bTransferFile >> header; if (header != "SL") errors++;
		bTransferFile >> header; if (header != "Rho") errors++;
		bTransferFile >> header; if (header != "StraightenPath") errors++;
		bTransferFile >> header; if (header != "SMtype") errors++;
		bTransferFile >> header; if (header != "SMconst") errors++;
		if (landtype == 0) {
			for (i = 0; i < maxNhab; i++) {
				colheader = "MortHab" + Int2Str(i + 1);
				bTransferFile >> header; if (header != colheader) morthaberrors++;
			}
		}
		break;
	} // end of CRW

	} // end of switch (transfer)
	// report any errors in headers, and if so, terminate validation
	if (errors > 0 || morthaberrors > 0 || costerrors > 0 || hrerrors > 0) {
		FormatError(filetype, errors + morthaberrors + costerrors);
		if (morthaberrors > 0) BatchError(filetype, -999, 333, "MortHab");
		if (costerrors > 0) BatchError(filetype, -999, 333, "CostHab");
		if (hrerrors > 0) BatchError(filetype, -999, 444, "Hr");
		return -111;
	}

	// Parse data lines
	int line = 1;
	simCheck current, prev;
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	bTransferFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {

		switch (gTransferType) {

		case 0: { // negative exponential dispersal kernel
			// read and validate columns relating to stage and sex-dependency and to IIV
			bTransferFile >> stagedep >> sexdep >> kerneltype >> distmort;
			bTransferFile >> indvar >> stage >> sex;
			current = CheckStageSex(filetype, line, simNb, prev, stagedep, sexdep, stage, sex, indvar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			// validate kernel type
			if (kerneltype < 0 || kerneltype > 1) {
				BatchError(filetype, line, 1, "KernelType"); errors++;
			}
			// validate mortality
			if (distmort < 0 || distmort > 1) {
				BatchError(filetype, line, 1, "DistMort"); errors++;
			}
			// read remaining columns of the current record
			bTransferFile >> meanDistI >> meanDistII >> ProbKernelI;
			bTransferFile >> mortProb >> slope >> inflPoint;


			if (!indvar) {
				if (meanDistI < resolution) {
					// NOTE - should also check whether emigration prob is constant and equal to 1
					//but checks across diffferent input files are not yet implemented
					BatchError(filetype, line, 2, "meanDistI", "Resolution"); errors++;
				}
				if (kerneltype != 0) {
					if (meanDistII < resolution) {
						// NOTE - DITTO
						BatchError(filetype, line, 2, "meanDistII", "Resolution"); errors++;
					}
					if (ProbKernelI <= 0.0 || ProbKernelI >= 1.0) {
						BatchError(filetype, line, 20, "ProbKernelI"); errors++;
					}
				}
			}

			if (stage == 0 && sex == 0) {
				if (distmort) { // distance-dependent mortality
					// WHAT CONDITIONS APPLY TO MORTALITY SLOPE AND INFLECTION POINT?
				}
				else { // constant mortality
					if (mortProb < 0.0 || mortProb >= 1.0) {
						BatchError(filetype, line, 20, "MortProb"); errors++;
					}
				}
			}

			break;
		} // end of negative exponential dispersal kernel

		case 1: { // SMS
			bTransferFile >> indvar;
			bTransferFile >> pr >> prMethod >> dp;
			bTransferFile >> memsize >> gb >> goaltype >> alphaDB >> betaDB;
			current = CheckStageSex(filetype, line, simNb, prev, 0, 0, 0, 0, 0, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			// validate SMS movement parameters
			if (pr < 1) {
				BatchError(filetype, line, 11, "PR"); errors++;
			}
			if (prMethod < 1 || prMethod > 3) {
				BatchError(filetype, line, 33, "PRmethod"); errors++;
			}
			if (!indvar && dp < 1.0) {
				BatchError(filetype, line, 11, "DP"); errors++;
			}
			if (memsize < 1 || memsize > 14) {
				BatchError(filetype, line, 0, "MemSize"); errors++;
				batchlog << "MemSize must be from 1 to 14" << endl;
			}
			if (!indvar && gb < 1.0) {
				BatchError(filetype, line, 11, "GB"); errors++;
			}
			if (goaltype < 0 || goaltype > 2) {
				BatchError(filetype, line, 2, "GoalType"); errors++;
			}
			if (!indvar && goaltype == 2) { // dispersal bias
				if (alphaDB <= 0.0) {
					BatchError(filetype, line, 10, "AlphaDB"); errors++;
				}
				if (betaDB <= 0.0) {
					BatchError(filetype, line, 10, "BetaDB"); errors++;
				}
			}
			bTransferFile >> inStraightenPath >> smtype >> smconst;
			if (inStraightenPath < 0 || inStraightenPath > 1) {
				BatchError(filetype, line, 1, "StraightenPath"); errors++;
			}
			if (landtype == 2) // habitat quality landscape 
			{ // must have constant mortality
				if (smtype != 0) {
					BatchError(filetype, line, 0, " "); errors++;
					batchlog << "SMtype must be 0 for LandType 2" << endl;
				}
			}
			else {
				if (smtype < 0 || smtype > 1) {
					BatchError(filetype, line, 1, "SMtype"); errors++;
				}
			}
			if (smtype == 0)
			{
				if (smconst < 0.0 || smconst >= 1.0) {
					BatchError(filetype, line, 20, "SMconst"); errors++;
				}
			}
			switch (landtype) {

			case 0: { // raster map with unique habitat codes
				//			batchlog << "for LandType = 0" << endl;
				for (i = 0; i < maxNhab; i++) {
					bTransferFile >> morthab;
					if (smtype == 1)
					{
						if (morthab < 0.0 || morthab >= 1.0) {
							colheader = "MortHab" + Int2Str(i + 1);
							BatchError(filetype, line, 20, colheader); errors++;
						}
					}
				}
				for (i = 0; i < maxNhab; i++) {
					bTransferFile >> costhab;
					if (gNameCostFile == "NULL") {
						if (costhab < 1) {
							colheader = "CostHab" + Int2Str(i + 1);
							BatchError(filetype, line, 11, colheader); errors++;
						}
					}
				}
				break;
			} // end of raster map with unique habitat codes

			case 2: { // raster map with habitat quality
				//			batchlog << "for LandType = 2" << endl;
				break;
			} // end of raster map with habitat quality

			case 9: { // artificial landscape
				//			batchlog << "for LandType = 9" << endl;
				bTransferFile >> morthab >> mortmatrix;
				bTransferFile >> costhab >> costmatrix;
				if (smtype) { // validate habitat-dependent mortality
					if (morthab < 0.0 || morthab >= 1.0) {
						BatchError(filetype, line, 20, "MortHabitat"); errors++;
					}
					if (mortmatrix < 0.0 || mortmatrix >= 1.0) {
						BatchError(filetype, line, 20, "MortMatrix"); errors++;
					}
				}
				if (costhab < 1) {
					BatchError(filetype, line, 11, "CostHabitat"); errors++;
				}
				if (costmatrix < 1) {
					BatchError(filetype, line, 11, "CostMatrix"); errors++;
				}
				break;
			} // end of artificial landscape

			} // end of switch (landtype)

			break;

		} // end of SMS

		case 2: { // CRW
			bTransferFile >> indvar >> SL >> rho >> inStraightenPath >> smtype >> smconst;
			current = CheckStageSex(filetype, line, simNb, prev, 0, 0, 0, 0, indvar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;

			if (SL <= 0.0) {
				BatchError(filetype, line, 10, "SL"); errors++;
			}
			if (rho <= 0.0 || rho >= 1.0) {
				BatchError(filetype, line, 20, "Rho"); errors++;
			}
			if (inStraightenPath < 0 || inStraightenPath > 1) {
				BatchError(filetype, line, 1, "StraightenPath"); errors++;
			}
			if (landtype == 0) { // real landscape with habitat types
				if (smtype < 0 || smtype > 1) {
					BatchError(filetype, line, 1, "SMtype"); errors++;
				}
				if (!smtype) {
					if (smconst < 0.0 || smconst >= 1.0) {
						BatchError(filetype, line, 20, "SMconst"); errors++;
					}
				}
				for (int i = 0; i < maxNhab; i++) {
					bTransferFile >> morthab;
					if (smtype) {
						if (morthab < 0.0 || morthab >= 1.0) {
							colheader = "MortHab" + Int2Str(i + 1);
							BatchError(filetype, line, 20, colheader); errors++;
						}
					}
				}
			}
			else { // real landscape with quality OR artificial landscape
				if (smtype != 0) {
					BatchError(filetype, line, 0, " "); errors++;
					batchlog << "SMtype must be 0 for LandType 2 or 9" << endl;
				}
				if (smconst < 0.0 || smtype >= 1.0) {
					BatchError(filetype, line, 20, "SMconst"); errors++;
				}
			}
			break;
		} // end of CRW

		} // end of switch (transfer)

		// read next simulation
		line++;
		simNb = -98765;
		bTransferFile >> simNb;
		if (bTransferFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (gTransferType == 0 // no. of lines checked for dispersal kernel transfer method only
		&& current.simLines != current.reqdSimLines) {
		BatchError(filetype, line, 0, " "); errors++;
		batchlog << msgnlines << current.simNb
			<< msgshldbe << current.reqdSimLines << endl;
	}
	if (!bTransferFile.eof()) {
		EOFerror(filetype);
		errors++;
	}
	costsfiles.clear();

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int ParseSettleFile(void)
{
	string header;
	int simNb, stagedep, sexdep, stage, sex, settletype;
	int densdep, indvar, findmate, minSteps, maxSteps, maxStepsYear;
	float s0, alphaS, betaS;
	int errors = 0;
	int simuls = 0;
	string filetype = "SettlementFile";

	// Parse header line;
	bSettlementFile >> header; if (header != "Simulation") errors++;
	bSettlementFile >> header; if (header != "StageDep") errors++;
	bSettlementFile >> header; if (header != "SexDep") errors++;
	bSettlementFile >> header; if (header != "Stage") errors++;
	bSettlementFile >> header; if (header != "Sex") errors++;
	if (gTransferType == 0)
	{ // dispersal kernel
		bSettlementFile >> header; if (header != "SettleType") errors++;
		bSettlementFile >> header; if (header != "FindMate") errors++;
	}
	else { // movement method
		bSettlementFile >> header; if (header != "DensDep") errors++;
		bSettlementFile >> header; if (header != "IndVar") errors++;
		bSettlementFile >> header; if (header != "FindMate") errors++;
		bSettlementFile >> header; if (header != "MinSteps") errors++;
		bSettlementFile >> header; if (header != "MaxSteps") errors++;
		bSettlementFile >> header; if (header != "MaxStepsYear") errors++;
		bSettlementFile >> header; if (header != "S0") errors++;
		bSettlementFile >> header; if (header != "AlphaS") errors++;
		bSettlementFile >> header; if (header != "BetaS") errors++;
	}
	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	simCheck current, prev;
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	bSettlementFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {
		if (gTransferType == 0)
		{ // dispersal kernel
			// read and validate columns relating to stage and sex-dependency (NB no IIV here)
			bSettlementFile >> stagedep >> sexdep >> stage >> sex >> settletype >> findmate;
			current = CheckStageSex(filetype, line, simNb, prev, stagedep, sexdep, stage, sex, 0, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			if (settletype < 0 || settletype > 3) {
				BatchError(filetype, line, 3, "SettleType"); errors++;
			}
			if (!stagestruct && (settletype == 1 || settletype == 3)) {
				BatchError(filetype, line, 0, " "); errors++;
				batchlog << "Invalid SettleType for a non-stage-structured population" << endl;
			}
			if (gNbSexesDisp > 1) {
				if (findmate < 0 || findmate > 1) {
					BatchError(filetype, line, 1, "FindMate"); errors++;
				}
			}
		}
		else { // movement method
			// read and validate columns relating to stage and sex-dependency (IIV psossible)
			bSettlementFile >> stagedep >> sexdep >> stage >> sex >> densdep >> indvar >> findmate;
			current = CheckStageSex(filetype, line, simNb, prev, stagedep, sexdep, stage, sex, indvar, true, false);
			if (current.isNewSim) simuls++;
			errors += current.errors;
			prev = current;
			if (densdep < 0 || densdep > 1) {
				BatchError(filetype, line, 1, "DensDep"); errors++;
			}
			if (densdep == 0) {
				if (indvar != 0) {
					BatchError(filetype, line, 0, " "); errors++;
					batchlog << "IndVar must be 0 if DensDep is 0" << endl;
				}
			}
			if (reproductn != 0 && gNbSexesDisp > 1) {
				if (findmate < 0 || findmate > 1) {
					BatchError(filetype, line, 1, "FindMate"); errors++;
				}
			}
			bSettlementFile >> minSteps >> maxSteps >> maxStepsYear;
			if (stage == 0 && sex == 0) {
				if (minSteps < 0) {
					BatchError(filetype, line, 19, "MinSteps"); errors++;
				}
				if (maxSteps < 0) {
					BatchError(filetype, line, 19, "MaxSteps"); errors++;
				}
			}
			if (maxStepsYear < 0) {
				BatchError(filetype, line, 19, "MaxStepsYear"); errors++;
			}
			bSettlementFile >> s0 >> alphaS >> betaS;

			if (densdep == 1) {

				if (s0 <= 0.0 || s0 > 1.0) {
					BatchError(filetype, line, 20, "S0"); errors++;
				}
				// NOTE: alphaS and betaS can take any value
			}
		}
		// read next simulation
		line++;
		simNb = -98765;
		bSettlementFile >> simNb;
		if (bSettlementFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (current.simLines != current.reqdSimLines) {
		BatchError(filetype, line, 0, " "); errors++;
		batchlog << msgnlines << current.simNb
			<< msgshldbe << current.reqdSimLines << endl;
	}
	if (!bSettlementFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int ParseTraitsFile(string indir)
{
	string header, colheader;
	int simNb;
	string filename, TraitType, Sex, initialDistribution, initialParameters,
		dominanceDistribution, dominanceParameters, isInherited, mutationDistribution, mutationParameters, positions, NbrOfPositions, expressionType, mutationRate;
	int errors = 0;
	int simuls = 0;
	vector <string> archfiles;
	const string filetype = "TraitsFile";

	// Parse header line;
	bTraitsFile >> header; if (header != "Simulation") errors++;
	bTraitsFile >> header; if (header != "TraitType") errors++;
	bTraitsFile >> header; if (header != "Sex") errors++;
	bTraitsFile >> header; if (header != "Positions") errors++;
	bTraitsFile >> header; if (header != "NbrOfPositions") errors++;
	bTraitsFile >> header; if (header != "ExpressionType") errors++;
	bTraitsFile >> header; if (header != "InitialDistribution") errors++;
	bTraitsFile >> header; if (header != "InitialParameters") errors++;
	bTraitsFile >> header; if (header != "DominanceDistribution") errors++;
	bTraitsFile >> header; if (header != "DominanceParameters") errors++;
	bTraitsFile >> header; if (header != "IsInherited") errors++;
	bTraitsFile >> header; if (header != "MutationDistribution") errors++;
	bTraitsFile >> header; if (header != "MutationParameters") errors++;
	bTraitsFile >> header; if (header != "MutationRate") errors++;

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	simCheck current, prev;
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	bTraitsFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	while (simNb != -98765) {
		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		bTraitsFile >> TraitType >> Sex >> positions >> NbrOfPositions >> expressionType >> initialDistribution >> initialParameters
			>> dominanceDistribution >> dominanceParameters >> isInherited >> mutationDistribution >> mutationParameters
			>> mutationRate;

		current = CheckStageSex(filetype, line, simNb, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) simuls++;
		errors += current.errors;
		prev = current;

		// validate parameters
		if ((isInherited == "true" || isInherited == "True" || isInherited == "TRUE") 
			&& (stof(mutationRate) < 0.0 || stof(mutationRate) > 1.0)) {
			BatchError(filetype, line, 20, "mutationRate"); errors++;
		}

		// read next simulation
		line++;
		simNb = -98765;
		bTraitsFile >> simNb;
		if (bTraitsFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (!(current.simLines >= current.reqdSimLines)) {
		BatchError(filetype, line, 0, " "); errors++;
		batchlog << msgnlines << current.simNb
			<< msgshldbe << current.reqdSimLines << endl;
	}
	if (!bTraitsFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------

int ParseGeneticsFile(string indir) {

	string header, colheader, tfName, ftype2;
	int i, simNb, err, NbrPatchesToSample, nIndividualsToSample;
	string filename, ChromosomeEnds, TraitsFile, PatchList, Stages,
		OutputNeutralStatistics, OutputPerLocusWCFstat, OutputPairwiseFst;
	int GenomeSize, OutputInterval;
	float RecombinationRate;
	bool traitsParsed = false;
	int errors = 0;
	int simuls = 0;
	vector <string> archfiles;
	string filetype = "GeneticsFile";

	// Parse header line;
	bGeneticsFile >> header; if (header != "Simulation") errors++;
	bGeneticsFile >> header; if (header != "GenomeSize") errors++;
	bGeneticsFile >> header; if (header != "ChromosomeEnds") errors++;
	bGeneticsFile >> header; if (header != "RecombinationRate") errors++;
	bGeneticsFile >> header; if (header != "OutputNeutralStatistics") errors++;
	bGeneticsFile >> header; if (header != "OutputPerLocusWCFstat") errors++;
	bGeneticsFile >> header; if (header != "OutputPairwiseFst") errors++;
	bGeneticsFile >> header; if (header != "OutputInterval") errors++;
	bGeneticsFile >> header; if (header != "PatchList") errors++;
	bGeneticsFile >> header; if (header != "NbrPatchesToSample") errors++;
	bGeneticsFile >> header; if (header != "nIndividualsToSample") errors++;
	bGeneticsFile >> header; if (header != "Stages") errors++;
	bGeneticsFile >> header; if (header != "TraitsFile") errors++;

	if (errors > 0) {
		FormatError(filetype, errors);
		return -111;
	}

	// Parse data lines
	int line = 1;
	simCheck current, prev;
	simNb = -98765;
	prev.simNb = -999;
	prev.simLines = prev.reqdSimLines = 0;
	bGeneticsFile >> simNb;
	// first simulation number must match first one in parameterFile
	if (simNb != gFirstSimNb) {
		BatchError(filetype, line, 111, "Simulation"); errors++;
	}
	current.simNb = 0; //dummy line to prevent warning message in VisualStudio 2019
	while (simNb != -98765) {
		// read and validate columns relating to stage and sex-dependency (NB no IIV here)
		bGeneticsFile >> GenomeSize >> ChromosomeEnds >> RecombinationRate >> OutputNeutralStatistics >>
			OutputPerLocusWCFstat >> OutputPairwiseFst >> OutputInterval >> PatchList >> NbrPatchesToSample
			>> nIndividualsToSample >> Stages >> TraitsFile;

		current = CheckStageSex(filetype, line, simNb, prev, 0, 0, 0, 0, 0, true, false);
		if (current.isNewSim) simuls++;
		errors += current.errors;
		prev = current;

		// validate parameters


		if (GenomeSize < 0) {
			BatchError(filetype, line, 10, "GenomeSize"); errors++;
		}

		if (TraitsFile == "NULL") {
			batchlog << "*** " << TraitsFile << " is compulsory for genetic models" << endl;
			errors++;
		}
		else {
			//if (!traitsParsed) { //only parse the first traits file for now, could have multiple different traits files if we want 
			tfName = indir + TraitsFile;
			ftype2 = "Traits file";
			batchlog << "Checking " << ftype2 << " " << tfName << endl;
			bTraitsFile.open(tfName.c_str());
			if (bTraitsFile.is_open()) {
				err = ParseTraitsFile(indir);
				if (err >= 0) FileHeadersOK(ftype2); else errors++;
				bTraitsFile.close();
			}
			else {
				OpenError(ftype2, tfName); errors++;
			}
			if (bTraitsFile.is_open()) bTraitsFile.close();
			bTraitsFile.clear();
		}

		// read next simulation
		line++;
		simNb = -98765;
		bGeneticsFile >> simNb;
		if (bGeneticsFile.eof()) simNb = -98765;
	} // end of while loop
	// check for correct number of lines for previous simulation
	if (current.simLines != current.reqdSimLines) {
		BatchError(filetype, line, 0, " "); errors++;
		batchlog << msgnlines << current.simNb
			<< msgshldbe << current.reqdSimLines << endl;
	}
	if (!bGeneticsFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int ParseInitFile(string indir)
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
			colheader = "PropStage" + Int2Str(i);
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
			batchlog << "SeedType must be 0 for an artificial landscape"
				<< endl;
		}
		if (!speciesdist && seedtype == 1) {
			BatchError(filetype, line, 0, " "); errors++;
			batchlog << "SeedType may not be 1 if there is no initial species distribution map"
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
				batchlog << "NCells may not be greater than the area specified (i.e. "
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
				batchlog << ftype2 << " is compulsory for SeedType 2" << endl;
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
					batchlog << "Checking " << ftype2 << " " << fname << endl;
					bInitIndsFile.open(fname.c_str());
					if (bInitIndsFile.is_open()) {
						err = ParseInitIndsFile();
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
				batchlog << ftype2 << " must be NULL for SeedType "
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
					colheader = "PropStage" + Int2Str(i);
					BatchError(filetype, line, 20, colheader); errors++;
				}
			}
			if (seedtype != 2 && (cumprop < 0.99999 || cumprop > 1.00001)) {
				BatchError(filetype, line, 0, " "); errors++;
				batchlog << "Initial proportions must sum to 1.0" << endl;
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
		batchlog << msgnlines << current.simNb
			<< msgshldbe << current.reqdSimLines << endl;
	}
	if (!bInitFile.eof()) {
		EOFerror(filetype);
		errors++;
	}

	if (errors > 0) return -111;
	else return simuls;

}

//---------------------------------------------------------------------------
int ParseInitIndsFile(void) {
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
			batchlog << "Species must be 0" << endl;
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
simCheck CheckStageSex(string filetype, int line, int simNb, simCheck prev,
	int stagedep, int sexdep, int stage, int sex, int indvar,
	bool checklines, bool stgdepindvarok)
{
	simCheck current;
	current.errors = 0;
	int iii;

	// has there been a change of simulation number?;
	if (simNb == prev.simNb) { // no
		current.isNewSim = false; current.simLines = prev.simLines + 1;
	}
	else { // yes
		// check for valid simulation number
		current.isNewSim = true; current.simLines = 1;
		if (line > 1 && simNb != prev.simNb + 1) {
			BatchError(filetype, line, 222, " "); current.errors++;
		}
		// check for correct number of lines for previous simulation
		if (checklines && !(prev.simLines >= prev.reqdSimLines)) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "No. of lines for previous Simulation " << prev.simNb
				<< msgshldbe << prev.reqdSimLines << endl;
		}
	}
	current.simNb = simNb;

	// validate stagedep
	if (stagestruct) {
		if (stagedep != 0 && stagedep != 1) {
			BatchError(filetype, line, 1, "StageDep"); current.errors++;
			stagedep = 1; // to calculate required number of lines
		}
	}
	else {
		if (stagedep != 0) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "StageDep must be 0 for non-stage-structured model" << endl;
			stagedep = 0; // to calculate required number of lines
		}
	}
	// validate sexdep
	if (gNbSexesDisp == 2) {
		if (sexdep != 0 && sexdep != 1) {
			BatchError(filetype, line, 1, "SexDep"); current.errors++;
			sexdep = 1; // to calculate required number of lines
		}
	}
	else {
		if (sexdep != 0) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "SexDep must be 0 for asexual model" << endl;
			sexdep = 0; // to calculate required number of lines
		}
	}
	if (current.isNewSim) { // set required number of lines
		if (stagedep) {
			if (sexdep) current.reqdSimLines = stages * gNbSexesDisp;
			else current.reqdSimLines = stages;
		}
		else {
			if (sexdep) current.reqdSimLines = gNbSexesDisp;
			else current.reqdSimLines = 1;
		}
	}
	else current.reqdSimLines = prev.reqdSimLines;

	// validate stage
	if (stagedep) { // there must be 1 or 2 lines for each stage
		if (sexdep) { // there must be 2 lines for each stage
			if (current.simLines % 2) iii = (current.simLines + 1) / 2; else  iii = current.simLines / 2;
			if (stage != iii - 1) {
				BatchError(filetype, line, 0, " "); current.errors++;
				batchlog << "Stages must be sequentially numbered from 0" << endl;
			}
		}
		else { // there must be 1 line for each stage
			if (stage != current.simLines - 1) {
				BatchError(filetype, line, 0, " "); current.errors++;
				batchlog << "Stages must be sequentially numbered from 0" << endl;
			}
		}
	}
	else { // no stage-dependent emigration
		if (stage != 0) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "Stage must be 0 for non-stage-structured model" << endl;
		}
	}
	// validate sex
	if (sexdep) {
		if (sex != (current.simLines + 1) % 2) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "Sex must be alternately 0 and 1 if SexDep is 1" << endl;
		}
	}
	else {
		if (sex != 0) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "Sex must be 0 if SexDep is 0" << endl;
		}
	}

	// validate indvar
	if (stagedep && !stgdepindvarok) {
		if (indvar != 0) {
			BatchError(filetype, line, 0, " "); current.errors++;
			batchlog << "IndVar must be 0 if stage-dependent" << endl;
		}
	}
	else {
		if (indvar < 0 || indvar > 1) {
			BatchError(filetype, line, 1, "IndVar"); current.errors++;
		}
	}
	return current;
}

// Functions to handle and report error conditions

void BatchError(string filename, int line, int option, string fieldname)
{
	if (line == -999) { // message does not cite line number
		batchlog << "*** Error in " << filename << ": ";
	}
	else {
		batchlog << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		batchlog << fieldname << " must be 0 or 1";
		break;
	case 2:
		batchlog << fieldname << " must be 0, 1 or 2";
		break;
	case 3:
		batchlog << fieldname << " must be 0, 1, 2 or 3";
		break;
	case 4:
		batchlog << fieldname << " must be from 0 to 4";
		break;
	case 5:
		batchlog << fieldname << " must be from 0 to 5";
		break;
	case 6:
		batchlog << fieldname << " must be from 0 to 6";
		break;
	case 7:
		batchlog << fieldname << " must be from 0 to 7";
		break;
	case 10:
		batchlog << fieldname << " must be greater than zero";
		break;
	case 11:
		batchlog << fieldname << " must be 1 or more";
		break;
	case 12:
		batchlog << fieldname << " must be 2 or more";
		break;
	case 13:
		batchlog << fieldname << " must be 3 or more";
		break;
	case 18:
		batchlog << fieldname << " must be greater than 1.0";
		break;
	case 19:
		batchlog << fieldname << " must be 0 or more";
		break;
	case 20:
		batchlog << fieldname << " must be between 0 and 1";
		break;
	case 21:
		batchlog << fieldname << " must be greater than 1";
		break;
	case 33:
		batchlog << fieldname << " must be 1, 2 or 3";
		break;
	case 44:
		batchlog << fieldname << " must be from 1 to 4";
		break;
	case 55:
		batchlog << fieldname << " must be from 1 to 5";
		break;
	case 66:
		batchlog << fieldname << " must be from 1 to 6";
		break;
	case 100:
		batchlog << fieldname << " must be between 0 and 100";
		break;
	case 111:
		batchlog << fieldname << " must match the first Simulation in ParameterFile";
		break;
	case 222:
		batchlog << "Simulation numbers must be sequential integers";
		break;
	case 333:
		batchlog << "No. of " << fieldname << " columns must equal max. no. of habitats ("
			<< maxNhab << ") and be sequentially numbered starting from 1";
		break;
	case 444:
		batchlog << "No. of " << fieldname << " columns must be one fewer than no. of stages, i.e. "
			<< stages - 1 << ", and be sequentially numbered starting from 1";
		break;
	case 555:
		batchlog << "No. of " << fieldname << " columns must equal no. of stages, i.e. "
			<< stages << ", and be sequentially numbered starting from 0";
		break;
	case 666:
		batchlog << fieldname << " must be a unique positive integer";
		break;
	default:
		batchlog << "*** Unspecified error regarding parameter " << fieldname;
	}
	if (option != 0) batchlog << endl;
}

void BatchError(string filename, int line, int option, string fieldname, string fieldname2)
{
	if (line == -999) { // message does not cite line number
		batchlog << "*** Error in " << filename << ": ";
	}
	else {
		batchlog << "*** Error in " << filename << " at line " << line << ": ";
	}
	switch (option) {
	case 0:
		break;
	case 1:
		batchlog << fieldname << " must be greater than " << fieldname2;
		break;
	case 2:
		batchlog << fieldname << " must be greater than or equal to " << fieldname2;
		break;
	case 3:
		batchlog << fieldname << " must be less than or equal to " << fieldname2;
		break;
	case 4:
		batchlog << fieldname << " must be less than " << fieldname2;
		break;
	default:
		batchlog << "*** Unspecified error regarding parameters " << fieldname
			<< " and " << fieldname2;
	}
	if (option != 0) batchlog << endl;
}

void CtrlFormatError(void)
{
	cout << "Format error in Control file" << endl;
	batchlog << endl << "***" << endl << "*** Format error in Control file:"
		<< msgcase << " and file names" << msgmatch
		<< endl
		<< "***" << endl;
}

void ArchFormatError(void)
{
	batchlog << "*** Format error in ArchFile:" << msgcase << msgmatch << endl;
}

void FormatError(string filename, int errors)
{
	batchlog << "*** Format error in header line of ";
	if (errors == 0) {
		batchlog << filename << endl;
	}
	else {
		batchlog << filename << ": " << errors << " error";
		if (errors > 1) batchlog << "s";
		batchlog << " detected" << endl;
	}
}

void OpenError(string ftype, string fname)
{
	batchlog << "*** Unable to open " << ftype << " " << fname << endl;
}

void EOFerror(string filename)
{
	batchlog << "*** Failed to read to EOF in " << filename << endl;
}

void FileOK(string ftype, int n, int option)
{
	batchlog << ftype << " OK: total no. of ";
	switch (option) {
	case 0:
		batchlog << "simulations = ";
		break;
	case 1:
		batchlog << "landscapes = ";
		break;
	case 2:
		batchlog << "parameters = ";
		break;
	default:
		batchlog << "PROBLEMS = ";
	}
	batchlog << n << endl;
}

void FileHeadersOK(string filename)
{
	batchlog << filename << " OK" << endl;
}

void SimulnCountError(string filename)
{
	batchlog << "*** No. of simulations in " << filename
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
		string dummy; // no longer necessary to read no. of habitats from landFile
		landfile >> ppLand.landNum >> dummy >> name_landscape >> name_patch;
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

int readGeneticsFile(int simulationN, Landscape* pLandscape) {

	string indir = paramsSim->getDir(1);
	ifstream inFile(geneticsFile.c_str());
	bool outputWCFstat, outputPerLocusWCFstat, outputPairwiseFst, outputGeneticInterval;

	//not ideal to reset these in here 
	pSpecies->resetGeneticParameters();

	if (inFile.is_open()) {
		//read first header line
		string headerLine, line, value;
		std::getline(inFile, headerLine);

		while (std::getline(inFile, line)) {
			stringstream ss(line);
			vector<string> parameters;
			while (std::getline(ss, value, '	'))
			{
				parameters.push_back(value);
			}

			if (stoi(parameters[0]) == simulationN) {

				int genomeSize = stoi(parameters[1]);

				outputWCFstat = (parameters[4] == "true");
				outputPerLocusWCFstat = (parameters[5] == "true");
				outputPairwiseFst = (parameters[6] == "true");
				outputGeneticInterval = stoi(parameters[7]);
				set<int> patchList;

				string nSampleCellsFst; //number of patches to sample for neutral markers in cell based landscape, not used in patch based landscape

				string patches = parameters[8];
				string n = parameters[9];

				if (pLandscape->getLandParams().patchModel) {// patch-based
					const vector<int> existingPatches = pLandscape->getTruePatchNums();
					patchList = convertStringToPatches(patches, stoi(n), existingPatches);
				}
				else { // cell-based
					if (patches == "all") nSampleCellsFst = "all";
					else if (patches == "random") nSampleCellsFst = n;
					else throw logic_error("Genetics File - ERROR: PatchList must be either 'all' or 'random' for cell-based landscapes.");
				}
				const int nbStages = pSpecies->getStageParams().nStages;
				set<int> stagesToSampleFrom = convertStringToStages(parameters[11], nbStages);

				pSpecies->setGeneticParameters(convertStringToChromosomeEnds(parameters[2], genomeSize), genomeSize, stof(parameters[3]),
					patchList, parameters[10], stagesToSampleFrom, nSampleCellsFst);

				paramsSim->setGeneticSim(outputWCFstat, outputPerLocusWCFstat, outputPairwiseFst, outputGeneticInterval);

				traitsFile = indir + parameters[12];
			}
		}
		inFile.close();
		inFile.clear();
	}
	return 0; //this is for error reporting, need to do error input checks in this function 
}

int readTraitsFile(int simulationN) {

	pSpecies->clearTraitTable();

	ifstream inFile(traitsFile.c_str());

	if (inFile.is_open()) {
		//read first header line
		string headerLine, line, value;
		std::getline(inFile, headerLine);

		while (std::getline(inFile, line)) {
			stringstream ss(line);
			vector<string> parameters;
			while (std::getline(ss, value, '	'))
			{
				parameters.push_back(value);
			}

			if (stoi(parameters[0]) == simulationN)
				setUpTrait(parameters);
			//create trait with parameters 

		}
		inFile.close();
		inFile.clear();
	}
	return 0;
}

void setUpTrait(vector<string> parameters) {
	SpeciesTrait* trait = new SpeciesTrait(parameters, pSpecies);
	TraitType type = trait->stringToTraitType(parameters[1], stringToSex(parameters[2]));
	pSpecies->addTrait(type, *trait);
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
		n = sstruct.nStages * maxNbSexes;

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
	int Nlines, simulationNb, gFirstSimNb = 0, inStage, inSex, inEmigstage;
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
			gFirstSimNb = simulationNb;
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

		if (simulationNb != gFirstSimNb) { // serious problem
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
				pSpecies->setEmigTraits(inStage, inSex, emigrationTraits);
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
				pSpecies->setEmigTraits(0, inSex, emigrationTraits);

			}
		}
		else { // !emig.sexDep
			if (emig.stgDep) {
				if (emig.densDep) {
					emigrationTraits.d0 = inD0; 
					emigrationTraits.alpha = inAlpha; 
					emigrationTraits.beta = inBeta;
					pSpecies->setEmigTraits(inStage, 0, emigrationTraits);
				}
				else {
					emigrationTraits.d0 = inEp; 
					emigrationTraits.alpha = emigrationTraits.beta = 0.0;
					pSpecies->setEmigTraits(inStage, 0, emigrationTraits);
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
				pSpecies->setEmigTraits(0, 0, emigrationTraits);
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
			gFirstSimNb = simNb;
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
		if (simNb != gFirstSimNb) { // serious problem
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
			pSpecies->setKernTraits(0, 0, kernParams, paramsLand.resol);
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
			pSpecies->setKernTraits(0, inSex, kernParams, paramsLand.resol);

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
			pSpecies->setKernTraits(inStage, 0, kernParams, paramsLand.resol);
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
			pSpecies->setKernTraits(inStage, inSex, kernParams, paramsLand.resol);
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
	int gFirstSimNb = 0; // bad, this is a global
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
	pSpecies->setMovtTraits(move);
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
	pSpecies->setMovtTraits(move);
	return error;
}

//---------------------------------------------------------------------------
int ReadSettlement(int option)
{
	int Nlines, simNb, gFirstSimNb = 0, inStageDep, inSexDep, inStage, inSex;
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
			gFirstSimNb = simNb;
			sett.stgDep = (inStageDep == 1);
			sett.sexDep = (inSexDep == 1);
			sett.indVar = (inIndVar == 1) && trfr.usesMovtProc; // no ind var for kernels
			pSpecies->setSettle(sett);

			// update no.of lines according to known stage- and sex-dependency
			Nlines = sett.sexDep ? gNbSexesDisp : 1;
			if (sett.stgDep) Nlines *= sstruct.nStages;
		}

		if (simNb != gFirstSimNb) { // serious problem
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
				pSpecies->setSettTraits(stageToSet, sexToSet, settleDD);
			}

			if (!sett.stgDep) {
				if (!sett.sexDep) {
					if (dem.stageStruct) { // model is structured - also set parameters for all stages
						for (int stg = 1; stg < sstruct.nStages; stg++) {
							pSpecies->setSettRules(stg, 0, srules);
							pSpecies->setSteps(stg, 0, ssteps);
							pSpecies->setSettTraits(stg, 0, settleDD);
							if (hasMales) { // model is sexual - also set parameters for males
								pSpecies->setSettRules(stg, 1, srules);
								pSpecies->setSteps(stg, 1, ssteps);
								if (srules.densDep && !sett.indVar) 
									pSpecies->setSettTraits(stg, 1, settleDD);
							}
						}
					}
					else {
						if (hasMales) { // model is sexual - also set parameters for males
							pSpecies->setSettRules(0, 1, srules);
							pSpecies->setSteps(0, 1, ssteps);
							if (srules.densDep) {
								pSpecies->setSettTraits(0, 1, settleDD);
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
								pSpecies->setSettTraits(stg, sexToSet, settleDD);
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
							pSpecies->setSettTraits(stageToSet, 1, settleDD);
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
			totalProps += propStage;
			paramsInit->setProp(stg, propStage);
		}
		if (totalProps != 1.0)
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

#if RSDEBUG
	DEBUGLOG << endl;
	DEBUGLOG << "RunBatch(): nSimuls=" << nSimuls << " nLandscapes=" << nLandscapes << endl;
	DEBUGLOG << "RunBatch(): landtype=" << landtype << " maxNhab=" << maxNhab << endl;
#endif

	t0 = (int)time(0);

	//int batch_line = 0;

	string name = paramsSim->getDir(2) + "Batch" + Int2Str(sim.batchNum) + "_RS_log.csv";
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
#if RSDEBUG
	rsLog << "WARNING,***** RSDEBUG mode is active *****,,," << endl;
#endif
	rsLog << "RANDOM SEED," << RS_random_seed << ",,," << endl;

	// Open landscape batch file and read header record
	if (ReadLandFile(0)) {
		cout << endl << "Error opening landFile - aborting batch run" << endl;
		return;
	}

	for (int j = 0; j < nLandscapes; j++) {
#if RSDEBUG
		DEBUGLOG << endl;
#endif
		// create new landscape
		if (pLandscape != NULL) delete pLandscape;
		pLandscape = new Landscape;
		bool landOK = true;

		t00 = (int)time(0);
		land_nr = ReadLandFile(1, pLandscape);
		if (land_nr <= 0) { // error condition
			string msg = "Error code " + Int2Str(-land_nr)
				+ " returned from reading LandFile - aborting batch run";
			cout << endl << msg << endl;
			ReadLandFile(9); // close the landscape file
			return;
		}

#if RSDEBUG
		DEBUGLOG << endl << "RunBatch(): j=" << j << " land_nr=" << land_nr
			<< " landtype=" << landtype;
		if (landtype != 9)
			DEBUGLOG << " name_landscape=" << name_landscape
			<< " name_patch=" << name_patch
			<< " name_costfile=" << gNameCostFile
			<< " name_sp_dist=" << name_sp_dist;
		DEBUGLOG << endl;
#endif
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
#if RSDEBUG
			landParams tempLand = pLandscape->getLandParams();
			DEBUGLOG << "RunBatch(): j=" << j
				<< " land_nr=" << land_nr
				<< " landcode=" << landcode
				<< " nHab=" << tempLand.nHab
				<< endl;
#endif

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
#if RSDEBUG
			DEBUGLOG << "RunBatch(): j=" << j
				<< " spDist=" << paramsLand.spDist
				<< endl;
#endif

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

				read_error = readGeneticsFile(i + 1, pLandscape); //simulations numbered >= 1 not 0
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}

				read_error = readTraitsFile(i + 1); //simulations numbered >= 1 not 0
				if (read_error) {
					rsLog << msgsim << sim.simulation << msgerr << read_error << msgabt << endl;
					params_ok = false;
				}

				if (params_ok) {
					simParams sim = paramsSim->getSim();

#if RSDEBUG
					DEBUGLOG << endl << "RunBatch(): i=" << i
						<< " simulation=" << sim.simulation << " landFile=" << landFile
						<< " outRange=" << sim.outRange << " outIntRange=" << sim.outIntRange
						<< endl;
#endif

					cout << endl << "Running simulation nr. " << Int2Str(sim.simulation)
						<< " on landscape no. " << Int2Str(land_nr) << endl;

					// for batch processing, include landscape number in parameter file name
					OutParameters(pLandscape);

					RunModel(pLandscape, i);
#if RSDEBUG
					//DEBUGLOG << endl << "RunBatch(): real landscape, i = " << i
					//	<< " simulation = " << sim.simulation << " landFile = " << landFile
					//	<< endl;
#endif

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
			if (stagestruct) ReadStageStructure(9);
			ReadEmigration(9);
			ReadTransferFile(9, pLandscape);
			ReadSettlement(9);
			ReadInitialisation(9, pLandscape);

			//		if (landtype != 9) 
			if (pLandscape != NULL)
			{
				delete pLandscape; pLandscape = NULL;
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


