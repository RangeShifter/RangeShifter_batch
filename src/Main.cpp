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

 /*------------------------------------------------------------------------------

 RangeShifter v2.0 Main

 Entry level function for BATCH MODE version

 For full details of RangeShifter, please see:
 Bocedi G., Palmer S.C.F., Pe’er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species’ responses to environmental changes.
 Methods in Ecology and Evolution, 5, 388-396. doi: 10.1111/2041-210X.12162

 Author: Steve Palmer, University of Aberdeen

 ------------------------------------------------------------------------------*/

#include <string>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <cassert>

using namespace std;

#include "./RScore/Parameters.h"
#include "./RScore/Landscape.h"
#include "./RScore/Species.h"
#include "./RScore/SubCommunity.h"
#include "./BatchMode.h"

#if LINUX_CLUSTER || R_CMD
#include <unistd.h>
#else
#include <tchar.h>
#include <io.h>
#include <direct.h>
#endif

paramGrad* paramsGrad;		// pointer to environmental gradient parameters
paramStoch* paramsStoch;	// pointer to environmental stochasticity parameters
paramInit* paramsInit;		// pointer to initialisation parameters
paramSim* paramsSim;		// pointer to simulation parameters

Species* pSpecies;		// pointer to species
Community* pComm;		// pointer to community
RSrandom* pRandom;		// pointer to random number routines

//---------------------------------------------------------------------------
#if LINUX_CLUSTER || RS_RCPP
int main(int argc, char* argv[])
#else
int _tmain(int argc, _TCHAR* argv[])
#endif
{

#ifdef NDEBUG
	cout << "RangeShifter Release Mode" << endl;
#else
	cout << "RangeShifter Debug Mode" << endl;
#endif

	int t0, t1;
	t0 = (int)time(0);

	// set up parameter objects
	paramsGrad = new paramGrad;
	paramsStoch = new paramStoch;
	paramsInit = new paramInit;
	paramsSim = new paramSim;

	// set up working directory and control file name
	string pathToControlFile;
#if LINUX_CLUSTER || RS_RCPP
	if (argc > 1) {
		// full path name of directory passed as a parameter
		paramsSim->setDir(argv[1]);
		if (argc > 2) {
			// control file number also passed as a parameter
			int i = atoi(argv[2]);
			pathToControlFile = paramsSim->getDir(0) + "Inputs/CONTROL" + to_string(i) + ".txt";
		}
		else {
			// default name is CONTROL.txt
			pathToControlFile = paramsSim->getDir(0) + "Inputs/CONTROL.txt";
		}
	}
	else {
		// use current directory - get name from first (automatic) parameter
		string nameS = argv[0];
		string path = argv[0];
		unsigned int loc = nameS.find("/", 0);
		while (loc < 999999) {
			nameS = nameS.substr(loc + 1);
			loc = nameS.find("/", 0);
		}
		path = path.substr(0, path.length() - nameS.length());
		paramsSim->setDir(path);
		// control file name is forced to be CONTROL.txt
		pathToControlFile = paramsSim->getDir(0) + "Inputs/CONTROL.txt";
	}
#else
	if (__argc > 1) {
		// full path name of directory passed as a parameter
		paramsSim->setDir(__argv[1]);
		if (__argc > 2) {
			// control file name also passed as a parameter
			pathToControlFile = paramsSim->getDir(0) + "Inputs\\" + __argv[2];
		}
		else {
			// default name is CONTROL.txt
			pathToControlFile = paramsSim->getDir(0) + "Inputs\\CONTROL.txt";
		}
	}
	else {
		// Get the current directory. 
		char* buffer = _getcwd(NULL, 0);
		string dir = buffer;
		free(buffer);
		dir = dir + "\\"; //Current directory path
		paramsSim->setDir(dir);
		// control file name is forced to be CONTROL.txt
		pathToControlFile = paramsSim->getDir(0) + "Inputs\\CONTROL.txt";
	}
#endif
#ifndef NDEBUG
	cout << endl << "Working directory: " << paramsSim->getDir(0) << endl;
	cout << endl << "Inputs folder:     " << paramsSim->getDir(1) << endl;
	cout << endl << "Control file:      " << pathToControlFile << endl << endl;
#endif

	bool errorfolder = CheckDirectory();
	if (errorfolder) {
		cout << endl << "***** Invalid working directory: " << paramsSim->getDir(0)
			<< endl << endl;
		cout << "***** Working directory must contain Inputs, Outputs and Output_Maps folders"
			<< endl << endl;
		cout << "*****" << endl;
		cout << "***** Simulation ABORTED" << endl;
		cout << "*****" << endl;
		return 666;
	}

	// set up species
	pSpecies = new Species;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	transferRules trfr = pSpecies->getTransferRules();

	string indir = paramsSim->getDir(1);
	string outdir = paramsSim->getDir(2);
	batchfiles b = ParseControlAndCheckInputFiles(pathToControlFile, indir, outdir);

	if (b.ok) {
		pSpecies = new Species(
			b.reproductn,
			b.repseasons,
			b.stagestruct == 1, // int to bool
			b.stages,
			b.transfer == 1,
			b.transfer
		);
		cout << endl << "Batch input files OK" << endl;

		simParams sim = paramsSim->getSim();
		sim.batchMode = true;
		sim.batchNum = b.batchNum;
		paramsSim->setSim(sim);
	}
	else {
		cout << endl << "Error in parsing batch input files - see BatchLog file for details" << endl;
	}

	// set up random number class
#if RS_RCPP
#ifndef NDEBUG
	pRandom = new RSrandom(666);
#else
	pRandom = new RSrandom(-1);  // need to be replaced with parameter from control file
#endif
#else
	pRandom = new RSrandom();
#endif

	if (b.ok) {
		try
		{
			RunBatch(b.nSimuls, b.nLandscapes);
		}
		catch (const std::exception& e)
		{
			cerr << endl << "Error: " << e.what() << endl;
		}
	}

	delete pRandom;
	delete paramsGrad;
	delete paramsStoch;
	delete paramsInit;
	delete paramsSim;
	delete pSpecies;

	t1 = (int)time(0);
	cout << endl << "***** Elapsed time " << t1 - t0 << " seconds" << endl << endl;
	cout << "*****" << endl;
	cout << "***** Simulation completed." << endl;
	cout << "*****" << endl;

	return 0;
}

//---------------------------------------------------------------------------

// Dummy functions corresponding to those used in GUI version

/* Batch mode of v2.0 currently has no facility to save maps (unless initiated from GUI).
To do so, we would need a form of bit map which is portable across platforms
and operating systems, rather than the Embarcadero VCL classes.
Does such exist?
*/

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

