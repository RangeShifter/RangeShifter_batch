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
 Bocedi G., Palmer S.C.F., Pe�er G., Heikkinen R.K., Matsinos Y.G., Watts K.
 and Travis J.M.J. (2014). RangeShifter: a platform for modelling spatial
 eco-evolutionary dynamics and species� responses to environmental changes.
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

string getProjectDir(const vector<string>& mainArgs) {
	
	string pathToProjectDir;
	if (mainArgs.size() > 1) {
		// full path name of directory passed as a parameter
		pathToProjectDir = mainArgs[1];
	}
	else {
		// use current directory - get name from first (automatic) parameter
		string nameS = mainArgs[0];
		pathToProjectDir = mainArgs[0];
		unsigned int loc = nameS.find("/", 0);
		while (loc < 999999) {
			nameS = nameS.substr(loc + 1);
			loc = nameS.find("/", 0);
		}
		pathToProjectDir = pathToProjectDir.substr(0, pathToProjectDir.length() - nameS.length());
	}
	return pathToProjectDir;
}

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

	// set up working directory and control file name
	vector<string> args(argc);
	args[0] = argv[0];
	if (argc > 1) args[1] = argv[1];

	string pathToProjectDir = getProjectDir(args);
	paramsSim = new paramSim(pathToProjectDir);

	string pathToControlFile = paramsSim->getDir(0) +
		"Inputs/CONTROL" + (argc > 2 ? argv[2] : "") + ".txt";

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

	// Set up species
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
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

