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
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
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

#if RSDEBUG
void run_batch_unit_tests() {
	cout << "******* Unit test output for batch interface *******" << endl;
	// call tests here
	cout << endl << "************************" << endl;
}
#endif // RSDEBUG

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

#ifdef _OPENMP
	cout << "OpenMP parallelisation enabled with up to " << omp_get_max_threads() << " threads." << endl;
#endif //_OPENMP

#if RSDEBUG
	assert(0.1 > 0.0); // assert does run correctly
	run_batch_unit_tests();
#else
	// assert does not run in Release mode
	assert(1 == 2);
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
		"Inputs/" + (argc > 2 ? argv[2] : "CONTROL.txt");

#ifndef NDEBUG
	cout << endl << "Working directory: " << paramsSim->getDir(0) << endl;
	cout << endl << "Inputs folder:     " << paramsSim->getDir(1) << endl;
	cout << endl << "Control file:      " << pathToControlFile << endl << endl;
#else 
if (__argc > 1) {
	if (!CheckDirectory(pathToProjectDir)) return 1;
}
DEBUGLOG << endl << "Normal(2.5,0.35):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Normal(2.5,0.35) << " ";
	}
	DEBUGLOG << endl;
}
DEBUGLOG << endl << "Normal(-564.7,123.4):" << endl;
for (int i = 0; i < 5; i++) {
	for (int j = 0; j < 10; j++) {
		DEBUGLOG << pRandom->Normal(-564.7,123.4) << " ";
	}
	DEBUGLOG << endl;
}

*/

/*
DEBUGLOG.close();
DEBUGLOG.clear();

cout << "*****" << endl;
cout << "***** Simulation completed - enter any number to terminate program" << endl;
cout << "*****" << endl;
cin >> i;

return 0;
*/


// set up species
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
			b.transfer
		);
		cout << endl << "Batch input files OK" << endl;
	}
	else {
		cout << endl << "Error in parsing batch input files - see BatchLog file for details" << endl;
	}

	t1 = (int)time(0);
	cout << endl << "***** Elapsed time " << t1 - t0 << " seconds" << endl << endl;
	cout << "*****" << endl;
	cout << "***** Simulation completed." << endl;
	cout << "*****" << endl;
#else
	pRandom = new RSrandom();
#endif

	delete paramsGrad;
	delete paramsStoch;
	delete paramsInit;
	delete paramsSim;
	delete pSpecies;

t1 = (int)time(0);
cout << endl << "***** Elapsed time " << t1-t0 << " seconds" << endl << endl;

cout << "*****" << endl;
cout << "***** Simulation completed - enter any number to terminate program" << endl;
cout << "*****" << endl;
cin >> i;

	return 0;
}

//---------------------------------------------------------------------------

// Dummy functions corresponding to those used in GUI version

/* Batch mode of v2.0 currently has no facility to save maps (unless initiated from GUI).
To do so, we would need a form of bit map which is portable across platforms
and operating systems, rather than the Embarcadero VCL classes.
Does such exist?
*/

traitCanvas SetupTraitCanvas(void) {
traitCanvas tcanv;
for (int i = 0; i < NTRAITS; i++) { tcanv.pcanvas[i] = 0; }
return tcanv;
}

void Landscape::setLandMap(void) { }
void Landscape::drawLandscape(int rep,int yr,int landnum) { }
void Community::viewOccSuit(int year,double mn,double se) { }
void Community::draw(int rep,int yr,int gen,int landNum) { }

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

