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

#include "Model.h"

ofstream outPar;
using namespace std::chrono;
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#if RS_RCPP && !R_CMD
Rcpp::List RunModel(Landscape* pLandscape, int seqsim, speciesMap_t allSpecies)
#else
int RunModel(Landscape* pLandscape, int seqsim, speciesMap_t allSpecies)
#endif
{
	int yr, totalInds;
	bool filesOK;

	landParams ppLand = pLandscape->getLandParams();
	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	transferRules trfr = pSpecies->getTransferRules();
	simParams sim = paramsSim->getSim();

	bool anyUsesGradient = false, anySavesVisits = false;
	for (auto& [sp, pSpecies] : allSpecies) {
		if (pSpecies->usesGradient()) anyUsesGradient = true;
		if (pSpecies->savesVisits()) anySavesVisits = true;
	}
	bool hasMultipleReplicates = sim.reps > 1;

	if (!ppLand.generated) {
		pComm = new Community(pLandscape, allSpecies);
		// Allocate patches, sample patches and set landscape limits
		pLandscape->initialise(allSpecies, ppLand, init);
	}

#if RS_RCPP && !R_CMD
	Rcpp::List list_outPop;
#endif

	// Loop through replicates
	for (int rep = 0; rep < sim.reps; rep++) {

		std::cout << "Running replicate " << rep << " / " << sim.reps - 1 << endl;

#if RS_RCPP && !R_CMD
		Rcpp::Rcout << endl << "starting replicate " << rep << endl;
#endif

		if (anySavesVisits && !ppLand.generated) {
			pLandscape->resetVisits();
		}
		if (sim.fixReplicateSeed) {
			pRandom->fixNewSeed(rep);
		}

		int iPatchChg = 0; // track outside year loop to reset between replicates
		int iCostChg = 0;
		
		// Create and select sampled patches for artifical landscapes
		if (ppLand.generated) { // then need to initialise for every replicate
			if (pComm != nullptr) delete pComm;
			pComm = new Community(pLandscape, allSpecies);
			pLandscape->resetLand();
			// Generate patches, sample patches and set landscape limits
			pLandscape->initialise(allSpecies, ppLand, init);
		}

		if (rep == 0) {
			if (!pComm->openOutputFiles(hasMultipleReplicates, ppLand.landNum)) {
				// abort if any file fails to open
#if RS_RCPP && !R_CMD
				return Rcpp::List::create(Rcpp::Named("Errors") = 666);
#else
				return 666;
#endif
			}
		}
		
		if (env.usesStoch && !env.stochIsLocal) {
			// create time series in case of global environmental stochasticity
			pLandscape->setGlobalStoch(sim.years + 1);
		}

		if (anyUsesGradient) pLandscape->drawGradientDev();

		if (ppLand.usesPatches) {
			for (auto& [sp, pSpecies] : allSpecies) 
				if (pSpecies->doesOutputConnect()) {
					pLandscape->createConnectMatrix(sp);
					pLandscape->outConnectHeaders(sp);
				}
		}

		// Dynamic landscape control
		bool updateland = false;
		int chgNb = 0; // landscape change index
		landChange landChg; 
		if (ppLand.dynamic) {
			landChg = pLandscape->getLandChange(0); // get first change year
		}

		// set up populations in the community
		pLandscape->updateCarryingCapacity(allSpecies, 0, 0);
		pComm->initialise(allSpecies, -1);

#if BATCH && RS_RCPP && !R_CMD
		Rcpp::Rcout << "RunModel(): completed initialisation " << endl;
#endif

		for (auto& [sp, pSpecies] : allSpecies) {

			// open a new individuals file for each replicate
			if (pSpecies->doesOutputInds())
				pComm->outIndsHeaders(sp, rep, ppLand.landNum, ppLand.usesPatches);

			if (pSpecies->doesOutputGeneValues()) {
				bool geneOutFileHasOpened = pComm->openOutGenesFile(sp, pSpecies->isDiploid(), ppLand.landNum, rep);
				if (!geneOutFileHasOpened) throw logic_error("Output gene value file could not be initialised.");
			}

			// open a new genetics file for each replicate for per locus and pairwise stats
			if (pSpecies->doesOutputWeirCockerham())
				pComm->openPerLocusFstFile(pSpecies, pLandscape, ppLand.landNum, rep);
			if (pSpecies->doesOutputWeirHill())
				pComm->openPairwiseFstFile(pSpecies, pLandscape, ppLand.landNum, rep);
		}
		
#if RS_RCPP
		// open a new movement paths file for each replicate
		if (sim.outPaths)
			pLandscape->outPathsHeaders(rep, 0);
#endif

		int coutYrFreq = sim.years < 30 ? 1 :
			sim.years < 300 ? 10 :
			sim.years < 3000 ? 100 : 
			sim.years < 30000 ? 1000 : 10000;

		// Years loop
		for (yr = 0; yr < sim.years; yr++) {

#if RS_RCPP && !R_CMD
			Rcpp::checkUserInterrupt();
#endif
			bool mustUpdateK = false;
			if (yr % coutYrFreq == 0) {
#if RS_RCPP && !R_CMD
				Rcpp::Rcout << "Starting year " << yr << "..." << endl;
#else
				std::cout << "Starting year " << yr << endl;
#endif
			}
			if (init.seedType == 0 && init.freeType < 2) {
				// apply any range restrictions
				if (yr == init.initFrzYr) {
					// release initial frozen range - reset landscape to its full extent
					pLandscape->resetLandLimits();
					mustUpdateK = true;
				}
				else if (init.restrictRange && yr > init.initFrzYr) {
					if (yr < init.finalFrzYr) {
						if ((yr - init.initFrzYr) % init.restrictFreq == 0) {
							// apply dynamic range restriction
							commStats s = pComm->getStats();
							int minY = s.maxY - init.restrictRows;
							if (minY < 0) minY = 0;
							pLandscape->setLandLimits(ppLand.minX, minY, ppLand.maxX, ppLand.maxY);
							mustUpdateK = true;
						}
					}
					else if (yr == init.finalFrzYr) {
						// apply final range restriction
						commStats s = pComm->getStats();
						pLandscape->setLandLimits(ppLand.minX, s.minY, ppLand.maxX, s.maxY);
						mustUpdateK = true;
					}
				}
			}
			// environmental gradient, stochasticity & local extinction
			// or dynamic landscape
			updateland = false;
			
			// Environmental stochasticity
			if (env.usesStoch && env.stochIsLocal) {
				pLandscape->updateLocalStoch();
				mustUpdateK = true;
			}

			// Environmental gradient
			for (auto& [sp, pSpecies] : allSpecies) {
				if (pSpecies->usesGradient()) {
					if (pSpecies->isGradientShifting(yr)) {
						pSpecies->incrementGradOptY();
					}
					mustUpdateK = true;
					pLandscape->updateEnvGradient(sp);
				}
			}
			
			// Dynamic landscape
			if (ppLand.dynamic && yr == landChg.chgyear) {
				chgNb = landChg.chgnum;
				updateland = mustUpdateK = true;

				if (ppLand.usesPatches) {
					iPatchChg = pLandscape->applyPatchChanges(chgNb, iPatchChg);
					// index used after years loop to reset between replicates
				}
				if (landChg.costfile != "none") {
					pLandscape->applyCostChanges(chgNb, iCostChg);
				}
				if (chgNb < pLandscape->numLandChanges()) { // get next change
					landChg = pLandscape->getLandChange(chgNb);
				}
				else {
					landChg.chgyear = 9999999;
				}
			}

			if (mustUpdateK) {
				pLandscape->updateCarryingCapacity(allSpecies, yr, chgNb);
			}

			if (ppLand.usesPatches) pLandscape->resetConnectMatrix();

			if (ppLand.dynamic && updateland) {
				if (trfr.usesMovtProc && trfr.moveType == 1) { // SMS
					if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
				}
				// apply effects of landscape change to species present in changed patches
				pComm->scanUnsuitablePatches();
				pComm->dispersal(chgNb, yr);
			}

			if (init.restrictRange) {
				// remove any population from region removed from restricted range
				if (yr > init.initFrzYr && yr < init.finalFrzYr) {
					if ((yr - init.initFrzYr) % init.restrictFreq == 0) {
						pComm->scanUnsuitablePatches();
					}
				}
			}

			if (init.seedType == 2) {
				// add any new initial individuals for the current year
				pComm->initialise(allSpecies, yr);
			}

			// Generation loop
			for (int gen = 0; gen < dem.repSeasons; gen++) {
				
				// Output and pop. visualisation before reproduction
				pComm->traitAndOccOutput(rep, yr, gen);

				// Non-structured pops: range and population output *before* reproductrion
				if (!dem.stageStruct) pComm->popAndRangeOutput(rep, yr, gen);

#if RS_RCPP && !R_CMD
				if (sim.ReturnPopRaster && sim.outPop && yr >= sim.outStartPop && yr % sim.outIntPop == 0) {
					list_outPop.push_back(pComm->addYearToPopList(rep, yr), "rep" + std::to_string(rep) + "_year" + std::to_string(yr));
				}
#endif
				if (gen == 0 && !ppLand.usesPatches) {
					// Local extinction applied before reproduction 
					// so nb juveniles can be reported
					if (env.usesLocalExt) pComm->applyRandLocExt(env.locExtProb);
					if (anyUsesGradient) pComm->applyLocalExtGrad();
				}

				// Reproduction
				pComm->reproduction(yr);

				if (dem.stageStruct && sstruct.survival == 0) { // at reproduction
					// Draw survival + devlpt for adults only
					pComm->drawSurvivalDevlpt(false, true, true, true);
				}

				// Stage-structured pops: range + pop output *after* reproductrion
				if (dem.stageStruct) pComm->popAndRangeOutput(rep, yr, gen);

				// Dispersal
				pComm->emigration();
				pComm->dispersal(chgNb, yr);

				// Draw survival and development
				bool drawJuvs = true;
				bool drawAdults = !dem.stageStruct 
					|| sstruct.survival != 0; // else already resolved for adults
				bool drawDevlpt = true;
				bool drawSurvival = !dem.stageStruct 
					|| sstruct.survival != 2; // else resolved at end of year
				pComm->drawSurvivalDevlpt(drawJuvs, drawAdults, drawDevlpt, drawSurvival);

				for (auto& [sp, pSpecies] : allSpecies) { // could subset ahead

					// Output Individuals
					if (pSpecies->isIndOutputYear(yr)) {
						pComm->outInds(sp, rep, yr, gen);
					}

					// Output Genetics
					if (pSpecies->isGeneticOutputYear(yr)) {
						if (pSpecies->getSamplingOption() == "random_occupied" 
							|| pSpecies->getSamplingOption() == "all")
							// then must re-sample every year
							pLandscape->samplePatches(pSpecies);
						pComm->sampleIndividuals(sp);
						if (pSpecies->doesOutputGeneValues())
							pComm->outputGeneValues(sp, yr, gen);
						if (pSpecies->doesOutputWeirCockerham() || pSpecies->doesOutputWeirHill())
							pComm->outNeutralGenetics(sp, rep, yr, gen);
					}
				}

				pComm->applySurvivalDevlpt();

			} // end of the generation loop

			if (dem.stageStruct) {
				if (sstruct.survival == 2) {
					// Draw survival for all stages
					pComm->drawSurvivalDevlpt(true, true, false, true);
					pComm->applySurvivalDevlpt();
				}
				// Apply age
				pComm->ageIncrement();
				for (auto& [sp, pSpecies] : allSpecies) {
					if (pSpecies->isIndOutputYear(yr))
						// list any individuals dying having reached maximum age
						pComm->outInds(sp, rep, yr, -1);
				}
				pComm->applySurvivalDevlpt(); // delete any such individuals
			}

			// Stop if community has gone extinct
			if (pComm->totalInds() <= 0) {
				std::cout << "All populations went extinct." << endl;
				yr++; 
				break; 
			}

			// Connectivity Matrix
			if (ppLand.usesPatches) {
				for (auto& [sp, pSpecies] : allSpecies) {
					if (pSpecies->isConnectOutputYear(yr))
						pLandscape->outConnect(sp, rep, yr);
				}
			}
		} // end of the years loop

		// Final summary output
		pComm->traitAndOccOutput(rep, yr, 0);
		pComm->popAndRangeOutput(rep, yr, 0);

		pComm->resetPopns();

		// Reset the gradient optimum
		for (auto& [sp, pSpecies] : allSpecies) {
			if (pSpecies->usesGradient())
				pSpecies->resetOptY();
		}
		pLandscape->resetLandLimits();
		const int lastChange = 666666;
		if (ppLand.usesPatches && ppLand.dynamic && iPatchChg > 0) {
			// apply any patch changes to reset landscape to original configuration
			// (provided that at least one has already occurred)
			pLandscape->applyPatchChanges(lastChange, iPatchChg);
		}
		if (ppLand.dynamic) {
			transferRules trfr = pSpecies->getTransferRules();
			if (trfr.usesMovtProc && trfr.moveType == 1) { // SMS
				if (iCostChg > 0) {
					// apply any cost changes to reset landscape to original configuration
					// (provided that at least one has already occurred)
					pLandscape->applyCostChanges(lastChange, iCostChg);
				}
				if (!trfr.costMap) pLandscape->resetCosts(); // in case habitats have changed
			}
		}

		if (ppLand.usesPatches) {
			for (auto& [sp, pSpecies] : allSpecies) {
				if (pSpecies->doesOutputConnect())
					pLandscape->resetConnectMatrix();
			}
		}
		
		pComm->closeYearlyOutputFiles();
		
		for (auto& [sp, pSpecies] : allSpecies) {
			if (pSpecies->savesVisits())
				pLandscape->outVisits(sp, rep, ppLand.landNum);
		}
		if (anySavesVisits) pLandscape->resetVisits();

#if RS_RCPP
		if (sim.outPaths)
			pLandscape->outPathsHeaders(rep, -999);
#endif

	} // end of the replicates loop

	if (ppLand.usesPatches) {
		for (auto& [sp, pSpecies] : allSpecies) {
			if (pSpecies->doesOutputConnect()) {
				pLandscape->deleteConnectMatrix(sp);
				pLandscape->closeConnectOfs(sp);
			}
		}
	}

	if (hasMultipleReplicates) {
		for (auto& [sp, pSpecies] : allSpecies) {
			if (pSpecies->doesOutputOccup()) {
				pComm->outOccupancy(pSpecies);
				pComm->outOccSuit(pSpecies);
			}
		}
	}

	pComm->closeGlobalOutputFiles(hasMultipleReplicates);
	pComm->closeYearlyOutputFiles(); // might still be open if the simulation was stopped by the user

	delete pComm; 
	pComm = nullptr;

#if RS_RCPP && !R_CMD
	return list_outPop;
#else
	return 0;
#endif

}

#if LINUX_CLUSTER || RS_RCPP
// Check whether a specified directory path exists
bool is_directory(const char* pathname) {
	struct stat info;
	if (stat(pathname, &info) != 0) return false; // path does not exist
	if (S_ISDIR(info.st_mode)) return true;
	return false;
}
#endif

//---------------------------------------------------------------------------
bool CheckDirectory(const string& pathToProjDir)
{
	bool errorfolder = false;

	string subfolder;

	subfolder = pathToProjDir + "Inputs";
	const char* inputs = subfolder.c_str();
	if (!is_directory(inputs)) errorfolder = true;
	subfolder = pathToProjDir + "Outputs";
	const char* outputs = subfolder.c_str();
	if (!is_directory(outputs)) errorfolder = true;
	subfolder = pathToProjDir + "Output_Maps";
	const char* outputmaps = subfolder.c_str();
	if (!is_directory(outputmaps)) errorfolder = true;

	if (errorfolder) {
		cout << endl << "***** Invalid working directory: " << pathToProjDir
			<< endl << endl;
		cout << "***** Working directory must contain Inputs, Outputs and Output_Maps folders"
			<< endl << endl;
		cout << "*****" << endl;
		cout << "***** Simulation ABORTED" << endl;
		cout << "*****" << endl;
		return false;
	}
	else return true;
}

//---------------------------------------------------------------------------
void OutParameters(Landscape* pLandscape, speciesMap_t allSpecies) {
	double k;
	int nsexes, nstages;

	landParams ppLand = pLandscape->getLandParams();
	genLandParams ppGenLand = pLandscape->getGenLandParams();
	envStochParams env = paramsStoch->getStoch();
	settleRules srules;
	settleSteps ssteps;
	settleTraits settleDD;
	simParams sim = paramsSim->getSim();
	//envGradParams grad = pSpecies->getEnvGradient();

	string name;
	if (sim.batchMode)
		name = paramsSim->getDir(2)
		+ "Batch" + to_string(sim.batchNum) + "_"
		+ "Sim" + to_string(sim.simulation)
		+ "_Land" + to_string(ppLand.landNum) + "_Parameters.txt";
	else
		name = paramsSim->getDir(2) + "Sim" + to_string(sim.simulation) + "_Parameters.txt";
	outPar.open(name.c_str());

	outPar << "RangeShifter 2.0 ";

	outPar << endl;

	outPar << "================ ";

	outPar << "   =====================";
	outPar << endl << endl;

	outPar << "BATCH MODE \t";
	if (sim.batchMode) outPar << "yes" << endl; 
	else outPar << "no" << endl;
	outPar << "SEED \t" << pRandom->getSeed() << endl;
	outPar << "REPLICATES \t" << sim.reps << endl;
	outPar << "YEARS \t" << sim.years << endl;
	if (ppLand.usesPatches) {
		outPar << "PATCH-BASED MODEL" << endl;
		outPar << "No. PATCHES: \n" << pLandscape->allPatchCount() - 1 << endl;
	}
	else
		outPar << "CELL-BASED MODEL" << endl;
	outPar << "BOUNDARIES \t";
	if (sim.absorbing) outPar << "absorbing" << endl;
	else outPar << "reflective" << endl;
	outPar << endl;

	outPar << "LANDSCAPE:\t";
	if (ppLand.generated) {
		outPar << "artificially generated map" << endl;
		outPar << "TYPE: \t";
		if (ppGenLand.continuous) outPar << "continuous \t";
		else outPar << "discrete \t";
		if (ppGenLand.fractal) outPar << "fractal";
		else outPar << "random";
		outPar << endl << "PROPORTION OF SUITABLE HABITAT (p)\t" << ppGenLand.propSuit << endl;
		if (ppGenLand.fractal) outPar << "HURST EXPONENT\t" << ppGenLand.hurst << endl;
	}
	else {
		outPar << "imported map" << endl;
		outPar << "TYPE: \t";
		switch (ppLand.rasterType) {
		case 0:
			outPar << "habitat codes" << endl;
			break;
		case 1:
			outPar << "habitat % cover" << endl;
			break;
		case 2:
			outPar << "habitat quality" << endl;
			break;
		}
		outPar << "FILE NAME: ";
#if RS_RCPP
		if (ppLand.dynamic) {
			outPar << gHabMapName << endl;
		}
		else {
			outPar << gHabMapName << endl;
		}
		if (ppLand.usesPatches) {
			outPar << "PATCH FILE: " << gPatchMapName << endl;
		}
		if (trfr.costMap) {
			outPar << "COSTS FILE: " << name_costfile << endl;
		}
#else
		if (sim.batchMode) outPar << " (see batch file) " << landFile << endl;
#endif
		outPar << "No. HABITATS:\t" << ppLand.nHab << endl;
	}
	outPar << "RESOLUTION (m): \t" << ppLand.resol << endl;
	outPar << "DIMENSIONS:  X " << ppLand.dimX << "  Y " << ppLand.dimY << endl;
	outPar << "AVAILABLE:   min.X " << ppLand.minX << " min.Y " << ppLand.minY
		<< "  max.X " << ppLand.maxX << " max.Y " << ppLand.maxY << endl;
	if (!ppLand.generated && ppLand.dynamic) {
		landChange chg;
		outPar << "DYNAMIC LANDSCAPE: " << endl;
		int nchanges = pLandscape->numLandChanges();
		for (int i = 0; i < nchanges; i++) {
			chg = pLandscape->getLandChange(i);
			outPar << "Change no. " << chg.chgnum << " in year " << chg.chgyear << endl;
			outPar << "Landscape: " << chg.habfile << endl;
			if (ppLand.usesPatches) {
				outPar << "Patches  : " << chg.pchfile << endl;
			}
			if (chg.costfile != "none") {
				outPar << "Costs    : " << chg.costfile << endl;
			}
		}
	}
	outPar << endl << "SPECIES DISTRIBUTION LOADED: \t";
	if (ppLand.useSpDist)
	{
		outPar << "yes" << endl;
		outPar << "RESOLUTION (m)\t" << ppLand.spResol << endl;
		outPar << "FILE NAME: ";
#if !RS_RCPP
		if (sim.batchMode) outPar << " (see batch file) " << landFile << endl;
#else
		outPar << gSpDistFileName << endl;
#endif
	}
	else outPar << "no" << endl;

	// Initialisation

	initParams init = paramsInit->getInit();
	outPar << endl << "INITIALISATION CONDITIONS:" << endl;
	switch (init.seedType) {
	case 0:
		outPar << "Free initialisation: \t";
		switch (init.freeType) {
		case 0:
			outPar << "Random \t";
			outPar << "No. of cells/patches: " << init.nSeedPatches << endl;
			break;
		case 1:
			outPar << "all suitable cells/patches" << endl;
			break;
		case 2:
			outPar << "manually selected cells/patches" << endl;
			break;
		}
		break;
	case 1:
		outPar << "From species distribution: \t" << endl;
		switch (init.spDistType) {
		case 0:
			outPar << "all presence cells/patches" << endl;
			break;
		case 1:
			outPar << "some random presence cells/patches" << endl;
			break;
		case 2:
			outPar << "all cells/patches within selected distribution cells" << endl;
			break;
		}
		break;
	case 2:
		outPar << "From initial individuals file: " << paramsSim->getDir(1) + init.indsFile << endl;
		break;
	case 3:
		outPar << "From file" << endl;
		break;
	}
	if (init.seedType != 2) {
		outPar << "INITIAL NO. OF INDIVIDUALS: \t";
		switch (init.initDens) {
		case 0:
			outPar << "at carrying capacity" << endl;
			break;
		case 1:
			outPar << "at half carrying capacity" << endl;
			break;
		case 2:
			if (ppLand.usesPatches) {
				outPar << init.indsHa << " individuals per ha" << endl;
			}
			else {
				outPar << init.indsCell << " individuals per cell" << endl;
			}
			break;
		}
		if (dem.stageStruct) {
			outPar << "INITIAL STAGE PROPORTIONS:" << endl;
			for (int i = 1; i < sstruct.nStages; i++) {
				outPar << "stage " << i << ": " << paramsInit->getProp(i) << " \t";
			}
			outPar << endl;
			outPar << "Initial age distribution: ";
			switch (init.initAge) {
			case 0:
				outPar << "lowest possible age";
				break;
			case 1:
				outPar << "randomised";
				break;
			case 2:
				outPar << "quasi-equilibrium";
				break;
			}
			outPar << endl;
		}
		outPar << "GEOGRAPHICAL CONSTRAINTS (cell numbers): " << endl;
		outPar << "min X: " << init.minSeedX << " max X: " << init.maxSeedX << endl;
		outPar << "min Y: " << init.minSeedY << " max Y: " << init.maxSeedY << endl;

		if (init.seedType == 0 && init.freeType < 2) {
			if (init.initFrzYr > 0) {
				outPar << "Freeze initial range until year " << init.initFrzYr << endl;
			}
			if (init.restrictRange) {
				outPar << "Restrict range to northern " << init.restrictRows
					<< " rows every " << init.restrictFreq << " years" << endl;
				if (init.finalFrzYr < sim.years) {
					outPar << "Freeze range at year " << init.finalFrzYr << endl;
				}
			}
		}
	}

	outPar << endl << "SPECIES PARAMETERS" << endl;

	for (auto& [sp, pSpecies] : allSpecies) {

		outPar << endl << "SPECIES " << to_string(sp) << endl;

		demogrParams dem = pSpecies->getDemogrParams();
		stageParams sstruct = pSpecies->getStageParams();
		emigRules emig = pSpecies->getEmigRules();
		transferRules trfr = pSpecies->getTransferRules();
		settleType sett = pSpecies->getSettle();


		outPar << endl << "ENVIRONMENTAL GRADIENT:\t ";
		if (pSpecies->usesGradient()) {
			envGradParams grad = pSpecies->getEnvGradient();
			switch (grad.gradType) {
			case 1:
				if (dem.stageStruct) outPar << "Density dependence strength (1/b)" << endl;
				else outPar << "Carrying capacity (K)" << endl;
				break;
			case 2:
				if (dem.stageStruct) outPar << "Fecundity" << endl;
				else outPar << "Intrinsic growth rate (r)" << endl;
				break;
			case 3:
				outPar << "Local extinction probability" << endl;
				break;
			default:
				outPar << "ERROR ERROR ERROR" << endl;
				;
			}
			outPar << "G:\t\t " << grad.gradIncr << endl;
			outPar << "optimum Y:\t " << grad.optY << endl;
			outPar << "f:\t\t " << grad.factor << endl;
			if (grad.gradType == 3) outPar << "Local extinction prob. at optimum:\t "
				<< grad.extProbOpt << endl;
			outPar << "GRADIENT SHIFTING:\t ";
			if (grad.doesShift)
			{
				outPar << "yes" << endl;
				outPar << "SHIFTING RATE  (rows/year):\t " << grad.shift_rate << endl;
				outPar << "SHIFTING START (year):\t\t " << grad.shiftBegin << endl;
				outPar << "SHIFTING STOP  (year):\t\t " << grad.shiftStop << endl;
			}
			else   outPar << "no" << endl;
		}
		else outPar << "no";
		outPar << endl;
		outPar << "ENVIRONMENTAL STOCHASTICITY:\t";
		if (env.usesStoch) {
			outPar << "yes" << endl;
			outPar << "TYPE\t in ";
			if (dem.stageStruct) {
				if (env.inK) outPar << "1/b" << endl;
				else outPar << "fecundity" << endl;
			}
			else {
				if (env.inK) outPar << "K" << endl;
				else outPar << "R" << endl;
			}
			outPar << "SPATIAL AUTOCORRELATION\t ";
			if (env.stochIsLocal) outPar << "local" << endl;
			else outPar << "global" << endl;
			outPar << "TEMPORAL AUTOCORRELATION (ac)\t" << env.ac << endl;
			outPar << "AMPLITUDE (std)\t" << env.std << endl;
			if (dem.stageStruct) {
				if (env.inK) {
					outPar << "MIN. 1/b\t" << pSpecies->getMinMax(0)
						* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
					outPar << "MAX. 1/b\t" << pSpecies->getMinMax(1)
						* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
				}
				else {
					outPar << "MIN. fecundity\t" << pSpecies->getMinMax(0) << endl;
					outPar << "MAX. fecundity\t" << pSpecies->getMinMax(1) << endl;
				}
			}
			else {
				if (env.inK) {
					outPar << "MIN. K\t" << pSpecies->getMinMax(0)
						* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
					outPar << "MAX. K\t" << pSpecies->getMinMax(1)
						* (10000.0 / (float)(ppLand.resol * ppLand.resol)) << endl;
				}
				else {
					outPar << "MIN. r\t" << pSpecies->getMinMax(0) << endl;
					outPar << "MAX. r\t" << pSpecies->getMinMax(1) << endl;
				}
			}
		}
		else outPar << "no" << endl;
		outPar << "LOCAL EXTINCTION PROBABILITY:\t";
		if (env.usesLocalExt) outPar << env.locExtProb << endl;
		else outPar << "0.0" << endl;

		outPar << "REPRODUCTION:" << endl;
		outPar << "TYPE: ";
		switch (dem.repType) {
		case 0:
			outPar << "Asexual / Only female model" << endl;
			break;
		case 1:
			outPar << "Sexual model (simple)";
			outPar << endl;
			outPar << "PROP. of MALES\t" << dem.propMales << endl;
			break;
		case 2:
			outPar << "Sexual model (explicit mating system)" << endl;
			outPar << "PROP. of MALES\t" << dem.propMales << endl;
			outPar << "MAX. HAREM SIZE (h)\t" << dem.harem << endl;
			break;
		}
		outPar << "STAGE STRUCTURE:\t";
		if (dem.stageStruct) {
			outPar << "yes" << endl;
			outPar << "PROBABILITY OF REPRODUCING IN SUBSEQUENT SEASONS\t" << sstruct.probRep << endl;
			outPar << "No. OF REP. SEASONS BEFORE SUBSEQUENT REPRODUCTIONS\t" << sstruct.repInterval << endl;
			if (!ppLand.generated && ppLand.dynamic) {
				outPar << "ACTION AFTER POPULATION DESTRUCTION: all individuals ";
				if (sstruct.disperseOnLoss) outPar << "disperse" << endl;
				else outPar << "die" << endl;
			}
			outPar << "No. STAGES\t" << sstruct.nStages << endl;
			outPar << "MAX. AGE\t" << sstruct.maxAge << endl;
			// no sex-specific demographic parameters
			if (dem.repType != 2) {
				outPar << "MIN. AGES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getMinAge(i, 0) << "\tyears" << endl;
				}
				outPar << "FECUNDITIES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getFec(i, 0) << endl;
				}
				outPar << "DEVELOPMENT PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getDev(i, 0) << endl;
				}
				outPar << "SURVIVAL PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "stage\t" << i << ":\t" << pSpecies->getSurv(i, 0) << endl;
				}
			}
			// sex-specific demographic parameters
			else {
				outPar << "MIN. AGES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males " << i << ":\t" << pSpecies->getMinAge(i, 1) << " years;\t";
					outPar << "females " << i << ":\t" << pSpecies->getMinAge(i, 0) << " years" << endl;
				}
				outPar << "FECUNDITIES:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males   " << i << ":\t" << pSpecies->getFec(i, 1) << endl;
					outPar << "females " << i << ":\t" << pSpecies->getFec(i, 0) << endl;
				}
				outPar << "DEVELOPMENT PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males   " << i << ":\t" << pSpecies->getDev(i, 1) << endl;
					outPar << "females " << i << ":\t" << pSpecies->getDev(i, 0) << endl;
				}
				outPar << "SURVIVAL PROB.:" << endl;
				for (int i = 0; i < sstruct.nStages; i++) {
					outPar << "males   " << i << ":\t" << pSpecies->getSurv(i, 1) << endl;
					outPar << "females " << i << ":\t" << pSpecies->getSurv(i, 0) << endl;
				}
			}

			outPar << "SCHEDULING OF SURVIVAL: ";
			switch (sstruct.survival) {
			case 0:
				outPar << "At reproduction" << endl;
				break;
			case 1:
				outPar << "Between reproductive events" << endl;
				break;
			case 2:
				outPar << "Annually" << endl;
				break;
			}

			int mSize; // index for weights matrices
			if (dem.repType == 2) mSize = sstruct.nStages * gMaxNbSexes;
			else mSize = sstruct.nStages;

			outPar << "DENSITY-DEPENDENCE IN FECUNDITY:\t";
			if (sstruct.fecDens) {
				outPar << "yes" << endl;
				if (sstruct.fecStageDens) {
					outPar << "STAGE'S WEIGHTS:" << endl;
					for (int i = 0; i < mSize; i++) {
						if (dem.repType == 2) {
							outPar << "stage " << i / gMaxNbSexes << " ";
							if (i % gMaxNbSexes == 0) outPar << "males  : \t";
							else outPar << "females: \t";
						}
						else outPar << "stage " << i << ": \t";
						for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtFec(j, i) << "\t";
						outPar << endl;
					}
				}
				else outPar << "not stage-dependent" << endl;
			}
			else outPar << "no" << endl;

			densDepParams ddparams = pSpecies->getDensDep();

			outPar << "DENSITY-DEPENDENCE IN DEVELOPMENT:\t";
			if (sstruct.devDens) {
				outPar << "yes - coefficient: " << ddparams.devCoeff << endl;
				if (sstruct.devStageDens) {
					outPar << "STAGE'S WEIGHTS:" << endl;
					for (int i = 0; i < mSize; i++) {
						if (dem.repType == 2) {
							outPar << "stage " << i / gMaxNbSexes << " ";
							if (i % gMaxNbSexes == 0) outPar << "males  : \t";
							else outPar << "females: \t";
						}
						else outPar << "stage " << i << ": \t";
						for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtDev(j, i) << "\t";
						outPar << endl;
					}
				}
				else outPar << "not stage-dependent" << endl;
			}
			else outPar << "no" << endl;

			outPar << "DENSITY-DEPENDENCE IN SURVIVAL:\t\t";
			if (sstruct.survDens) {
				outPar << "yes - coefficient: " << ddparams.survCoeff << endl;
				if (sstruct.survStageDens) {
					outPar << "STAGE'S WEIGHTS:" << endl;
					for (int i = 0; i < mSize; i++) {
						if (dem.repType == 2) {
							outPar << "stage " << i / gMaxNbSexes << " ";
							if (i % gMaxNbSexes == 0) outPar << "males  : \t";
							else outPar << "females: \t";
						}
						else outPar << "stage " << i << ": \t";
						for (int j = 0; j < mSize; j++) outPar << pSpecies->getDDwtSurv(j, i) << "\t";
						outPar << endl;
					}
				}
				else outPar << "not stage-dependent" << endl;
			}
			else outPar << "no" << endl;
		} // end of if (dem.stageStruct)
		else { // not stage-strutured
			outPar << "no" << endl;
			outPar << "Rmax\t" << dem.lambda << endl;
			outPar << "bc\t" << dem.bc << endl;
		}

		if (dem.stageStruct) {
			outPar << endl << "HABITAT SPECIFIC 1/b:" << endl;
		}
		else {
			outPar << endl << "CARRYING CAPACITIES:" << endl;
		}
		int nhab = ppLand.nHab;
		if (ppLand.generated) {
			if (ppGenLand.continuous) nhab = 1;
		}
		for (int i = 0; i < nhab; i++) {
			k = pSpecies->getHabK(i) * (10000.0 / (float)(ppLand.resol * ppLand.resol));
			if (!ppLand.generated && ppLand.rasterType == 0) { // imported & habitat codes
				outPar << "Habitat " << pLandscape->getHabCode(i) << ": \t";
			}
			else {
				outPar << "Habitat " << i << ": ";
			}
			if (dem.stageStruct) outPar << "1/b ";
			else outPar << "K ";
			outPar << k << endl;
		}
		
		outPar << "REPRODUCTIVE SEASONS / YEAR\t" << dem.repSeasons << endl;
		
		emigTraits ep0, ep1;
		string sexdept = "SEX-DEPENDENT:   ";
		string stgdept = "STAGE-DEPENDENT: ";
		string indvar = "INDIVIDUAL VARIABILITY: ";
		string emigstage = "EMIGRATION STAGE: ";

		outPar << endl << "DISPERSAL - EMIGRATION:\t";
		if (emig.densDep) {
			outPar << "density-dependent" << endl;
			if (emig.sexDep) {
				outPar << sexdept << "yes" << endl;
				if (emig.stgDep) {
					outPar << stgdept << "yes" << endl;
					outPar << indvar << "no" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage " << i << ":" << endl;
						ep0 = pSpecies->getSpEmigTraits(i, 0);
						ep1 = pSpecies->getSpEmigTraits(i, 1);
						outPar << "D0:    females " << ep0.d0 << "  males " << ep1.d0 << endl;
						outPar << "alpha: females " << ep0.alpha << "  males " << ep1.alpha << endl;
						outPar << "beta:  females " << ep0.beta << "  males " << ep1.beta << endl;
					}
				}
				else { // !emig.stgDep
					outPar << stgdept << "no" << endl;
					ep0 = pSpecies->getSpEmigTraits(0, 0);
					ep1 = pSpecies->getSpEmigTraits(0, 1);
					outPar << "D0:    females " << ep0.d0 << "  males " << ep1.d0 << endl;
					outPar << "alpha: females " << ep0.alpha << "  males " << ep1.alpha << endl;
					outPar << "beta:  females " << ep0.beta << "  males " << ep1.beta << endl;
				}
			}
			else { // !emig.sexDep
				outPar << sexdept << "no" << endl;
				if (emig.stgDep) {
					outPar << stgdept << "yes" << endl;
					outPar << indvar << "no" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						ep0 = pSpecies->getSpEmigTraits(i, 0);
						outPar << "stage " << i << ": \t" << "D0: " << ep0.d0;
						outPar << " \talpha: " << ep0.alpha << " \tbeta: " << ep0.beta << endl;
					}
				}
				else { // !emig.stgDep
					outPar << stgdept << "no" << endl;
					ep0 = pSpecies->getSpEmigTraits(0, 0);
					outPar << "D0:    " << ep0.d0 << endl;
					outPar << "alpha: " << ep0.alpha << endl;
					outPar << "beta:  " << ep0.beta << endl;
				}
			}
		}
		else { // not density-dependent
			string initprob = "INITIAL EMIGRATION PROB. ";
			outPar << "density-independent" << endl;
			if (!trfr.usesMovtProc) { // transfer by kernel
				outPar << "USE FULL KERNEL TO DETERMINE EMIGRATION: ";
				if (pSpecies->useFullKernel()) outPar << "yes";
				else outPar << "no";
				outPar << endl;
			}

			if (emig.sexDep) {
				outPar << sexdept << "yes" << endl;
				if (emig.stgDep) {
					outPar << stgdept << "yes" << endl;
					outPar << indvar << "no" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage " << i << ": \t" << "EMIGRATION PROB.: \tfemales "
							<< pSpecies->getSpEmigD0(i, 0) << " \tmales " << pSpecies->getSpEmigD0(i, 1) << endl;
					}
				}
				else { // !emig.stgDep
					outPar << stgdept << "no" << endl;
					outPar << "EMIGRATION PROB.: \tfemales " << pSpecies->getSpEmigD0(0, 0)
						<< "\t males " << pSpecies->getSpEmigD0(0, 1) << endl;
				}
			}
			else { // !emig.sexDep
				outPar << sexdept << "no" << endl;
				if (emig.stgDep) {
					outPar << stgdept << "yes" << endl;
					outPar << indvar << "no" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage " << i << ": \t" << "EMIGRATION PROB.: "
							<< pSpecies->getSpEmigD0(i, 0) << endl;
					}
				}
				else { // !emig.stgDep
					outPar << stgdept << "no" << endl;
					outPar << "EMIGRATION PROB.:\t" << pSpecies->getSpEmigD0(0, 0) << endl;
				}
			}
		}

		// Transfer

		outPar << endl << "DISPERSAL - TRANSFER: \t";

		if (trfr.usesMovtProc) {
			bool straightenPath;
			if (trfr.moveType == 1) { // SMS
				trfrSMSTraits move = pSpecies->getSpSMSTraits();
				straightenPath = move.straightenPath;
				if (trfr.costMap) {
					outPar << "SMS\tcosts from imported cost map" << endl;
				}
				else {
					outPar << "SMS\tcosts:" << endl;
					if (!ppLand.generated && ppLand.rasterType == 0) {
						for (int i = 0; i < ppLand.nHab; i++)
							outPar << "\thab. " << pLandscape->getHabCode(i) << "\t"
							<< pSpecies->getHabCost(i) << endl;
					}
					else {
						for (int i = 0; i < ppLand.nHab; i++)
							outPar << "\thab. " << i << "\t"
							<< pSpecies->getHabCost(i) << endl;
					}
				}
				string pr = "PERCEPTUAL RANGE";
				outPar << pr << ":        " << move.pr << endl;
				outPar << pr << " METHOD: " << move.prMethod << endl;
				if (!trfr.indVar) outPar << "DIRECTIONAL PERSISTENCE: " << move.dp << endl;
				outPar << "MEMORY SIZE: " << move.memSize << endl;
				outPar << "GOAL TYPE:   " << move.goalType << endl;
				if (!trfr.indVar) {
					if (move.goalType == 2) { //  dispersal bias
						outPar << "GOAL BIAS:   " << move.gb << endl;
						outPar << "ALPHA DB:    " << move.alphaDB << endl;
						outPar << "BETA DB:     " << move.betaDB << endl;
					}
				}
				outPar << indvar << "no " << endl;
			}
			else { // CRW
				trfrCRWTraits move = pSpecies->getSpCRWTraits();
				straightenPath = move.straightenPath;
				outPar << "CRW" << endl;
				string lgth = "STEP LENGTH (m) ";
				string corr = "STEP CORRELATION";
				outPar << lgth << ": " << move.stepLength << endl;
				outPar << corr << ": " << move.rho << endl;
			}
			outPar << "STRAIGHTEN PATH AFTER DECISION NOT TO SETTLE: ";
			if (straightenPath) outPar << "yes" << endl;
			else outPar << "no" << endl;
			outPar << "STEP MORTALITY:\t" << endl;
			if (trfr.habMort)
			{
				outPar << "habitat dependent:\t" << endl;
				if (!ppLand.generated && ppLand.rasterType == 0) {
					for (int i = 0; i < ppLand.nHab; i++)
						outPar << "\thab. " << pLandscape->getHabCode(i) << "\t"
						<< pSpecies->getHabMort(i) << endl;
				}
				else {
					for (int i = 0; i < ppLand.nHab; i++)
						outPar << "\thab. " << i << "\t"
						<< pSpecies->getHabMort(i) << endl;
				}
			}
			else
			{
				trfrCRWTraits move = pSpecies->getSpCRWTraits();
				outPar << "constant " << move.stepMort << endl;
			}
		} // end of movement process
		else { // kernel
			string meandist = "MEAN DISTANCE";
			string probkern = "PROB. KERNEL I";
			trfrKernelParams kern0, kern1;
			outPar << "dispersal kernel" << endl << "TYPE: \t";
			if (trfr.twinKern) outPar << "double ";
			outPar << "negative exponential" << endl;

			if (trfr.sexDep) {
				outPar << sexdept << "yes" << endl;
				if (trfr.stgDep) {
					outPar << stgdept << "yes" << endl;
					outPar << indvar << "no" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						outPar << "stage " << i << ":" << endl;
						kern0 = pSpecies->getSpKernTraits(i, 0);
						kern1 = pSpecies->getSpKernTraits(i, 1);
						outPar << meandist << " I: \tfemales " << kern0.meanDist1 << " \tmales " << kern1.meanDist1 << endl;
						if (trfr.twinKern)
						{
							outPar << meandist << " II: \tfemales " << kern0.meanDist2 << " \tmales " << kern1.meanDist2 << endl;
							outPar << probkern << ": \tfemales " << kern0.probKern1 << " \tmales " << kern1.probKern1 << endl;
						}
					}
				}
				else { // !trfr.stgDep
					outPar << stgdept << "no" << endl;
					kern0 = pSpecies->getSpKernTraits(0, 0);
					kern1 = pSpecies->getSpKernTraits(0, 1);
					outPar << meandist << " I: \tfemales " << kern0.meanDist1 << " \tmales " << kern1.meanDist1 << endl;
					if (trfr.twinKern)
					{
						outPar << meandist << " II: \tfemales " << kern0.meanDist2 << " \tmales " << kern1.meanDist2 << endl;
						outPar << probkern << ": \tfemales " << kern0.probKern1 << " \tmales " << kern1.probKern1 << endl;
					}
				}
			}
			else { // !trfr.sexDep
				outPar << sexdept << "no" << endl;
				if (trfr.stgDep) {
					outPar << stgdept << "yes" << endl;
					outPar << indvar << "no" << endl;
					for (int i = 0; i < sstruct.nStages; i++) {
						kern0 = pSpecies->getSpKernTraits(i, 0);
						outPar << "stage " << i << ": \t" << meandist << " I: " << kern0.meanDist1;
						if (trfr.twinKern)
						{
							outPar << " \t" << meandist << " II: " << kern0.meanDist2;
							outPar << " \t" << probkern << ": " << kern0.probKern1;
						}
						outPar << endl;
					}
				}
				else { // !trfr.stgDep
					outPar << stgdept << "no" << endl;
					kern0 = pSpecies->getSpKernTraits(0, 0);
					outPar << meandist << " I: \t" << kern0.meanDist1 << endl;
					if (trfr.twinKern)
					{
						outPar << meandist << " II: \t" << kern0.meanDist2 << endl;
						outPar << probkern << ": \t" << kern0.probKern1 << endl;
					}
				}
			}

			outPar << "DISPERSAL MORTALITY:   ";
			trfrMortParams mort = pSpecies->getMortParams();
			if (trfr.distMort) {
				outPar << "distance-dependent" << endl;
				outPar << "SLOPE: " << mort.mortAlpha << " \tINFLECTION POINT: " << mort.mortBeta << endl;
			}
			else {
				outPar << "constant" << endl << "MORTALITY PROBABILITY: " << mort.fixedMort << endl;
			}
		} // end of kernel transfer

		// Settlement

		outPar << endl << "DISPERSAL - SETTLEMENT:" << endl;

		if (trfr.usesMovtProc) {
			string plusmating = "+ mating requirements";

			if (sett.sexDep) {
				nsexes = 2;
				outPar << sexdept << "yes" << endl;
				if (sett.stgDep) {
					nstages = sstruct.nStages;
					outPar << stgdept << "yes" << endl;
					for (int i = 0; i < nstages; i++) {
						if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
						for (int sx = 0; sx < nsexes; sx++) {
							if (sx == 0) outPar << "FEMALES:" << endl;
							else outPar << "MALES:" << endl;
							ssteps = pSpecies->getSteps(i, sx);

							outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
							outPar << "MAX. No. OF STEPS:\t ";
							if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
							else outPar << ssteps.maxSteps << endl;
						}
					}
				}
				else { // !sett.stgDep
					nstages = 1;
					outPar << stgdept << "no" << endl;
					for (int sx = 0; sx < nsexes; sx++) {
						if (sx == 0) outPar << "FEMALES:" << endl;
						else outPar << "MALES:" << endl;
						ssteps = pSpecies->getSteps(0, sx);

						outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
						outPar << "MAX. No. OF STEPS:\t ";
						if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
						else outPar << ssteps.maxSteps << endl;
					}
				}
			}
			else { // !sett.sexDep
				nsexes = 1;
				outPar << sexdept << "no" << endl;
				if (sett.stgDep) {
					nstages = sstruct.nStages;
					outPar << stgdept << "yes" << endl;
					for (int i = 0; i < nstages; i++) {
						if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
						ssteps = pSpecies->getSteps(i, 0);

						outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
						outPar << "MAX. No. OF STEPS:\t ";
						if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
						else outPar << ssteps.maxSteps << endl;
					}
				}
				else { // !sett.stgDep
					nstages = 1;
					outPar << stgdept << "no" << endl;
					ssteps = pSpecies->getSteps(0, 0);

					outPar << "MIN. No. OF STEPS:\t " << ssteps.minSteps << endl;
					outPar << "MAX. No. OF STEPS:\t ";
					if (ssteps.maxSteps == 99999999) outPar << "not applied" << endl;
					else outPar << ssteps.maxSteps << endl;
				}
			}
			for (int sx = 0; sx < nsexes; sx++) {
				if (sett.sexDep) {
					if (sx == 0) outPar << "FEMALES:" << endl;
					else outPar << "MALES:" << endl;
				}
				outPar << "SETTLE IF: ";
				for (int i = 0; i < nstages; i++) {
					if (dem.stageStruct && nstages > 1) outPar << "stage " << i << ": " << endl;
					outPar << "find a suitable cell/patch ";
					srules = pSpecies->getSettRules(i, sx);
					if (srules.densDep) {
						settleDD = pSpecies->getSpSettTraits(i, sx);
						outPar << "+ density dependence ";
						if (srules.findMate) outPar << plusmating;
						outPar << endl;
						if (!sett.indVar) {
							outPar << "S0: " << settleDD.s0 << "  AlphaS: " << settleDD.alpha
								<< "  BetaS: " << settleDD.beta << endl;
						}
					}
					else {
						if (srules.findMate) outPar << plusmating << endl;
						else outPar << "(not the natal one)" << endl;
					}
					if (dem.stageStruct) {
						ssteps = pSpecies->getSteps(i, sx);
						outPar << "MAX. No. OF STEPS/YEAR:\t ";
						if (ssteps.maxStepsYr == 99999999) outPar << "not applied" << endl;
						else outPar << ssteps.maxStepsYr << endl;
					}
				}
			}
		}
		else { // kernel-based transfer
			string notsuit = "IF THE ARRIVAL CELL/PATCH IS UNSUITABLE: ";
			string rchoose = " randomly choose a suitable neighb. cell/patch or ";
			string matereq = "MATING REQUIREMENTS: ";
			if (sett.sexDep) {
				nsexes = 2;
				outPar << sexdept << "yes" << endl;
				if (sett.stgDep) {
					nstages = sstruct.nStages;
					outPar << stgdept << "yes" << endl;
					outPar << notsuit << endl;
				}
				else {
					nstages = 1;
					outPar << stgdept << "no" << endl;
				}
			}
			else {
				nsexes = 1;
				outPar << sexdept << "no" << endl;
				if (sett.stgDep) {
					nstages = sstruct.nStages;
					outPar << stgdept << "yes" << endl;
					outPar << notsuit << endl;
				}
				else {
					nstages = 1;
					outPar << stgdept << "no" << endl;
					outPar << notsuit;
				}
			}
			for (int i = 0; i < nstages; i++) {
				if (sett.stgDep) {
					outPar << "stage " << i << ":" << endl;
				}
				for (int sx = 0; sx < nsexes; sx++) {
					if (sett.sexDep) {
						if (sx == 0) outPar << "FEMALES: ";
						else outPar << "MALES:   ";
						if (!sett.stgDep) outPar << notsuit;
					}
					srules = pSpecies->getSettRules(i, sx);
					if (srules.goToNeighbourLocn) {
						outPar << rchoose;
						if (srules.wait) outPar << "wait" << endl;
						else outPar << "die" << endl;
					}
					else {
						if (srules.wait) outPar << "wait" << endl;
						else outPar << "die" << endl;
					}
					outPar << matereq;
					if (srules.findMate) outPar << "yes" << endl;
					else outPar << "no" << endl;
				}
			}
		}

		// Genetics
		outPar << endl << "GENETICS:" << endl;
		set<TraitType> traitList = pSpecies->getTraitTypes();

		if (pSpecies->isDiploid()) outPar << "DIPLOID" << endl; else outPar << "HAPLOID" << endl;
		outPar << "Genome size: " << pSpecies->getGenomeSize() << endl;
		outPar << "Chromosome breaks : ";

		for (auto end : pSpecies->getChromosomeEnds())
			outPar << end << " ";
		outPar << endl;
		outPar << "Recombination rate: " << pSpecies->getRecombinationRate() << endl;
		outPar << "Traits modelled:  " << endl;
		for (auto trait : traitList)
			outPar << trait << endl;

		outPar << endl << "OUTPUTS:" << endl;
		outputParams out = pSpecies->getOutputParams();
		if (pSpecies->doesOutputRange()) {
		outPar << "Range - every " << out.outIntRange << " year";
		if (out.outIntRange > 1) outPar << "s";
			outPar << endl;
		}
		if (pSpecies->doesOutputOccup()) {
			outPar << "Occupancy - every " << pSpecies->getOutOccInt() << " year";
			if (pSpecies->getOutOccInt() > 1) outPar << "s";
			outPar << endl;
		}
		if (pSpecies->doesOutputPop()) {
			outPar << "Populations - every " << out.outIntPop << " year";
			if (out.outIntPop > 1) outPar << "s";
			if (out.outStartPop > 0) outPar << " starting year " << out.outStartPop;
			outPar << endl;
		}
		if (pSpecies->doesOutputInds()) {
			outPar << "Individuals - every " << out.outIntInd << " year";
			if (out.outIntInd > 1) outPar << "s";
			if (out.outStartInd > 0) outPar << " starting year " << out.outStartInd;
			outPar << endl;
		}
		if (pSpecies->doesOutputWeirCockerham() 
			|| pSpecies->doesOutputWeirHill()) {
			outPar << "Neutral genetics - every " << out.outputGeneticInterval << " year";
			if (out.outputGeneticInterval > 1) outPar << "s";
			if (out.outputWeirHill) outPar << " outputting pairwise patch fst";
			if (out.outputWeirCockerham) outPar << " outputting per locus fst ";
			outPar << endl;
		}

		if (pSpecies->doesOutputTraitCell()) {
			outPar << "Traits per ";
			if (ppLand.usesPatches) outPar << "patch"; else outPar << "cell";
			outPar << " - every " << out.outIntTraitCell << " year";
			if (out.outIntTraitCell > 1) outPar << "s";
			if (out.outStartTraitCell > 0) outPar << " starting year " << out.outStartTraitCell;
			outPar << endl;
		}
		if (pSpecies->doesOutputTraitRows()) {
			outPar << "Traits per row - every " << out.outIntTraitRow << " year";
			if (out.outIntTraitRow > 1) outPar << "s";
			if (out.outStartTraitRow > 0) outPar << " starting year " << out.outStartTraitRow;
			outPar << endl;
		}
		if (pSpecies->doesOutputConnect()) {
			outPar << "Connectivity matrix - every " << out.outIntConn << " year";
			if (out.outIntConn > 1) outPar << "s";
			if (out.outStartConn > 0) outPar << " starting year " << out.outStartConn;
			outPar << endl;
		}
#if RS_RCPP
		if (pSpecies->doesOutputPaths()) {
			outPar << "SMS paths - every " << out.outIntPaths << " year";
			if (out.outIntPaths > 1) outPar << "s";
			if (out.outStartPaths > 0) outPar << " starting year " << out.outStartPaths;
			outPar << endl;
		}
#endif

		if (trfr.usesMovtProc && trfr.moveType == 1) {
			outPar << "SMS HEAT MAPS: ";
			if (out.saveVisits) outPar << "yes" << endl;
			else outPar << "no" << endl;
		}
	}
	
	outPar.close(); 
	outPar.clear();
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
