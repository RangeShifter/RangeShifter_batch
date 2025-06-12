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

#include "Community.h"

//---------------------------------------------------------------------------

Community::Community(Landscape* pLand, speciesMap_t allSpecies) {
	pLandscape = pLand;
	speciesMap = allSpecies;
	indIx = 0;
	resetActiveSpecies(); // all species are active by default

	// Populate species maps
	for (auto& [sp, pSpecies] : allSpecies) {
		allPopns.emplace(sp, vector<Population*>());
		neutralStatsMaps.emplace(sp, nullptr);
		// Output file streams
		if (pSpecies->doesOutputInds()) outIndsOfs.emplace(sp, ofstream());
		if (pSpecies->doesOutputPop()) outPopOfs.emplace(sp, ofstream());
		if (pSpecies->doesOutputRange()) outRangeOfs.emplace(sp, ofstream());
		if (pSpecies->doesOutputGeneValues()) ofsGenes.emplace(sp, ofstream());
		if (pSpecies->doesOutputWeirHill()) outPairwiseFstOfs.emplace(sp, ofstream());
		if (pSpecies->doesOutputWeirCockerham()) {
			outWCFstatOfs.emplace(sp, ofstream());
			outPerLocusFstat.emplace(sp, ofstream());
		}
		if (pSpecies->doesOutputOccup()) {
			outSuitOfs.emplace(sp, ofstream());
		}
		if (pSpecies->doesOutputTraitRows())
			outTraitsRows.emplace(sp, ofstream());
		if (pSpecies->doesOutputTraitCell())
			outTraitsOfs.emplace(sp, ofstream());
	}
}

Community::~Community() {
	for (auto& [speciesID, popns] : allPopns) {
		for (int i = 0; i < popns.size(); i++) {
			delete popns[i];
		}
		popns.clear();
	}
	for (auto& [speciesID, mtxPop] : matrixPops) {
		delete mtxPop;
	}
}

void Community::initialise(Species* pSpecies, int year) {
	int npatches, ndistcells, patchnum, candidatePatch = 0;
	locn distloc;
	patchData pch;
	patchLimits limits = patchLimits();
	set<int> selectedPatches;
	set<int> suitablePatches;
	Patch* pPatch;
	Cell* pCell;
	landParams ppLand = pLandscape->getLandParams();
	initParams init = pSpecies->getInitParams();
	species_id sp = pSpecies->getID();
	int spratio;
	if (init.seedType == 1)
		spratio = pLandscape->getSpDistResol(sp) / ppLand.resol;

	// Initialise (empty) matrix populations
	matrixPops.emplace(sp,
		new Population(pSpecies, pLandscape->findPatch(sp, 0), 0, ppLand.resol)
	);

	switch (init.seedType) {

	case 0:	// free initialisation

		switch (init.freeType) {

		case 0:	// random
			// determine no. of patches / cells within the specified initialisation limits
			// and record their corresponding sub-communities in a list
			// parallel list records which have been selected
			npatches = pLandscape->getPatchCount(sp);
			limits.xMin = init.minSeedX;
			limits.xMax = init.maxSeedX;
			limits.yMin = init.minSeedY;
			limits.yMax = init.maxSeedY;

			for (int i = 0; i < npatches; i++) {
				pch = pLandscape->getPatchData(sp, i);
				patchnum = pch.pPatch->getPatchNum();
				if (pch.pPatch->withinLimits(limits)) {
					if (ppLand.usesPatches) {
						if (patchnum != 0) {
							suitablePatches.insert(patchnum);
						}
					}
					else { // cell-based model - is cell(patch) suitable
						if (pch.pPatch->isSuitable()) {
							suitablePatches.insert(patchnum);
						}
					}
				}
			}

			// select specified no. of patches/cells at random
			sample(
				suitablePatches.begin(),
				suitablePatches.end(),
				inserter(selectedPatches, selectedPatches.begin()),
				init.nSeedPatches,
				pRandom->getRNG()
			);
			break;

		case 1:	// all suitable patches/cells
			npatches = pLandscape->getPatchCount(sp);
			limits.xMin = init.minSeedX;
			limits.xMax = init.maxSeedX;
			limits.yMin = init.minSeedY;
			limits.yMax = init.maxSeedY;

			for (int i = 0; i < npatches; i++) {
				pch = pLandscape->getPatchData(sp, i);
				if (pch.pPatch->withinLimits(limits)) {
					patchnum = pch.pPatch->getPatchNum();
					if (!pch.pPatch->isMatrix() && pch.pPatch->isSuitable()) {
						selectedPatches.insert(patchnum);
					}
				}
			}

			break;

		} // end of switch (init.freeType)

		for (auto pchNum : selectedPatches) {
			Patch* pPatch = pLandscape->findPatch(sp, pchNum);
			// Determine size of initial population
			int nInds = pPatch->getInitNbInds(init, ppLand.usesPatches, ppLand.resol);
			if (nInds > 0) {
				Population* pPop = new Population(pSpecies, pPatch, nInds, ppLand.resol);
				allPopns.at(sp).push_back(pPop); // add new population to community list
			}
		}
		break;

	case 1:	// from species distribution
		// initialise from loaded species distribution
		switch (init.spDistType) {
		case 0: // all presence cells
			pLandscape->setDistribution(sp, 0); // activate all patches
			break;
		case 1: // some randomly selected presence cells
			pLandscape->setDistribution(sp, init.nSpDistPatches); // activate random patches
			break;
		}

		ndistcells = pLandscape->distCellCount(0);
		for (int i = 0; i < ndistcells; i++) {
			distloc = pLandscape->getSelectedDistnCell(0, i);
			if (distloc.x >= 0) { // distribution cell is selected
				// process each landscape cell within the distribution cell

				for (int x = 0; x < spratio; x++) {
					for (int y = 0; y < spratio; y++) {
						pCell = pLandscape->findCell(distloc.x * spratio + x, distloc.y * spratio + y);
						if (pCell != nullptr) { // not a no-data cell
							pPatch = pCell->getPatch(sp);
							if (pPatch != nullptr) {
								if (!pPatch->isMatrix()) { // not the matrix patch
									selectedPatches.insert(pPatch->getPatchNum());
								}
							}
						}
					}
				}

			}
		}

		for (auto pchNum : selectedPatches) {
			Patch* pPatch = pLandscape->findPatch(sp, pchNum);
			// Determine size of initial population
			int nInds = pPatch->getInitNbInds(init, ppLand.usesPatches, ppLand.resol);
			if (nInds > 0) {
				Population* pPop = new Population(pSpecies, pPatch, nInds, ppLand.resol);
				allPopns.at(pSpecies->getID()).push_back(pPop); // add new population to community list
			}
		}
		break;

	case 2:	// initial individuals in specified patches/cells

		if (year == 0) indIx = 0; // reset index
		initInd iind = initInd();
		iind.year = 0;
		int ninds = pSpecies->getNbInitInds();
		while (indIx < ninds && iind.year <= year) {
			iind = pSpecies->getInitInd(indIx);
			while (iind.year == year && iind.speciesID == sp) {
				if (ppLand.usesPatches) {
					pPatch = pLandscape->findPatch(sp, iind.patchID);
					if (pPatch != nullptr) { // exists
						if (pPatch->isSuitable()) {
							initialInd(pLandscape, pSpecies, pPatch, pPatch->getRandomCell(), indIx);
						}
					}
				}
				else { // cell-based model
					pCell = pLandscape->findCell(iind.x, iind.y);
					if (pCell != nullptr) {
						pPatch = pCell->getPatch(sp);
						if (pPatch != nullptr) {
							if (pPatch->isSuitable()) {
								initialInd(pLandscape, pSpecies, pPatch, pCell, indIx);
							}
						}
					}
				}
				indIx++;
				if (indIx < ninds) {
					iind = pSpecies->getInitInd(indIx);
				}
				else {
					iind.year = 99999999;
				}
			}
		}
		break;
	} // end of switch (init.seedType)

}

Species* Community::findSpecies(species_id id) {
	if (auto search = speciesMap.find(id); search != speciesMap.end()) {
		return search->second;
	}
	else throw logic_error("Species " + to_string(id) + " couldn't be found.");
}

void Community::resetPopns() {

	for (auto& [sp, mtxPop] : matrixPops)
		mtxPop->getPatch()->resetPop();

	for (auto& [sp, popns] : allPopns) {
		for (auto pop : popns) {
			pop->getPatch()->resetPop();
		}
		popns.clear();
	}

	// reset the individual ids to start from zero
	Individual::indCounter = 0;
}

void Community::resetActiveSpecies() {
	activeSpecies.clear();
	for (auto& sp : views::keys(speciesMap))
		activeSpecies.insert(sp);
}

void Community::disableInactiveSpecies(int gen) {
	set<species_id> actSpCopy = activeSpecies; 
	// Since we're erasing elements we need to iterate over a copy

	for (auto& sp : actSpCopy) {
		int nbSeasons = speciesMap.at(sp)->getDemogrParams().repSeasons;
		if (gen >= nbSeasons) activeSpecies.erase(sp);
	}
}

void Community::applyRandLocExt(const float& probExt) {
	for (auto& sp : activeSpecies) {
		const float probExt = speciesMap.at(sp)->getLocalExtProb();
		if (probExt == 0.0) continue;
		for (auto pop : allPopns.at(sp)) {
			if (pRandom->Bernoulli(probExt))
				pop->extirpate();
		}
	}
}

void Community::applyLocalExtGrad() {
	for (auto& sp : activeSpecies) {
		for (auto pop : allPopns.at(sp)) {
			pop->applyLocalExtGrad();
		}
	}
}

void Community::scanUnsuitablePatches(Species* pSpecies) {
	for (auto pPop : allPopns.at(pSpecies->getID())) {
		float localK = pPop->getPatch()->getK();
		if (localK <= 0.0) { // patch in dynamic landscape has become unsuitable
			if (pSpecies->getDemogrParams().stageStruct) {
				if (pSpecies->getStageParams().disperseOnLoss)
					pPop->allEmigrate();
				else pPop->extirpate();
			}
			else { // non-stage-structured species is destroyed
				pPop->extirpate();
			}
		}
	}
}

void Community::reproduction(int yr)
{
	float eps = 0.0; // epsilon for environmental stochasticity
	landParams land = pLandscape->getLandParams();
	
	for (auto& sp : activeSpecies) {
		for (auto pop : allPopns.at(sp)) {
			Patch* pPatch = pop->getPatch();
			float localK = pPatch->getK();
			if (localK > 0.0) {
				pop->reproduction(localK, land.resol);
				pop->fledge();
			}
		}
	}
}

void Community::emigration()
{
	for (auto& sp : activeSpecies) {
		for (auto pop : allPopns.at(sp)) {
			pop->emigration(pop->getPatch()->getK());
		}
	}
}

void Community::dispersal(short landIx, short nextseason)
{
	simParams sim = paramsSim->getSim();

	// initiate dispersal - all emigrants leave their natal community and join matrix community
	for (auto& sp : activeSpecies) {
		for (auto pop : allPopns.at(sp)) {

			for (int j = 0; j < pop->getStats().nInds; j++) {
				disperser disp = pop->extractDisperser(j);
				if (disp.isDispersing) { // disperser - has already been removed from natal population
					short spID = pop->getSpecies()->getID();
					auto it = matrixPops.find(spID);
					if (it != matrixPops.end())
						it->second->recruit(disp.pInd); // add to matrix population
					else throw runtime_error("");
				}
			}
			// remove pointers to emigrants
			pop->clean();
		}
	}

	// dispersal is undertaken by all individuals now in the matrix patch
	// (even if not physically in the matrix)
	int ndispersers = 0;
	do {
		// Reset possible settlers for all patches before transfer
		for (auto& sp : activeSpecies) {
			matrixPops.at(sp)->getPatch()->resetPossSettlers();
		}
		for (auto& sp : activeSpecies) {
			for (auto pop : allPopns.at(sp)) {
				pop->getPatch()->resetPossSettlers();
			}
		}

		// Transfer takes place in the matrix
		for (auto& sp : activeSpecies) {
			ndispersers = matrixPops.at(sp)->transfer(pLandscape, landIx, nextseason);
		}
		completeDispersal(pLandscape);

	} while (ndispersers > 0);
}

// Remove emigrants from patch 0 (matrix) and transfer to the population
// in which their destination co-ordinates fall
// This function is executed for matrix patch populations only
void Community::completeDispersal(Landscape* pLandscape)
{
	Population* pPop;
	Patch* pPrevPatch;
	Patch* pNewPatch;
	Cell* pPrevCell;

	for (auto& sp : activeSpecies) {

		Population* mtxPop = matrixPops.at(sp);
		int popsize = mtxPop->getNInds();
		for (int j = 0; j < popsize; j++) {

			disperser settler = mtxPop->extractSettler(j);
			if (settler.isSettling) {
				// settler - has already been removed from matrix population
				// find new patch
				pNewPatch = settler.pCell->getPatch(sp);

				// find population within the patch (if there is one)
				pPop = pNewPatch->getPop();

				if (pPop == nullptr) { // settler is the first in a previously uninhabited patch
					// create a new population in the corresponding sub-community
					pPop = new Population(mtxPop->getSpecies(), pNewPatch, 0, pLandscape->getLandParams().resol);
					allPopns.at(sp).push_back(pPop); // add new pop to community list
				}

				pPop->recruit(settler.pInd);

				if (mtxPop->getSpecies()->doesOutputConnect()) { // increment connectivity totals
					int newpatch = pNewPatch->getSeqNum();
					pPrevCell = settler.pInd->getLocn(0); // previous cell
					Patch* pPatch = pPrevCell->getPatch(sp);
					if (pPatch != nullptr) {
						pPrevPatch = pPatch;
						int prevpatch = pPrevPatch->getSeqNum();
						pLandscape->incrConnectMatrix(sp, prevpatch, newpatch);
					}
				}
			}
		}
		// remove pointers in the matrix popn to settlers
		mtxPop->clean();
	}
}

// initialise a specified individual
void Community::initialInd(Landscape* pLandscape, Species* pSpecies,
	Patch* pPatch, Cell* pCell, int ix)
{
	demogrParams dem = pSpecies->getDemogrParams();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	short stg, age, repInt;

	Population* pPop = pPatch->getPop();
	// create new population if not already in existence
	if (pPop == nullptr) {
		pPop = new Population(pSpecies, pPatch, 0, pLandscape->getLandParams().resol);
		pPatch->setPop(pPop);
		allPopns.at(pSpecies->getID()).push_back(pPop);
	}

	// create new individual
	initInd iind = pSpecies->getInitInd(ix);
	if (dem.stageStruct) {
		stg = iind.stage;
		age = iind.age;
		repInt = pSpecies->getStageParams().repInterval;
	}
	else {
		age = stg = 1;
		repInt = 0;
	}
	float probmale = (dem.repType != 0 && iind.sex == 1) ? 1.0 : 0.0;

	Individual* pInd = new Individual(pSpecies, pCell, pPatch, stg, age, repInt, probmale, trfr.usesMovtProc, trfr.moveType);

	// add new individual to the population
	pPop->recruit(pInd);

	if (pSpecies->getNTraits() > 0) {
		// individual variation - set up genetics
		pInd->setUpGenes(pLandscape->getLandData().resol);
	}
}

void Community::drawSurvivalDevlpt(const int phase)
{	
	for (auto& sp: activeSpecies) {
		
		bool hasStages = speciesMap.at(sp)->stageStructured();
		short survOption = speciesMap.at(sp)->getStageParams().survival;

		switch (phase) {
		case 0: { // After reproduction, before dispersal
			if (hasStages && survOption == phase) {
				// Survival + developments adults
				matrixPops.at(sp)->drawSurvivalDevlpt(false, true, true, true);
				for (auto pop : allPopns.at(sp)) {
					pop->drawSurvivalDevlpt(false, true, true, true);
				}
			}
			break;
		}
		case 1: { // After dispersal
			bool resolveJuvs = true;
			bool resolveAdults = !(hasStages && survOption == 0); // else already resolved after reproduction
			bool resolveDev = true;
			bool resolveSurv = !(hasStages && survOption == 2); // else resolved yearly

			matrixPops.at(sp)->drawSurvivalDevlpt(resolveJuvs, resolveAdults, resolveDev, resolveSurv);
			for (auto pop : allPopns.at(sp)) {
				pop->drawSurvivalDevlpt(resolveJuvs, resolveAdults, resolveDev, resolveSurv);
			}
			break;
		}
		case 2: { // End of year
			if (hasStages && survOption == phase) {
				// Survival juveniles + adults
				matrixPops.at(sp)->drawSurvivalDevlpt(true, true, false, true);
				for (auto pop : allPopns.at(sp)) {
					pop->drawSurvivalDevlpt(true, true, false, true);
				}
			}
			break;
		}
		}
	}
}

void Community::applySurvivalDevlpt() {
	for (auto& sp : activeSpecies) {
		matrixPops.at(sp)->applySurvivalDevlpt();
		for (auto pop : allPopns.at(sp))
			pop->applySurvivalDevlpt();
	}
}

void Community::ageIncrement() {
	for (auto& [spId, mtxPop] : matrixPops) {
		mtxPop->ageIncrement();
	}
	for (auto& [sp, popns] : allPopns) {
		for (auto pop : popns)
			pop->ageIncrement();
	}
}

// Calculate total no. of individuals of all species
int Community::totalInds() {
	int total = 0;
	for (auto& [spId, mtxPop] : matrixPops) {
		total += mtxPop->getStats().nInds;
	}
	for (auto& [sp, popns] : allPopns) {
		for (auto pop : popns) {
			total += pop->getStats().nInds;
		}
	}
	return total;
}

//---------------------------------------------------------------------------
void Community::createOccupancy(species_id sp, int nbOutputRows, int nbReps) {

	pLandscape->createOccupancy(sp, nbOutputRows);
	
	// Initialise array for occupancy of suitable cells/patches
	vector<vector<int>> occupancyMap;
	for (int i = 0; i < nbOutputRows; i++) {
		occupancyMap.push_back(vector<int>(nbReps, 0));
	}
	occupancyMaps.emplace(sp, occupancyMap);
}

void Community::updateOccupancy(species_id sp, int yr, int rep) {

	int whichRow = yr / speciesMap.at(sp)->getOutOccInt();
	matrixPops.at(sp)->getPatch()->updateOccupancy(whichRow);
	for (auto pop : allPopns.at(sp)) {
		pop->getPatch()->updateOccupancy(whichRow);
	}
	commStats s = getStats(sp);
	occupancyMaps.at(sp)[whichRow][rep] = trunc(s.occupied / static_cast<double>(s.suitable));
}

//---------------------------------------------------------------------------
// Count no. of sub-communities (suitable patches) and those occupied (non-zero populations)
// Determine range margins
commStats Community::getStats(species_id sp)
{
	commStats s = commStats();
	landParams ppLand = pLandscape->getLandParams();
	s.suitable = s.occupied = 0;
	s.minX = ppLand.maxX;
	s.minY = ppLand.maxY;
	s.maxX = s.maxY = 0;
	float localK;
	popStats patchPop;

	// Count individuals for the matrix
	s.ninds = 0;
	s.nnonjuvs = 0;
	s.ninds += matrixPops.at(sp)->getStats().nInds;
	s.nnonjuvs += matrixPops.at(sp)->getStats().nNonJuvs;

	for (auto pPop : allPopns.at(sp)) {

		patchPop = pPop->getStats();
		s.ninds += patchPop.nInds;
		s.nnonjuvs += patchPop.nNonJuvs;

		if (patchPop.pPatch != nullptr) {

			if (patchPop.pPatch->isSuitable() > 0.0) s.suitable++;
			if (patchPop.nInds > 0 && patchPop.breeding) {
				s.occupied++;
				patchLimits pchlim = patchPop.pPatch->getLimits();
				if (pchlim.xMin < s.minX) s.minX = pchlim.xMin;
				if (pchlim.xMax > s.maxX) s.maxX = pchlim.xMax;
				if (pchlim.yMin < s.minY) s.minY = pchlim.yMin;
				if (pchlim.yMax > s.maxY) s.maxY = pchlim.yMax;
			}
		}
	}
	return s;
}

//---------------------------------------------------------------------------

// Functions to control production of output files

// For outputs and population visualisations pre-reproduction
void Community::popAndRangeOutput(int rep, int yr, int gen) {

	for (auto& sp : activeSpecies) {

		Species* pSpecies = speciesMap.at(sp);
		if (pSpecies->isRangeOutputYear(yr))
			outRange(sp, rep, yr, gen);
		if (pSpecies->isPopOutputYear(yr))
			outPop(sp, rep, yr, gen);
	}
}

// Open population file and write header record
bool Community::outPopHeaders(Species* pSpecies) {

	landParams land = pLandscape->getLandParams();
	simParams sim = paramsSim->getSim();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	string name = paramsSim->getDir(2)
		+ (sim.batchMode ? "Batch" + to_string(sim.batchNum) + "_" : "")
		+ "Sim" + to_string(sim.simulation) + "_Land" + to_string(land.landNum) + 
		+"_Species" + to_string(pSpecies->getID()) + "_Pop.txt";

	ofstream& popOfs = outPopOfs.at(pSpecies->getID());

	popOfs.open(name.c_str());

	popOfs << "Rep\tYear\tRepSeason";
	if (land.usesPatches) popOfs << "\tPatchID\tNcells";
	else popOfs << "\tx\ty";

	// determine whether environmental data need be written for populations
	bool writeEnv = pSpecies->usesGradient() || paramsStoch->envStoch();
	if (writeEnv) popOfs << "\tEpsilon\tGradient\tLocal_K";

	popOfs << "\tNInd";
	if (dem.stageStruct) {
		if (dem.repType == 0) {
			for (int i = 1; i < sstruct.nStages; i++) popOfs << "\tNInd_stage" << i;
			popOfs << "\tNJuvs";
		}
		else {
			for (int i = 1; i < sstruct.nStages; i++)
				popOfs << "\tNfemales_stage" << i << "\tNmales_stage" << i;
			popOfs << "\tNJuvFemales\tNJuvMales";
		}
	}
	else {
		if (dem.repType != 0) popOfs << "\tNfemales\tNmales";
	}
	popOfs << endl;
	return popOfs.is_open();
}

bool Community::closePopOfs(species_id sp) {
	if (outPopOfs.at(sp).is_open())
		outPopOfs.at(sp).close();
	outPopOfs.at(sp).clear();
	return true;
}

// Write records to population file
void Community::outPop(species_id sp, int rep, int yr, int gen) {
	Species* pSpecies = speciesMap.at(sp);
	landParams land = pLandscape->getLandParams();
	envStochParams env = paramsStoch->getStoch();
	bool writeEnv = pSpecies->usesGradient() || env.usesStoch;
	bool gradK = pSpecies->usesGradient() && pSpecies->getEnvGradient().gradType == 1;

	float eps = 0.0;
	if (env.usesStoch && !env.stochIsLocal) {
		eps = pLandscape->getGlobalStoch(yr);
	}

	// generate output for each population (patch x species) in the community
	if (matrixPops.at(sp)->totalPop() > 0) {
		matrixPops.at(sp)->outPopulation(outPopOfs.at(sp), rep, yr, gen, env.stochIsLocal, eps, land.usesPatches, writeEnv, gradK);
	}
	for (auto pop : allPopns.at(sp)) {
		if (pop->getPatch()->isSuitable() || pop->totalPop() > 0) {
			pop->outPopulation(outPopOfs.at(sp), rep, yr, gen, env.stochIsLocal, eps, land.usesPatches, writeEnv, gradK);
		}
	}
}

void Community::indsAndGeneticsOutput(int rep, int yr, int gen) {

	for (auto& sp : activeSpecies) {

		Species* pSpecies = speciesMap.at(sp);

		// Output Individuals
		if (pSpecies->isIndOutputYear(yr))
			outInds(sp, rep, yr, gen);

		// Output Genetics
		if (pSpecies->isGeneticOutputYear(yr)) {
			if (pSpecies->getSamplingOption() == "random_occupied"
				|| pSpecies->getSamplingOption() == "all")
				// then must re-sample every year
				pLandscape->samplePatches(pSpecies);

			sampleIndividuals(sp);

			if (pSpecies->doesOutputGeneValues())
				outputGeneValues(sp, yr, gen);
			if (pSpecies->doesOutputWeirCockerham() 
				|| pSpecies->doesOutputWeirHill())
				outNeutralGenetics(sp, rep, yr, gen);
		}
	}
}

// Open individuals file and write header record
void Community::outIndsHeaders(species_id sp, int rep, int landNr, bool usesPatches)
{
	string name;
	Species* pSpecies = speciesMap.at(sp);
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();
	bool hasGenLoad = pSpecies->getNbGenLoadTraits() > 0;

	name = paramsSim->getDir(2)
		+ (sim.batchMode ? "Batch" + to_string(sim.batchNum) + "_" : "")
		+ "Sim" + to_string(sim.simulation)
		+ "_Land" + to_string(landNr) +
		"_Rep" + to_string(rep) +
		"_Species" + to_string(sp) +
		"_Inds.txt";

	ofstream& indsOfs = outIndsOfs.at(sp);

	indsOfs.open(name.c_str());
	indsOfs << "Rep\tYear\tRepSeason\tIndID\tStatus";
	if (usesPatches) indsOfs << "\tNatal_patch\tPatchID";
	else indsOfs << "\tNatal_X\tNatal_Y\tX\tY";
	if (dem.repType != 0) indsOfs << "\tSex";
	if (dem.stageStruct) indsOfs << "\tAge\tStage";
	if (hasGenLoad) indsOfs << "\tProbViable";
	if (emig.indVar) {
		if (emig.densDep) indsOfs << "\tD0\tAlpha\tBeta";
		else indsOfs << "\tEP";
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 1) { // SMS
				indsOfs << "\tDP\tGB\tAlphaDB\tBetaDB";
			}
			if (trfr.moveType == 2) { // CRW
				indsOfs << "\tStepLength\tRho";
			}
		}
		else { // kernel
			indsOfs << "\tMeanDistI";
			if (trfr.twinKern) indsOfs << "\tMeanDistII\tPKernelI";
		}
	}
	if (sett.indVar) {
		indsOfs << "\tS0\tAlphaS\tBetaS";
	}
	indsOfs << "\tDistMoved";
#ifndef NDEBUG
	indsOfs << "\tNsteps";
#else
	if (trfr.usesMovtProc) indsOfs << "\tNsteps";
#endif
	indsOfs << endl;
}

void Community::closeOutIndsOfs(species_id sp) {
	if (outIndsOfs.at(sp).is_open()) {
		outIndsOfs.at(sp).close();
		outIndsOfs.at(sp).clear();
	}
	return;
}

// Write records to individuals file
void Community::outInds(species_id sp, int rep, int yr, int gen) {
	matrixPops.at(sp)->outIndividual(outIndsOfs.at(sp), pLandscape, rep, yr, gen);
	for (Population* pop : allPopns.at(sp)) // all sub-communities
		pop->outIndividual(outIndsOfs.at(sp), pLandscape, rep, yr, gen);
}

bool Community::closeRangeOfs(species_id sp) {
	if (outRangeOfs.at(sp).is_open())
		outRangeOfs.at(sp).close();
	outRangeOfs.at(sp).clear();
	return true;
}

// Open range file and write header record
bool Community::outRangeHeaders(species_id sp, int landNr)
{
	string name;
	landParams ppLand = pLandscape->getLandParams();
	envStochParams env = paramsStoch->getStoch();
	simParams sim = paramsSim->getSim();

	Species* pSpecies = speciesMap.at(sp);
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + to_string(sim.batchNum) +
			+"_Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr)
			+ "_Species" + to_string(sp)
			+ "_Range.txt";
	}
	else {
		name = paramsSim->getDir(2)
			+ "Sim" + to_string(sim.simulation)
			+ "_Species" + to_string(sp)
			+ "_Range.txt";
	}

	ofstream& rangeOfs = outRangeOfs.at(sp);
	rangeOfs.open(name.c_str());
	rangeOfs << "Rep\tYear\tRepSeason";
	if (env.usesStoch && !env.stochIsLocal) rangeOfs << "\tEpsilon";

	rangeOfs << "\tNInds";
	if (dem.stageStruct) {
		for (int i = 1; i < sstruct.nStages; i++) rangeOfs << "\tNInd_stage" << i;
		rangeOfs << "\tNJuvs";
	}
	if (ppLand.usesPatches) rangeOfs << "\tNOccupPatches";
	else rangeOfs << "\tNOccupCells";
	rangeOfs << "\tOccup/Suit\tmin_X\tmax_X\tmin_Y\tmax_Y";

	if (emig.indVar) {
		if (emig.sexDep) {
			if (emig.densDep) {
				rangeOfs << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
				rangeOfs << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
				rangeOfs << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
			}
			else {
				rangeOfs << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
			}
		}
		else {
			if (emig.densDep) {
				rangeOfs << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
				rangeOfs << "\tmeanBeta\tstdBeta";
			}
			else {
				rangeOfs << "\tmeanEP\tstdEP";
			}
		}
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 1) {
				rangeOfs << "\tmeanDP\tstdDP\tmeanGB\tstdGB";
				rangeOfs << "\tmeanAlphaDB\tstdAlphaDB\tmeanBetaDB\tstdBetaDB";
			}
			if (trfr.moveType == 2) {
				rangeOfs << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
			}
		}
		else {
			if (trfr.sexDep) {
				rangeOfs << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
				if (trfr.twinKern)
					rangeOfs << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
			}
			else {
				rangeOfs << "\tmean_distI\tstd_distI";
				if (trfr.twinKern)
					rangeOfs << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
			}
		}
	}
	if (sett.indVar) {
		if (sett.sexDep) {
			rangeOfs << "\tF_meanS0\tF_stdS0\tM_meanS0\tM_stdS0";
			rangeOfs << "\tF_meanAlphaS\tF_stdAlphaS\tM_meanAlphaS\tM_stdAlphaS";
			rangeOfs << "\tF_meanBetaS\tF_stdBetaS\tM_meanBetaS\tM_stdBetaS";

		}
		else {
			rangeOfs << "\tmeanS0\tstdS0";
			rangeOfs << "\tmeanAlphaS\tstdAlphaS";
			rangeOfs << "\tmeanBetaS\tstdBetaS";
		}
	}
	rangeOfs << endl;
	return rangeOfs.is_open();
}

// Write record to range file
void Community::outRange(species_id sp, int rep, int yr, int gen)
{
	landParams ppLand = pLandscape->getLandParams();
	envStochParams env = paramsStoch->getStoch();

	Species* pSpecies = speciesMap.at(sp);
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	ofstream& rangeOfs = outRangeOfs.at(sp);
	rangeOfs << rep << "\t" << yr << "\t" << gen;
	if (env.usesStoch && !env.stochIsLocal) // write global environmental stochasticity
		rangeOfs << "\t" << pLandscape->getGlobalStoch(yr);

	commStats s = getStats(sp);

	if (dem.stageStruct) {
		rangeOfs << "\t" << s.nnonjuvs;
		int stagepop;
		// all non-juvenile stages
		for (int stg = 1; stg < sstruct.nStages; stg++) {
			stagepop = 0;
			for (auto& [spId, mtxPop] : matrixPops) {
				stagepop += mtxPop->stagePop(stg);
			}
			for (auto pop : allPopns.at(sp)) {
				stagepop += pop->stagePop(stg);
			}
			rangeOfs << "\t" << stagepop;
		}
		// juveniles born in current reproductive season
		stagepop = 0;
		for (auto& [spId, mtxPop] : matrixPops) {
			stagepop += mtxPop->stagePop(0);
		}
		for (auto pop : allPopns.at(sp)) {
			stagepop += pop->stagePop(0);
		}
		rangeOfs << "\t" << stagepop;
	}
	else { // non-structured species
		rangeOfs << "\t" << s.ninds;
	}

	float occsuit = 0.0;
	if (s.suitable > 0) occsuit = (float)s.occupied / (float)s.suitable;
	rangeOfs << "\t" << s.occupied << "\t" << occsuit;
	// RANGE MINIMA AND MAXIMA NEED TO BECOME A PROPERTY OF THE SPECIES
	if (s.ninds > 0) {
		landOrigin originVal = pLandscape->getOrigin();
		rangeOfs << "\t" << (float)s.minX * (float)ppLand.resol + originVal.minEast
			<< "\t" << (float)(s.maxX + 1) * (float)ppLand.resol + originVal.minEast
			<< "\t" << (float)s.minY * (float)ppLand.resol + originVal.minNorth
			<< "\t" << (float)(s.maxY + 1) * (float)ppLand.resol + originVal.minNorth;
	}
	else rangeOfs << "\t0\t0\t0\t0";

	if (emig.indVar || trfr.indVar || sett.indVar) { // output trait means
		traitsums ts = traitsums();
		traitsums popTraits;
		int ngenes, popsize;
		for (auto& [spId, mtxPop] : matrixPops) {
			popTraits = mtxPop->outTraits(outTraitsOfs.at(sp), false);
			for (int j = 0; j < gMaxNbSexes; j++) {
				ts.ninds[j] += popTraits.ninds[j];
				ts.sumD0[j] += popTraits.sumD0[j];     ts.ssqD0[j] += popTraits.ssqD0[j];
				ts.sumAlpha[j] += popTraits.sumAlpha[j];  ts.ssqAlpha[j] += popTraits.ssqAlpha[j];
				ts.sumBeta[j] += popTraits.sumBeta[j];   ts.ssqBeta[j] += popTraits.ssqBeta[j];
				ts.sumDist1[j] += popTraits.sumDist1[j];  ts.ssqDist1[j] += popTraits.ssqDist1[j];
				ts.sumDist2[j] += popTraits.sumDist2[j];  ts.ssqDist2[j] += popTraits.ssqDist2[j];
				ts.sumProp1[j] += popTraits.sumProp1[j];  ts.ssqProp1[j] += popTraits.ssqProp1[j];
				ts.sumDP[j] += popTraits.sumDP[j];     ts.ssqDP[j] += popTraits.ssqDP[j];
				ts.sumGB[j] += popTraits.sumGB[j];     ts.ssqGB[j] += popTraits.ssqGB[j];
				ts.sumAlphaDB[j] += popTraits.sumAlphaDB[j]; ts.ssqAlphaDB[j] += popTraits.ssqAlphaDB[j];
				ts.sumBetaDB[j] += popTraits.sumBetaDB[j];  ts.ssqBetaDB[j] += popTraits.ssqBetaDB[j];
				ts.sumStepL[j] += popTraits.sumStepL[j];  ts.ssqStepL[j] += popTraits.ssqStepL[j];
				ts.sumRho[j] += popTraits.sumRho[j];    ts.ssqRho[j] += popTraits.ssqRho[j];
				ts.sumS0[j] += popTraits.sumS0[j];     ts.ssqS0[j] += popTraits.ssqS0[j];
				ts.sumAlphaS[j] += popTraits.sumAlphaS[j]; ts.ssqAlphaS[j] += popTraits.ssqAlphaS[j];
				ts.sumBetaS[j] += popTraits.sumBetaS[j];  ts.ssqBetaS[j] += popTraits.ssqBetaS[j];
			}
		}
		int npops = static_cast<int>(allPopns.at(sp).size());
		for (int i = 0; i < npops; i++) {
			popTraits = allPopns.at(sp)[i]->outTraits(outTraitsOfs.at(sp), false);
			for (int j = 0; j < gMaxNbSexes; j++) {
				ts.ninds[j] += popTraits.ninds[j];
				ts.sumD0[j] += popTraits.sumD0[j];     ts.ssqD0[j] += popTraits.ssqD0[j];
				ts.sumAlpha[j] += popTraits.sumAlpha[j];  ts.ssqAlpha[j] += popTraits.ssqAlpha[j];
				ts.sumBeta[j] += popTraits.sumBeta[j];   ts.ssqBeta[j] += popTraits.ssqBeta[j];
				ts.sumDist1[j] += popTraits.sumDist1[j];  ts.ssqDist1[j] += popTraits.ssqDist1[j];
				ts.sumDist2[j] += popTraits.sumDist2[j];  ts.ssqDist2[j] += popTraits.ssqDist2[j];
				ts.sumProp1[j] += popTraits.sumProp1[j];  ts.ssqProp1[j] += popTraits.ssqProp1[j];
				ts.sumDP[j] += popTraits.sumDP[j];     ts.ssqDP[j] += popTraits.ssqDP[j];
				ts.sumGB[j] += popTraits.sumGB[j];     ts.ssqGB[j] += popTraits.ssqGB[j];
				ts.sumAlphaDB[j] += popTraits.sumAlphaDB[j]; ts.ssqAlphaDB[j] += popTraits.ssqAlphaDB[j];
				ts.sumBetaDB[j] += popTraits.sumBetaDB[j];  ts.ssqBetaDB[j] += popTraits.ssqBetaDB[j];
				ts.sumStepL[j] += popTraits.sumStepL[j];  ts.ssqStepL[j] += popTraits.ssqStepL[j];
				ts.sumRho[j] += popTraits.sumRho[j];    ts.ssqRho[j] += popTraits.ssqRho[j];
				ts.sumS0[j] += popTraits.sumS0[j];     ts.ssqS0[j] += popTraits.ssqS0[j];
				ts.sumAlphaS[j] += popTraits.sumAlphaS[j]; ts.ssqAlphaS[j] += popTraits.ssqAlphaS[j];
				ts.sumBetaS[j] += popTraits.sumBetaS[j];  ts.ssqBetaS[j] += popTraits.ssqBetaS[j];
			}
		}

		if (emig.indVar) {
			if (emig.sexDep) { // must be a sexual species
				ngenes = 2;
			}
			else {
				if (dem.repType == 0) { // asexual reproduction
					ngenes = 1;
				}
				else { // sexual reproduction
					ngenes = 1;
				}
			}
			double mnD0[2], mnAlpha[2], mnBeta[2], sdD0[2], sdAlpha[2], sdBeta[2];
			for (int g = 0; g < ngenes; g++) {
				mnD0[g] = mnAlpha[g] = mnBeta[g] = sdD0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ngenes == 2) popsize = ts.ninds[g];
				else popsize = ts.ninds[0] + ts.ninds[1];
				if (popsize > 0) {
					mnD0[g] = ts.sumD0[g] / (double)popsize;
					mnAlpha[g] = ts.sumAlpha[g] / (double)popsize;
					mnBeta[g] = ts.sumBeta[g] / (double)popsize;
					if (popsize > 1) {
						sdD0[g] = ts.ssqD0[g] / (double)popsize - mnD0[g] * mnD0[g];
						if (sdD0[g] > 0.0) sdD0[g] = sqrt(sdD0[g]); else sdD0[g] = 0.0;
						sdAlpha[g] = ts.ssqAlpha[g] / (double)popsize - mnAlpha[g] * mnAlpha[g];
						if (sdAlpha[g] > 0.0) sdAlpha[g] = sqrt(sdAlpha[g]); else sdAlpha[g] = 0.0;
						sdBeta[g] = ts.ssqBeta[g] / (double)popsize - mnBeta[g] * mnBeta[g];
						if (sdBeta[g] > 0.0) sdBeta[g] = sqrt(sdBeta[g]); else sdBeta[g] = 0.0;
					}
					else {
						sdD0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
					}
				}
			}
			if (emig.sexDep) {
				rangeOfs << "\t" << mnD0[0] << "\t" << sdD0[0];
				rangeOfs << "\t" << mnD0[1] << "\t" << sdD0[1];
				if (emig.densDep) {
					rangeOfs << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
					rangeOfs << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
					rangeOfs << "\t" << mnBeta[0] << "\t" << sdBeta[0];
					rangeOfs << "\t" << mnBeta[1] << "\t" << sdBeta[1];
				}
			}
			else { // sex-independent
				rangeOfs << "\t" << mnD0[0] << "\t" << sdD0[0];
				if (emig.densDep) {
					rangeOfs << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
					rangeOfs << "\t" << mnBeta[0] << "\t" << sdBeta[0];
				}
			}
		}

		if (trfr.indVar) {
			if (trfr.usesMovtProc) {
				// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				ngenes = 1;
			}
			else {
				if (trfr.sexDep) { // must be a sexual species
					ngenes = 2;
				}
				else {
					ngenes = 1;
				}
			}
			double mnDist1[2], mnDist2[2], mnProp1[2], mnStepL[2], mnRho[2];
			double sdDist1[2], sdDist2[2], sdProp1[2], sdStepL[2], sdRho[2];
			double mnDP[2], mnGB[2], mnAlphaDB[2], mnBetaDB[2];
			double sdDP[2], sdGB[2], sdAlphaDB[2], sdBetaDB[2];
			for (int g = 0; g < ngenes; g++) {
				mnDist1[g] = mnDist2[g] = mnProp1[g] = mnStepL[g] = mnRho[g] = 0.0;
				sdDist1[g] = sdDist2[g] = sdProp1[g] = sdStepL[g] = sdRho[g] = 0.0;
				mnDP[g] = mnGB[g] = mnAlphaDB[g] = mnBetaDB[g] = 0.0;
				sdDP[g] = sdGB[g] = sdAlphaDB[g] = sdBetaDB[g] = 0.0;
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ngenes == 2) popsize = ts.ninds[g];
				else popsize = ts.ninds[0] + ts.ninds[1];
				if (popsize > 0) {
					mnDist1[g] = ts.sumDist1[g] / (double)popsize;
					mnDist2[g] = ts.sumDist2[g] / (double)popsize;
					mnProp1[g] = ts.sumProp1[g] / (double)popsize;
					mnStepL[g] = ts.sumStepL[g] / (double)popsize;
					mnRho[g] = ts.sumRho[g] / (double)popsize;
					mnDP[g] = ts.sumDP[g] / (double)popsize;
					mnGB[g] = ts.sumGB[g] / (double)popsize;
					mnAlphaDB[g] = ts.sumAlphaDB[g] / (double)popsize;
					mnBetaDB[g] = ts.sumBetaDB[g] / (double)popsize;
					if (popsize > 1) {
						sdDist1[g] = ts.ssqDist1[g] / (double)popsize - mnDist1[g] * mnDist1[g];
						if (sdDist1[g] > 0.0) sdDist1[g] = sqrt(sdDist1[g]); else sdDist1[g] = 0.0;
						sdDist2[g] = ts.ssqDist2[g] / (double)popsize - mnDist2[g] * mnDist2[g];
						if (sdDist2[g] > 0.0) sdDist2[g] = sqrt(sdDist2[g]); else sdDist2[g] = 0.0;
						sdProp1[g] = ts.ssqProp1[g] / (double)popsize - mnProp1[g] * mnProp1[g];
						if (sdProp1[g] > 0.0) sdProp1[g] = sqrt(sdProp1[g]); else sdProp1[g] = 0.0;
						sdStepL[g] = ts.ssqStepL[g] / (double)popsize - mnStepL[g] * mnStepL[g];
						if (sdStepL[g] > 0.0) sdStepL[g] = sqrt(sdStepL[g]); else sdStepL[g] = 0.0;
						sdRho[g] = ts.ssqRho[g] / (double)popsize - mnRho[g] * mnRho[g];
						if (sdRho[g] > 0.0) sdRho[g] = sqrt(sdRho[g]); else sdRho[g] = 0.0;
						sdDP[g] = ts.ssqDP[g] / (double)popsize - mnDP[g] * mnDP[g];
						if (sdDP[g] > 0.0) sdDP[g] = sqrt(sdDP[g]); else sdDP[g] = 0.0;
						sdGB[g] = ts.ssqGB[g] / (double)popsize - mnGB[g] * mnGB[g];
						if (sdGB[g] > 0.0) sdGB[g] = sqrt(sdGB[g]); else sdGB[g] = 0.0;
						sdAlphaDB[g] = ts.ssqAlphaDB[g] / (double)popsize - mnAlphaDB[g] * mnAlphaDB[g];
						if (sdAlphaDB[g] > 0.0) sdAlphaDB[g] = sqrt(sdAlphaDB[g]); else sdAlphaDB[g] = 0.0;
						sdBetaDB[g] = ts.ssqBetaDB[g] / (double)popsize - mnBetaDB[g] * mnBetaDB[g];
						if (sdBetaDB[g] > 0.0) sdBetaDB[g] = sqrt(sdBetaDB[g]); else sdBetaDB[g] = 0.0;
					}
				}
			}
			if (trfr.usesMovtProc) {
				if (trfr.moveType == 1) {
					rangeOfs << "\t" << mnDP[0] << "\t" << sdDP[0];
					rangeOfs << "\t" << mnGB[0] << "\t" << sdGB[0];
					rangeOfs << "\t" << mnAlphaDB[0] << "\t" << sdAlphaDB[0];
					rangeOfs << "\t" << mnBetaDB[0] << "\t" << sdBetaDB[0];
				}
				if (trfr.moveType == 2) {
					rangeOfs << "\t" << mnStepL[0] << "\t" << sdStepL[0];
					rangeOfs << "\t" << mnRho[0] << "\t" << sdRho[0];
				}
			}
			else {
				if (trfr.sexDep) {
					rangeOfs << "\t" << mnDist1[0] << "\t" << sdDist1[0];
					rangeOfs << "\t" << mnDist1[1] << "\t" << sdDist1[1];
					if (trfr.twinKern)
					{
						rangeOfs << "\t" << mnDist2[0] << "\t" << sdDist2[0];
						rangeOfs << "\t" << mnDist2[1] << "\t" << sdDist2[1];
						rangeOfs << "\t" << mnProp1[0] << "\t" << sdProp1[0];
						rangeOfs << "\t" << mnProp1[1] << "\t" << sdProp1[1];
					}
				}
				else { // sex-independent
					rangeOfs << "\t" << mnDist1[0] << "\t" << sdDist1[0];
					if (trfr.twinKern)
					{
						rangeOfs << "\t" << mnDist2[0] << "\t" << sdDist2[0];
						rangeOfs << "\t" << mnProp1[0] << "\t" << sdProp1[0];
					}
				}
			}
		}

		if (sett.indVar) {
			if (sett.sexDep) { // must be a sexual species
				ngenes = 2;
			}
			else {
				if (dem.repType == 0) { // asexual reproduction
					ngenes = 1;
				}
				else { // sexual reproduction
					ngenes = 1;
				}
			}
			// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
			double mnS0[2], mnAlpha[2], mnBeta[2], sdS0[2], sdAlpha[2], sdBeta[2];
			for (int g = 0; g < ngenes; g++) {
				mnS0[g] = mnAlpha[g] = mnBeta[g] = sdS0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ngenes == 2) popsize = ts.ninds[g];
				else popsize = ts.ninds[0] + ts.ninds[1];
				if (popsize > 0) {
					mnS0[g] = ts.sumS0[g] / (double)popsize;
					mnAlpha[g] = ts.sumAlphaS[g] / (double)popsize;
					mnBeta[g] = ts.sumBetaS[g] / (double)popsize;
					if (popsize > 1) {
						sdS0[g] = ts.ssqS0[g] / (double)popsize - mnS0[g] * mnS0[g];
						if (sdS0[g] > 0.0) sdS0[g] = sqrt(sdS0[g]); else sdS0[g] = 0.0;
						sdAlpha[g] = ts.ssqAlphaS[g] / (double)popsize - mnAlpha[g] * mnAlpha[g];
						if (sdAlpha[g] > 0.0) sdAlpha[g] = sqrt(sdAlpha[g]); else sdAlpha[g] = 0.0;
						sdBeta[g] = ts.ssqBetaS[g] / (double)popsize - mnBeta[g] * mnBeta[g];
						if (sdBeta[g] > 0.0) sdBeta[g] = sqrt(sdBeta[g]); else sdBeta[g] = 0.0;
					}
					else {
						sdS0[g] = sdAlpha[g] = sdBeta[g] = 0.0;
					}
				}
			}
			if (sett.sexDep) {
				rangeOfs << "\t" << mnS0[0] << "\t" << sdS0[0];
				rangeOfs << "\t" << mnS0[1] << "\t" << sdS0[1];
				rangeOfs << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
				rangeOfs << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
				rangeOfs << "\t" << mnBeta[0] << "\t" << sdBeta[0];
				rangeOfs << "\t" << mnBeta[1] << "\t" << sdBeta[1];
			}
			else {
				rangeOfs << "\t" << mnS0[0] << "\t" << sdS0[0];
				rangeOfs << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
				rangeOfs << "\t" << mnBeta[0] << "\t" << sdBeta[0];
			}
		}

	}

	rangeOfs << endl;
}

void Community::closeOccSuitOfs(species_id sp) {
	if (outSuitOfs.at(sp).is_open()) outSuitOfs.at(sp).close();
	outSuitOfs.at(sp).clear();
}

// Open occupancy file, write header record and set up occupancy array
bool Community::outOccSuitHeaders(Species* pSpecies)
{
	simParams sim = paramsSim->getSim();
	landParams ppLand = pLandscape->getLandParams();
	species_id sp = pSpecies->getID();

	string name = paramsSim->getDir(2);
	if (sim.batchMode) {
		name += "Batch" + to_string(sim.batchNum) + "_";
		name += "Sim" + to_string(sim.simulation) + "_Land" + to_string(ppLand.landNum);
	}
	else
		name += "Sim" + to_string(sim.simulation);
	name += "_Species" + to_string(pSpecies->getID()) + "_Occupancy_Stats.txt";
	
	ofstream& suitOfs = outSuitOfs.at(sp);
	suitOfs.open(name.c_str());
	suitOfs << "Year\tMean_OccupSuit\tStd_error" << endl;

	return suitOfs.is_open();
}

void Community::outOccSuit(Species* pSpecies) {
	double sum, ss, mean, sd, se;
	simParams sim = paramsSim->getSim();
	int occInt = pSpecies->getOutOccInt();
	species_id sp = pSpecies->getID();
	ofstream& suitOfs = outSuitOfs.at(sp);

	for (int i = 0; i < (sim.years / occInt) + 1; i++) {

		sum = ss = 0.0;
		for (int rep = 0; rep < sim.reps; rep++) {
			int occ = occupancyMaps.at(sp)[i][rep];
			sum += occ;
			ss += occ * occ;
		}
		mean = sum / (double)sim.reps;
		sd = (ss - (sum * sum / (double)sim.reps)) / (double)(sim.reps - 1);
		if (sd > 0.0) sd = sqrt(sd);
		else sd = 0.0;
		se = sd / sqrt((double)(sim.reps));

		suitOfs << i * occInt << "\t" << mean << "\t" << se << endl;
	}
}

bool Community::closeOutTraitOfs(species_id sp) {
	if (outTraitsOfs.at(sp).is_open()) outTraitsOfs.at(sp).close();
	outTraitsOfs.at(sp).clear();
	return true;
}

// Open traits file and write header record
bool Community::outTraitsHeaders(species_id sp, Landscape* pLandscape, int landNr)
{
	landParams land = pLandscape->getLandParams();
	string name;
	Species* pSpecies = speciesMap.at(sp);
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();
	demogrParams dem = pSpecies->getDemogrParams();
	bool hasGenLoad = pSpecies->getNbGenLoadTraits() > 0;

	string DirOut = paramsSim->getDir(2);
	if (sim.batchMode) {
		name = DirOut
			+ "Batch" + to_string(sim.batchNum) + "_"
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr)
			+ "_Species" + to_string(sp)
			+ (land.usesPatches ? "_TraitsXpatch.txt" : "_TraitsXcell.txt");
	}
	else {
		name = DirOut + "Sim" + to_string(sim.simulation)
			+ "_Species" + to_string(sp)
			+ (land.usesPatches ? "_TraitsXpatch.txt" : "_TraitsXcell.txt");
	}

	ofstream& traitsOfs = outTraitsOfs.at(sp);
	traitsOfs.open(name.c_str());

	traitsOfs << "Rep\tYear\tRepSeason";
	if (land.usesPatches) traitsOfs << "\tPatchID";
	else traitsOfs << "\tx\ty";

	if (emig.indVar) {
		if (emig.sexDep) {
			if (emig.densDep) {
				traitsOfs << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
				traitsOfs << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
				traitsOfs << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
			}
			else {
				traitsOfs << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
			}
		}
		else {
			if (emig.densDep) {
				traitsOfs << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
				traitsOfs << "\tmeanBeta\tstdBeta";
			}
			else {
				traitsOfs << "\tmeanEP\tstdEP";
			}
		}
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 1) {
				traitsOfs << "\tmeanDP\tstdDP\tmeanGB\tstdGB";
				traitsOfs << "\tmeanAlphaDB\tstdAlphaDB\tmeanBetaDB\tstdBetaDB";
			}
			if (trfr.moveType == 2) {
				traitsOfs << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
			}
		}
		else {
			if (trfr.sexDep) {
				traitsOfs << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
				if (trfr.twinKern)
					traitsOfs << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
			}
			else {
				traitsOfs << "\tmean_distI\tstd_distI";
				if (trfr.twinKern)
					traitsOfs << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
			}
		}
	}
	if (sett.indVar) {
		if (sett.sexDep) {
			traitsOfs << "\tF_meanS0\tF_stdS0\tM_meanS0\tM_stdS0";
			traitsOfs << "\tF_meanAlphaS\tF_stdAlphaS\tM_meanAlphaS\tM_stdAlphaS";
			traitsOfs << "\tF_meanBetaS\tF_stdBetaS\tM_meanBetaS\tM_stdBetaS";
		}
		else {
			traitsOfs << "\tmeanS0\tstdS0";
			traitsOfs << "\tmeanAlphaS\tstdAlphaS";
			traitsOfs << "\tmeanBetaS\tstdBetaS";
		}
	}
	if (hasGenLoad) {
		if (dem.repType > 0) {
			traitsOfs << "\tF_meanGenFitness\tF_stdGenFitness\tM_meanGenFitness\tM_stdGenFitness";
		}
		else {
			traitsOfs << "\tmeanGenFitness\tstdGenFitness";
		}
	}

	traitsOfs << endl;

	return traitsOfs.is_open();
}

// Write records to traits file
/* NOTE: for summary traits by rows, which is permissible for a cell-based landscape
only, this function relies on the fact that populations are created in the same
sequence as patches, which is in ascending order of x nested within descending
order of y
*/
void Community::outTraits(species_id sp, int rep, int yr, int gen)
{
	landParams land = pLandscape->getLandParams();
	traitsums* ts = 0;
	traitsums popTraits;

	Species* pSpecies = speciesMap.at(sp);
	const bool mustOutputTraitRows = pSpecies->isTraitRowOutYear(yr);
	const bool mustOutputTraitCells = pSpecies->doesOutputTraitCell()
		&& yr % pSpecies->getOutTrCellInt() == 0;

	if (mustOutputTraitRows) {
		// create array of traits means, etc., one for each row
		ts = new traitsums[land.dimY];
	}

	if (pSpecies->isTraitCellOutYear(yr)
		|| mustOutputTraitRows) {

		// Generate output for each population in the community
		if (mustOutputTraitCells) {
			matrixPops.at(sp)->outputTraitPatchInfo(outTraitsOfs.at(sp), rep, yr, gen, land.usesPatches);
		}
		if (mustOutputTraitRows) {
			popTraits = matrixPops.at(sp)->outTraits(outTraitsOfs.at(sp), mustOutputTraitCells);
			int y = matrixPops.at(sp)->getPatch()->getCellLocn(0).y;
			for (int s = 0; s < gMaxNbSexes; s++) {
				ts[y].ninds[s] += popTraits.ninds[s];
				ts[y].sumD0[s] += popTraits.sumD0[s];
				ts[y].ssqD0[s] += popTraits.ssqD0[s];
				ts[y].sumAlpha[s] += popTraits.sumAlpha[s];
				ts[y].ssqAlpha[s] += popTraits.ssqAlpha[s];
				ts[y].sumBeta[s] += popTraits.sumBeta[s];
				ts[y].ssqBeta[s] += popTraits.ssqBeta[s];
				ts[y].sumDist1[s] += popTraits.sumDist1[s];
				ts[y].ssqDist1[s] += popTraits.ssqDist1[s];
				ts[y].sumDist2[s] += popTraits.sumDist2[s];
				ts[y].ssqDist2[s] += popTraits.ssqDist2[s];
				ts[y].sumProp1[s] += popTraits.sumProp1[s];
				ts[y].ssqProp1[s] += popTraits.ssqProp1[s];
				ts[y].sumStepL[s] += popTraits.sumStepL[s];
				ts[y].ssqStepL[s] += popTraits.ssqStepL[s];
				ts[y].sumRho[s] += popTraits.sumRho[s];
				ts[y].ssqRho[s] += popTraits.ssqRho[s];
				ts[y].sumS0[s] += popTraits.sumS0[s];
				ts[y].ssqS0[s] += popTraits.ssqS0[s];
				ts[y].sumAlphaS[s] += popTraits.sumAlphaS[s];
				ts[y].ssqAlphaS[s] += popTraits.ssqAlphaS[s];
				ts[y].sumBetaS[s] += popTraits.sumBetaS[s];
				ts[y].ssqBetaS[s] += popTraits.ssqBetaS[s];
				ts[y].sumGeneticFitness[s] += popTraits.sumGeneticFitness[s];
				ts[y].ssqGeneticFitness[s] += popTraits.ssqGeneticFitness[s];
			}
		}
		for (auto pop : allPopns.at(sp)) {
			if (mustOutputTraitCells) {
				pop->outputTraitPatchInfo(outTraitsOfs.at(sp), rep, yr, gen, land.usesPatches);
			}
			popTraits = pop->outTraits(outTraitsOfs.at(sp), mustOutputTraitCells);
			int y = pop->getPatch()->getCellLocn(0).y;
			if (mustOutputTraitRows) {
				for (int s = 0; s < gMaxNbSexes; s++) {
					ts[y].ninds[s] += popTraits.ninds[s];
					ts[y].sumD0[s] += popTraits.sumD0[s];
					ts[y].ssqD0[s] += popTraits.ssqD0[s];
					ts[y].sumAlpha[s] += popTraits.sumAlpha[s];
					ts[y].ssqAlpha[s] += popTraits.ssqAlpha[s];
					ts[y].sumBeta[s] += popTraits.sumBeta[s];
					ts[y].ssqBeta[s] += popTraits.ssqBeta[s];
					ts[y].sumDist1[s] += popTraits.sumDist1[s];
					ts[y].ssqDist1[s] += popTraits.ssqDist1[s];
					ts[y].sumDist2[s] += popTraits.sumDist2[s];
					ts[y].ssqDist2[s] += popTraits.ssqDist2[s];
					ts[y].sumProp1[s] += popTraits.sumProp1[s];
					ts[y].ssqProp1[s] += popTraits.ssqProp1[s];
					ts[y].sumStepL[s] += popTraits.sumStepL[s];
					ts[y].ssqStepL[s] += popTraits.ssqStepL[s];
					ts[y].sumRho[s] += popTraits.sumRho[s];
					ts[y].ssqRho[s] += popTraits.ssqRho[s];
					ts[y].sumS0[s] += popTraits.sumS0[s];
					ts[y].ssqS0[s] += popTraits.ssqS0[s];
					ts[y].sumAlphaS[s] += popTraits.sumAlphaS[s];
					ts[y].ssqAlphaS[s] += popTraits.ssqAlphaS[s];
					ts[y].sumBetaS[s] += popTraits.sumBetaS[s];
					ts[y].ssqBetaS[s] += popTraits.ssqBetaS[s];
					ts[y].sumGeneticFitness[s] += popTraits.sumGeneticFitness[s];
					ts[y].ssqGeneticFitness[s] += popTraits.ssqGeneticFitness[s];
				}
			}
		}

		if (allPopns.at(sp).size() > 0 && mustOutputTraitRows) {
			for (int y = 0; y < land.dimY; y++) {
				if ((ts[y].ninds[0] + ts[y].ninds[1]) > 0) {
					writeTraitsRows(sp, rep, yr, gen, y, ts[y]);
				}
			}
		}
	}
	if (ts != 0) {
		delete[] ts;
		ts = 0;
	}
}

// Write records to trait rows file
void Community::writeTraitsRows(species_id sp, int rep, int yr, int gen, int y,
	traitsums ts)
{
	Species* pSpecies = speciesMap.at(sp);
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	bool hasGenLoad = pSpecies->getNbGenLoadTraits() > 0;
	double mn, sd;

	ofstream& traitRowsOfs = outTraitsRows.at(sp);

	// calculate population size in case one phase is sex-dependent and the other is not
	// (in which case numbers of individuals are recorded by sex)
	int popsize = ts.ninds[0] + ts.ninds[1];
	traitRowsOfs << rep << "\t" << yr << "\t" << gen
		<< "\t" << y;
	if ((emig.indVar && emig.sexDep) || (trfr.indVar && trfr.sexDep))
		traitRowsOfs << "\t" << ts.ninds[0] << "\t" << ts.ninds[1];
	else
		traitRowsOfs << "\t" << popsize;

	if (emig.indVar) {
		if (emig.sexDep) {
			if (ts.ninds[0] > 0) mn = ts.sumD0[0] / (double)ts.ninds[0]; else mn = 0.0;
			if (ts.ninds[0] > 1) sd = ts.ssqD0[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			traitRowsOfs << "\t" << mn << "\t" << sd;
			if (ts.ninds[1] > 0) mn = ts.sumD0[1] / (double)ts.ninds[1]; else mn = 0.0;
			if (ts.ninds[1] > 1) sd = ts.ssqD0[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			traitRowsOfs << "\t" << mn << "\t" << sd;
			if (emig.densDep) {
				if (ts.ninds[0] > 0) mn = ts.sumAlpha[0] / (double)ts.ninds[0]; else mn = 0.0;
				if (ts.ninds[0] > 1) sd = ts.ssqAlpha[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (ts.ninds[1] > 0) mn = ts.sumAlpha[1] / (double)ts.ninds[1]; else mn = 0.0;
				if (ts.ninds[1] > 1) sd = ts.ssqAlpha[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (ts.ninds[0] > 0) mn = ts.sumBeta[0] / (double)ts.ninds[0]; else mn = 0.0;
				if (ts.ninds[0] > 1) sd = ts.ssqBeta[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (ts.ninds[1] > 0) mn = ts.sumBeta[1] / (double)ts.ninds[1]; else mn = 0.0;
				if (ts.ninds[1] > 1) sd = ts.ssqBeta[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
			}
		}
		else { // no sex dependence in emigration
			if (popsize > 0) mn = ts.sumD0[0] / (double)popsize; else mn = 0.0;
			if (popsize > 1) sd = ts.ssqD0[0] / (double)popsize - mn * mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			traitRowsOfs << "\t" << mn << "\t" << sd;
			if (emig.densDep) {
				if (popsize > 0) mn = ts.sumAlpha[0] / (double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqAlpha[0] / (double)popsize - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (popsize > 0) mn = ts.sumBeta[0] / (double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqBeta[0] / (double)popsize - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
			}
		}
	}

	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 2) { // CRW
				// NB - CURRENTLY CANNOT BE SEX-DEPENDENT...
				if (popsize > 0) mn = ts.sumStepL[0] / (double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqStepL[0] / (double)popsize - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (popsize > 0) mn = ts.sumRho[0] / (double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqRho[0] / (double)popsize - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
			}
		}
		else { // dispersal kernel
			if (trfr.sexDep) {
				if (ts.ninds[0] > 0) mn = ts.sumDist1[0] / (double)ts.ninds[0]; else mn = 0.0;
				if (ts.ninds[0] > 1) sd = ts.ssqDist1[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (ts.ninds[1] > 0) mn = ts.sumDist1[1] / (double)ts.ninds[1]; else mn = 0.0;
				if (ts.ninds[1] > 1) sd = ts.ssqDist1[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (trfr.twinKern)
				{
					if (ts.ninds[0] > 0) mn = ts.sumDist2[0] / (double)ts.ninds[0]; else mn = 0.0;
					if (ts.ninds[0] > 1) sd = ts.ssqDist2[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
					if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
					traitRowsOfs << "\t" << mn << "\t" << sd;
					if (ts.ninds[1] > 0) mn = ts.sumDist2[1] / (double)ts.ninds[1]; else mn = 0.0;
					if (ts.ninds[1] > 1) sd = ts.ssqDist2[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
					if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
					traitRowsOfs << "\t" << mn << "\t" << sd;
					if (ts.ninds[0] > 0) mn = ts.sumProp1[0] / (double)ts.ninds[0]; else mn = 0.0;
					if (ts.ninds[0] > 1) sd = ts.ssqProp1[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
					if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
					traitRowsOfs << "\t" << mn << "\t" << sd;
					if (ts.ninds[1] > 0) mn = ts.sumProp1[1] / (double)ts.ninds[1]; else mn = 0.0;
					if (ts.ninds[1] > 1) sd = ts.ssqProp1[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
					if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
					traitRowsOfs << "\t" << mn << "\t" << sd;
				}
			}
			else { // sex-independent
				if (popsize > 0) mn = ts.sumDist1[0] / (double)popsize; else mn = 0.0;
				if (popsize > 1) sd = ts.ssqDist1[0] / (double)popsize - mn * mn; else sd = 0.0;
				if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
				traitRowsOfs << "\t" << mn << "\t" << sd;
				if (trfr.twinKern)
				{
					if (popsize > 0) mn = ts.sumDist2[0] / (double)popsize; else mn = 0.0;
					if (popsize > 1) sd = ts.ssqDist2[0] / (double)popsize - mn * mn; else sd = 0.0;
					if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
					traitRowsOfs << "\t" << mn << "\t" << sd;
					if (popsize > 0) mn = ts.sumProp1[0] / (double)popsize; else mn = 0.0;
					if (popsize > 1) sd = ts.ssqProp1[0] / (double)popsize - mn * mn; else sd = 0.0;
					if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
					traitRowsOfs << "\t" << mn << "\t" << sd;
				}
			}
		}
	}

	if (sett.indVar) {
		if (popsize > 0) mn = ts.sumS0[0] / (double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqS0[0] / (double)popsize - mn * mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		traitRowsOfs << "\t" << mn << "\t" << sd;
		if (popsize > 0) mn = ts.sumAlphaS[0] / (double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqAlphaS[0] / (double)popsize - mn * mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		traitRowsOfs << "\t" << mn << "\t" << sd;
		if (popsize > 0) mn = ts.sumBetaS[0] / (double)popsize; else mn = 0.0;
		if (popsize > 1) sd = ts.ssqBetaS[0] / (double)popsize - mn * mn; else sd = 0.0;
		if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
		traitRowsOfs << "\t" << mn << "\t" << sd;
	}

	if (hasGenLoad) {
		if (gMaxNbSexes > 1) {
			if (ts.ninds[0] > 0) mn = ts.sumGeneticFitness[0] / (double)ts.ninds[0]; else mn = 0.0;
			if (ts.ninds[0] > 1) sd = ts.ssqGeneticFitness[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			traitRowsOfs << "\t" << mn << "\t" << sd;
			if (ts.ninds[1] > 0) mn = ts.sumGeneticFitness[1] / (double)ts.ninds[1]; else mn = 0.0;
			if (ts.ninds[1] > 1) sd = ts.ssqGeneticFitness[1] / (double)ts.ninds[1] - mn * mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			traitRowsOfs << "\t" << mn << "\t" << sd;
		}
		else {
			if (ts.ninds[0] > 0) mn = ts.sumGeneticFitness[0] / (double)ts.ninds[0]; else mn = 0.0;
			if (ts.ninds[0] > 1) sd = ts.ssqGeneticFitness[0] / (double)ts.ninds[0] - mn * mn; else sd = 0.0;
			if (sd > 0.0) sd = sqrt(sd); else sd = 0.0;
			traitRowsOfs << "\t" << mn << "\t" << sd;
		}
	}
	traitRowsOfs << endl;
}

bool Community::closeTraitRows(species_id sp) {
	if (outTraitsRows.at(sp).is_open()) outTraitsRows.at(sp).close();
	outTraitsRows.at(sp).clear();
	return true;
}

// Open trait rows file and write header record
bool Community::outTraitsRowsHeaders(species_id sp, int landNr) {

	string name;
	Species* pSpecies = speciesMap.at(sp);
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();
	bool hasGenLoad = pSpecies->getNbGenLoadTraits() > 0;

	string DirOut = paramsSim->getDir(2);
	if (sim.batchMode) {
		name = DirOut
			+ "Batch" + to_string(sim.batchNum) + "_"
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr)
			+ "_Species" + to_string(sp);
		+"_TraitsXrow.txt";
	}
	else {
		name = DirOut
			+ "Sim" + to_string(sim.simulation)
			+ "_Species" + to_string(sp);
		+"_TraitsXrow.txt";
	}

	ofstream& traitRowsOfs = outTraitsRows.at(sp);

	traitRowsOfs.open(name.c_str());

	traitRowsOfs << "Rep\tYear\tRepSeason\ty";
	if ((emig.indVar && emig.sexDep) || (trfr.indVar && trfr.sexDep))
		traitRowsOfs << "\tN_females\tN_males";
	else
		traitRowsOfs << "\tN";

	if (emig.indVar) {
		if (emig.sexDep) {
			if (emig.densDep) {
				traitRowsOfs << "\tF_meanD0\tF_stdD0\tM_meanD0\tM_stdD0";
				traitRowsOfs << "\tF_meanAlpha\tF_stdAlpha\tM_meanAlpha\tM_stdAlpha";
				traitRowsOfs << "\tF_meanBeta\tF_stdBeta\tM_meanBeta\tM_stdBeta";
			}
			else {
				traitRowsOfs << "\tF_meanEP\tF_stdEP\tM_meanEP\tM_stdEP";
			}
		}
		else {
			if (emig.densDep) {
				traitRowsOfs << "\tmeanD0\tstdD0\tmeanAlpha\tstdAlpha";
				traitRowsOfs << "\tmeanBeta\tstdBeta";
			}
			else {
				traitRowsOfs << "\tmeanEP\tstdEP";
			}
		}
	}
	if (trfr.indVar) {
		if (trfr.usesMovtProc) {
			if (trfr.moveType == 2) {
				traitRowsOfs << "\tmeanStepLength\tstdStepLength\tmeanRho\tstdRho";
			}
		}
		else { // dispersal kernel
			if (trfr.sexDep) {
				traitRowsOfs << "\tF_mean_distI\tF_std_distI\tM_mean_distI\tM_std_distI";
				if (trfr.twinKern)
					traitRowsOfs << "\tF_mean_distII\tF_std_distII\tM_mean_distII\tM_std_distII"
					<< "\tF_meanPfirstKernel\tF_stdPfirstKernel"
					<< "\tM_meanPfirstKernel\tM_stdPfirstKernel";
			}
			else {
				traitRowsOfs << "\tmean_distI\tstd_distI";
				if (trfr.twinKern)
					traitRowsOfs << "\tmean_distII\tstd_distII\tmeanPfirstKernel\tstdPfirstKernel";
			}
		}
	}

	if (sett.indVar) {
		traitRowsOfs << "\tmeanS0\tstdS0";
		traitRowsOfs << "\tmeanAlphaS\tstdAlphaS";
		traitRowsOfs << "\tmeanBetaS\tstdBetaS";
	}

	if (hasGenLoad) {
		if (gMaxNbSexes > 1) {
			traitRowsOfs << "\tF_meanProbViable\tF_stdProbViable\tM_meanProbViable\tM_stdProbViable";
		}
		else
			traitRowsOfs << "\tmeanProbViable\tstdProbViable";
	}

	traitRowsOfs << endl;

	return traitRowsOfs.is_open();

}

#if RS_RCPP && !R_CMD
Rcpp::IntegerMatrix Community::addYearToPopList(int rep, int yr) {  // TODO: define new simparams to control start and interval of output

	landParams ppLand = pLandscape->getLandParams();
	Rcpp::IntegerMatrix pop_map_year(ppLand.dimY, ppLand.dimX);
	Patch* pPatch = nullptr;
	Population* pPop = nullptr;

	for (int y = 0; y < ppLand.dimY; y++) {
		for (int x = 0; x < ppLand.dimX; x++) {
			Cell* pCell = pLandscape->findCell(x, y);
			if (pCell == nullptr) { // no-data cell
				pop_map_year(ppLand.dimY - 1 - y, x) = NA_INTEGER;
			}
			else {
				pPatch = pCell->getPatch();
				if (pPatch == nullptr) { // matrix cell
					pop_map_year(ppLand.dimY - 1 - y, x) = 0;
				}
				else {
					pPop = pPatch->getPop();
					if (pPop == nullptr) { // check if population exists
						pop_map_year(ppLand.dimY - 1 - y, x) = 0;
					}
					else {
						pop_map_year(ppLand.dimY - 1 - y, x) = pPop->getPopStats().nInds; // use indices like this because matrix gets transposed upon casting it into a raster on R-level
					}
				}
			}
		}
	}
	return pop_map_year;
}
#endif

bool Community::closeOutGenesOfs(species_id sp) {
	if (ofsGenes.at(sp).is_open()) {
		ofsGenes.at(sp).close();
		ofsGenes.at(sp).clear();
	}
	return true;
}

bool Community::openOutGenesFile(species_id sp, const bool& isDiploid, const int landNr, const int rep)
{
	string name;
	simParams sim = paramsSim->getSim();

	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + to_string(sim.batchNum)
			+ "_Sim" + to_string(sim.simulation) 
			+ "_Land" + to_string(landNr) 
			+ "_Rep" + to_string(rep)
			+ "_Species" + to_string(sp)
			+ "_geneValues.txt";
	}
	else {
		name = paramsSim->getDir(2)
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr)
			+ "_Rep" + to_string(rep)
			+ "_Species" + to_string(sp)
			+ "_geneValues.txt";
	}

	ofsGenes.at(sp).open(name.c_str());
	ofsGenes.at(sp) << "Year\tGeneration\tIndID\ttraitType\tlocusPosition"
		<< "\talleleValueA\tdomCoefA";
	if (isDiploid) ofsGenes.at(sp) << "\talleleValueB\tdomCoefB";
	ofsGenes.at(sp) << endl;

	return ofsGenes.at(sp).is_open();
}

void Community::outputGeneValues(species_id sp, const int& year, const int& gen) {

	if (!ofsGenes.at(sp).is_open())
		throw runtime_error("Could not open output gene values file.");

	const set<int> patchList = speciesMap.at(sp)->getSamplePatches();
	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(sp, patchId);
		if (patch == nullptr) {
			throw runtime_error("Sampled patch does not exist");
		}
		const auto pPop = patch->getPop();
		if (pPop != nullptr) {
			pPop->outputGeneValues(ofsGenes.at(sp), year, gen);
		}
	}
}

// ----------------------------------------------------------------------------------------
// Sample individuals from sample patches
// ----------------------------------------------------------------------------------------

void Community::sampleIndividuals(species_id sp) {

	Species* pSpecies = speciesMap.at(sp);
	const set<int> patchList = pSpecies->getSamplePatches();
	string nbIndsToSample = pSpecies->getNIndsToSample();
	const set<int> stagesToSampleFrom = pSpecies->getStagesToSample();

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(sp, patchId);
		if (patch == nullptr) {
			throw runtime_error("Can't sample individuals: patch" + to_string(patchId) + "doesn't exist.");
		}
		auto pPop = patch->getPop();
		if (pPop != nullptr) {
			pPop->sampleIndsWithoutReplacement(nbIndsToSample, stagesToSampleFrom);
		}
	}
}

// ----------------------------------------------------------------------------------------
// Open population level Fstat output file
// ----------------------------------------------------------------------------------------

bool Community::closeNeutralOutputOfs(species_id sp) {
	if (outWCFstatOfs.at(sp).is_open())
		outWCFstatOfs.at(sp).close();
	outWCFstatOfs.at(sp).clear();
	return true;
}

bool Community::openNeutralOutputFile(species_id sp, int landNr)
{
	string name;
	simParams sim = paramsSim->getSim();

	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + to_string(sim.batchNum)
			+ "_Sim" + to_string(sim.simulation) 
			+ "_Land" + to_string(landNr) 
			+ "_Species" + to_string(sp)
			+ "_neutralGenetics.txt";
	}
	else {
		paramsSim->getDir(2)
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr)
			+ "_Species" + to_string(sp)
			+ "_neutralGenetics.txt";
	}
	outWCFstatOfs.at(sp).open(name.c_str());
	outWCFstatOfs.at(sp) << "Rep\tYear\tRepSeason\tnExtantPatches\tnIndividuals\tFstWC\tFisWC\tFitWC\tFstWH\tmeanAllelePerLocus\tmeanAllelePerLocusPatches\tmeanFixedLoci\tmeanFixedLociPatches\tmeanObHeterozygosity";
	outWCFstatOfs.at(sp) << endl;

	return outWCFstatOfs.at(sp).is_open();
}

// ----------------------------------------------------------------------------------------
// open per locus WC fstat using MS approach, this will output MS calculated FIS, FIT, FST
// in general population neutral genetics output file
// ----------------------------------------------------------------------------------------

bool Community::closePerLocusFstFile(species_id sp) {
	if (outPerLocusFstat.at(sp).is_open())
		outPerLocusFstat.at(sp).close();
	outPerLocusFstat.at(sp).clear();
	return true;
}

bool Community::openPerLocusFstFile(Species* pSpecies, Landscape* pLandscape, const int landNr, const int rep)
{
	species_id sp = pSpecies->getID();
	set<int> patchList = pSpecies->getSamplePatches();
	if (patchList.size() == 0) {
		// list of patches is not known yet and may change every generation,
		// e.g. for randomOccupied sampling option
		// instead, header patch numbers range from 1 to nb of sampled patches
		for (int i = 0; i < pSpecies->getNbPatchesToSample(); i++) {
			patchList.emplace(i + 1);
		}
	}

	string name;
	simParams sim = paramsSim->getSim();

	name = paramsSim->getDir(2)
		+ (sim.batchMode ? "Batch" + to_string(sim.batchNum) + "_" : "")
		+ "Sim" + to_string(sim.simulation)
		+ "_Land" + to_string(landNr)
		+ "_Rep" + to_string(rep)
		+ "_Species" + to_string(pSpecies->getID()) +
		+"_perLocusNeutralGenetics.txt";

	outPerLocusFstat.at(sp).open(name.c_str());
	outPerLocusFstat.at(sp) << "Year\tRepSeason\tLocus\tFst\tFis\tFit\tHet";
	for (int patchId : patchList) {
		outPerLocusFstat.at(sp) << "\tpatch_" + to_string(patchId) + "_Het";
	}
	outPerLocusFstat.at(sp) << endl;

	return outPerLocusFstat.at(sp).is_open();
}

// ----------------------------------------------------------------------------------------
// open pairwise fst file
// ----------------------------------------------------------------------------------------

bool Community::closePairwiseFstFile(species_id sp) {
	if (outPairwiseFstOfs.at(sp).is_open())
		outPairwiseFstOfs.at(sp).close();
	outPairwiseFstOfs.at(sp).clear();
	return true;
}

bool Community::openPairwiseFstFile(Species* pSpecies, Landscape* pLandscape, const int landNr, const int rep) {

	species_id sp = pSpecies->getID();
	const set<int> patchList = pSpecies->getSamplePatches();
	string name;
	simParams sim = paramsSim->getSim();

	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + to_string(sim.batchNum)
			+ "_Sim" + to_string(sim.simulation) 
			+ "_Land" + to_string(landNr) 
			+ "_Rep" + to_string(rep)
			+ "_Species" + to_string(sp) +
			+ "_pairwisePatchNeutralGenetics.txt";
	}
	else {
		name = paramsSim->getDir(2)
			+ "Sim" + to_string(sim.simulation)
			+ "_Land" + to_string(landNr)
			+ "_Rep" + to_string(rep)
			+ "_Species" + to_string(sp) +
			+"_pairwisePatchNeutralGenetics.txt";
	}
	outPairwiseFstOfs.at(sp).open(name.c_str());
	outPairwiseFstOfs.at(sp) << "Year\tRepSeason\tpatchA\tpatchB\tFst";
	outPairwiseFstOfs.at(sp) << endl;

	return outPairwiseFstOfs.at(sp).is_open();
}

// ----------------------------------------------------------------------------------------
// Write population level FST results file
// ----------------------------------------------------------------------------------------

void Community::writeNeutralOutputFile(const species_id& sp, int rep, int yr, int gen, bool outWeirCockerham, bool outWeirHill) {

	outWCFstatOfs.at(sp) << rep << "\t" << yr << "\t" << gen << "\t";
	outWCFstatOfs.at(sp) << neutralStatsMaps.at(sp)->getNbPopulatedSampledPatches()
		<< "\t" << neutralStatsMaps.at(sp)->getTotalNbSampledInds() << "\t";

	if (outWeirCockerham) {
		outWCFstatOfs.at(sp) << neutralStatsMaps.at(sp)->getFstWC() << "\t"
			<< neutralStatsMaps.at(sp)->getFisWC() << "\t"
			<< neutralStatsMaps.at(sp)->getFitWC() << "\t";
	}
	else outWCFstatOfs.at(sp) << "N/A" << "\t" << "N/A" << "\t" << "N/A" << "\t";

	if (outWeirHill) outWCFstatOfs.at(sp) << neutralStatsMaps.at(sp)->getWeightedFst() << "\t";
	else outWCFstatOfs.at(sp) << "N/A" << "\t";

	outWCFstatOfs.at(sp) << neutralStatsMaps.at(sp)->getMeanNbAllPerLocus() << "\t"
		<< neutralStatsMaps.at(sp)->getMeanNbAllPerLocusPerPatch() << "\t"
		<< neutralStatsMaps.at(sp)->getTotalFixdAlleles() << "\t"
		<< neutralStatsMaps.at(sp)->getMeanFixdAllelesPerPatch() << "\t"
		<< neutralStatsMaps.at(sp)->getHo();

	outWCFstatOfs.at(sp) << endl;
}

// ----------------------------------------------------------------------------------------
// Write per locus FST results file
// ----------------------------------------------------------------------------------------

void Community::writePerLocusFstatFile(Species* pSpecies, const int yr, const int gen, const  int nAlleles, const int nLoci, set<int> const& patchList)
{
	const species_id sp = pSpecies->getID();
	const set<int> positions = pSpecies->getSpTrait(NEUTRAL)->getGenePositions();
	int thisLocus = 0;
	for (int position : positions) {

		outPerLocusFstat.at(sp) << yr << "\t"
			<< gen << "\t"
			<< position << "\t";
		outPerLocusFstat.at(sp) << neutralStatsMaps.at(sp)->getPerLocusFst(thisLocus) << "\t"
			<< neutralStatsMaps.at(sp)->getPerLocusFis(thisLocus) << "\t"
			<< neutralStatsMaps.at(sp)->getPerLocusFit(thisLocus) << "\t"
			<< neutralStatsMaps.at(sp)->getPerLocusHo(thisLocus);

		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(sp, patchId);
			const auto pPop = patch->getPop();
			int popSize = 0;
			int het = 0;
			if (pPop != nullptr) {
				popSize = pPop->sampleSize();
				if (popSize == 0) {
					outPerLocusFstat.at(sp) << "\t" << "N/A";
				}
				else {
					for (int a = 0; a < nAlleles; ++a) {
						het += static_cast<int>(pPop->getHeteroTally(thisLocus, a));
					}
					outPerLocusFstat.at(sp) << "\t"
						<< het / (2.0 * popSize);
				}
			}
			else {
				outPerLocusFstat.at(sp) << "\t" << "N/A";
			}
		}
		++thisLocus;
		outPerLocusFstat.at(sp) << endl;
	}
}


// ----------------------------------------------------------------------------------------
// Write pairwise FST results file
// ----------------------------------------------------------------------------------------
void Community::writePairwiseFstFile(Species* pSpecies, const int yr, const int gen, const  int nAlleles, const int nLoci, set<int> const& patchList) {

	const species_id sp = pSpecies->getID();
	// within patch fst (diagonal of matrix)
	int i = 0;
	for (int patchId : patchList) {
		outPairwiseFstOfs.at(sp) << yr << "\t" << gen << "\t";
		outPairwiseFstOfs.at(sp) << patchId << "\t" << patchId << "\t"
			<< neutralStatsMaps.at(sp)->getPairwiseFst(i, i)
			<< endl;
		++i;
	}

	// between patch fst
	i = 0;
	for (int patchIdA : patchList | std::views::take(patchList.size() - 1)) {
		int j = i + 1;
		for (int patchIdB : patchList | std::views::drop(j)) {
			outPairwiseFstOfs.at(sp) << yr << "\t" << gen << "\t";
			outPairwiseFstOfs.at(sp) << patchIdA << "\t" << patchIdB << "\t"
				<< neutralStatsMaps.at(sp)->getPairwiseFst(i, j)
				<< endl;
			++j;
		}
		++i;
	}
}


// ----------------------------------------------------------------------------------------
// Output and calculate neutral statistics
// ----------------------------------------------------------------------------------------
void Community::outNeutralGenetics(species_id sp, int rep, int yr, int gen) {

	Species* pSpecies = speciesMap.at(sp);
	const int maxNbNeutralAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const int nLoci = (int)pSpecies->getNPositionsForTrait(NEUTRAL);
	const set<int> patchList = pSpecies->getSamplePatches();
	int nInds = 0, nbPops = 0;

	for (int patchId : patchList) {
		const auto pPatch = pLandscape->findPatch(sp, patchId);
		if (pPatch == nullptr) {
			throw runtime_error("Sampled patch does not exist");
		}
		const auto pPop = pPatch->getPop();
		if (pPop != nullptr) { // empty patches do not contribute
			nInds += pPop->sampleSize();
			nbPops++;
		}
	}

	if (neutralStatsMaps.at(sp) == nullptr)
		neutralStatsMaps.at(sp) = make_unique<NeutralStatsManager>(patchList.size(), nLoci);

	neutralStatsMaps.at(sp)->updateAllNeutralTables(pSpecies, pLandscape, patchList);
	neutralStatsMaps.at(sp)->calculateHo(patchList, nInds, nLoci, pSpecies, pLandscape);
	neutralStatsMaps.at(sp)->calculatePerLocusHo(patchList, nInds, nLoci, pSpecies, pLandscape);
	neutralStatsMaps.at(sp)->calcAllelicDiversityMetrics(patchList, nInds, pSpecies, pLandscape);

	bool outWeirCockerham = pSpecies->doesOutputWeirCockerham();
	if (outWeirCockerham) {
		neutralStatsMaps.at(sp)->calculateFstatWC(patchList, nInds, nLoci, maxNbNeutralAlleles, pSpecies, pLandscape);
	}
	bool outWeirHill = pSpecies->doesOutputWeirHill();
	if (outWeirHill) {
		neutralStatsMaps.at(sp)->calcPairwiseWeightedFst(patchList, nInds, nLoci, pSpecies, pLandscape);
	}

	writeNeutralOutputFile(sp, rep, yr, gen, outWeirCockerham, outWeirHill);

	if (outWeirCockerham) {
		writePerLocusFstatFile(pSpecies, yr, gen, maxNbNeutralAlleles, nLoci, patchList);
	}
	if (outWeirHill) {
		writePairwiseFstFile(pSpecies, yr, gen, maxNbNeutralAlleles, nLoci, patchList);
	}
}

//---------------------------------------------------------------------------
//For outputs and population visualisations pre-reproduction
void Community::traitAndOccOutput(int rep, int yr, int gen) {

	for (auto& sp : activeSpecies) {
		Species* pSpecies = speciesMap.at(sp);
		// trait outputs and visualisation
		if (pSpecies->isTraitCellOutYear(yr)
			|| pSpecies->isTraitRowOutYear(yr))
			outTraits(sp, rep, yr, gen);
		if (gen == 0 && pSpecies->isOccOutputYear(yr))
			updateOccupancy(sp, yr, rep);
	}
}


bool Community::openOutputFiles(bool hasMultipleReplicates, const int landNum) {

	bool filesOK = true;

	// Open output files
	for (auto& [sp, pSpecies] : speciesMap) {
		if (pSpecies->doesOutputRange()) { // open Range file
			if (!outRangeHeaders(sp, landNum)) {
				filesOK = false;
			}
		}
		if (pSpecies->doesOutputOccup() && hasMultipleReplicates) {
			if (!outOccSuitHeaders(pSpecies))
				filesOK = false;
			if (!pLandscape->outOccupancyHeaders(pSpecies))
					filesOK = false;
		}

		if (pSpecies->doesOutputPop())
			if (!outPopHeaders(pSpecies))
				filesOK = false;

		if (pSpecies->doesOutputTraitCell())
			if (!outTraitsHeaders(sp, pLandscape, landNum)) {
				filesOK = false;
			}
		if (pSpecies->doesOutputTraitRows())
			if (!outTraitsRowsHeaders(sp, landNum)) {
				filesOK = false;
			}
		if (pSpecies->doesOutputWeirCockerham()
			|| pSpecies->doesOutputWeirHill()) { // open neutral genetics file
			if (!openNeutralOutputFile(sp, landNum)) {
				filesOK = false;
			}
		}
	}

	if (!filesOK) {
		closeGlobalOutputFiles(hasMultipleReplicates);
	}

	return filesOK;
}

void Community::closeGlobalOutputFiles(bool hasMultipleReplicates) {

	for (auto& [sp, pSpecies] : speciesMap) {

		if (pSpecies->doesOutputTraitRows())
			closeTraitRows(sp);
		if (pSpecies->doesOutputTraitCell()) closeOutTraitOfs(sp);

		if (pSpecies->doesOutputRange()) closeRangeOfs(sp);
		if (pSpecies->doesOutputPop()) closePopOfs(sp);

		if (pSpecies->doesOutputWeirCockerham()
			|| pSpecies->doesOutputWeirCockerham())
			closeNeutralOutputOfs(sp);

		if (hasMultipleReplicates && pSpecies->doesOutputOccup()) {
			closeOccSuitOfs(sp);
			pLandscape->closeOccupancyOfs(sp);
		}
	}
}

void Community::closeYearlyOutputFiles() {

	for (auto& [sp, pSpecies] : speciesMap) {
		if (pSpecies->doesOutputInds()) closeOutIndsOfs(sp);
		if (pSpecies->doesOutputGeneValues()) closeOutGenesOfs(sp);
		if (pSpecies->doesOutputWeirCockerham()) closePerLocusFstFile(sp);
		if (pSpecies->doesOutputWeirHill()) closePairwiseFstFile(sp);
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
