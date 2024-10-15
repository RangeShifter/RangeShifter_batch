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

#include "SubCommunity.h"
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

SubCommunity::SubCommunity(Patch* pPch, int subCommId) {
	subCommNum = subCommId;
	pPatch = pPch;
	// record the new sub-community no. in the patch
	pPatch->setSubComm(this);
	initialSubComm = false;
	occupancy = 0;
}

SubCommunity::~SubCommunity() {
	pPatch->setSubComm(0);
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		delete popns[i];
	}
	popns.clear();
	if (occupancy != 0) delete[] occupancy;
}

int SubCommunity::getNum(void) { return subCommNum; }

Patch* SubCommunity::getPatch(void) { return pPatch; }

locn SubCommunity::getLocn(void) {
	locn loc = pPatch->getCellLocn(0);
	return loc;
}

void SubCommunity::setInitial(bool b) { initialSubComm = b; }

void SubCommunity::initialise(Landscape* pLandscape, Species* pSpecies)
{
	int ncells;
	landParams ppLand = pLandscape->getLandParams();
	initParams init = paramsInit->getInit();

	// determine size of initial population
	int nInds = 0;
	if (subCommNum == 0 // matrix patch
		|| !initialSubComm)   		// not in initial region or distribution
		nInds = 0;
	else {
		float k = pPatch->getK();
		if (k > 0.0) { // patch is currently suitable for this species
			switch (init.initDens) {
			case 0: // at carrying capacity
				nInds = (int)k;
				break;
			case 1: // at half carrying capacity
				nInds = (int)(k / 2.0);
				break;
			case 2: // specified no. per cell or density
				ncells = pPatch->getNCells();
				if (ppLand.patchModel) {
					nInds = (int)(init.indsHa * (float)(ncells * ppLand.resol * ppLand.resol) / 10000.0);
				}
				else {
					nInds = init.indsCell * ncells;
				}
				break;
			}
		}
		else nInds = 0;
	}

	// create new population only if it is non-zero or the matrix popn
	if (subCommNum == 0 || nInds > 0) {
		newPopn(pLandscape, pSpecies, pPatch, nInds);
	}

}

// initialise a specified individual
void SubCommunity::initialInd(Landscape* pLandscape, Species* pSpecies,
	Patch* pPatch, Cell* pCell, int ix)
{
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	short stg, age, repInt;
	Individual* pInd;
	float probmale;

	// create new population if not already in existence
	int npopns = (int)popns.size();
	if (npopns < 1) {
		newPopn(pLandscape, pSpecies, pPatch, 0);
	}

	// create new individual
	initInd iind = paramsInit->getInitInd(ix);
	if (dem.stageStruct) {
		stg = iind.stage; 
		age = iind.age;   
		repInt = sstruct.repInterval;
	}
	else {
		age = stg = 1;
		repInt = 0;
	}
	probmale = (dem.repType != 0 && iind.sex == 1) ? 1.0 : 0.0;
	
	pInd = new Individual(pCell, pPatch, stg, age, repInt, probmale, trfr.usesMovtProc, trfr.moveType);

	// add new individual to the population
	popns[0]->recruit(pInd);

	if (pSpecies->getNTraits() > 0)
	{
		// individual variation - set up genetics
		pInd->setUpGenes(pSpecies, pLandscape->getLandData().resol);
	}

}

// Create a new population, and return its address
Population* SubCommunity::newPopn(Landscape* pLandscape, Species* pSpecies,
	Patch* pPatch, int nInds)
{
	landParams land = pLandscape->getLandParams();
	int npopns = (int)popns.size();
	popns.push_back(new Population(pSpecies, pPatch, nInds, land.resol));
	return popns[npopns];
}

void SubCommunity::resetPossSettlers(void) {
	if (subCommNum == 0) return; // not applicable in the matrix
	pPatch->resetPossSettlers();
}

// Extirpate all populations according to
// option 0 - random local extinction probability
// option 1 - local extinction probability gradient
// NB only applied for cell-based model
void SubCommunity::localExtinction(int option) {
	double pExtinct = 0.0;
	if (option == 0) {
		envStochParams env = paramsStoch->getStoch();
		if (env.localExt) pExtinct = env.locExtProb;
	}
	else {
		envGradParams grad = paramsGrad->getGradient();
		Cell* pCell = pPatch->getRandomCell(); // get only cell in the patch
		// extinction prob is complement of cell gradient value plus any non-zero prob at the optimum
		pExtinct = 1.0 - pCell->getEnvVal() + grad.extProbOpt;
		if (pExtinct > 1.0) pExtinct = 1.0;
	}
	if (pRandom->Bernoulli(pExtinct)) {
		int npops = (int)popns.size();
		for (int i = 0; i < npops; i++) { // all populations
			popns[i]->extirpate();
		}
	}
}

// Action in event of patch becoming unsuitable owing to landscape change
void SubCommunity::patchChange(void) {
	if (subCommNum == 0) return; // no reproduction in the matrix
	Species* pSpecies;
	float localK = 0.0;
	int npops = (int)popns.size();
	// THE FOLLOWING MAY BE MORE EFFICIENT WHILST THERE IS ONLY ONE SPECIES ...
	if (npops < 1) return;
	localK = pPatch->getK();
	if (localK <= 0.0) { // patch in dynamic landscape has become unsuitable
		for (int i = 0; i < npops; i++) { // all populations
			pSpecies = popns[i]->getSpecies();
			demogrParams dem = pSpecies->getDemogrParams();
			if (dem.stageStruct) {
				stageParams sstruct = pSpecies->getStageParams();
				if (sstruct.disperseOnLoss) popns[i]->allEmigrate();
				else popns[i]->extirpate();
			}
			else { // non-stage-structured species is destroyed
				popns[i]->extirpate();
			}
		}
	}
}

void SubCommunity::reproduction(int resol, float epsGlobal, short rasterType, bool patchModel)
{
	if (subCommNum == 0) return; // no reproduction in the matrix
	float localK, envval;
	Cell* pCell;
	envGradParams grad = paramsGrad->getGradient();
	envStochParams env = paramsStoch->getStoch();

	int npops = popns.size();
	if (npops < 1) return;

	localK = pPatch->getK();
	if (localK > 0.0) {
		if (patchModel) {
			envval = 1.0; // environmental gradient is currently not applied for patch-based model
		}
		else { // cell-based model
			if (grad.gradient && grad.gradType == 2)
			{ // gradient in fecundity
				Cell* pCell = pPatch->getRandomCell(); // locate the only cell in the patch
				envval = pCell->getEnvVal();
			}
			else envval = 1.0;
		}
		if (env.stoch && !env.inK) { // stochasticity in fecundity
			if (env.local) {
				if (!patchModel) { // only permitted for cell-based model
					pCell = pPatch->getRandomCell();
					if (pCell != 0) envval += pCell->getEps();
				}
			}
			else { // global stochasticity
				envval += epsGlobal;
			}
		}
		for (int i = 0; i < npops; i++) { // all populations
			popns[i]->reproduction(localK, envval, resol);
			popns[i]->fledge();
		}
	}
}

void SubCommunity::emigration(void)
{
	if (subCommNum == 0) return; // no emigration from the matrix
	int npops = static_cast<int>(popns.size());
	if (npops < 1) return;
	float localK = pPatch->getK();
	// NOTE that even if K is zero, it could have been >0 in previous time-step, and there
	// might be emigrants if there is non-juvenile emigration
	for (int i = 0; i < npops; i++) { // all populations
		popns[i]->emigration(localK);
	}
}

// Remove emigrants from their natal patch and add to patch 0 (matrix)
void SubCommunity::initiateDispersal(SubCommunity* matrix) {

	if (subCommNum == 0) return; // no dispersal initiation in the matrix
	
	popStats pop;
	disperser disp;

	int npops = popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		pop = popns[i]->getStats();
		for (int j = 0; j < pop.nInds; j++) {
			disp = popns[i]->extractDisperser(j);
			if (disp.yes) { // disperser - has already been removed from natal population
				// add to matrix population
				matrix->recruit(disp.pInd, pop.pSpecies);
			}
		}
		// remove pointers to emigrants
		popns[i]->clean();
	}
}

//---------------------------------------------------------------------------

// Remove emigrants from patch 0 (matrix) and transfer to sub-community
// in which their destination co-ordinates fall
// This function is executed for the matrix patch only

void SubCommunity::completeDispersal(Landscape* pLandscape, bool connect)
{
	int popsize;
	disperser settler;
	Species* pSpecies;
	Population* pPop;
	Patch* pPrevPatch;
	Patch* pNewPatch;
	Cell* pPrevCell;
	SubCommunity* pSubComm;

	int npops = popns.size();
	for (int i = 0; i < npops; i++) { // all populations

		pSpecies = popns[i]->getSpecies();
		popsize = popns[i]->getNInds();

		for (int j = 0; j < popsize; j++) {

			bool settled;
			settler = popns[i]->extractSettler(j);
			settled = settler.yes;
			if (settled) {
			// settler - has already been removed from matrix population
			// find new patch
				pNewPatch = settler.pCell->getPatch();
				// find population within the patch (if there is one)
				pPop = pNewPatch->getPopn(pSpecies);
				if (pPop == 0) { // settler is the first in a previously uninhabited patch
					// create a new population in the corresponding sub-community
					pSubComm = (SubCommunity*)pNewPatch->getSubComm();
					pPop = pSubComm->newPopn(pLandscape, pSpecies, pNewPatch, 0);
				}
				pPop->recruit(settler.pInd);
				if (connect) { // increment connectivity totals
					int newpatch = pNewPatch->getSeqNum();
					pPrevCell = settler.pInd->getLocn(0); // previous cell
					Patch* patch = pPrevCell->getPatch();
					if (patch != 0) {
						pPrevPatch = patch;
						int prevpatch = pPrevPatch->getSeqNum();
						pLandscape->incrConnectMatrix(prevpatch, newpatch);
					}
				}
			}
		}
		// remove pointers in the matrix popn to settlers
		popns[i]->clean();
	}

}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


