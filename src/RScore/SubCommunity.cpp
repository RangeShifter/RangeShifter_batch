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
}

SubCommunity::~SubCommunity() {
	pPatch->setSubComm(0);
	int npops = (int)popns.size();
	for (int i = 0; i < npops; i++) { // all populations
		delete popns[i];
	}
	popns.clear();
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
					pPop = new Population(pSpecies, pNewPatch, 0, pLandscape->getLandParams().resol);
					
					// TMP = below should be for new subcomm, not this one
					popns.push_back(pPop);
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


