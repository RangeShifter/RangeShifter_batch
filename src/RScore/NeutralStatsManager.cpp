#include "NeutralStatsManager.h"
#include "Population.h"

/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2026 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Roslyn Henry, Théo Pannetier, Jette Wolff, Damaris Zurell
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
 * File Created by Roslyn Henry March 2023. Code adapted from NEMO (https://nemo2.sourceforge.io/)
 --------------------------------------------------------------------------*/

 // ----------------------------------------------------------------------------------------
 // Constructor
 // ----------------------------------------------------------------------------------------
NeutralStatsManager::NeutralStatsManager(const int& nbSampledPatches, const int nLoci) {
	this->pairwiseFstMatrix = PatchMatrix(nbSampledPatches, nbSampledPatches);
	commNeutralCountTables.reserve(nLoci); //don't have to be pointers, not shared or moved

	perLocusFst = perLocusFis = perLocusFit = vector<double>(nLoci, 0.0);
}

// ----------------------------------------------------------------------------------------
// Populate population and community-level NEUTRAL count tables
// Update allele occurrence and heterozygosity counts, and allele frequencies
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::updateAllNeutralTables(Species* pSpecies, Landscape* pLandscape, set<int> const& patchList) {

	const int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	const int maxNbNeutralAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;

	// Create / Update community-level NEUTRAL counts table
	if (!commNeutralCountTables.empty()) {
		resetCommNeutralTables();
	}
	else { // populate the tables with default values
		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			NeutralCountsTable newNeutralTbl = NeutralCountsTable(maxNbNeutralAlleles);
			commNeutralCountTables.push_back(newNeutralTbl);
		}
	}

	int nbSampledInds = 0;
	int patchAlleleCount;
	// Update counts for each population
	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
		const auto pPop = patch->getPop();
		if (pPop != nullptr) {
			// Update this population's NEUTRAL counts tables
			pPop->updatePopNeutralTables();
			nbSampledInds += pPop->sampleSize();
		}
		// Add population-level counts to community-level counts 
		for (int thisLocus = 0; thisLocus < nLoci; thisLocus++) {
			for (int allele = 0; allele < maxNbNeutralAlleles; allele++) {

				if (pPop != nullptr) {
					patchAlleleCount = pPop->getAlleleTally(thisLocus, allele);
				}
				else {
					patchAlleleCount = 0;
				}
				commNeutralCountTables[thisLocus].incrementTallyBy(patchAlleleCount, allele);
			}
		}
	}

	// Update community-level frequencies
	std::for_each(commNeutralCountTables.begin(),
		commNeutralCountTables.end(),
		[&](NeutralCountsTable& v) -> void {
			v.setFrequencies(nbSampledInds * ploidy);
		});
}

// ----------------------------------------------------------------------------------------
// Reset allele tables in NeutralTable structs
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::resetCommNeutralTables() {
	for (auto& entry : commNeutralCountTables) {
		entry.reset();
	}
}

// ----------------------------------------------------------------------------------------
//	Calculate allelic diversity metrics
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calcAllelicDiversityMetrics(set<int> const& patchList, const int nInds, Species* pSpecies, Landscape* pLandscape)
{
	int i, j;
	const int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	const int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;
	unsigned int nbPops = 0;
	int nbAllelesInPatch = 0;
	double meanAllelicDivInPatch = 0;
	bool alleleExistsInPop = 0;

	vector<vector<bool>> alleleExistsInCommTable(nLoci);
	for (i = 0; i < nLoci; ++i) {
		alleleExistsInCommTable[i] = vector<bool>(nAlleles, false);
	}

	// Compute mean nb alleles per locus per patch
	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
		const auto pPop = patch->getPop();
		if (pPop != nullptr) {
			if (pPop->sampleSize() > 0) {
				nbPops++;
				nbAllelesInPatch = 0;
				for (i = 0; i < nLoci; ++i)
					for (j = 0; j < nAlleles; ++j) {
						alleleExistsInPop = pPop->getAlleleTally(i, j) != 0;
						nbAllelesInPatch += alleleExistsInPop;
						alleleExistsInCommTable[i][j] = alleleExistsInCommTable[i][j] || alleleExistsInPop;
					}
				// add mean nb of alleles per locus for Patch k to the pop mean
				meanAllelicDivInPatch += static_cast<double>(nbAllelesInPatch) / nLoci;
			}
		}
	}
	meanNbAllelesPerLocusPerPatch = nbPops > 0 ? meanAllelicDivInPatch / nbPops : 0;

	// Compute mean nb alleles per locus
	meanNbAllelesPerLocus = 0;
	for (i = 0; i < nLoci; ++i)
		for (j = 0; j < nAlleles; ++j)
			meanNbAllelesPerLocus += alleleExistsInCommTable[i][j];
	meanNbAllelesPerLocus /= nLoci;
	
	// Compute number of fixed loci per patch
	// mean number of loci that are fixed at pop level per pop
	meanNbFixedLociPerPatch = 0;
	if (nbPops > 0) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
			const auto pPop = patch->getPop();
			if (pPop != nullptr) {
				for (i = 0; i < nLoci; ++i)
					for (j = 0; j < nAlleles; ++j)
						meanNbFixedLociPerPatch += pPop->getAlleleFrequency(i, j) == 1;
			}
		}
		meanNbFixedLociPerPatch /= nbPops;
	}

	// Compute number of fixed loci
	meanFixedLoci = 0;
	for (i = 0; i < nLoci; ++i)
		for (j = 0; j < nAlleles; ++j)
			meanFixedLoci += commNeutralCountTables[i].getFrequency(j) == 1;
}

// ----------------------------------------------------------------------------------------
// Calculate Ho per Nei and Chesser
// Average (observed) heterozygosity per individual
// Sum (nb of heterozygote loci) across individuals / nb individuals / nb loci
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	int nbHetero = 0;

	if (nbInds != 0 && pSpecies->isDiploid()) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
			const auto pPop = patch->getPop();
			if (pPop != nullptr) nbHetero += pPop->countHeterozygoteLoci();
		}
		ho = static_cast<double>(nbHetero) / (nbInds * nbrLoci);
	}
	else ho = 0.0;
}

// ----------------------------------------------------------------------------------------
// Calculate Hs per Nei and Chesser
// Average expected population-level heterozygosity per locus per population
// currently not used but may be useful
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHs(set<int> const& patchList, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	double hs = 0;
	int nPatches = 0;

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
		const auto pPop = patch->getPop();
		if (pPop->sampleSize() > 0) {
			nPatches++;
			hs += pPop->computeHs();
		}
	}
	hs = (nPatches != 0 ? hs / (nbrLoci * nPatches) : 0.0);
}

// ----------------------------------------------------------------------------------------
// Calculate Ht per Nei and Chesser
// Average expected community-level heterozygosity per locus
// Currently not used but may be useful
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateHt(Species* pSpecies, Landscape* pLandscape, const int nLoci, const int nAlleles) {

	double ht = 0;
	int nPatches = 0;
	vector<double>locihet(nLoci, 1);
	double freq;

	for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
		for (int allele = 0; allele < nAlleles; ++allele) {
			freq = commNeutralCountTables[thisLocus].getFrequency(allele);
			freq *= freq; //squared frequencies
			locihet[thisLocus] -= freq;  //1 - sum of p2 = expected heterozygosity
		}
		ht += locihet[thisLocus];
	}
	ht = ht / nLoci;
}

// ----------------------------------------------------------------------------------------
// Calculate Ho per locus as per Nei and Chesser
// Observed proportion of heterozygote individuals for each locus
// Sum (nb of heterozygote individuals) / nb individuals for each locus
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculatePerLocusHo(set<int> const& patchList, const int nbInds, const int nbrLoci, Species* pSpecies, Landscape* pLandscape) {

	vector<int> nbHeterosInComm(nbrLoci, 0);
	vector<int> nbHeterosInPop(nbrLoci);

	if (pSpecies->isDiploid()) {
		for (int patchId : patchList) {
			const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
			const auto pPop = patch->getPop();
			if (pPop != nullptr) {
				if (pPop->sampleSize() > 0) {
					nbHeterosInPop = pPop->countNbHeterozygotesEachLocus();
					// Add counts to community total
					transform(nbHeterosInComm.begin(), nbHeterosInComm.end(), nbHeterosInPop.begin(),
						nbHeterosInComm.begin(), plus<int>());
				}
			}
		}
	}

	perLocusHo = vector<double>(nbrLoci, 0);
	if (nbInds != 0) {
		for (int i = 0; i < nbHeterosInComm.size(); i++) {
			perLocusHo[i] = static_cast<double>(nbHeterosInComm[i]) / nbInds;
		}
	}
}

// ----------------------------------------------------------------------------------------
// Fstat Weir & Cockerham
// ----------------------------------------------------------------------------------------
void NeutralStatsManager::calculateGlobalFst(set<int> const& patchList, const int nbSampledIndsInComm, const int nLoci, const int nAlleles, 
	Species* pSpecies, Landscape* pLandscape, bool isPairwise) {

	double inverseNtotal;
	double sumWeights = 0;
	double nBar, nC, inverseNbar;
	unsigned int nbPops = 0;
	const int totalSampleSize = nbSampledIndsInComm; // r * n_bar
	const double ploidy = pSpecies->isDiploid() ? 2.0 : 1.0;

	// Reset per-locus vectors between generations
	perLocusFst = perLocusFis = perLocusFit = vector<double>(nLoci, 0.0);

	for (int patchId : patchList) {
		const auto patch = pLandscape->findPatch(pSpecies->getID(), patchId);
		const auto pPop = patch->getPop();
		if (pPop != nullptr) {
			int popSampleSize = pPop->sampleSize(); // n_i
			if (popSampleSize > 0) {
				nbPops++;
				sumWeights += static_cast<double>(popSampleSize * popSampleSize) / totalSampleSize; // sum(n_i^2/rn_bar)
			}
		}
	}

	nbExtantPops = nbPops; // r
	totalNbSampledInds = nbSampledIndsInComm;

	if (nbPops > 1) {

		// Calculate F stats
		nBar = static_cast<double>(totalSampleSize) / nbPops; // average sample size, cannot be less than 1
		nC = (totalSampleSize - sumWeights) / (nbPops - 1);
		double nBarMinusOne = (nBar == 1.0) ? 1.0 : nBar - 1.0; // avoid / 0 if exactly 1 ind per pop
		inverseNbar = 1.0 / nBarMinusOne;
		inverseNtotal = 1.0 / totalSampleSize;

		double var, intermediateTerm;
		double s2, pBar, hBar;
		int pairwiseAlleleCount;
		double s2Denom = (nbPops - 1) * nBar;
		double rTerm = static_cast<double>(nbPops - 1) / nbPops;
		double hBarFactor = (2 * nBar - 1) / (4 * nBar);

		double numFst = 0.0, numFis = 0.0, numFit = 0.0;
		double denomFst = 0.0, denomFis = 0.0, denomFit = 0.0;

		for (int l = 0; l < nLoci; ++l) {

			// Sums of a_u, b_u, c_u for all alleles u at locus l
			double a_l = 0, b_l = 0, c_l = 0;

			for (int u = 0; u < nAlleles; ++u) {

				pBar = s2 = hBar = 0;
				pairwiseAlleleCount = 0;

				//if global wc approach use this
				if (!isPairwise) {
					pBar = commNeutralCountTables[l].getFrequency(u);
				}
				//else calculate total frequencies just in pair of patches for pairwise 
				else {
					for (int patchId : patchList) {
						const auto pPatch = pLandscape->findPatch(pSpecies->getID(), patchId);
						const auto pPop = pPatch->getPop();
						double patchLocusAlleleTally = pPop->getAlleleTally(l, u);
						pairwiseAlleleCount += patchLocusAlleleTally;
					}
					const double denomAlleleCount = totalSampleSize * ploidy;
					pBar = pairwiseAlleleCount / denomAlleleCount;

				}
			
				for (int patchId : patchList) {
					const auto pPatch = pLandscape->findPatch(pSpecies->getID(), patchId);
					const auto pPop = pPatch->getPop();
					if (pPop != 0) {
						var = pPop->getAlleleFrequency(l, u) - pBar;
						var *= var;
						s2 += var * pPop->sampleSize();
						hBar += pPop->getHeteroTally(l, u); // n_i * h_i

					}
				} //end for pop

				s2 /= s2Denom;
				hBar /= static_cast<float>(totalSampleSize); // / (r * n_bar)

				intermediateTerm = pBar * (1 - pBar) - rTerm * s2;
				a_l += s2 - inverseNbar * (intermediateTerm - 0.25 * hBar);
				b_l += intermediateTerm - hBarFactor * hBar;
				c_l += hBar;

			} // end for allele 

			a_l *= nBar / nC;
			b_l *= nBar / nBarMinusOne;
			c_l *= 0.5;

			perLocusFst[l] = a_l / (a_l + b_l + c_l);
			perLocusFis[l] = (b_l + c_l) == 0.0 ? 0.0 : b_l / (b_l + c_l);
			perLocusFit[l] = a_l + b_l / (a_l + b_l + c_l);

			numFst += a_l;
			numFis += b_l;
			numFit += a_l + b_l;
			denomFst += a_l + b_l + c_l;
			denomFis += b_l + c_l;
			
		} // end for locus

		denomFit = denomFst; // same quantity

		fst = numFst / denomFst; // theta hat in eq. 1 in WC 1984
		fis = (denomFis == 0.0) ? 0.0 : numFis / denomFis; // f hat
		fit = numFit / denomFit; // F hat
	}
	else { // zero or one sampled pops, cannot compute F-stats
		fst = 0.0;
		fis = 0.0;
		fit = 0.0;
	}
}

////////New pairwise function ///////////


void NeutralStatsManager::calculatePairwiseFst(set<int> const& patchList, const int nLoci, const int nAlleles, Species* pSpecies, Landscape* pLandscape) {

	// Needs to be in vector to iterate over, copy preserves order
	vector<int> patchVect;
	copy(patchList.begin(), patchList.end(), std::back_inserter(patchVect));

	const int nPatches = static_cast<int>(patchVect.size());

	// Initialise 
	pairwiseFstMatrix = PatchMatrix(nPatches, nPatches);

	// Reset table
	pairwiseFstMatrix.setAll(0.0); // or nanf("NULL")?


	for (int i = 0; i < nPatches - 1; ++i)
	{
		const auto patchA = pLandscape->findPatch(pSpecies->getID(), patchVect[i]);
		const auto pPopA = patchA->getPop();

		// Skip if patch A has no individuals
		if (pPopA->sampleSize() == 0)
			continue;

		for (int j = i + 1; j < nPatches; ++j)
		{
			const auto patchB = pLandscape->findPatch(pSpecies->getID(), patchVect[j]);
			const auto pPopB = patchB->getPop();

			// Skip if patch B has no individuals
			if (pPopB->sampleSize() == 0)
				continue;

			// Build pair-of-patches list
			set<int> pairofPatchesList;
			pairofPatchesList.insert(patchVect[i]);
			pairofPatchesList.insert(patchVect[j]);

			// Total individuals in the pair
			int nbSampledIndsInPair = pPopA->sampleSize() + pPopB->sampleSize();

			//NB this overwrites global fst variable so be careful!!
			calculateGlobalFst(pairofPatchesList, nbSampledIndsInPair,
				 nLoci, nAlleles,pSpecies, pLandscape, true);

			pairwiseFstMatrix.set(i, j, fst);
			// else remain 0
		}
	}
	
}
