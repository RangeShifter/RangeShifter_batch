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

#include "Population.h"
#include "Patch.h"

#include <algorithm>
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
Population::Population(Species* pSp, Patch* pPch, int ninds, int resol)
{
	int n, nindivs, age = 0, minage, maxage, nAges = 0;
	int cumtotal = 0;
	float probmale;
	double ageprob, ageprobsum;
	std::vector <double> ageProb; // for quasi-equilibrium initial age distribution
	Cell* pCell;

	if (ninds > 0) {
		inds.reserve(ninds);
		newborns.reserve(ninds);
	}

	pSpecies = pSp;
	pPatch = pPch;

	// record the new population in the patch
	pPatch->setPop(this);

	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	initParams init = pSpecies->getInitParams();

	// determine no. of stages and sexes of species to initialise
	if (dem.stageStruct) {
		nStages = sstruct.nStages;
	}
	else // non-structured population has 2 stages, but user only ever sees stage 1
		nStages = 2;
	if (dem.repType == 0) { nSexes = 1; probmale = 0.0; }
	else { nSexes = 2; probmale = dem.propMales; }

	// set up population sub-totals
	for (int stg = 0; stg < gMaxNbStages; stg++) {
		for (int sex = 0; sex < gMaxNbSexes; sex++) {
			nInds[stg][sex] = 0;
		}
	}

	// set up local copy of minimum age table
	short minAge[gMaxNbStages][gMaxNbSexes];
	for (int stg = 0; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			if (dem.stageStruct) {
				int sexIx = dem.repType == 1 ? 0 : sex;
				minAge[stg][sex] = pSpecies->getMinAge(stg, sexIx);
			}
			else { // non-structured population
				minAge[stg][sex] = 0;
			}
		}
	}

	fecInitdEffects = fecRecdEffects = fecResDepEffects
		= devInitdEffects = devRecdEffects = devResDepEffects
		= survInitdEffects = survRecdEffects = survResDepEffects
		= vector(nStages, 0.0);

	// individuals of new population must be >= stage 1
	for (int stg = 1; stg < nStages; stg++) {
		if (dem.stageStruct) { // allocate to stages according to initialisation conditions
			// final stage is treated separately to ensure that correct total
			// no. of individuals is created
			if (stg == nStages - 1) {
				n = ninds - cumtotal;
			}
			else {
				n = (int)(ninds * pSpecies->getProp(stg) + 0.5);
				cumtotal += n;
			}
		}
		else { // non-structured - all individuals go into stage 1
			n = ninds;
		}
			// establish initial age distribution
		minage = maxage = stg;
		if (dem.stageStruct) {
			// allow for stage-dependent minimum ages (use whichever sex is greater)
			if (minAge[stg][0] > 0 && minage < minAge[stg][0]) minage = minAge[stg][0];
			if (nSexes == 2 && minAge[stg][1] > 0 && minage < minAge[stg][1]) minage = minAge[stg][1];
			// allow for specified age distribution
			if (init.initAge != 0) { // not lowest age
				if (stg == nStages - 1) maxage = sstruct.maxAge; // final stage
				else { // all other stages - use female max age, as sex of individuals is not predetermined
					maxage = minAge[stg + 1][0] - 1;
				}
				if (maxage < minage) maxage = minage;
				nAges = maxage - minage + 1;
				if (init.initAge == 2) { // quasi-equilibrium distribution
					double psurv = (double)pSpecies->getSurv(stg, 0); // use female survival for the stage
					ageProb.clear();
					ageprobsum = 0.0;
					ageprob = 1.0;
					for (int i = 0; i < nAges; i++) {
						ageProb.push_back(ageprob); ageprobsum += ageprob; ageprob *= psurv;
					}
					for (int i = 0; i < nAges; i++) {
						ageProb[i] /= ageprobsum;
						if (i > 0) ageProb[i] += ageProb[i - 1]; // to give cumulative probability
					}
				}
			}
		}
	// create individuals
		int sex;
		nindivs = (int)inds.size();
		for (int i = 0; i < n; i++) {
			pCell = pPatch->getRandomCell();
			if (dem.stageStruct) {
				switch (init.initAge) {
				case 0: // lowest possible age
					age = minage;
					break;
				case 1: // randomised
					if (maxage > minage) age = pRandom->IRandom(minage, maxage);
					else age = minage;
					break;
				case 2: // quasi-equilibrium
					if (nAges > 1) {
						double rrr = pRandom->Random();
						int ageclass = 0;
						while (rrr > ageProb[ageclass]) ageclass++;
						age = minage + ageclass;
					}
					else age = minage;
					break;
				}
			}
			else age = stg;

			Individual* pInd = new Individual(pSpecies, pCell, pPatch, stg, age, sstruct.repInterval,
				probmale, trfr.usesMovtProc, trfr.moveType);
			inds.push_back(pInd);

			if (pSpecies->getNTraits() > 0) {
				// individual variation - set up genetics
				pInd->setUpGenes(resol);
			}
			nInds[stg][pInd->getSex()]++;
		}
	}
}

Population::~Population() {
	int ninds = static_cast<int>(inds.size());
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != nullptr) delete inds[i];
	}
	inds.clear();
	int njuvs = static_cast<int>(newborns.size());
	for (int i = 0; i < njuvs; i++) {
		if (newborns[i] != nullptr) delete newborns[i];
	}
	newborns.clear();
}

traitsums Population::getIndTraitsSums() {
	int g;
	traitsums ts = traitsums();

	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	for (auto& ind : inds) {

		int sex = ind->getSex();
		ts.ninds[sex] += 1;

		// emigration traits
		emigTraits e = ind->getIndEmigTraits();
		ts.sumD0[sex] += e.d0;    
		ts.ssqD0[sex] += e.d0 * e.d0;
		ts.sumAlpha[sex] += e.alpha; 
		ts.ssqAlpha[sex] += e.alpha * e.alpha;
		ts.sumBeta[sex] += e.beta;  
		ts.ssqBeta[sex] += e.beta * e.beta;

		// transfer traits
		if (trfr.usesMovtProc) {

			switch (trfr.moveType) {

			case 1: { // SMS
				trfrSMSTraits sms = ind->getIndSMSTraits();
				ts.sumDP[sex] += sms.dp; 
				ts.ssqDP[sex] += sms.dp * sms.dp;
				ts.sumGB[sex] += sms.gb;
				ts.ssqGB[sex] += sms.gb * sms.gb;
				ts.sumAlphaDB[sex] += sms.alphaDB;
				ts.ssqAlphaDB[sex] += sms.alphaDB * sms.alphaDB;
				ts.sumBetaDB[sex] += sms.betaDB; 
				ts.ssqBetaDB[sex] += sms.betaDB * sms.betaDB;
				break;
			}
			case 2: {
				trfrCRWTraits c = ind->getIndCRWTraits();
				ts.sumStepL[sex] += c.stepLength;
				ts.ssqStepL[sex] += c.stepLength * c.stepLength;
				ts.sumRho[sex] += c.rho;       
				ts.ssqRho[sex] += c.rho * c.rho;
				break;
			}
			default:
				throw runtime_error("usesMoveProcess is ON but moveType is neither 1 (SMS) or 2 (CRW).");
				break;
			}
		}
		else {
			trfrKernelParams k = ind->getIndKernTraits();
			ts.sumDist1[sex] += k.meanDist1; 
			ts.ssqDist1[sex] += k.meanDist1 * k.meanDist1;
			ts.sumDist2[sex] += k.meanDist2;
			ts.ssqDist2[sex] += k.meanDist2 * k.meanDist2;
			ts.sumProp1[sex] += k.probKern1; 
			ts.ssqProp1[sex] += k.probKern1 * k.probKern1;
		}
		// settlement traits
		settleTraits s = ind->getIndSettTraits();
		ts.sumS0[sex] += s.s0;     
		ts.ssqS0[sex] += s.s0 * s.s0;
		ts.sumAlphaS[sex] += s.alpha; 
		ts.ssqAlphaS[sex] += s.alpha * s.alpha;
		ts.sumBetaS[sex] += s.beta;   
		ts.ssqBetaS[sex] += s.beta * s.beta;

		double fitness = ind->getGeneticFitness();
		ts.sumGeneticFitness[sex] += fitness;
		ts.ssqGeneticFitness[sex] += fitness * fitness;
	}
	return ts;
}

// ----------------------------------------------------------------------------------------
// reset allele table
// ----------------------------------------------------------------------------------------
void Population::resetPopNeutralTables() {
	for (auto& entry : popNeutralCountTables) {
		entry.reset();
	}
}

// ----------------------------------------------------------------------------------------
// Populate population-level NEUTRAL count tables
// Update allele occurrence and heterozygosity counts, and allele frequencies
// ----------------------------------------------------------------------------------------
void Population::updatePopNeutralTables() {

	const int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	const int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	const auto& positions = pSpecies->getSpTrait(NEUTRAL)->getGenePositions();
	const int ploidy = pSpecies->isDiploid() ? 2 : 1;

	// Create /reset empty tables
	if (popNeutralCountTables.size() != 0)
		resetPopNeutralTables();
	else {
		popNeutralCountTables.reserve(nLoci);

		for (int l = 0; l < nLoci; l++) {
			popNeutralCountTables.push_back(NeutralCountsTable(nAlleles));
		}
	}

	// Fill tallies for each locus
	for (Individual* individual : sampledInds) {

		const auto trait = individual->getTrait(NEUTRAL);
		int whichLocus = 0;
		for (auto position : positions) {

			int alleleOnChromA = (int)trait->getAlleleValueAtLocus(0, position);
			popNeutralCountTables[whichLocus].incrementTally(alleleOnChromA);

			if (ploidy == 2) { // second allele and heterozygosity
				int alleleOnChromB = (int)trait->getAlleleValueAtLocus(1, position);
				popNeutralCountTables[whichLocus].incrementTally(alleleOnChromB);

				bool isHetero = alleleOnChromA != alleleOnChromB;
				if (isHetero) {
					popNeutralCountTables[whichLocus].incrementHeteroTally(alleleOnChromA);
					popNeutralCountTables[whichLocus].incrementHeteroTally(alleleOnChromB);
				}
			}
			whichLocus++;
		}
	}

	// Fill frequencies
	if (sampledInds.size() > 0) {
		std::for_each(
			popNeutralCountTables.begin(),
			popNeutralCountTables.end(),
			[&](NeutralCountsTable& thisLocus) -> void {
				thisLocus.setFrequencies(static_cast<int>(sampledInds.size()) * ploidy);
			});
	}
}

double Population::getAlleleFrequency(int thisLocus, int whichAllele) {
	return popNeutralCountTables[thisLocus].getFrequency(whichAllele);
}

int Population::getAlleleTally(int thisLocus, int whichAllele) {
	return popNeutralCountTables[thisLocus].getTally(whichAllele);
}

int Population::getHeteroTally(int thisLocus, int whichAllele) {
	return popNeutralCountTables[thisLocus].getHeteroTally(whichAllele);
}

// ----------------------------------------------------------------------------------------
// Count number of heterozygotes loci in sampled individuals
// ----------------------------------------------------------------------------------------
int Population::countHeterozygoteLoci() {
	int nbHetero = 0;
	if (pSpecies->isDiploid()) {
		for (Individual* ind : sampledInds) {
			const NeutralTrait* trait = (NeutralTrait*)(ind->getTrait(NEUTRAL));
			nbHetero += trait->countHeterozygoteLoci();
		}
	}
	return nbHetero;
}

// ----------------------------------------------------------------------------------------
// Count number of heterozygotes among sampled individuals for each locus
// ----------------------------------------------------------------------------------------
vector<int> Population::countNbHeterozygotesEachLocus() {
	const auto& positions = pSpecies->getSpTrait(NEUTRAL)->getGenePositions();
	vector<int> hetero(positions.size(), 0);

	if (pSpecies->isDiploid()) {
		for (Individual* ind : sampledInds) {
			const NeutralTrait* trait = (NeutralTrait*)ind->getTrait(NEUTRAL);
			int counter = 0;
			for (auto position : positions) {
				hetero[counter] += trait->isHeterozygoteAtLocus(position);
				counter++;
			}
		}
	}
	return hetero;
}

// ----------------------------------------------------------------------------------------
//	Compute the expected heterozygosity for population
// ----------------------------------------------------------------------------------------
double Population::computeHs() {
	int nLoci = pSpecies->getNPositionsForTrait(NEUTRAL);
	int nAlleles = pSpecies->getSpTrait(NEUTRAL)->getNbNeutralAlleles();
	double hs = 0;
	double freq;
	vector<double> locihet(nLoci, 1);

	if (sampledInds.size() > 0) {
		for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
			for (int allele = 0; allele < nAlleles; ++allele) {
				freq = getAlleleFrequency(thisLocus, allele);
				freq *= freq; //squared frequencies (expected _homozygosity)
				locihet[thisLocus] -= freq; // 1 - sum of p2 = expected heterozygosity
			}
			hs += locihet[thisLocus];
		}
	}
	return hs;
}

popStats Population::getStats()
{
	popStats p = popStats();
	int ninds;
	float fec;
	bool breeders[2] = { false, false };
	demogrParams dem = pSpecies->getDemogrParams();
	p.pSpecies = pSpecies;
	p.pPatch = pPatch;
	p.speciesID = pSpecies->getID();
	p.nInds = (int)inds.size();
	p.nNonJuvs = p.nAdults = 0;
	p.breeding = false;
	for (int stg = 1; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			ninds = nInds[stg][sex];
			p.nNonJuvs += ninds;
			if (ninds > 0) {
				if (pSpecies->stageStructured()) {
					if (dem.repType == 2) fec = pSpecies->getFec(stg, sex);
					else fec = pSpecies->getFec(stg, 0);
					if (fec > 0.0) { breeders[sex] = true; p.nAdults += ninds; }
				}
				else breeders[sex] = true;
			}
		}
	}
	// is there a breeding population present?
	if (nSexes == 1) {
		p.breeding = breeders[0];
	}
	else {
		if (breeders[0] && breeders[1]) p.breeding = true;
	}
	return p;
}

Species* Population::getSpecies() { return pSpecies; }

int Population::getNbInds() const {
	return inds.size();
}

int Population::getNbInds(int stg) const {
	int t = 0;
	if (stg < 0 || stg >= nStages) throw runtime_error("Attempt to get nb individuals for stage " + to_string(stg) + ", no such stage.");
	for (int sex = 0; sex < nSexes; sex++) {
		t += nInds[stg][sex];
	}
	return t;
}

int Population::getNbInds(int stg, int sex) const {
	if (stg < 0 || stg >= nStages) throw runtime_error("Attempt to get nb individuals for stage " + to_string(stg) + ", no such stage.");
	return nInds[stg][sex];
}

//---------------------------------------------------------------------------
// Extirpate all populations according to
// local extinction probability gradient
// NB only applied for cell-based model
void Population::applyLocalExtGrad() {

	envGradParams grad = pSpecies->getEnvGradient();
	if (grad.gradType != 3) return;

	// extinction prob is complement of cell gradient value plus any non-zero prob at the optimum
	double pExtinct = min(1.0, 1.0 - pPatch->getGradVal() + grad.extProbOpt);
	
	if (pRandom->Bernoulli(pExtinct)) {
		extirpate();
	}
}

// Remove all Individuals
void Population::extirpate() {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != nullptr) delete inds[i];
	}
	inds.clear();
	int njuvs = (int)newborns.size();
	for (int i = 0; i < njuvs; i++) {
		if (newborns[i] != nullptr) delete newborns[i];
	}
	newborns.clear();
	for (int sex = 0; sex < nSexes; sex++) {
		for (int stg = 0; stg < nStages; stg++) {
			nInds[stg][sex] = 0;
		}
	}
}

//---------------------------------------------------------------------------
// Produce new individuals and hold them in the newborns vector
void Population::reproduction(const float localK, const int resol)
{
	if (inds.size() == 0) return;

	int stage, sex, nbOffspring, nmales, nfemales;
	Cell* pCell;
	double exptdNbOffspring;
	bool skipBreeding;

	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	int nsexes = dem.repType == 0 ? 1 : 2;
	// if sexual, must also calculate male fecundity
	// to decide which males reproduce (those with fec > 0)

	// Base fecundity
	float fec[gMaxNbStages][gMaxNbSexes];
	for (int stg = 1; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				fec[stg][sex] = pSpecies->getFec(stg, dem.repType == 0 ? 0 : sex);
			}
			else { // non-structured population
				fec[stg][sex] = dem.lambda;
			}
		}
	}

	// Environmental effects
	for (int stg = 1; stg < sstruct.nStages; stg++) {

		if (fec[stg][0] <= 0.0) continue;

		// Environmental gradient + stochasticty
		fec[stg][0] *= pPatch->getGradVal();

		// Limits if applicable
		if (env.usesStoch && !env.inK) {
			float minFec = pSpecies->getMinMax(0);
			float maxFec = pSpecies->getMinMax(1);
			fec[stg][0] = min(maxFec, max(minFec, fec[stg][0]));
		}
	}

	if (dem.stageStruct) {

		// Received contributions from interspecific interactions
		for (int stg = 1; stg < nStages; stg++) {
			fec[stg][0] += fecRecdEffects[stg];
		}

		// Intraspecific density-dependence
		for (int stg = 1; stg < nStages; stg++) {

			if (fec[stg][0] <= 0.0 || !sstruct.fecDens) continue;

			float densDepEffect = 0.0;
			if (sstruct.fecStageDens) { // stage-specific density dependence
				// NOTE: matrix entries represent effect of ROW on COLUMN 
				// AND males precede females
				float weight = 0.0;
				for (int effstg = 0; effstg < nStages; effstg++) {
					for (int effsex = 0; effsex < nSexes; effsex++) {
						if (dem.repType == 2) {
							int sexColIndex = effsex == 0 ? 1 : 0;
							weight = pSpecies->getDDwtFec(2 * stg + 1, 2 * effstg + sexColIndex);
						}
						else weight = pSpecies->getDDwtFec(stg, effstg);

						densDepEffect += (float)nInds[effstg][effsex] * weight;
					}
				}
			}
			else densDepEffect = static_cast<float>(getNbInds());

			if (localK > 0.0) fec[stg][0] *= exp(-densDepEffect / localK);
		}

		// Other interaction effects
		for (int stg = 1; stg < nStages; stg++) {
			// Contribution from resource-dependent interactions
			fec[stg][0] *= exp(fecResDepEffects[stg]);

			// Contributions from initiated interspecific interactions
			fec[stg][0] += fecInitdEffects[stg];
		}
	}
	else { // non-structured - set fecundity for adult females only
		if (localK > 0.0) {
			if (dem.repType != 0) // sexual model
				fec[1][0] *= 2.0; // cf manual
			fec[1][0] /= (1.0f + fabs(dem.lambda - 1.0f) * pow((static_cast<float>(inds.size()) / localK), dem.bc));
		}
	}

	double propBreed;
	Individual* newJuv;
	Individual* pFather = nullptr;
	std::vector <Individual*> fathers;

	switch (dem.repType) {

	case 0: // asexual model
		for (auto& pInd : inds) {

			// Skip juveniles, males and dispersing or dead females
			if (!pInd->isBreedingFem()) continue;
			
			if (dem.stageStruct) // check for reproduction cooldown (fallow)
				if (!pInd->breedsThisSeason(sstruct))
					continue;

			exptdNbOffspring = fec[pInd->getStats().stage][0];
			nbOffspring = exptdNbOffspring <= 0.0 ? 0 : pRandom->Poisson(exptdNbOffspring);
			pCell = pPatch->getRandomCell();

			for (int j = 0; j < nbOffspring; j++) {

				newJuv = new Individual(pSpecies, pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);
				if (pSpecies->getNTraits() > 0)
					newJuv->inheritTraits(pInd, resol);
				if (!newJuv->isViable())
					delete newJuv; 
				else {
					newborns.push_back(newJuv);
					nInds[0][0]++;
				}
			}
		}
		break;

	case 1: // simple sexual model
	case 2: // complex sexual model
		// count breeding females and males
		// add breeding males to list of potential fathers
		nfemales = nmales = 0;
		for (auto& pInd : inds) {
			indStats ind = pInd->getStats();
			if (ind.sex == 0 && fec[ind.stage][0] > 0.0) nfemales++;
			if (ind.sex == 1 && fec[ind.stage][1] > 0.0) {
				fathers.push_back(pInd);
				nmales++;
			}
		}
		if (nfemales == 0 || nmales == 0)
			break; // no reproduction

		if (dem.repType == 2) { // complex sexual model
			propBreed = min(1.0, (2.0 * dem.harem * nmales)
				/ (nfemales + dem.harem * nmales));
		}

		for (auto& pInd : inds) {
			stage = pInd->getStats().stage;
			if (!pInd->isBreedingFem() || fec[stage][0] <= 0.0) continue;

			if (dem.stageStruct) // check for reproduction cooldown (fallow)
				if (!pInd->breedsThisSeason(sstruct))
					continue;

			// NOTE: FOR COMPLEX SEXUAL MODEL, NO. OF FEMALES *ACTUALLY* BREEDING DOES NOT
			// NECESSARILY EQUAL THE EXPECTED NO. FROM EQN. 7 IN THE MANUAL...
			if (dem.repType == 2)
				if (!pRandom->Bernoulli(propBreed)) continue;

			exptdNbOffspring = fec[stage][0];
			nbOffspring = exptdNbOffspring > 0.0 ? pRandom->Poisson(exptdNbOffspring) : 0;
			if (nbOffspring <= 0) continue;

			// Select father randomly among eligible individuals
			int whichFather = nmales > 1 ? pRandom->IRandom(0, nmales - 1) : 0;
			pFather = fathers[whichFather];
			pCell = pPatch->getRandomCell();

			for (int j = 0; j < nbOffspring; j++) {

				newJuv = new Individual(pSpecies, pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);
				if (pSpecies->getNTraits() > 0)
					newJuv->inheritTraits(pInd, pFather, resol);
				if (!newJuv->isViable())
					delete newJuv;
				else {
					newborns.push_back(newJuv);
					sex = newJuv->getSex();
					nInds[0][sex]++;
				}
			}
		}
		break;
	} // end of switch (dem.repType)
}

// Following reproduction of ALL species, add juveniles to the population prior to dispersal
void Population::fledge()
{
	demogrParams dem = pSpecies->getDemogrParams();

	if (!dem.stageStruct) { // all adults die
		int ninds = inds.size();
		for (int i = 0; i < ninds; i++)
			delete inds[i];
		inds.clear();
		for (int sex = 0; sex < nSexes; sex++)
			nInds[1][sex] = 0; // set count of adults to zero
	}
	inds = std::move(newborns);
	newborns.clear();
}

Individual* Population::sampleInd() const {
	int index = pRandom->IRandom(0, static_cast<int>(inds.size() - 1));
	return inds[index];
}

void Population::sampleIndsWithoutReplacement(string strNbToSample, const set<int>& sampleStages) {

	if (sampledInds.size() > 0) {
		sampledInds.clear();
	}
	auto rng = pRandom->getRNG();
	vector<Individual*> stagedInds;

	// Stage individuals in eligible stages
	for (int stage : sampleStages) {
		vector<Individual*> toAdd = getIndividualsInStage(stage);
		stagedInds.insert(stagedInds.begin(), toAdd.begin(), toAdd.end());
	}

	if (strNbToSample == "all") {
		// Sample all individuals in selected stages
		sampledInds = stagedInds;
	}
	else { // random
		int nbToSample = stoi(strNbToSample);
		if (stagedInds.size() <= nbToSample) {
			// Sample all individuals in selected stages
			sampledInds = stagedInds;
		}
		else {
			// Sample n individuals across selected stages
			sample(stagedInds.begin(), stagedInds.end(), std::back_inserter(sampledInds), nbToSample, rng);
		}
	}
}

int Population::sampleSize() const {
	return static_cast<int>(sampledInds.size());
}

vector<Individual*> Population::getIndividualsInStage(int stage) {
	vector<Individual*> indsInStage;
	for (auto ind : inds) {
		if (ind->getStats().stage == stage)
			indsInStage.push_back(ind);
	}
	return indsInStage;
}

// Determine which individuals will disperse
void Population::emigration(float localK)
{
	int nsexes;
	double disp, pbDisp, NK;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	emigTraits eparams;
	indStats ind;

	// to avoid division by zero, assume carrying capacity is at least one individual
	// localK can be zero if there is a moving gradient or stochasticity in K
	if (localK < 1.0) localK = 1.0;
	NK = static_cast<float>(getNbInds()) / localK;

	// set up local copy of emigration probability table
	// used when there is no individual variability
	// NB - IT IS DOUBTFUL THIS CONTRIBUTES ANY SUBSTANTIAL TIME SAVING
	if (dem.repType == 0) nsexes = 1; 
	else nsexes = 2;
	double pbEmig[gMaxNbStages][gMaxNbSexes];

	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {

			int stgId = emig.stgDep ? stg : 0;
			int sexId = emig.sexDep ? sex : 0;

			if (emig.indVar) pbEmig[stg][sex] = 0.0;
			else {
				if (emig.densDep) {
					eparams = pSpecies->getSpEmigTraits(stgId, sexId);
					pbEmig[stg][sex] = eparams.d0 / (1.0 + exp(-(NK - eparams.beta) * eparams.alpha));
				}
				else {
					pbEmig[stg][sex] = pSpecies->getSpEmigD0(stgId, sexId);
				}
			}
		}
	}

	int ninds = inds.size();
	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();
		int stgId = emig.stgDep ? ind.stage : 0;
		int sexId = emig.sexDep ? ind.sex : 0;
		if (ind.status == initial) {
			if (emig.indVar) { // individual variability in emigration
				if (dem.stageStruct && ind.stage != emig.emigStage) {
					// emigration may not occur
					pbDisp = 0.0;
				}
				else { // non-structured or individual is in emigration stage
					eparams = inds[i]->getIndEmigTraits();
					if (emig.densDep) { // density-dependent
						NK = (float)getNbInds() / localK;
						pbDisp = eparams.d0 / (1.0 + exp(-(NK - eparams.beta) * eparams.alpha));
					}
					else { // density-independent
						pbDisp = pbEmig[0][sexId] + eparams.d0;
					}
				}
			} // end of individual variability
			else { // no individual variability
				pbDisp = pbEmig[stgId][sexId];
			}

			disp = pRandom->Bernoulli(pbDisp);

			if (disp == 1) { // emigrant
				inds[i]->setStatus(dispersing);
			}
		} // end of if (ind.status == initial) condition
	} // end of for loop
}

// All individuals emigrate after patch destruction
void Population::allEmigrate() {
	for (auto ind : inds) {
		ind->setStatus(dispersing);
	}
}

// If an Individual has been identified as an emigrant, remove it from the Population
disperser Population::extractDisperser(int ix) {
	disperser d = disperser();
	indStats ind = inds[ix]->getStats();
	if (ind.status == 1) { // emigrant
		d.pInd = inds[ix]; 
		d.isDispersing = true;
		inds[ix] = nullptr;
		nInds[ind.stage][ind.sex]--;
	}
	else {
		d.pInd = nullptr; 
		d.isDispersing = false;
	}
	return d;
}

// Remove emigrants from their natal patch and add to a map of vectors
void Population::recruitDispersers(std::vector<Individual*>& disperserPool) {
	
	for (auto& pInd : inds) {
		if (pInd->getStatus() == dispersing) {
			disperserPool.push_back(std::move(pInd));
		}
	}
	clean();
}

// For an individual identified as being in the matrix population:
// if it is a settler, return its new location and remove it from the current population
// otherwise, leave it in the matrix population for possible reporting before deletion
disperser Population::extractSettler(int ix) {
	
	disperser d = disperser();
	indStats ind = inds[ix]->getStats();
	Cell* pCell = inds[ix]->getCurrCell();
	d.pInd = inds[ix];  
	d.pCell = pCell;
	d.isSettling = false;
	if (ind.status == settled || ind.status == settledNeighbour) {
		d.isSettling = true;
		inds[ix] = nullptr;
		nInds[ind.stage][ind.sex]--;
	}
	return d;
}

// Add a specified individual to the new/current dispersal group
void Population::recruit(Individual* pInd) {
	indStats ind = pInd->getStats();
	nInds[ind.stage][ind.sex]++;
#ifdef _OPENMP
	const std::lock_guard<std::mutex> lock(inds_mutex);
#endif // _OPENMP
	inds.push_back(pInd);
}

//---------------------------------------------------------------------------

// Add all individuals in the population to the disperser pool
void Population::disperseMatrix(std::vector<Individual*>& dispPool) 
{
	dispPool = move(inds);
}

// Transfer is run for populations in the matrix only
int Population::resolveTransfer(vector<Individual*>& dispPool, Landscape* pLandscape, short landIx)
{
	int nbDispersers = 0;
	bool isDispersing;
	short oppositeSex;
	bool mateOK, densdepOK;
	Patch* patch;
	Population* pPop;
	int patchnum;
	double localK, settProb;
	int density;
	Patch* pPatch = nullptr;
	Cell* pCell = nullptr;
	indStats ind;
	locn newloc = locn();
	locn neighbourLoc = locn();

	landData ppLand = pLandscape->getLandData();
	short reptype = pSpecies->getRepType();
	transferRules trfr = pSpecies->getTransferRules();
	settleType settletype = pSpecies->getSettle();
	settleRules sett;
	settleTraits settDensDep;
	settlePatch settle;
	simParams sim = paramsSim->getSim();

	for (auto& pInd : dispPool) {

		if (trfr.usesMovtProc) {
			// Resolve a single movement step
			isDispersing = pInd->moveStep(pLandscape, landIx, sim.absorbing);
		}
		else {
			// Resolve the (only) movement step
			isDispersing = pInd->moveKernel(pLandscape, sim.absorbing);
		}
		if (isDispersing) nbDispersers++;

		// Record potential settlers to each patch
		if (isDispersing
			&& reptype > 0 // always settle if asexual 
			&& pInd->getStatus() == waitSettlement // disperser has found a patch
			) {
			pCell = pInd->getCurrCell();
			pPatch = pCell->getPatch(pSpecies->getID());
			if (pPatch != nullptr) { // not no-data area
				pPatch->incrPossSettler(pInd->getSex());
			}
		}
	}

	for (auto& pInd : inds) {
		ind = pInd->getStats();

		short stgId = settletype.stgDep ? ind.stage : 0;
		short sexId = settletype.sexDep ? ind.sex : 0;
		sett = pSpecies->getSettRules(stgId, sexId);

		// Resolve candidate settlers
		if (ind.status == waitSettlement) { // awaiting settlement
			pCell = pInd->getCurrCell();
			if (pCell == nullptr) {
				// this condition can occur in a patch-based model at the time of a dynamic landscape
				// change when there is a range restriction in place, since a patch can straddle the
				// range restriction and an individual forced to disperse upon patch removal could
				// start its trajectory beyond the boundary of the restrictyed range - such a model is 
				// not good practice, but the condition must be handled by killing the individual conceerned
				ind.status = diedInTransfer;
			}
			else {
				// Check for reqt of potential mate present, if applicable
				oppositeSex = ind.sex == 0 ? 1 : 0;
				mateOK = sett.findMate ? isMatePresent(pCell, oppositeSex) : true;
				
				densdepOK = false;
				settle = pInd->getSettPatch();

				if (sett.densDep) {
					pPatch = pCell->getPatch(pSpecies->getID());
					if (pPatch != nullptr) { // not no-data area
						
						if (settle.settleStatus == 0 // not yet resolved
							|| settle.pSettPatch != pPatch) {
							// note: second condition allows for having moved from one patch to another
							// adjacent one
						
							// Resolve settlement density-dependence
							localK = pPatch->getK();
							pPop = pPatch->getPop();

							// Get local density
							if (pPop == nullptr) { // empty patch
								density = 0.0;
							} else {
								density = pPop->getNbInds();
							}
							if (localK > 0.0) {

								if (settletype.indVar) settDensDep = pInd->getIndSettTraits();
								else settDensDep = pSpecies->getSpSettTraits(ind.stage, ind.sex);

								settProb = settDensDep.s0 / 
									(1.0 + exp(-(density / localK - settDensDep.beta) * settDensDep.alpha));
								
								if (pRandom->Bernoulli(settProb)) {
									densdepOK = true;
									settle.settleStatus = 2; // can settle
								}
								else {
									settle.settleStatus = 1; // won't settle
								}
								settle.pSettPatch = pPatch;
							}
							pInd->setSettPatch(settle);
						}
						else if (settle.settleStatus == 2) { // previously allowed to settle
							densdepOK = true;
						}
					}
				} else { // no density-dependent settlement
					densdepOK = true;
					settle.settleStatus = 2;
					settle.pSettPatch = pPatch;
					pInd->setSettPatch(settle);
				}

				// Update individual status
				if (mateOK && densdepOK) { // can recruit to patch
					ind.status = settled;
					nbDispersers--;
				} else { // does not recruit
					if (trfr.usesMovtProc) {
						ind.status = dispersing; // continue dispersing, 
						// unless maximum steps has been exceeded
						pathSteps steps = pInd->getSteps();
						settleSteps settsteps = pSpecies->getSteps(ind.stage, ind.sex);
						if (steps.year >= settsteps.maxStepsYr) {
							ind.status = waitNextDispersal;
						}
						if (steps.total >= settsteps.maxSteps) {
							ind.status = diedInTransfer;
						}
					}
					else { // dispersal kernel
						ind.status = sett.wait ? waitNextDispersal : diedInTransfer;
						nbDispersers--;
					}
				}
			}

			pInd->setStatus(ind.status);
		}

#if RS_RCPP
		// write each individuals current movement step and status to paths file
		if (trfr.usesMovtProc && sim.outPaths) {
			if (nextseason >= sim.outStartPaths && nextseason % sim.outIntPaths == 0) {
				inds[i]->outMovePath(nextseason);
			}
		}
#endif

		// Determine whether move to neighbouring cell is possible
		if (!trfr.usesMovtProc // kernel only
			&& sett.goToNeighbourLocn 
			&& (ind.status == waitNextDispersal || ind.status == diedInTransfer)) {

			pCell = pInd->getCurrCell();
			newloc = pCell->getLocn();
			vector <Cell*> neighbourCells;

			// Find which adjacent cells are suitable, if any
			for (int dx = -1; dx < 2; dx++) {
				for (int dy = -1; dy < 2; dy++) {

					if (dx == 0 && dy == 0) break; // skip the current cell

					neighbourLoc.x = newloc.x + dx;
					neighbourLoc.y = newloc.y + dy;

					bool neighbourWithinLandscape = neighbourLoc.x >= 0 
						&& neighbourLoc.x <= ppLand.maxX
						&& neighbourLoc.y >= 0 
						&& neighbourLoc.y <= ppLand.maxY;
					if (!neighbourWithinLandscape) break; // unsuitable

					pCell = pLandscape->findCell(neighbourLoc.x, neighbourLoc.y);
					if (pCell != nullptr) { // not no-data area
						pPatch = pCell->getPatch(pSpecies->getID());
						if (pPatch != nullptr) { // not no-data area

							// Check whether patch is suitable
							if (!pPatch->isMatrix()
								&& pPatch != pInd->getNatalPatch()// or natal patch
								&& pPatch->isSuitable()) {

								// Check mate reqt if applicable
								if (sett.findMate) {
									oppositeSex = ind.sex == 0 ? 1 : 0;
									if (isMatePresent(pCell, oppositeSex))
										neighbourCells.push_back(pCell);
								}
								else neighbourCells.push_back(pCell);
							}
						}
					}
				}
			}
			// Draw a suitable adjacent cell, if any
			if (neighbourCells.size() > 0) {
				vector<Cell*> destCell;
				sample(neighbourCells.begin(), neighbourCells.end(), std::back_inserter(destCell), 1, pRandom->getRNG());
				pInd->moveTo(destCell[0]);
			}
		}

	} // end of individuals loop

	return nbDispersers;
}

// Determine whether there is a potential mate present in a patch which a potential
// settler has reached
bool Population::isMatePresent(Cell* pCell, short othersex)
{
	Patch* pPatch;
	Population* pNewPopn;

	pPatch = pCell->getPatch(pSpecies->getID());
	if (pPatch != nullptr) {
		if (!pPatch->isMatrix() && pPatch->isSuitable()) {

			pNewPopn = pPatch->getPop();
			if (pNewPopn != nullptr) {
				for (int stg = 0; stg < nStages; stg++) {
					if (pNewPopn->nInds[stg][othersex] > 0)
						return true; // one is enough
				}
			}
			// If empty, check for incoming settlers
			if (pPatch->getPossSettlers(othersex) > 0)
				return true;
		}
	}
	return false; // no mates? :(
}

// Add specified individuals to the population
void Population::recruitMany(std::vector<Individual*>& recruits) {
	if (recruits.empty()) return;
	for (Individual* pInd : recruits) {
		indStats ind = pInd->getStats();
		nInds[ind.stage][ind.sex]++;
	}
#ifdef _OPENMP
	const std::lock_guard<std::mutex> lock(inds_mutex);
#endif // _OPENMP
	inds.insert(inds.end(), recruits.begin(), recruits.end());
	recruits.clear();
}

//---------------------------------------------------------------------------
// Determine survival and development and record in individual's status code
// Changes are NOT applied to the Population at this stage
void Population::drawSurvivalDevlpt(bool resolveJuvs, bool resolveAdults, bool resolveDev, bool resolveSurv)
{
	densDepParams ddparams = pSpecies->getDensDep();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	
	double localK = pPatch->getK();

	// get current population size
	int ninds = inds.size();
	if (ninds == 0) return;

	// set up local copies of species development and survival tables
	int nsexes = dem.repType == 0 ? 1 : 2;
	float dev[gMaxNbStages][gMaxNbSexes];
	float surv[gMaxNbStages][gMaxNbSexes];
	short minAge[gMaxNbStages][gMaxNbSexes];

	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {

			if (dem.stageStruct) {
				int sexId = dem.repType == 1 ? 0 : sex; // if simple sexual, both sexes use female parameters
				dev[stg][sex] = resolveDev ? pSpecies->getDev(stg, sexId) : 0.0;
				surv[stg][sex] = resolveSurv ? pSpecies->getSurv(stg, sexId) : 1.0;
				minAge[stg][sex] = pSpecies->getMinAge(stg, sexId);
			}
			else { // non-structured population
				// all juveniles survive and develop, all adults die
				dev[stg][sex] = stg == 0 ? 1.0 : 0.0;
				surv[stg][sex] = stg == 0 ? 1.0 : 0.0;
				minAge[stg][sex] = 0;
			}
		}
	}

	if (dem.stageStruct) {
		
		// if (dev)
		// += dev recd
		// +=

		// if resolve dev and devDD and stage > 0
		// dev DD


	}

	if (dem.stageStruct) {

		for (int stg = 0; stg < nStages; stg++) {
			for (int sex = 0; sex < nsexes; sex++) {

				if (resolveDev) {

					// Contribution of received effects from interspecific interactions
					dev[stg][sex] += devRecdEffects[stg];

					// Calculate development density-dependence
					if (sstruct.devDens && stg > 0 && localK > 0.0) {
						// NB DD in development does NOT apply to juveniles,
						float density = 0.0;

						if (sstruct.devStageDens) { // stage-specific density dependence
							// NOTE: matrix entries represent effect of ROW on COLUMN 
							// AND males precede females
							float weight = 0.0;
							for (int effStg = 0; effStg < nStages; effStg++) {
								for (int effSex = 0; effSex < nSexes; effSex++) {
									if (dem.repType == 2) {
										int rowIncr = effSex == 0 ? 1 : 0;
										int colIncr = sex == 0 ? 1 : 0;
										weight = pSpecies->getDDwtDev(2 * stg + colIncr, 2 * effStg + rowIncr);
									}
									else weight = pSpecies->getDDwtDev(stg, effStg);
									density += nInds[effStg][effSex] * weight;
								}
							}
						}
						else density = getNbInds(); // no stage-dependence

						dev[stg][sex] *= exp(-(ddparams.devCoeff * density) / localK);
					}
					
					// Contribution from interspecific resource-dependent interactions
					dev[stg][sex] *= exp(devResDepEffects[stg]);

					// Contribution from initiated species interactions
					dev[stg][sex] += devInitdEffects[stg];

					// Ensure development probability remains bounded between 0 and 1
					dev[stg][sex] = max(0.0f, min(1.0f, dev[stg][sex]));
				}
				
				if (resolveSurv) {
					// Contribution of received effects from interspecific interactions
					surv[stg][sex] += survRecdEffects[stg];

					// Calculate survival density-dependence
					if (sstruct.survDens) {
						float density = 0.0;

						if (sstruct.survStageDens) { // stage-specific density dependence
							// NOTE: matrix entries represent effect of ROW on COLUMN 
							// AND males precede females
							float weight = 0.0;
							for (int effStg = 0; effStg < nStages; effStg++) {
								for (int effSex = 0; effSex < nSexes; effSex++) {
									if (dem.repType == 2) {
										int rowIncr = effSex == 0 ? 1 : 0;
										int colIncr = sex == 0 ? 1 : 0;
										weight = pSpecies->getDDwtSurv(2 * stg + colIncr, 2 * effStg + rowIncr);
									}
									else weight = pSpecies->getDDwtSurv(stg, effStg);
									density += nInds[effStg][effSex] * weight;
								}
							}
						}
						else density = getNbInds();

						surv[stg][sex] *= exp(-(ddparams.survCoeff * density) / localK);
					}
					// Contribution from interspecific resource-dependent interactions
					surv[stg][sex] *= exp(survResDepEffects[stg]);

					// Contribution from initiated species interactions
					surv[stg][sex] += survInitdEffects[stg];

					// Ensure survival probability remains bounded between 0 and 1
					surv[stg][sex] = max(0.0f, min(1.0f, surv[stg][sex]));
				}
			} // sex loop
		} // stage loop

	} // if stage structure

	// Draw survival and development
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();

		if ((ind.stage == 0 && resolveJuvs) 
			|| (ind.stage > 0 && resolveAdults)
			&& isAlive(ind.status)
			){ 

			// Does the individual survive?
			double probSurvives = surv[ind.stage][ind.sex];
			if (!pRandom->Bernoulli(probSurvives)) {
				inds[i]->setStatus(diedDemogrMort); // doomed to die
			}
			else { 
				// Does the individual develop?
				double probDevelops = dev[ind.stage][ind.sex];
				if (ind.stage < nStages - 1 // not final stage
					&& ind.age >= minAge[ind.stage + 1][ind.sex] // old enough
					) {
					if (pRandom->Bernoulli(probDevelops)) {
						inds[i]->setToDevelop();
					}
				}
			}
		}
	}
}

// Apply survival changes to the population
void Population::applySurvivalDevlpt()
{
	int ninds = inds.size();
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();

		if (!isAlive(ind.status)) {
			delete inds[i];
			inds[i] = nullptr;
			nInds[ind.stage][ind.sex]--;
		}
		else {
			if (ind.isDeveloping) { // develops to next stage
				nInds[ind.stage][ind.sex]--;
				inds[i]->develop();
				nInds[ind.stage + 1][ind.sex]++;
			}
		}
	}
	clean();
}

void Population::ageIncrement(void) {
	int ninds = (int)inds.size();
	stageParams sstruct = pSpecies->getStageParams();
	for (int i = 0; i < ninds; i++) {
		inds[i]->ageIncrement(sstruct.maxAge);
	}
}

void Population::resolveResMedtdInteractions() {

	int nbStg = pSpecies->getStageParams().nStages;

	for (int stg = 0; stg < nbStg; stg++) {

		const auto& allResDepInteractions = pSpecies->getAllResDepInteractions(stg);

		for (auto& [partnerSpStg, interaction] : allResDepInteractions) {

			// Find all populations of target species that are in contact with this one
			const auto& patchesInContact = pPatch->getOverlappingPatches(partnerSpStg.first);
			for (auto& [pContactPatch, overlap] : patchesInContact) {

				auto pTargetPop = pContactPatch->getPop();
				if (pTargetPop == nullptr) continue; // empty patch

				// Get abundance scaled down by the % of overlap between the two patches
				double partnerAbundance = pTargetPop->getNbInds(partnerSpStg.second);
				partnerAbundance *= overlap;

				for (auto& [whichProcess, alpha] : interaction.alphas) {
					
					switch (whichProcess)
					{
					case FEC:
						this->fecResDepEffects[stg] += alpha * partnerAbundance;
						break;
					case DEV:
						this->devResDepEffects[stg] += alpha * partnerAbundance;
						break;
					case SURV:
						this->survResDepEffects[stg] += alpha * partnerAbundance;
						break;
					default:
						break;
					}
				}
			}
		}
	}
}

void Population::resolveInitiatedInteractions() {

	int nbStg = pSpecies->getStageParams().nStages;

	for (int stg = 0; stg < nbStg; stg++) {

		double initiatorAbundance = getNbInds(stg);
		if (initiatorAbundance == 0.0) continue; // no eligible individuals for interaction

		map<pair<Population*, int>, double> interactionRates; // C_i, one entry per target population and stage
		double totalIntrctRate = 0.0; // sum_k (h_k * C_k)
		double totalPreference = 0.0; // sum_k (pi_k * N_k)

		// Initiator to recipient abundance ratios
		map<pair<species_id, int>, double> ratios;

		// Loop through all initiated interactions involving this process and stage
		const auto& allInitdInteractions = pSpecies->getAllInitdInteractions(stg);

		for (auto& [targetSpStg, interaction] : allInitdInteractions) {

			const species_id tgtSp = targetSpStg.first;
			const int tgtStg = targetSpStg.second;
			double ratio = 0.0, denominator = 0.0;

			// Find all populations of target species that are in contact with this one
			const auto& patchesInContact = pPatch->getOverlappingPatches(tgtSp);
			for (auto& [pContactPatch, overlap] : patchesInContact) {

				auto pTargetPop = pContactPatch->getPop();
				if (pTargetPop == nullptr || overlap == 0.0) continue; // no interaction

				// Get abundances scaled down by the % of overlap between the two patches
				double targetAbundance = pTargetPop->getNbInds(tgtStg);
				if (targetAbundance == 0.0) continue; // no individuals of that stage here
				targetAbundance *= overlap;
				double effctvInitrAbundance = initiatorAbundance * overlap;
				ratio += effctvInitrAbundance;
				denominator += targetAbundance;

				// Calculate interaction rate terms
				double intrctRate = interaction.attackRate
					* pow(targetAbundance, interaction.hullCoeff); // a_i * N_i^h
				double interference = interaction.interfIntercept // omega_i + N_p^q
					+ pow(effctvInitrAbundance, interaction.interfExponent);
				if (interference != 0.0) intrctRate /= interference;

				if (interaction.usesRelPref) {
					double targetPreference = interaction.relPreference * targetAbundance; // pi_i * N_i
					totalPreference += targetPreference; // sum_k (pi_k * N_k)
				}

				totalIntrctRate += interaction.handlingTime * intrctRate; // h_i * C_i

				interactionRates.emplace(make_pair(pTargetPop, tgtStg), intrctRate);

			}

			ratio /= denominator;
			ratios.emplace(targetSpStg, ratio);
		}

		// Re-scale the preference terms once we know their sum
		if (totalPreference > 0.0) {
			for (auto& [targetStgPop, intrctRate] : interactionRates) {
				intrctRate /= totalPreference; // divide by sum_k (pi_k * c_k)
			}
			totalIntrctRate /= totalPreference; // divide by sum_k(pi_k * c_k)
		}
		
		// Finish the calculation of the functional reponse
		// and increment the corresponding sum of effects
		for (auto& [targetStgPop, intrctRate] : interactionRates) {

			double funcResp = intrctRate / (1 + totalIntrctRate);

			const species_id tgtSp = targetStgPop.first->getSpecies()->getID();
			const int tgtStg = targetStgPop.second;
			const initdIntrctParams intrct = allInitdInteractions.at(make_pair(tgtSp, tgtStg));

			// Add interaction effect on initiator to each relevant process
			for (auto& [whichProcess, beta] : intrct.betas) {
				switch (whichProcess)
				{
				case FEC:
					this->fecInitdEffects[stg] += beta * funcResp;
					break;
				case DEV:
					this->devInitdEffects[stg] += beta * funcResp;
					break;
				case SURV:
					this->survInitdEffects[stg] += beta * funcResp;
					break;
				default:
					break;
				}
			}
			
			// Add interaction effects on recipient 
			auto initiatorSpStg = make_pair(pSpecies->getID(), stg);
			targetStgPop.first->addRecvdIntrctEffect(tgtStg, initiatorSpStg, ratios, funcResp);
		}

	} // stage loop
}

void Population::addRecvdIntrctEffect(const int& stg, const pair<species_id, int>& initiatorSpStg, 
	const map<pair<species_id, int>, double>& ratios, const double& funcResp) {
	
	const recdIntrctParams intrct = pSpecies->getAllRecdInteractions(stg).at(initiatorSpStg);
	const auto tgtSpStg = make_pair(pSpecies->getID(), stg);
	double ratio = ratios.at(tgtSpStg);

	for (auto& [whichProcess, delta] : intrct.deltas) {
		switch (whichProcess)
		{
		case FEC:
			fecRecdEffects[stg] += delta * funcResp * ratio;
			break;
		case DEV:
			devRecdEffects[stg] += delta * funcResp * ratio;
			break;
		case SURV:
			survRecdEffects[stg] += delta * funcResp * ratio;
			break;
		default:
			break;
		}
	}
}

void Population::resetIntrctEffects() {
	for (int stg = 0; stg < pSpecies->getStageParams().nStages; stg++) {
		fecInitdEffects[stg] = 0.0;
		fecRecdEffects[stg] = 0.0;
		fecResDepEffects[stg] = 0.0;
		devInitdEffects[stg] = 0.0;
		devRecdEffects[stg] = 0.0;
		devResDepEffects[stg] = 0.0;
		survInitdEffects[stg] = 0.0;
		survRecdEffects[stg] = 0.0;
		survResDepEffects[stg] = 0.0;
	}
}

//---------------------------------------------------------------------------
// Remove zero pointers to dead or dispersed individuals
void Population::clean()
{
	int ninds = inds.size();
	if (ninds > 0) {
		inds.erase(std::remove(inds.begin(), inds.end(), (Individual *)nullptr), inds.end());
#ifdef RS_RCPP || NDEBUG
		// do not randomise individuals in DEBUG mode, as the function uses rand()
		// and therefore the randomisation will differ between identical runs of RS
		shuffle(inds.begin(), inds.end(), pRandom->getRNG());
#endif
	}
}

//---------------------------------------------------------------------------
// Write record to population file
void Population::outPopulation(ofstream& outPopOfs, int rep, int yr, int gen, bool envLocal, float epsGlobal,
	bool usesPatches, bool writeEnv, bool gradK)
{
	Cell* pCell;
	demogrParams dem = pSpecies->getDemogrParams();

	popStats p;
	outPopOfs << rep << "\t" << yr << "\t" << gen;
	if (usesPatches) {
		outPopOfs << "\t" << pPatch->getPatchNum();
		outPopOfs << "\t" << pPatch->getNCells();
	}
	else {
		locn loc = pPatch->getCellLocn(0);
		outPopOfs << "\t" << loc.x << "\t" << loc.y;
	}

	float eps = 0.0;
	if (writeEnv) {
		if (envLocal) { // then override eps with local value
			Cell* pCell = pPatch->getRandomCell();
			if (pCell != nullptr) eps = pCell->getEps();
		}
		if (pPatch->isMatrix()) {
			outPopOfs << "\t0\t0\t0";
		}
		else {
			float k = pPatch->getK();
			outPopOfs << "\t" << eps << "\t" << pPatch->getGradVal() << "\t" << k;
		}
	}
	if (dem.stageStruct) {
		p = getStats();
		outPopOfs << "\t" << p.nNonJuvs;
		// non-juvenile stage totals from permanent array
		for (int stg = 1; stg < nStages; stg++) {
			for (int sex = 0; sex < nSexes; sex++) {
				outPopOfs << "\t" << nInds[stg][sex];
			}
		}
		// juveniles from permanent array
		for (int sex = 0; sex < nSexes; sex++) {
			outPopOfs << "\t" << nInds[0][sex];
		}
	}
	else { // non-structured population
		outPopOfs << "\t" << getNbInds();
		if (dem.repType != 0)
		{ // sexual model
			outPopOfs << "\t" << nInds[1][0] << "\t" << nInds[1][1];
		}
	}
	outPopOfs << endl;
}

//---------------------------------------------------------------------------
// Write records to individuals file
void Population::outIndividual(ofstream& outIndsOfs, Landscape* pLandscape, int rep, int yr, int gen)
{
	bool writeInd;
	pathSteps steps;
	Cell* pCell;
	landParams ppLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	species_id speciesID = pSpecies->getID();
	int patchNum = pPatch->getPatchNum();

	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();
		if (yr == -1) { // write all initialised individuals
			writeInd = true;
			outIndsOfs << rep << "\t" << yr << "\t" << dem.repSeasons - 1;
		}
		else {
			if (dem.stageStruct && gen < 0) { // write dying old age individuals only
				if (ind.status == diedOldAge) {
					writeInd = true;
					outIndsOfs << rep << "\t" << yr << "\t" << dem.repSeasons - 1;
				}
				else writeInd = false;
			}
			else {
				writeInd = true;
				outIndsOfs << rep << "\t" << yr << "\t" << gen;
			}
		}
		if (writeInd) {
			outIndsOfs << "\t" << inds[i]->getId();
			if (dem.stageStruct) outIndsOfs << "\t" << to_string(ind.status);
			else { // non-structured population
				outIndsOfs << "\t" << to_string(ind.status);
			}
			pCell = inds[i]->getCurrCell();
			locn loc;
			if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
			else loc = pCell->getLocn();
			pCell = inds[i]->getPrevCell();
			locn natalloc = pCell->getLocn();
			if (ppLand.usesPatches) {
				outIndsOfs << "\t" << inds[i]->getNatalPatch()->getPatchNum();
				if (loc.x == -1) outIndsOfs << "\t-1";
				else outIndsOfs << "\t" << patchNum;
			}
			else { // cell-based model
				outIndsOfs << "\t" << (float)natalloc.x << "\t" << natalloc.y;
				outIndsOfs << "\t" << (float)loc.x << "\t" << (float)loc.y;
			}
			if (dem.repType != 0) outIndsOfs << "\t" << ind.sex;
			if (dem.stageStruct) outIndsOfs << "\t" << ind.age << "\t" << ind.stage;

			if (pSpecies->getNbGenLoadTraits() > 0) outIndsOfs << "\t" << inds[i]->getGeneticFitness();
		
			if (emig.indVar) {
				emigTraits e = inds[i]->getIndEmigTraits();
				if (emig.densDep) {
					outIndsOfs << "\t" << e.d0 << "\t" << e.alpha << "\t" << e.beta;
				}
				else {
					outIndsOfs << "\t" << e.d0;
				}
			} // end of if (emig.indVar)
			if (trfr.indVar) {
				if (trfr.usesMovtProc) {
					if (trfr.moveType == 1) { // SMS
						trfrSMSTraits s = inds[i]->getIndSMSTraits();
						outIndsOfs << "\t" << s.dp << "\t" << s.gb;
						outIndsOfs << "\t" << s.alphaDB << "\t" << s.betaDB;
					} // end of SMS
					if (trfr.moveType == 2) { // CRW
						trfrCRWTraits c = inds[i]->getIndCRWTraits();
						outIndsOfs << "\t" << c.stepLength << "\t" << c.rho;
					} // end of CRW
				}
				else { // kernel
					trfrKernelParams k = inds[i]->getIndKernTraits();
					if (trfr.twinKern)
					{
						outIndsOfs << "\t" << k.meanDist1 << "\t" << k.meanDist2 << "\t" << k.probKern1;
					}
					else {
						outIndsOfs << "\t" << k.meanDist1;
					}
				}
			}

			if (sett.indVar) {
				settleTraits s = inds[i]->getIndSettTraits();
				outIndsOfs << "\t" << s.s0 << "\t" << s.alpha << "\t" << s.beta;
			}

			// distance moved (metres)
			if (loc.x == -1) outIndsOfs << "\t-1";
			else {
				float d = ppLand.resol * sqrt((float)((natalloc.x - loc.x) * (natalloc.x - loc.x)
					+ (natalloc.y - loc.y) * (natalloc.y - loc.y)));
				outIndsOfs << "\t" << d;
			}
#ifndef NDEBUG
			// ALWAYS WRITE NO. OF STEPS
			steps = inds[i]->getSteps();
			outIndsOfs << "\t" << steps.year;
#else
			if (trfr.usesMovtProc) {
				steps = inds[i]->getSteps();
				outIndsOfs << "\t" << steps.year;
			}
#endif
			outIndsOfs << endl;
		} // end of writeInd condition

	}
}

void Population::outputTraitPatchInfo(ofstream& outtraits, int rep, int yr, int gen, bool usesPatches)
{
	if (pPatch->isSuitable() && this->getNbInds() > 0) {
		outtraits << rep << "\t" << yr << "\t" << gen;
		if (usesPatches) {
			outtraits << "\t" << pPatch->getPatchNum();
		}
		else {
			locn loc = pPatch->getCellLocn(0);
			outtraits << "\t" << loc.x << "\t" << loc.y;
		}
	}
}

// Write records to traits file and return aggregated sums
traitsums Population::outTraits(ofstream& outtraits, const bool& writefile)
{
	int popsize, ploidy;
	simParams sim = paramsSim->getSim();
	traitsums ts, indTraitsSums;

	// generate output for each population within the sub-community (patch)
	// provided that the patch is suitable (i.e. non-zero carrying capacity)
	
	Species* pSpecies;
	if (pPatch->isSuitable() && this->getNbInds() > 0) {

		pSpecies = this->getSpecies();
		demogrParams dem = pSpecies->getDemogrParams();
		emigRules emig = pSpecies->getEmigRules();
		transferRules trfr = pSpecies->getTransferRules();
		settleType sett = pSpecies->getSettle();

		indTraitsSums = this->getIndTraitsSums();

		if (emig.indVar) {

			ploidy = emig.sexDep ? 2 : 1;
			vector<double> mnD0, mnAlpha, mnBeta, sdD0, sdAlpha, sdBeta;
			mnD0 = mnAlpha = mnBeta = sdD0 = sdAlpha = sdBeta = vector<double>(2, 0.0);

			for (int whichChr = 0; whichChr < ploidy; whichChr++) {

				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ploidy == 2) popsize = indTraitsSums.ninds[whichChr];
				else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];

				if (popsize > 0) {
					mnD0[whichChr] = indTraitsSums.sumD0[whichChr] / (double)popsize;
					mnAlpha[whichChr] = indTraitsSums.sumAlpha[whichChr] / (double)popsize;
					mnBeta[whichChr] = indTraitsSums.sumBeta[whichChr] / (double)popsize;
					if (popsize > 1) {
						sdD0[whichChr] = indTraitsSums.ssqD0[whichChr] / (double)popsize - mnD0[whichChr] * mnD0[whichChr];
						if (sdD0[whichChr] > 0.0) sdD0[whichChr] = sqrt(sdD0[whichChr]); else sdD0[whichChr] = 0.0;
						sdAlpha[whichChr] = indTraitsSums.ssqAlpha[whichChr] / (double)popsize - mnAlpha[whichChr] * mnAlpha[whichChr];
						if (sdAlpha[whichChr] > 0.0) sdAlpha[whichChr] = sqrt(sdAlpha[whichChr]); else sdAlpha[whichChr] = 0.0;
						sdBeta[whichChr] = indTraitsSums.ssqBeta[whichChr] / (double)popsize - mnBeta[whichChr] * mnBeta[whichChr];
						if (sdBeta[whichChr] > 0.0) sdBeta[whichChr] = sqrt(sdBeta[whichChr]); else sdBeta[whichChr] = 0.0;
					}
				}
			}
			if (writefile) {
				if (emig.sexDep) {
					outtraits << "\t" << mnD0[0] << "\t" << sdD0[0];
					outtraits << "\t" << mnD0[1] << "\t" << sdD0[1];
					if (emig.densDep) {
						outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
						outtraits << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
						outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
						outtraits << "\t" << mnBeta[1] << "\t" << sdBeta[1];
					}
				}
				else { // sex-independent
					outtraits << "\t" << mnD0[0] << "\t" << sdD0[0];
					if (emig.densDep) {
						outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
						outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
					}
				}
			}
		}

		if (trfr.indVar) {
			ploidy = !trfr.usesMovtProc && trfr.sexDep ? 2 : 1;

			vector<double> mnDist1, mnDist2, mnProp1, mnStepL, mnRho,
				sdDist1, sdDist2, sdProp1, sdStepL, sdRho,
				mnDP, mnGB, mnAlphaDB, mnBetaDB,
				sdDP, sdGB, sdAlphaDB, sdBetaDB;
			mnDist1 = mnDist2 = mnProp1 = mnStepL = mnRho =
				sdDist1 = sdDist2 = sdProp1 = sdStepL = sdRho =
				mnDP = mnGB = mnAlphaDB = mnBetaDB =
				sdDP = sdGB = sdAlphaDB = sdBetaDB = vector<double>(2, 0.0);

			for (int whichChr = 0; whichChr < ploidy; whichChr++) {
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ploidy == 2) popsize = indTraitsSums.ninds[whichChr];
				else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];

				if (popsize > 0) {
					mnDist1[whichChr] = indTraitsSums.sumDist1[whichChr] / (double)popsize;
					mnDist2[whichChr] = indTraitsSums.sumDist2[whichChr] / (double)popsize;
					mnProp1[whichChr] = indTraitsSums.sumProp1[whichChr] / (double)popsize;
					mnStepL[whichChr] = indTraitsSums.sumStepL[whichChr] / (double)popsize;
					mnRho[whichChr] = indTraitsSums.sumRho[whichChr] / (double)popsize;
					mnDP[whichChr] = indTraitsSums.sumDP[whichChr] / (double)popsize;
					mnGB[whichChr] = indTraitsSums.sumGB[whichChr] / (double)popsize;
					mnAlphaDB[whichChr] = indTraitsSums.sumAlphaDB[whichChr] / (double)popsize;
					mnBetaDB[whichChr] = indTraitsSums.sumBetaDB[whichChr] / (double)popsize;
					if (popsize > 1) {
						sdDist1[whichChr] = indTraitsSums.ssqDist1[whichChr] / (double)popsize - mnDist1[whichChr] * mnDist1[whichChr];
						if (sdDist1[whichChr] > 0.0) sdDist1[whichChr] = sqrt(sdDist1[whichChr]); else sdDist1[whichChr] = 0.0;
						sdDist2[whichChr] = indTraitsSums.ssqDist2[whichChr] / (double)popsize - mnDist2[whichChr] * mnDist2[whichChr];
						if (sdDist2[whichChr] > 0.0) sdDist2[whichChr] = sqrt(sdDist2[whichChr]); else sdDist2[whichChr] = 0.0;
						sdProp1[whichChr] = indTraitsSums.ssqProp1[whichChr] / (double)popsize - mnProp1[whichChr] * mnProp1[whichChr];
						if (sdProp1[whichChr] > 0.0) sdProp1[whichChr] = sqrt(sdProp1[whichChr]); else sdProp1[whichChr] = 0.0;
						sdStepL[whichChr] = indTraitsSums.ssqStepL[whichChr] / (double)popsize - mnStepL[whichChr] * mnStepL[whichChr];
						if (sdStepL[whichChr] > 0.0) sdStepL[whichChr] = sqrt(sdStepL[whichChr]); else sdStepL[whichChr] = 0.0;
						sdRho[whichChr] = indTraitsSums.ssqRho[whichChr] / (double)popsize - mnRho[whichChr] * mnRho[whichChr];
						if (sdRho[whichChr] > 0.0) sdRho[whichChr] = sqrt(sdRho[whichChr]); else sdRho[whichChr] = 0.0;
						sdDP[whichChr] = indTraitsSums.ssqDP[whichChr] / (double)popsize - mnDP[whichChr] * mnDP[whichChr];
						if (sdDP[whichChr] > 0.0) sdDP[whichChr] = sqrt(sdDP[whichChr]); else sdDP[whichChr] = 0.0;
						sdGB[whichChr] = indTraitsSums.ssqGB[whichChr] / (double)popsize - mnGB[whichChr] * mnGB[whichChr];
						if (sdGB[whichChr] > 0.0) sdGB[whichChr] = sqrt(sdGB[whichChr]); else sdGB[whichChr] = 0.0;
						sdAlphaDB[whichChr] = indTraitsSums.ssqAlphaDB[whichChr] / (double)popsize - mnAlphaDB[whichChr] * mnAlphaDB[whichChr];
						if (sdAlphaDB[whichChr] > 0.0) sdAlphaDB[whichChr] = sqrt(sdAlphaDB[whichChr]); else sdAlphaDB[whichChr] = 0.0;
						sdBetaDB[whichChr] = indTraitsSums.ssqBetaDB[whichChr] / (double)popsize - mnBetaDB[whichChr] * mnBetaDB[whichChr];
						if (sdBetaDB[whichChr] > 0.0) sdBetaDB[whichChr] = sqrt(sdBetaDB[whichChr]); else sdBetaDB[whichChr] = 0.0;
					}
				}
			}
			if (writefile) {
				if (trfr.usesMovtProc) {
					if (trfr.moveType == 1) {
						outtraits << "\t" << mnDP[0] << "\t" << sdDP[0];
						outtraits << "\t" << mnGB[0] << "\t" << sdGB[0];
						outtraits << "\t" << mnAlphaDB[0] << "\t" << sdAlphaDB[0];
						outtraits << "\t" << mnBetaDB[0] << "\t" << sdBetaDB[0];
					}
					if (trfr.moveType == 2) {
						outtraits << "\t" << mnStepL[0] << "\t" << sdStepL[0];
						outtraits << "\t" << mnRho[0] << "\t" << sdRho[0];
					}
				}
				else {
					if (trfr.sexDep) {
						outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
						outtraits << "\t" << mnDist1[1] << "\t" << sdDist1[1];
						if (trfr.twinKern)
						{
							outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
							outtraits << "\t" << mnDist2[1] << "\t" << sdDist2[1];
							outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
							outtraits << "\t" << mnProp1[1] << "\t" << sdProp1[1];
						}
					}
					else { // sex-independent
						outtraits << "\t" << mnDist1[0] << "\t" << sdDist1[0];
						if (trfr.twinKern)
						{
							outtraits << "\t" << mnDist2[0] << "\t" << sdDist2[0];
							outtraits << "\t" << mnProp1[0] << "\t" << sdProp1[0];
						}
					}
				}
			}
		}

		if (sett.indVar) {

			ploidy = sett.sexDep ? 2 : 1;

			// CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
			double mnS0[2], mnAlpha[2], mnBeta[2], sdS0[2], sdAlpha[2], sdBeta[2];

			for (int whichChr = 0; whichChr < ploidy; whichChr++) {

				mnS0[whichChr] = mnAlpha[whichChr] = mnBeta[whichChr] = sdS0[whichChr] = sdAlpha[whichChr] = sdBeta[whichChr] = 0.0;
				// individuals may have been counted by sex if there was
				// sex dependency in another dispersal phase
				if (ploidy == 2) popsize = indTraitsSums.ninds[whichChr];
				else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];

				if (popsize > 0) {

					mnS0[whichChr] = indTraitsSums.sumS0[whichChr] / (double)popsize;
					mnAlpha[whichChr] = indTraitsSums.sumAlphaS[whichChr] / (double)popsize;
					mnBeta[whichChr] = indTraitsSums.sumBetaS[whichChr] / (double)popsize;

					if (popsize > 1) {
						sdS0[whichChr] = indTraitsSums.ssqS0[whichChr] / (double)popsize - mnS0[whichChr] * mnS0[whichChr];
						if (sdS0[whichChr] > 0.0) sdS0[whichChr] = sqrt(sdS0[whichChr]); else sdS0[whichChr] = 0.0;
						sdAlpha[whichChr] = indTraitsSums.ssqAlphaS[whichChr] / (double)popsize - mnAlpha[whichChr] * mnAlpha[whichChr];
						if (sdAlpha[whichChr] > 0.0) sdAlpha[whichChr] = sqrt(sdAlpha[whichChr]); else sdAlpha[whichChr] = 0.0;
						sdBeta[whichChr] = indTraitsSums.ssqBetaS[whichChr] / (double)popsize - mnBeta[whichChr] * mnBeta[whichChr];
						if (sdBeta[whichChr] > 0.0) sdBeta[whichChr] = sqrt(sdBeta[whichChr]); else sdBeta[whichChr] = 0.0;
					}
					else {
						sdS0[whichChr] = sdAlpha[whichChr] = sdBeta[whichChr] = 0.0;
					}
				}
			}
			if (writefile) {
				if (sett.sexDep) {
					outtraits << "\t" << mnS0[0] << "\t" << sdS0[0];
					outtraits << "\t" << mnS0[1] << "\t" << sdS0[1];
					outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
					outtraits << "\t" << mnAlpha[1] << "\t" << sdAlpha[1];
					outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
					outtraits << "\t" << mnBeta[1] << "\t" << sdBeta[1];
				}
				else { // sex-independent
					outtraits << "\t" << mnS0[0] << "\t" << sdS0[0];
					outtraits << "\t" << mnAlpha[0] << "\t" << sdAlpha[0];
					outtraits << "\t" << mnBeta[0] << "\t" << sdBeta[0];
				}
			}
		}

		// Genetic load
		if (pSpecies->getNbGenLoadTraits() > 0) {

			ploidy = pSpecies->isDiploid() ? 2 : 1;
			double mnGenFitness[2], sdGenFitness[2];

			for (int whichChr = 0; whichChr < ploidy; whichChr++) {
				mnGenFitness[whichChr] = sdGenFitness[whichChr] = 0.0;

				if (ploidy == 2) popsize = indTraitsSums.ninds[whichChr];
				else popsize = indTraitsSums.ninds[0] + indTraitsSums.ninds[1];

				if (popsize > 0) {

					mnGenFitness[whichChr] = indTraitsSums.sumGeneticFitness[whichChr] / (double)popsize;
					if (popsize > 1) {
						sdGenFitness[whichChr] = indTraitsSums.ssqGeneticFitness[whichChr] / (double)popsize - mnGenFitness[whichChr] * mnGenFitness[whichChr];
						if (sdGenFitness[whichChr] > 0.0) sdGenFitness[whichChr] = sqrt(sdGenFitness[whichChr]); else sdGenFitness[whichChr] = 0.0;
					}
					else {
						sdGenFitness[whichChr] = 0.0;
					}
				}
			}

			if (writefile) {
				if (pSpecies->getDemogrParams().repType > 0) {
					outtraits << "\t" << mnGenFitness[0] << "\t" << sdGenFitness[0];
					outtraits << "\t" << mnGenFitness[1] << "\t" << sdGenFitness[1];
				}
				else { // sex-independent
					outtraits << "\t" << mnGenFitness[0] << "\t" << sdGenFitness[0];
				}
			}
		}

		if (writefile) outtraits << endl;

		for (int iSex = 0; iSex < gMaxNbSexes; iSex++) {
			ts.ninds[iSex] += indTraitsSums.ninds[iSex];
			ts.sumD0[iSex] += indTraitsSums.sumD0[iSex];
			ts.ssqD0[iSex] += indTraitsSums.ssqD0[iSex];
			ts.sumAlpha[iSex] += indTraitsSums.sumAlpha[iSex];
			ts.ssqAlpha[iSex] += indTraitsSums.ssqAlpha[iSex];
			ts.sumBeta[iSex] += indTraitsSums.sumBeta[iSex];
			ts.ssqBeta[iSex] += indTraitsSums.ssqBeta[iSex];
			ts.sumDist1[iSex] += indTraitsSums.sumDist1[iSex];
			ts.ssqDist1[iSex] += indTraitsSums.ssqDist1[iSex];
			ts.sumDist2[iSex] += indTraitsSums.sumDist2[iSex];
			ts.ssqDist2[iSex] += indTraitsSums.ssqDist2[iSex];
			ts.sumProp1[iSex] += indTraitsSums.sumProp1[iSex];
			ts.ssqProp1[iSex] += indTraitsSums.ssqProp1[iSex];
			ts.sumDP[iSex] += indTraitsSums.sumDP[iSex];
			ts.ssqDP[iSex] += indTraitsSums.ssqDP[iSex];
			ts.sumGB[iSex] += indTraitsSums.sumGB[iSex];
			ts.ssqGB[iSex] += indTraitsSums.ssqGB[iSex];
			ts.sumAlphaDB[iSex] += indTraitsSums.sumAlphaDB[iSex];
			ts.ssqAlphaDB[iSex] += indTraitsSums.ssqAlphaDB[iSex];
			ts.sumBetaDB[iSex] += indTraitsSums.sumBetaDB[iSex];
			ts.ssqBetaDB[iSex] += indTraitsSums.ssqBetaDB[iSex];
			ts.sumStepL[iSex] += indTraitsSums.sumStepL[iSex];
			ts.ssqStepL[iSex] += indTraitsSums.ssqStepL[iSex];
			ts.sumRho[iSex] += indTraitsSums.sumRho[iSex];
			ts.ssqRho[iSex] += indTraitsSums.ssqRho[iSex];
			ts.sumS0[iSex] += indTraitsSums.sumS0[iSex];
			ts.ssqS0[iSex] += indTraitsSums.ssqS0[iSex];
			ts.sumAlphaS[iSex] += indTraitsSums.sumAlphaS[iSex];
			ts.ssqAlphaS[iSex] += indTraitsSums.ssqAlphaS[iSex];
			ts.sumBetaS[iSex] += indTraitsSums.sumBetaS[iSex];
			ts.ssqBetaS[iSex] += indTraitsSums.ssqBetaS[iSex];
			ts.sumGeneticFitness[iSex] += indTraitsSums.sumGeneticFitness[iSex];
			ts.ssqGeneticFitness[iSex] += indTraitsSums.ssqGeneticFitness[iSex];
		}
	}
	return ts;
}

void Population::outputGeneValues(ofstream& ofsGenes, const int& yr, const int& gen) const {
	
	const bool isDiploid = pSpecies->isDiploid();
	int indID;
	float alleleOnChromA, alleleOnChromB;
	float domCoefA, domCoefB;

	// Subset traits that are selected to be output
	set<TraitType> traitTypes = pSpecies->getTraitTypes();
	set<TraitType> outputTraitTypes;
	for (auto trType : traitTypes) {
		if (pSpecies->getSpTrait(trType)->isOutput())
			outputTraitTypes.insert(trType);
	}

	// Fetch map to positions for each trait
	// Presumably faster than fetching for every individual
	map<TraitType, set<int>> allGenePositions;
	for (auto trType : outputTraitTypes) {
		set<int> traitPositions = pSpecies->getSpTrait(trType)->getGenePositions();
		allGenePositions.insert(make_pair(trType, traitPositions));
	}

	set<int> positions;
	for (Individual* ind : sampledInds) {
		indID = ind->getId();
		for (auto trType : outputTraitTypes) {
			positions = allGenePositions[trType];
			auto indTrait = ind->getTrait(trType);
			for (auto pos : positions) {
				alleleOnChromA = indTrait->getAlleleValueAtLocus(0, pos);
				if (trType == GENETIC_LOAD1 || trType == GENETIC_LOAD2 || trType == GENETIC_LOAD3 || trType == GENETIC_LOAD4 || trType == GENETIC_LOAD5) {
					domCoefA = indTrait->getDomCoefAtLocus(0, pos);
				}
				else {
					domCoefA = 0.0;
				}
				ofsGenes << yr << '\t' << gen << '\t' << indID << '\t' << to_string(trType) << '\t' << pos << '\t' << alleleOnChromA << '\t' << domCoefA;
				if (isDiploid) {
					alleleOnChromB = indTrait->getAlleleValueAtLocus(1, pos);
					if (trType == GENETIC_LOAD1 || trType == GENETIC_LOAD2 || trType == GENETIC_LOAD3 || trType == GENETIC_LOAD4 || trType == GENETIC_LOAD5) {
						domCoefB = indTrait->getDomCoefAtLocus(1, pos);
					}
					else {
						domCoefB = 0.0;
					}
					ofsGenes << '\t' << alleleOnChromB << '\t' << domCoefB;
				}
				ofsGenes << endl;
			}
		}
	}
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
