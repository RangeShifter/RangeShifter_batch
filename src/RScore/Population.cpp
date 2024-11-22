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

//---------------------------------------------------------------------------

Population::Population() {
	nSexes = nStages = 0;
	pPatch = NULL;
	pSpecies = NULL;
	return;
}

Population::Population(Species* pSp, Patch* pPch, int ninds, int resol)
{
	// constructor for a Population of a specified size

	int n, nindivs, age = 0, minage, maxage, nAges = 0;
	int cumtotal = 0;
	float probmale;
	double ageprob, ageprobsum;
	std::vector<double> ageProb; // for quasi-equilibrium initial age distribution

	if (ninds > 0) {
		inds.reserve(ninds);
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
	initParams init = paramsInit->getInit();

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
				if (dem.repType == 1) { // simple sexual model
					// both sexes use minimum ages recorded for females
					minAge[stg][sex] = pSpecies->getMinAge(stg, 0);
				}
				else {
					minAge[stg][sex] = pSpecies->getMinAge(stg, sex);
				}
			}
			else { // non-structured population
				minAge[stg][sex] = 0;
			}
		}
	}

	// individuals of new population must be >= stage 1
	for (int stg = 1; stg < nStages; stg++) {
		if (dem.stageStruct) { // allocate to stages according to initialisation conditions
			// final stage is treated separately to ensure that correct total
			// no. of individuals is created
			if (stg == nStages - 1) {
				n = ninds - cumtotal;
			}
			else {
				n = (int)(ninds * paramsInit->getProp(stg) + 0.5);
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
			Cell* pCell = pPatch->getRandomCell();
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

			inds.push_back(DBG_NEW Individual(pCell, pPatch, stg, age, sstruct.repInterval,
				probmale, trfr.usesMovtProc, trfr.moveType));

			sex = inds[nindivs + i]->getSex();
			if (pSpecies->getNTraits() > 0) {
				// individual variation - set up genetics
				inds[nindivs + i]->setUpGenes(pSpecies, resol);
			}
			nInds[stg][sex]++;
		}
	}
}

Population::~Population() {
	int ninds = inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != nullptr) delete inds[i];
	}
	inds.clear();
	juvs.clear();
	sampledInds.clear();
}

traitsums Population::getIndTraitsSums(Species* pSpecies) {
	int g;
	traitsums ts = traitsums();

	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	int ninds = (int)inds.size();
	for (int iInd = 0; iInd < ninds; iInd++) {
		int sex = inds[iInd]->getSex();
		if (emig.sexDep || trfr.sexDep || sett.sexDep) 
			g = sex; 
		else g = 0;
		ts.ninds[g] += 1;

		// emigration traits
		emigTraits e = inds[iInd]->getIndEmigTraits();
		if (emig.sexDep) g = sex; 
		else g = 0;
		ts.sumD0[g] += e.d0;    
		ts.ssqD0[g] += e.d0 * e.d0;
		ts.sumAlpha[g] += e.alpha; 
		ts.ssqAlpha[g] += e.alpha * e.alpha;
		ts.sumBeta[g] += e.beta;  
		ts.ssqBeta[g] += e.beta * e.beta;

		// transfer traits
		if (trfr.usesMovtProc) {

			switch (trfr.moveType) {

			case 1: // SMS
			{
				trfrSMSTraits sms = inds[iInd]->getIndSMSTraits();
				g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				ts.sumDP[g] += sms.dp; 
				ts.ssqDP[g] += sms.dp * sms.dp;
				ts.sumGB[g] += sms.gb;
				ts.ssqGB[g] += sms.gb * sms.gb;
				ts.sumAlphaDB[g] += sms.alphaDB;
				ts.ssqAlphaDB[g] += sms.alphaDB * sms.alphaDB;
				ts.sumBetaDB[g] += sms.betaDB; 
				ts.ssqBetaDB[g] += sms.betaDB * sms.betaDB;
				break;
			}
			case 2:
			{
				trfrCRWTraits c = inds[iInd]->getIndCRWTraits();
				g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				ts.sumStepL[g] += c.stepLength;
				ts.ssqStepL[g] += c.stepLength * c.stepLength;
				ts.sumRho[g] += c.rho;       
				ts.ssqRho[g] += c.rho * c.rho;
				break;
			}
			default:
				throw runtime_error("usesMoveProcess is ON but moveType is neither 1 (SMS) or 2 (CRW).");
				break;
			}
		}
		else {
			trfrKernelParams k = inds[iInd]->getIndKernTraits();
			if (trfr.sexDep) g = sex; 
			else g = 0;
			ts.sumDist1[g] += k.meanDist1; 
			ts.ssqDist1[g] += k.meanDist1 * k.meanDist1;
			ts.sumDist2[g] += k.meanDist2;
			ts.ssqDist2[g] += k.meanDist2 * k.meanDist2;
			ts.sumProp1[g] += k.probKern1; 
			ts.ssqProp1[g] += k.probKern1 * k.probKern1;
		}
		// settlement traits
		settleTraits s = inds[iInd]->getIndSettTraits();
		if (sett.sexDep) g = sex; 
		else g = 0;
		//	g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
		ts.sumS0[g] += s.s0;     
		ts.ssqS0[g] += s.s0 * s.s0;
		ts.sumAlphaS[g] += s.alpha; 
		ts.ssqAlphaS[g] += s.alpha * s.alpha;
		ts.sumBetaS[g] += s.beta;   
		ts.ssqBetaS[g] += s.beta * s.beta;

		if (gMaxNbSexes > 1) g = sex; 
		else g = 0;

		ts.sumGeneticFitness[g] += inds[iInd]->getGeneticFitness();
		ts.ssqGeneticFitness[g] += inds[iInd]->getGeneticFitness() * inds[iInd]->getGeneticFitness();
	}
	return ts;
}

int Population::getNInds(void) { return (int)inds.size(); }

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

popStats Population::getStats(void)
{
	popStats p = popStats();
	int ninds;
	float fec;
	bool breeders[2] = { false, false };
	demogrParams dem = pSpecies->getDemogrParams();
	p.pSpecies = pSpecies;
	p.pPatch = pPatch;
	p.spNum = pSpecies->getSpNum();
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

Species* Population::getSpecies(void) { return pSpecies; }

int Population::totalPop() {
	int t = 0;
	for (int stg = 0; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			t += nInds[stg][sex];
		}
	}
	return t;
}

int Population::stagePop(int stg) {
	int t = 0;
	if (stg < 0 || stg >= nStages) return t;
	for (int sex = 0; sex < nSexes; sex++) {
		t += nInds[stg][sex];
	}
	return t;
}

//---------------------------------------------------------------------------
// Extirpate all populations according to
// option 0 - random local extinction probability
// option 1 - local extinction probability gradient
// NB only applied for cell-based model
void Population::localExtinction(int option) {
	double pExtinct = 0.0;

	if (option == 0) {
		envStochParams env = paramsStoch->getStoch();
		if (env.localExt) pExtinct = env.locExtProb;
	}
	else {
		Cell* pCell = pPatch->getRandomCell(); // get only cell in the patch
		// extinction prob is complement of cell gradient value plus any non-zero prob at the optimum
		pExtinct = 1.0 - pCell->getEnvVal() + paramsGrad->getGradient().extProbOpt;
		if (pExtinct > 1.0) pExtinct = 1.0;
	}

	if (pRandom->Bernoulli(pExtinct)) {
		extirpate();
	}
}

// Remove all Individuals
void Population::extirpate() {
	int ninds = inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != nullptr) delete inds[i];
	}
	inds.clear();
	juvs.clear();
	for (int sex = 0; sex < nSexes; sex++) {
		for (int stg = 0; stg < nStages; stg++) {
			nInds[stg][sex] = 0;
		}
	}
}

//---------------------------------------------------------------------------
// Produce juveniles and hold them in the juvs vector
void Population::reproduction(const float localK, const float envval, const int resol)
{
	// get population size at start of reproduction
	int ninds = inds.size();
	if (ninds == 0) return;

	int nsexes, stage, sex, njuvs, nmales, nfemales;
	Cell* pCell;
	indStats ind;
	double expected;
	bool skipbreeding;

	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	if (dem.repType == 0)
		nsexes = 1; 
	else nsexes = 2;


// set up local copy of species fecundity table
	float fec[gMaxNbStages][gMaxNbSexes];
	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use fecundity recorded for females
					fec[stg][sex] = pSpecies->getFec(stg, 0);
				}
				else
					fec[stg][sex] = pSpecies->getFec(stg, sex);
			}
			else { // non-structured population
				if (stg == 1) fec[stg][sex] = dem.lambda; // adults
				else fec[stg][sex] = 0.0; // juveniles
			}
		}
	}

	if (dem.stageStruct) {
	// apply environmental effects and density dependence
	// to all non-zero female non-juvenile stages
		for (int stg = 1; stg < nStages; stg++) {
			if (fec[stg][0] > 0.0) {
				// apply any effect of environmental gradient and/or stochasticty
				fec[stg][0] *= envval;
				if (env.stoch && !env.inK) {
					// fecundity (at low density) is constrained to lie between limits specified
					// for the species
					float limit;
					limit = pSpecies->getMinMax(0);
					if (fec[stg][0] < limit) fec[stg][0] = limit;
					limit = pSpecies->getMinMax(1);
					if (fec[stg][0] > limit) fec[stg][0] = limit;
				}
				if (sstruct.fecDens) { // apply density dependence
					float effect = 0.0;
					if (sstruct.fecStageDens) { // stage-specific density dependence
						// NOTE: matrix entries represent effect of ROW on COLUMN 
						// AND males precede females
						float weight = 0.0;
						for (int effstg = 0; effstg < nStages; effstg++) {
							for (int effsex = 0; effsex < nSexes; effsex++) {
								if (dem.repType == 2) {
									if (effsex == 0) weight = pSpecies->getDDwtFec(2 * stg + 1, 2 * effstg + 1);
									else weight = pSpecies->getDDwtFec(2 * stg + 1, 2 * effstg);
								}
								else {
									weight = pSpecies->getDDwtFec(stg, effstg);
								}
								effect += (float)nInds[effstg][effsex] * weight;
							}
						}
					}
					else // not stage-specific
						effect = (float)totalPop();
					if (localK > 0.0) fec[stg][0] *= exp(-effect / localK);
				}
			}
		}
	}
	else { // non-structured - set fecundity for adult females only
		// apply any effect of environmental gradient and/or stochasticty
		fec[1][0] *= envval;
		if (env.stoch && !env.inK) {
			// fecundity (at low density) is constrained to lie between limits specified
			// for the species
			float limit;
			limit = pSpecies->getMinMax(0);
			if (fec[1][0] < limit) fec[1][0] = limit;
			limit = pSpecies->getMinMax(1);
			if (fec[1][0] > limit) fec[1][0] = limit;
		}
		// apply density dependence
		if (localK > 0.0) {
			if (dem.repType == 1 || dem.repType == 2) { // sexual model
				// apply factor of 2 (as in manual, eqn. 6)
				fec[1][0] *= 2.0;
			}
			fec[1][0] /= (1.0f + fabs(dem.lambda - 1.0f) * pow(((float)ninds / localK), dem.bc));
		}
	}

	double propBreed;
	Individual* father = nullptr;
	std::vector <Individual*> fathers;

	switch (dem.repType) {

	case 0: // asexual model
		for (int i = 0; i < ninds; i++) {
			stage = inds[i]->breedingFem();
			if (stage > 0) { // female of breeding age
				if (dem.stageStruct) {
					// determine whether she must miss current breeding attempt
					ind = inds[i]->getStats();
					if (ind.fallow >= sstruct.repInterval) {
						if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
						else skipbreeding = true;
					}
					else skipbreeding = true; // cannot breed this time
				}
				else skipbreeding = false; // not structured - always breed
				if (skipbreeding) {
					inds[i]->incFallow();
				}
				else { // attempt to breed
					inds[i]->resetFallow();
					expected = fec[stage][0];
					if (expected <= 0.0) njuvs = 0;
					else njuvs = pRandom->Poisson(expected);
					pCell = pPatch->getRandomCell();
					for (int j = 0; j < njuvs; j++) {

						unique_ptr<Individual> newJuv = make_unique<Individual>(pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);

						if (pSpecies->getNTraits() > 0) {
							newJuv->inheritTraits(pSpecies, inds[i], resol);
						}
						if (newJuv->isViable()) {
							juvs.push_back(move(newJuv));
							nInds[0][0]++;
						}
						// else (unviable) newJuv is deleted when going out of scope
					}
				}
			}
		}
		break;

	case 1: // simple sexual model
	case 2: // complex sexual model
		// count breeding females and males
		// add breeding males to list of potential fathers
		nfemales = nmales = 0;
		for (int i = 0; i < ninds; i++) {
			ind = inds[i]->getStats();
			if (ind.sex == 0 && fec[ind.stage][0] > 0.0) nfemales++;
			if (ind.sex == 1 && fec[ind.stage][1] > 0.0) {
				fathers.push_back(inds[i]);
				nmales++;
			}
		}
		if (nfemales > 0 && nmales > 0) { // population can breed
			if (dem.repType == 2) { // complex sexual model
				// calculate proportion of eligible females which breed
				propBreed = (2.0 * dem.harem * nmales) / (nfemales + dem.harem * nmales);
				if (propBreed > 1.0) propBreed = 1.0;
			}
			else propBreed = 1.0;
			for (int i = 0; i < ninds; i++) {
				stage = inds[i]->breedingFem();
				if (stage > 0 && fec[stage][0] > 0.0) { // (potential) breeding female
					if (dem.stageStruct) {
						// determine whether she must miss current breeding attempt
						ind = inds[i]->getStats();
						if (ind.fallow >= sstruct.repInterval) {
							if (pRandom->Bernoulli(sstruct.probRep)) skipbreeding = false;
							else skipbreeding = true;
						}
						else skipbreeding = true; // cannot breed this time
					}
					else skipbreeding = false; // not structured - always breed
					if (skipbreeding) {
						inds[i]->incFallow();
					}
					else { // attempt to breed
						inds[i]->resetFallow();
						// NOTE: FOR COMPLEX SEXUAL MODEL, NO. OF FEMALES *ACTUALLY* BREEDING DOES NOT
						// NECESSARILY EQUAL THE EXPECTED NO. FROM EQN. 7 IN THE MANUAL...
						if (pRandom->Bernoulli(propBreed)) {
							expected = fec[stage][0]; // breeds
						}
						else expected = 0.0; // fails to breed
						if (expected <= 0.0) njuvs = 0;
						else njuvs = pRandom->Poisson(expected);

						if (njuvs > 0) {
							// select father at random from breeding males ...
							int rrr = 0;
							if (nmales > 1) rrr = pRandom->IRandom(0, nmales - 1);
							father = fathers[rrr];
							pCell = pPatch->getRandomCell();
							for (int j = 0; j < njuvs; j++) {

								unique_ptr<Individual> newJuv = make_unique<Individual>(pCell, pPatch, 0, 0, 0, dem.propMales, trfr.usesMovtProc, trfr.moveType);

								if (pSpecies->getNTraits() > 0) {
									newJuv->inheritTraits(pSpecies, inds[i], father, resol);
								}

								if (newJuv->isViable()) {
									sex = newJuv->getSex();
									nInds[0][sex]++;
									juvs.push_back(move(newJuv));
								}
								// else (unviable) newJuv is deleted when going out of scope

							}
						}
					}
				}
			}
		}
		fathers.clear();
		break;

	} // end of switch (dem.repType)

// THIS MAY NOT BE CORRECT FOR MULTIPLE SPECIES IF THERE IS SOME FORM OF
// CROSS-SPECIES DENSITY-DEPENDENT FECUNDITY
}

// Following reproduction of ALL species, add juveniles to the population prior to dispersal
void Population::fledge()
{
	demogrParams dem = pSpecies->getDemogrParams();

	if (dem.stageStruct) { 
		// juveniles are moved to the individuals vector
		for (int j = 0; j < juvs.size(); j++) {
			inds.push_back(juvs[j].release());
		}
	}
	else { // all adults die and juveniles replace adults
		int ninds = inds.size();
		for (int i = 0; i < ninds; i++) {
			delete inds[i];
		}
		inds.clear();
		for (int sex = 0; sex < nSexes; sex++) {
			nInds[1][sex] = 0; // set count of adults to zero
		}
		for (int j = 0; j < juvs.size(); j++) {
			inds.push_back(juvs[j].release());
		}
	}
	juvs.clear();
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
	NK = static_cast<float>(totalPop()) / localK;

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
						NK = (float)totalPop() / localK;
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
		d.yes = true;
		inds[ix] = 0;
		nInds[ind.stage][ind.sex]--;
	}
	else {
		d.pInd = NULL; 
		d.yes = false;
	}
	return d;
}


// For an individual identified as being in the matrix population:
// if it is a settler, return its new location and remove it from the current population
// otherwise, leave it in the matrix population for possible reporting before deletion
disperser Population::extractSettler(int ix) {
	disperser d = disperser();
	Cell* pCell;

	indStats ind = inds[ix]->getStats();
	pCell = inds[ix]->getLocn(1);
	d.pInd = inds[ix];  
	d.pCell = pCell;
	d.yes = false;
	if (ind.status == settled || ind.status == settledNeighbour) {
		d.yes = true;
		inds[ix] = 0;
		nInds[ind.stage][ind.sex]--;
	}
	return d;
}

// Add a specified individual to the new/current dispersal group
void Population::recruit(Individual* pInd) {
	inds.push_back(pInd);
	indStats ind = pInd->getStats();
	nInds[ind.stage][ind.sex]++;
}

//---------------------------------------------------------------------------

// Transfer is run for populations in the matrix only
int Population::transfer(Landscape* pLandscape, short landIx, short nextseason)
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

	int ninds = inds.size();
	for (int i = 0; i < ninds; i++) {

		if (trfr.usesMovtProc) {
			// Resolve a single movement step
			isDispersing = inds[i]->moveStep(pLandscape, pSpecies, landIx, sim.absorbing);
		}
		else {
			// Resolve the (only) movement step
			isDispersing = inds[i]->moveKernel(pLandscape, pSpecies, sim.absorbing);
		}
		nbDispersers += isDispersing;

		// Record potential settlers to each patch
		if (isDispersing
			&& reptype > 0 // always settle if asexual 
			&& inds[i]->getStatus() == waitSettlement // disperser has found a patch
			) {
			pCell = inds[i]->getLocn(1);
			pPatch = pCell->getPatch();
			if (pPatch != nullptr) { // not no-data area
				pPatch->incrPossSettler(pSpecies, inds[i]->getSex());
			}
		}
	}

	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();

		short stgId = settletype.stgDep ? ind.stage : 0;
		short sexId = settletype.sexDep ? ind.sex : 0;
		sett = pSpecies->getSettRules(stgId, sexId);

		// Resolve candidate settlers
		if (ind.status == waitSettlement) { // awaiting settlement
			pCell = inds[i]->getLocn(1);
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
				settle = inds[i]->getSettPatch();

				if (sett.densDep) {
					pPatch = pCell->getPatch();
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
								density = pPop->totalPop();
							}
							if (localK > 0.0) {

								if (settletype.indVar) settDensDep = inds[i]->getIndSettTraits();
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
							inds[i]->setSettPatch(settle);
						}
						else if (settle.settleStatus == 2) { // previously allowed to settle
							densdepOK = true;
						}
					}
				} else { // no density-dependent settlement
					densdepOK = true;
					settle.settleStatus = 2;
					settle.pSettPatch = pPatch;
					inds[i]->setSettPatch(settle);
				}

				// Update individual status
				if (mateOK && densdepOK) { // can recruit to patch
					ind.status = settled;
					nbDispersers--;
				} else { // does not recruit
					if (trfr.usesMovtProc) {
						ind.status = dispersing; // continue dispersing, 
						// unless maximum steps has been exceeded
						pathSteps steps = inds[i]->getSteps();
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

			inds[i]->setStatus(ind.status);
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

			pCell = inds[i]->getLocn(1);
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
						pPatch = pCell->getPatch();
						if (pPatch != nullptr) { // not no-data area

							// Check whether patch is suitable
							if (pPatch->getPatchNum() > 0 // not the matrix
								&& pPatch != inds[i]->getNatalPatch()// or natal patch
								&& pPatch->getK() > 0.0) { // suitable

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
				inds[i]->moveTo(destCell[0]);
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

	pPatch = pCell->getPatch();
	if (pPatch != nullptr) {
		if (pPatch->getPatchNum() > 0 // not the matrix patch
			&& pPatch->getK() > 0.0) { // suitable
			 
				pNewPopn = pPatch->getPop();
				if (pNewPopn != nullptr) {
					for (int stg = 0; stg < nStages; stg++) {
						if (pNewPopn->nInds[stg][othersex] > 0) 
							return true; // one is enough
					}
				}
				// If empty, check for incoming settlers
				if (pPatch->getPossSettlers(pSpecies, othersex) > 0) 
					return true;
		}
	}
	return false; // no mates? :(
}

//---------------------------------------------------------------------------
// Determine survival and development and record in individual's status code
// Changes are NOT applied to the Population at this stage

// FOR MULTIPLE SPECIES, MAY NEED TO SEPARATE OUT THIS IDENTIFICATION STAGE,
// SO THAT IT CAN BE PERFORMED FOR ALL SPECIES BEFORE ANY UPDATING OF POPULATIONS

void Population::drawSurvivalDevlpt(bool resolveJuvs, bool resolveAdults, bool resolveDev, bool resolveSurv)
{
	// option0:	0 - stage 0 (juveniles) only
	//			1 - all stages
	//			2 - stage 1 and above (all non-juveniles)
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

	if (dem.stageStruct && localK > 0.0) {

		for (int stg = 0; stg < nStages; stg++) {
			for (int sex = 0; sex < nsexes; sex++) {

				// Calculate development density-dependence
				if (resolveDev && sstruct.devDens && stg > 0) {
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
								else {
									weight = pSpecies->getDDwtDev(stg, effStg);
								}
								density += nInds[effStg][effSex] * weight;

							}
						}
					}
					else density = totalPop(); // no stage-dependence

					dev[stg][sex] *= exp(-(ddparams.devCoeff * density) / localK);

				}

				// Calculate survival density-dependence
				if (resolveSurv && sstruct.survDens) {
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
								else {
									weight = pSpecies->getDDwtSurv(stg, effStg);
								}
								density += nInds[effStg][effSex] * weight;
							}
						}
					}
					else // not stage-specific
						density = totalPop();
					
					surv[stg][sex] *= exp(-(ddparams.survCoeff * density) / localK);
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

//---------------------------------------------------------------------------
// Remove zero pointers to dead or dispersed individuals
void Population::clean()
{
	int ninds = inds.size();
	if (ninds > 0) {
			// ALTERNATIVE METHOD: AVOIDS SLOW SORTING OF POPULATION
		std::vector <Individual*> survivors; // all surviving individuals
		for (int i = 0; i < ninds; i++) {
			if (inds[i] != NULL) {
				survivors.push_back(inds[i]);
			}
		}
		inds.clear();
		inds = survivors;

#if RS_RCPP || NDEBUG
		// do not randomise individuals in DEBUG mode, as the function uses rand()
		// and therefore the randomisation will differ between identical runs of RS
		shuffle(inds.begin(), inds.end(), pRandom->getRNG());
#endif
	}
}

//---------------------------------------------------------------------------
// Write record to population file
void Population::outPopulation(ofstream& outPopOfs, int rep, int yr, int gen, bool envLocal, float epsGlobal,
	bool patchModel, bool writeEnv, bool gradK)
{
	Cell* pCell;
	demogrParams dem = pSpecies->getDemogrParams();

	popStats p;
	outPopOfs << rep << "\t" << yr << "\t" << gen;
	if (patchModel) {
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
			if (pCell != 0) eps = pCell->getEps();
		}
		if (pPatch->getPatchNum() == 0) { // matrix
			outPopOfs << "\t0\t0\t0";
		}
		else {
			float k = pPatch->getK();
			float envval = 0.0;
			pCell = pPatch->getRandomCell();
			if (pCell != 0) envval = pCell->getEnvVal();
			outPopOfs << "\t" << eps << "\t" << envval << "\t" << k;
		}
	}
	outPopOfs << "\t" << pSpecies->getSpNum();
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
		outPopOfs << "\t" << totalPop();
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
	short spNum = pSpecies->getSpNum();
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
			outIndsOfs << "\t" << spNum << "\t" << inds[i]->getId();
			if (dem.stageStruct) outIndsOfs << "\t" << to_string(ind.status);
			else { // non-structured population
				outIndsOfs << "\t" << to_string(ind.status);
			}
			pCell = inds[i]->getLocn(1);
			locn loc = locn();
			if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
			else loc = pCell->getLocn();
			pCell = inds[i]->getLocn(0);
			locn natalloc = pCell->getLocn();
			if (ppLand.patchModel) {
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

void Population::outputTraitPatchInfo(ofstream& outtraits, int rep, int yr, int gen, bool patchModel)
{
	if (pPatch->getK() > 0.0 && this->getNInds() > 0) {
		outtraits << rep << "\t" << yr << "\t" << gen;
		if (patchModel) {
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
	float localK = pPatch->getK();
	if (localK > 0.0 && this->getNInds() > 0) {

		pSpecies = this->getSpecies();
		demogrParams dem = pSpecies->getDemogrParams();
		emigRules emig = pSpecies->getEmigRules();
		transferRules trfr = pSpecies->getTransferRules();
		settleType sett = pSpecies->getSettle();

		indTraitsSums = this->getIndTraitsSums(pSpecies);

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
