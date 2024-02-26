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
//---------------------------------------------------------------------------

ofstream outPop;
ofstream outInds;

//---------------------------------------------------------------------------

Population::Population(void) {
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
	std::vector <double> ageProb; // for quasi-equilibrium initial age distribution
	Cell* pCell;

	if (ninds > 0) {
		inds.reserve(ninds);
		juvs.reserve(ninds);
	}

	pSpecies = pSp;
	pPatch = pPch;
	// record the new population in the patch
	patchPopn pp = patchPopn();
	pp.pSp = (intptr)pSpecies; pp.pPop = (intptr)this;
	pPatch->addPopn(pp);
#if RSDEBUG
	//DEBUGLOG << "Population::Population(): this=" << this
	//	<< " added population to patch " << endl;
#endif

	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	//trfrSMSTraits sms = pSpecies->getSMSTraits();
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
	for (int stg = 0; stg < maxNbStages; stg++) {
		for (int sex = 0; sex < maxNbSexes; sex++) {
			nInds[stg][sex] = 0;
		}
	}

	// set up local copy of minimum age table
	short minAge[maxNbStages][maxNbSexes];
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
		//	for (int sex = 0; sex < nSexes; sex++) {
		//		if (n < nSexes) n = nSexes; // to ensure at least one individual of each age is created
		//		subPops.push_back(new SubPop(loc,stg,sex,n/nSexes));
		//	}
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
#if RSDEBUG
			// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
			inds.push_back(new Individual(pCell, pPatch, stg, age, sstruct.repInterval,
				probmale, true, trfr.moveType));
#else
			inds.push_back(new Individual(pCell, pPatch, stg, age, sstruct.repInterval,
				probmale, trfr.moveModel, trfr.moveType));
#endif
			sex = inds[nindivs + i]->getSex();
			if (pSpecies->getNTraits() > 0) {
				// individual variation - set up genetics
				inds[nindivs + i]->setUpGenes(pSpecies, resol);
			}
			nInds[stg][sex]++;
		}
	}
}

Population::~Population(void) {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != NULL) delete inds[i];
	}
	inds.clear();
	int njuvs = (int)juvs.size();
	for (int i = 0; i < njuvs; i++) {
		if (juvs[i] != NULL) delete juvs[i];
	}
	juvs.clear();
}

traitsums Population::getIndTraitsSums(Species* pSpecies) {
	int g;
	traitsums ts = traitsums();
	for (int sex = 0; sex < maxNbSexes; sex++) {
		ts.ninds[sex] = 0;
		ts.sumD0[sex] = ts.ssqD0[sex] = 0.0;
		ts.sumAlpha[sex] = ts.ssqAlpha[sex] = 0.0; 
		ts.sumBeta[sex] = ts.ssqBeta[sex] = 0.0;
		ts.sumDist1[sex] = ts.ssqDist1[sex] = 0.0;
		ts.sumDist2[sex] = ts.ssqDist2[sex] = 0.0;
		ts.sumProp1[sex] = ts.ssqProp1[sex] = 0.0;
		ts.sumDP[sex] = ts.ssqDP[sex] = 0.0;
		ts.sumGB[sex] = ts.ssqGB[sex] = 0.0;
		ts.sumAlphaDB[sex] = ts.ssqAlphaDB[sex] = 0.0;
		ts.sumBetaDB[sex] = ts.ssqBetaDB[sex] = 0.0;
		ts.sumStepL[sex] = ts.ssqStepL[sex] = 0.0; 
		ts.sumRho[sex] = ts.ssqRho[sex] = 0.0;
		ts.sumS0[sex] = ts.ssqS0[sex] = 0.0;
		ts.sumAlphaS[sex] = ts.ssqAlphaS[sex] = 0.0;
		ts.sumBetaS[sex] = ts.ssqBetaS[sex] = 0.0;
		ts.sumFitness[sex] = ts.ssqFitness[sex] = 0.0;
	}

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
		emigTraits e = inds[iInd]->getEmigTraits();
		if (emig.sexDep) g = sex; 
		else g = 0;
		ts.sumD0[g] += e.d0;    
		ts.ssqD0[g] += e.d0 * e.d0;
		ts.sumAlpha[g] += e.alpha; 
		ts.ssqAlpha[g] += e.alpha * e.alpha;
		ts.sumBeta[g] += e.beta;  
		ts.ssqBeta[g] += e.beta * e.beta;

		// transfer traits
		if (trfr.moveModel) {

			switch (trfr.moveType) {

			case 1: // SMS
			{
				trfrSMSTraits sms = inds[iInd]->getSMSTraits();
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
				trfrCRWTraits c = inds[iInd]->getCRWTraits();
				g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
				ts.sumStepL[g] += c.stepLength;
				ts.ssqStepL[g] += c.stepLength * c.stepLength;
				ts.sumRho[g] += c.rho;       
				ts.ssqRho[g] += c.rho * c.rho;
				break;
			}
			default:
				throw runtime_error("moveModel enabled but moveType is neither 1 (SMS) or 2 (CRW).");
				break;
			}
		}
		else {
			trfrKernelParams k = inds[iInd]->getKernTraits();
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
		settleTraits s = inds[iInd]->getSettTraits();
		if (sett.sexDep) g = sex; 
		else g = 0;
		//	g = 0; // CURRENTLY INDIVIDUAL VARIATION CANNOT BE SEX-DEPENDENT
		ts.sumS0[g] += s.s0;     
		ts.ssqS0[g] += s.s0 * s.s0;
		ts.sumAlphaS[g] += s.alpha; 
		ts.ssqAlphaS[g] += s.alpha * s.alpha;
		ts.sumBetaS[g] += s.beta;   
		ts.ssqBetaS[g] += s.beta * s.beta;

		if (maxNbSexes > 1) g = sex; 
		else g = 0;
		ts.sumFitness[g] += inds[iInd]->getFitness();
		ts.ssqFitness[g] += inds[iInd]->getFitness() * inds[iInd]->getFitness();
	}
	return ts;
}

int Population::getNInds(void) { return (int)inds.size(); }

// ----------------------------------------------------------------------------------------
// reset allele table
// ----------------------------------------------------------------------------------------
void Population::resetAlleleTable() {
	for (auto& entry : alleleTable) {
		entry.reset();
	}
}

// ----------------------------------------------------------------------------------------
//  allele frequency in population of sampled individuals 
// ----------------------------------------------------------------------------------------

void Population::updateAlleleTable() {

	const int nLoci = pSpecies->getNPositionsForTrait(SNP);
	const int nAlleles = (int)pSpecies->getSpTrait(SNP)->getMutationParameters().find(MAX)->second;
	const auto& positions = pSpecies->getSpTrait(SNP)->getPositions();

	if (alleleTable.size() != 0)
		resetAlleleTable();
	else {
		alleleTable.reserve(nLoci);

		for (int l = 0; l < nLoci; l++) {
			alleleTable.push_back(NeutralData(nAlleles));
		}
	}

	for (Individual* individual : sampledInds) {

		const auto trait = individual->getTrait(SNP);

		int lociCounter = 0;
		for (auto position : positions) {

			auto a = (int)trait->getSelectionCoefAtLoci(0, position);
			auto b = (int)trait->getSelectionCoefAtLoci(1, position);

			int isHetero = a != b;
			alleleTable[lociCounter].incrementHeteroBy(isHetero, a);
			alleleTable[lociCounter].incrementHeteroBy(isHetero, b);

			alleleTable[lociCounter].incrementCount(a);
			alleleTable[lociCounter].incrementCount(b);

			lociCounter++;
		}

	}

	if (sampledInds.size() > 0) {
		std::for_each(alleleTable.begin(),
			alleleTable.end(),
			[&](NeutralData& v) -> void {
				v.setFrequencies(static_cast<int>(sampledInds.size()) * 2);
				//v->divideHeteros(sampledInds.size()); //weir and cockerham doesn't need this division??
			});
	}
}

double Population::getAlleleFrequency(int locus, int allele) {
	return alleleTable[locus].getFrequency(allele);
}


int Population::getAlleleCount(int locus, int allele) {
	return alleleTable[locus].getCount(allele);
}

double Population::getHetero(int locus, int allele) {
	return alleleTable[locus].getHetero(allele);
}

// ----------------------------------------------------------------------------------------
// Count number of heterozygotes loci in sampled individuals
// ----------------------------------------------------------------------------------------

int Population::countHeterozygoteLoci() {
	int hetero = 0;
	for (Individual* ind : sampledInds) {
		const auto trait = ind->getTrait(SNP);
		hetero += trait->countHeterozygoteLoci();
	}
	return hetero;
}

// ----------------------------------------------------------------------------------------
// Count number of heterozygotes per loci loci in sampled individuals
// ----------------------------------------------------------------------------------------

vector<double> Population::countLociHeterozyotes() {
	const auto& positions = pSpecies->getSpTrait(SNP)->getPositions();
	vector<double> hetero(positions.size(), 0);

	for (Individual* ind : sampledInds) {
		const auto trait = ind->getTrait(SNP);
		int counter = 0;
		for (auto position : positions) {
			hetero[counter] += trait->isHeterozygoteAtLocus(position);
			counter++;
		}
	}
	return hetero;
}

// ----------------------------------------------------------------------------------------
//	compute the expected heterozygosity for population
// ----------------------------------------------------------------------------------------

double Population::computeHs() {
	int nLoci = pSpecies->getNPositionsForTrait(SNP);
	int nAlleles = (int)pSpecies->getSpTrait(SNP)->getInitialParameters().find(MAX)->second;
	double hs = 0;
	double freq;

	vector<double>locihet(nLoci, 1);

	if (sampledInds.size() > 0) {
		for (int thisLocus = 0; thisLocus < nLoci; ++thisLocus) {
			for (int allele = 0; allele < nAlleles; ++allele) {
				freq = getAlleleFrequency(thisLocus, allele);
				freq *= freq; //squared frequencies (expected _homozygosity)
				locihet[thisLocus] -= freq; //1 - sum of p2 = expected heterozygosity
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
#if RSDEBUG
	//DEBUGLOG << "Population::getStats(): this=" << this
	////	<< " p.pSpecies=" << p.pSpecies << " p.spNum=" << p.spNum
	//	<< " p.pPatch=" << p.pPatch << " patchNum=" << p.pPatch->getPatchNum()
	//	<< " nStages=" << nStages << " nSexes=" << nSexes << " p.nInds=" << p.nInds
	//	<< endl;
#endif
	for (int stg = 1; stg < nStages; stg++) {
		for (int sex = 0; sex < nSexes; sex++) {
			ninds = nInds[stg][sex];
			p.nNonJuvs += ninds;
#if RSDEBUG
			//DEBUGLOG << "Population::getStats(): this=" << this
			//	<< " stg=" << stg << " sex=" << sex
			//	<< " nInds[stg][sex]=" << nInds[stg][sex] << " p.nNonJuvs=" << p.nNonJuvs
			//	<< endl;
#endif
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
#if RSDEBUG
	//DEBUGLOG << "Population::getStats(): this=" << this
	//	<< " p.nInds=" << p.nInds << " p.nAdults=" << p.nAdults << " p.nNonJuvs=" << p.nNonJuvs
	//	<< " breeders[0]=" << breeders[0] << " breeders[1]=" << breeders[1]
	//	<< " p.breeding=" << p.breeding
	//	<< endl;
#endif
	return p;
}

Species* Population::getSpecies(void) { return pSpecies; }

int Population::totalPop(void) {
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
// Remove all Individuals
void Population::extirpate(void) {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		if (inds[i] != NULL) delete inds[i];
	}
	inds.clear();
	int njuvs = (int)juvs.size();
	for (int i = 0; i < njuvs; i++) {
		if (juvs[i] != NULL) delete juvs[i];
	}
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
	int ninds = (int)inds.size();
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//	<< " ninds=" << ninds
//	<< endl;
#endif // RSDEBUG 
	if (ninds == 0) return;

	int nsexes, stage, sex, njuvs, nj, nmales, nfemales;
	Cell* pCell;
	indStats ind;
	double expected;
	bool skipbreeding;

	//envGradParams grad = paramsGrad->getGradient();
	envStochParams env = paramsStoch->getStoch();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();

	if (dem.repType == 0) nsexes = 1; else nsexes = 2;

#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): this=" << this
//	<< " pSpecies=" << pSpecies
//	<< " localK=" << localK << " envval=" << envval << " resol=" << resol
//	<< " sstruct.nStages=" << sstruct.nStages << " nsexes=" << nsexes << " ninds=" << ninds
//	<< endl;
#endif

// set up local copy of species fecundity table
	float fec[maxNbStages][maxNbSexes];
	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use fecundity recorded for females
					fec[stg][sex] = pSpecies->getFec(stg, 0);
				}
				else
					fec[stg][sex] = pSpecies->getFec(stg, sex);
				//			if (sex == 0 && fec[stg][sex] > dem.lambda) dem.lambda = fec[stg][sex];
			}
			else { // non-structured population
				if (stg == 1) fec[stg][sex] = dem.lambda; // adults
				else fec[stg][sex] = 0.0; // juveniles
			}
#if RSDEBUG
			//if (ninds > 0) {
			//DEBUGLOG << "Population::reproduction(): fec[" << stg << "][" << sex << "] = " << fec[stg][sex]
			//	<< endl;
			//}
#endif
		}
	}

	if (dem.stageStruct) {
#if RSDEBUG
		//if (ninds > 0) {
		//	DEBUGLOG << "Population::reproduction(): ninds=" << ninds << " localK=" << localK
		//		<< " effect of density dependence:" << endl;
		//}
#endif
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
#if RSDEBUG
								//if (ninds > 0) {
								//	DEBUGLOG << " effstg=" << effstg << " effsex=" << effsex << " nInds=" << nInds[effstg][effsex];
								//	DEBUGLOG << " weight=" << weight << " effect=" << effect
								//		<< endl;
								//}
#endif
							}
						}
					}
					else // not stage-specific
						effect = (float)totalPop();
					if (localK > 0.0) fec[stg][0] *= exp(-effect / localK);
#if RSDEBUG
					//if (ninds > 0) {
					//	DEBUGLOG << " eff popn=" << effect << " exponential=" << exp(-effect/localK);
					//	DEBUGLOG << " fec[" << stg << "][0]=" << fec[stg][0]
					//		<< endl;
					//}
#endif
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
			//#if GOBYMODEL
			//		ddeffect[1] = (float)ninds/localK;
			//#else
			if (dem.repType == 1 || dem.repType == 2) { // sexual model
				// apply factor of 2 (as in manual, eqn. 6)
				fec[1][0] *= 2.0;
			}
			fec[1][0] /= (1.0f + fabs(dem.lambda - 1.0f) * pow(((float)ninds / localK), dem.bc));
			//#endif
		}
#if RSDEBUG
		//DEBUGLOG << "Population::reproduction(): dem.lambda=" << dem.lambda << " ninds=" << ninds
		//	<< " localK=" << localK << " dem.bc=" << dem.bc << " fec[1][0]=" << fec[1][0]
		//	<< endl;
#endif
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
					nj = (int)juvs.size();
					pCell = pPatch->getRandomCell();
					for (int j = 0; j < njuvs; j++) {

						Individual* newJuv;
#if RSDEBUG
						// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
						newJuv = new Individual(pCell, pPatch, 0, 0, 0, dem.propMales, true, trfr.moveType);
#else
						newJuv = new Individual(pCell, pPatch, 0, 0, 0, dem.propMales, trfr.moveModel, trfr.moveType);
#endif

						if (pSpecies->getNTraits() > 0) {
							newJuv->inheritTraits(pSpecies, inds[i], father, resol);
						}

						if (newJuv->getFitness() < pRandom->Random()) {
							delete newJuv;
						}
						else {
							juvs.push_back(newJuv);
							nInds[0][0]++;
						}
					}
				}
			}
		}
		break;

	case 1: // simple sexual model
	case 2: // complex sexual model
		// count breeding females and males
		// add breeding males to list of potential fathers
#if RSDEBUG
//DEBUGLOG << "Population::reproduction(): case 1:"
//	<< " fec[1][0]=" << fec[1][0] << " fec[1][1]=" << fec[1][1]
//	<< endl;
#endif
		nfemales = nmales = 0;
		for (int i = 0; i < ninds; i++) {
			ind = inds[i]->getStats();
			if (ind.sex == 0 && fec[ind.stage][0] > 0.0) nfemales++;
			if (ind.sex == 1 && fec[ind.stage][1] > 0.0) {
				fathers.push_back(inds[i]);
				nmales++;
#if RSDEBUG
				//DEBUGLOG << "Population::reproduction(): i=" << i << " nmales=" << nmales
				//	<< " inds[i]=" << inds[i] << endl;
#endif
			}
		}
#if RSDEBUG
		//DEBUGLOG << "Population::reproduction(): breeding nfemales=" << nfemales
		//	<< " breeding nmales=" << nmales << endl;
#endif
		if (nfemales > 0 && nmales > 0)
		{ // population can breed
			if (dem.repType == 2) { // complex sexual model
				// calculate proportion of eligible females which breed
				propBreed = (2.0 * dem.harem * nmales) / (nfemales + dem.harem * nmales);
				if (propBreed > 1.0) propBreed = 1.0;
#if RSDEBUG
				//DEBUGLOG << "Population::reproduction(): harem=" << dem.harem
				//	<< " nfemales=" << nfemales << " nmales=" << nmales << " propBreed=" << propBreed
				//	<< endl;
#endif
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
#if RSDEBUG
							//DEBUGLOG << "Population::reproduction(): THIS LINE SHOULD NOT APPEAR FOR GOBY MODEL"
							//	<< " expected=" << expected
							//	<< endl;
#endif
						}
						else expected = 0.0; // fails to breed
						if (expected <= 0.0) njuvs = 0;
						else njuvs = pRandom->Poisson(expected);
#if RSDEBUG
						//DEBUGLOG << "Population::reproduction():"
						//	<< " i " << i << " ID=" << inds[i]->getId() << " stage=" << stage
						//	<< " expected=" << expected << " njuvs=" << njuvs << endl;
#endif
						if (njuvs > 0)
						{
							nj = (int)juvs.size();
							// select father at random from breeding males ...
							int rrr = 0;
							if (nmales > 1) rrr = pRandom->IRandom(0, nmales - 1);
							father = fathers[rrr];
							pCell = pPatch->getRandomCell();
							for (int j = 0; j < njuvs; j++) {

								Individual* newJuv;
#if RSDEBUG
								// NOTE: CURRENTLY SETTING ALL INDIVIDUALS TO RECORD NO. OF STEPS ...
								newJuv = new Individual(pCell, pPatch, 0, 0, 0, dem.propMales, true, trfr.moveType);
#else
								newJuv = new Individual(pCell, pPatch, 0, 0, 0, dem.propMales, trfr.moveModel, trfr.moveType);
#endif
								if (pSpecies->getNTraits() > 0) {
									newJuv->inheritTraits(pSpecies, inds[i], father, resol);
								}
								if (newJuv->getFitness() < pRandom->Random()) {
									delete newJuv;
								}
								else {
									juvs.push_back(newJuv);
									sex = newJuv->getSex();
									nInds[0][sex]++;
								}
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
void Population::fledge(void)
{
#if RSDEBUG
	//DEBUGLOG << "Population::fledge(): this=" << this
	//	<< " ninds=" << (int)inds.size()
	//	<< " njuvs=" << (int)juvs.size()
	//	<< endl;
#endif
	demogrParams dem = pSpecies->getDemogrParams();

	if (dem.stageStruct) { // juveniles are added to the individuals vector
		inds.insert(inds.end(), juvs.begin(), juvs.end());
		//	nInds += nJuvs; nJuvs = 0;
	}
	else { // all adults die and juveniles replace adults
		int ninds = (int)inds.size();
		for (int i = 0; i < ninds; i++) {
			delete inds[i];
		}
		inds.clear();
		for (int sex = 0; sex < nSexes; sex++) {
			nInds[1][sex] = 0; // set count of adults to zero
		}
		inds = juvs;
	}
	juvs.clear();

}

Individual* Population::sampleInd() const {
	int index = pRandom->IRandom(0, static_cast<int>(inds.size() - 1));
	return inds[index];
}

void Population::sampleIndsWithoutReplacement(string n, const set<int>& sampleStages) {

	sampledInds.clear();
	auto rng = pRandom->getRNG();
	set<Individual*> stagedInds;

	for (int stage : sampleStages) {
		auto sInds = getIndividualsInStage(stage);
		stagedInds.insert(sInds.begin(), sInds.end());
	}

	if (n == "all" || stagedInds.size() < stoi(n))
		// Sample all individuals in selected stages
		sampledInds = stagedInds;
	else {
		vector<Individual*> out;
		// Sample n individuals across selected stages
		sample(stagedInds.begin(), stagedInds.end(), std::back_inserter(out), stoi(n), rng);
		std::copy(out.begin(), out.end(), std::inserter(sampledInds, sampledInds.end()));
	}
}

int Population::sampleSize() const {
	return static_cast<int>(sampledInds.size());
}

set<Individual*> Population::getIndividualsInStage(int stage) {
	set<Individual*> filteredInds;
	for (auto ind : inds) {
		if (ind->getStats().stage == stage)
			filteredInds.insert(ind);
	}
	return filteredInds;
}

// Determine which individuals will disperse
void Population::emigration(float localK)
{
	int nsexes;
	double disp, Pdisp, NK;
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	emigRules emig = pSpecies->getEmigRules();
	emigTraits eparams;
	indStats ind;

	// to avoid division by zero, assume carrying capacity is at least one individual
	// localK can be zero if there is a moving gradient or stochasticity in K
	if (localK < 1.0) localK = 1.0;
	NK = static_cast<float>(totalPop()) / localK;

	int ninds = static_cast<int>(inds.size());

	// set up local copy of emigration probability table
	// used when there is no individual variability
	// NB - IT IS DOUBTFUL THIS CONTRIBUTES ANY SUBSTANTIAL TIME SAVING
	if (dem.repType == 0) nsexes = 1; 
	else nsexes = 2;
	double Pemig[maxNbStages][maxNbSexes];

	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (emig.indVar) Pemig[stg][sex] = 0.0;
			else {
				if (emig.densDep) {
					if (emig.sexDep) {
						if (emig.stgDep) {
							eparams = pSpecies->getEmigTraits(stg, sex);
						}
						else {
							eparams = pSpecies->getEmigTraits(0, sex);
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							eparams = pSpecies->getEmigTraits(stg, 0);
						}
						else {
							eparams = pSpecies->getEmigTraits(0, 0);
						}
					}
					Pemig[stg][sex] = eparams.d0 / (1.0 + exp(-(NK - eparams.beta) * eparams.alpha));
				}
				else { // density-independent
					if (emig.sexDep) {
						if (emig.stgDep) {
							Pemig[stg][sex] = pSpecies->getEmigD0(stg, sex);
						}
						else { // !emig.stgDep
							Pemig[stg][sex] = pSpecies->getEmigD0(0, sex);
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							Pemig[stg][sex] = pSpecies->getEmigD0(stg, 0);
						}
						else { // !emig.stgDep
							Pemig[stg][sex] = pSpecies->getEmigD0(0, 0);
						}
					}
				}
			} // end of !emig.indVar
		}
	}

	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();
		if (ind.status < 1) {
			if (emig.indVar) { // individual variability in emigration
				if (dem.stageStruct && ind.stage != emig.emigStage) {
					// emigration may not occur
					Pdisp = 0.0;
				}
				else { // non-structured or individual is in emigration stage
					eparams = inds[i]->getEmigTraits();
					if (emig.densDep) { // density-dependent
						NK = (float)totalPop() / localK;
						Pdisp = eparams.d0 / (1.0 + exp(-(NK - eparams.beta) * eparams.alpha));
					}
					else { // density-independent
						if (emig.sexDep) {
							Pdisp = Pemig[0][ind.sex] + eparams.d0;
						}
						else {
							Pdisp = Pemig[0][0] + eparams.d0;
						}
					}
				}
			} // end of individual variability
			else { // no individual variability

				if (emig.densDep) {
					if (emig.sexDep) {
						if (emig.stgDep) {
							Pdisp = Pemig[ind.stage][ind.sex];
						}
						else {
							Pdisp = Pemig[0][ind.sex];
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							Pdisp = Pemig[ind.stage][0];
						}
						else {
							Pdisp = Pemig[0][0];
						}
					}
				}
				else { // density-independent
					if (emig.sexDep) {
						if (emig.stgDep) {
							Pdisp = Pemig[ind.stage][ind.sex];
						}
						else { // !emig.stgDep
							Pdisp = Pemig[0][ind.sex];
						}
					}
					else { // !emig.sexDep
						if (emig.stgDep) {
							Pdisp = Pemig[ind.stage][0];
						}
						else { // !emig.stgDep
							Pdisp = Pemig[0][0];
						}
					}
				}
			} // end of no individual variability

			disp = pRandom->Bernoulli(Pdisp);

			if (disp == 1) { // emigrant
				inds[i]->setStatus(1);
			}
		} // end of if (ind.status < 1) condition
	} // end of for loop
}

// All individuals emigrate after patch destruction
void Population::allEmigrate(void) {
	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		inds[i]->setStatus(1);
	}
}

// If an Individual has been identified as an emigrant, remove it from the Population
disperser Population::extractDisperser(int ix) {
	disperser d = disperser();
	indStats ind = inds[ix]->getStats();
	if (ind.status == 1) { // emigrant
		d.pInd = inds[ix]; d.yes = true;
		inds[ix] = 0;
		nInds[ind.stage][ind.sex]--;
	}
	else {
		d.pInd = NULL; d.yes = false;
	}
	return d;
}


// For an individual identified as being in the matrix population:
// if it is a settler, return its new location and remove it from the current population
// otherwise, leave it in the matrix population for possible reporting before deletion
disperser Population::extractSettler(int ix) {
	disperser d = disperser();
	Cell* pCell;
//Patch* pPatch;

	indStats ind = inds[ix]->getStats();
	pCell = inds[ix]->getLocn(1);
	d.pInd = inds[ix];  d.pCell = pCell; d.yes = false;
	if (ind.status == 4 || ind.status == 5) { // settled
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
#if RS_RCPP // included also SEASONAL
int Population::transfer(Landscape* pLandscape, short landIx, short nextseason)
#else
int Population::transfer(Landscape* pLandscape, short landIx)
#endif
{
	int ndispersers = 0;
	int disperser;
	short othersex;
	bool mateOK, densdepOK;
	intptr patch, popn;
	int patchnum;
	double localK, popsize, settprob;
	Patch* pPatch = 0;
	Cell* pCell = 0;
	indStats ind;
	Population* pNewPopn = 0;
	locn newloc = locn();
	locn nbrloc = locn();

	landData ppLand = pLandscape->getLandData();
	short reptype = pSpecies->getRepType();
	transferRules trfr = pSpecies->getTransferRules();
	settleType settletype = pSpecies->getSettle();
	settleRules sett;
	settleTraits settDD;
	settlePatch settle;
	simParams sim = paramsSim->getSim();
	// each individual takes one step
	// for dispersal by kernel, this should be the only step taken
	int ninds = (int)inds.size();
#if RSDEBUG
	//DEBUGLOG << "Population::transfer(): 0000: ninds = " << ninds
	//	<< " ndispersers = " << ndispersers << endl;
#endif
	for (int i = 0; i < ninds; i++) {
#if RSDEBUG
		//DEBUGLOG << "Population::transfer(): 1111: i = " << i << " ID = " << inds[i]->getId()
		//	<< endl;
#endif
		if (trfr.moveModel) {
#if RSDEBUG
			//pCell = inds[i]->getLocn(1);
			//locn loc = pCell->getLocn();
			//DEBUGLOG << "Population::transfer(): 1112: i = " << i << " ID = " << inds[i]->getId()
			//	<< " before:" << " x = " << loc.x << " y = " << loc.y
			//	<< endl;
#endif
			disperser = inds[i]->moveStep(pLandscape, pSpecies, landIx, sim.absorbing);
#if RSDEBUG
			//pCell = inds[i]->getLocn(1);
			//newloc = pCell->getLocn();
			//DEBUGLOG << "Population::transfer(): 1113: i = " << i << " ID = " << inds[i]->getId()
			//	<< " after: " << " x = " << newloc.x << " y = " << newloc.y
			//	<< endl;
#endif
		}
		else {
			disperser = inds[i]->moveKernel(pLandscape, pSpecies, reptype, sim.absorbing);
		}
		ndispersers += disperser;
		if (disperser) {
			if (reptype > 0)
			{ // sexual species - record as potential settler in new patch
				if (inds[i]->getStatus() == 2)
				{ // disperser has found a patch
					pCell = inds[i]->getLocn(1);
					patch = pCell->getPatch();
					if (patch != 0) { // not no-data area
						pPatch = (Patch*)patch;
						pPatch->incrPossSettler(pSpecies, inds[i]->getSex());
					}
				}
			}
		}
	}
#if RSDEBUG
	//DEBUGLOG << "Population::transfer(): 5555: ninds=" << ninds
	//	<< " ndispersers=" << ndispersers << endl;
#endif

	for (int i = 0; i < ninds; i++) {
		ind = inds[i]->getStats();
		if (ind.sex == 0) othersex = 1; else othersex = 0;
		if (settletype.stgDep) {
			if (settletype.sexDep) sett = pSpecies->getSettRules(ind.stage, ind.sex);
			else sett = pSpecies->getSettRules(ind.stage, 0);
		}
		else {
			if (settletype.sexDep) sett = pSpecies->getSettRules(0, ind.sex);
			else sett = pSpecies->getSettRules(0, 0);
		}
		if (ind.status == 2)
		{ // awaiting settlement
			pCell = inds[i]->getLocn(1);
			if (pCell == 0) {
				// this condition can occur in a patch-based model at the time of a dynamic landscape
				// change when there is a range restriction in place, since a patch can straddle the
				// range restriction and an individual forced to disperse upon patch removal could
				// start its trajectory beyond the boundary of the restrictyed range - such a model is 
				// not good practice, but the condition must be handled by killing the individual conceerned
				ind.status = 6;
			}
			else {

#if RSDEBUG
				//newloc = pCell->getLocn();
				//DEBUGLOG << "Population::transfer(): 6666: i=" << i << " ID=" << inds[i]->getId()
				//	<< " sex=" << ind.sex << " status=" << ind.status
				//	<< " pCell=" << pCell << " x=" << newloc.x << " y=" << newloc.y
				//	<< " findMate=" << sett.findMate
				////	<< " wait=" << sett.wait
				////	<< " go2nbrLocn=" << sett.go2nbrLocn
				//	<< endl;
#endif

				mateOK = false;
				if (sett.findMate) {
					// determine whether at least one individual of the opposite sex is present in the
					// new population
					if (matePresent(pCell, othersex)) mateOK = true;
#if RSDEBUG
					//DEBUGLOG << "Population::transfer(): 7777: othersex=" << othersex
					//	<< " this=" << this << " pNewPopn=" << pNewPopn << " popsize=" << popsize << " mateOK=" << mateOK
					//	<< endl;
#endif
				}
				else { // no requirement to find a mate
					mateOK = true;
				}
				densdepOK = false;
				settle = inds[i]->getSettPatch();
				if (sett.densDep)
				{
					patch = pCell->getPatch();
#if RSDEBUG
					//DEBUGLOG << "Population::transfer(): 8880: i=" << i << " patch=" << patch
					//	<< endl;
#endif
					if (patch != 0) { // not no-data area
						pPatch = (Patch*)patch;
						if (settle.settleStatus == 0
							|| settle.pSettPatch != pPatch)
							// note: second condition allows for having moved from one patch to another
							// adjacent one
						{
							//						inds[i]->resetPathOut(); // reset steps out of patch to zero
													// determine whether settlement occurs in the (new) patch
							localK = (double)pPatch->getK();
							popn = pPatch->getPopn((intptr)pSpecies);
#if RSDEBUG
							//DEBUGLOG << "Population::transfer(): 8881: i=" << i << " patchNum=" << pPatch->getPatchNum()
							//	<< " localK=" << localK << " popn=" << popn << endl;
#endif
							if (popn == 0) { // population has not been set up in the new patch
								popsize = 0.0;
							}
							else {
								pNewPopn = (Population*)popn;
								popsize = (double)pNewPopn->totalPop();
							}
							if (localK > 0.0) {
								// make settlement decision
								if (settletype.indVar) settDD = inds[i]->getSettTraits();
								else settDD = pSpecies->getSettTraits(ind.stage, ind.sex);
								settprob = settDD.s0 /
									(1.0 + exp(-(popsize / localK - (double)settDD.beta) * (double)settDD.alpha));
#if RSDEBUG
								//DEBUGLOG << "Population::transfer(): 8888: i=" << i << " ind.stage=" << ind.stage
								//	<< " this=" << this << " pNewPopn=" << pNewPopn << " popsize=" << popsize
								//	<< " localK=" << localK << " alpha=" << settDD.alpha << " beta=" << settDD.beta
								//	<< " settprob=" << settprob
								//	<< endl;
#endif
								if (pRandom->Bernoulli(settprob)) { // settlement allowed
									densdepOK = true;
									settle.settleStatus = 2;
								}
								else { // settlement procluded
									settle.settleStatus = 1;
								}
								settle.pSettPatch = pPatch;
							}
							inds[i]->setSettPatch(settle);
						}
						else {
							if (settle.settleStatus == 2) { // previously allowed to settle
								densdepOK = true;
							}
						}
					}
				}
				else { // no density-dependent settlement
					densdepOK = true;
					settle.settleStatus = 2;
					settle.pSettPatch = pPatch;
					inds[i]->setSettPatch(settle);
				}

				if (mateOK && densdepOK) { // can recruit to patch
					ind.status = 4;
					ndispersers--;
				}
				else { // does not recruit
					if (trfr.moveModel) {
						ind.status = 1; // continue dispersing, unless ...
						// ... maximum steps has been exceeded
						pathSteps steps = inds[i]->getSteps();
						settleSteps settsteps = pSpecies->getSteps(ind.stage, ind.sex);
						if (steps.year >= settsteps.maxStepsYr) {
							ind.status = 3; // waits until next year
						}
						if (steps.total >= settsteps.maxSteps) {
							ind.status = 6; // dies
						}
					}
					else { // dispersal kernel
						if (sett.wait) {
							ind.status = 3; // wait until next dispersal event
						}
						else {
							ind.status = 6; // (dies unless a neighbouring cell is suitable)
						}
						ndispersers--;
					}
				}
			}

			inds[i]->setStatus(ind.status);
		}
#if RS_RCPP
		// write each individuals current movement step and status to paths file
		if (trfr.moveModel && sim.outPaths) {
			if (nextseason >= sim.outStartPaths && nextseason % sim.outIntPaths == 0) {
				inds[i]->outMovePath(nextseason);
			}
		}
#endif
		if (!trfr.moveModel && sett.go2nbrLocn && (ind.status == 3 || ind.status == 6))
		{
			// for kernel-based transfer only ...
			// determine whether recruitment to a neighbouring cell is possible
#if RSDEBUG
//DEBUGLOG << "Population::transfer(): neighbour cell search: sett.go2nbrLocn = " << sett.go2nbrLocn
//	<< " ind.status = " << ind.status
//	<< endl;
#endif
			pCell = inds[i]->getLocn(1);
			newloc = pCell->getLocn();
			vector <Cell*> nbrlist;
			for (int dx = -1; dx < 2; dx++) {
				for (int dy = -1; dy < 2; dy++) {
					if (dx != 0 || dy != 0) { //cell is not the current cell
						nbrloc.x = newloc.x + dx; nbrloc.y = newloc.y + dy;
						if (nbrloc.x >= 0 && nbrloc.x <= ppLand.maxX
							&& nbrloc.y >= 0 && nbrloc.y <= ppLand.maxY) { // within landscape
							// add to list of potential neighbouring cells if suitable, etc.
							pCell = pLandscape->findCell(nbrloc.x, nbrloc.y);
							if (pCell != 0) { // not no-data area
								patch = pCell->getPatch();
								if (patch != 0) { // not no-data area
									pPatch = (Patch*)patch;
									patchnum = pPatch->getPatchNum();
									if (patchnum > 0 && pPatch != inds[i]->getNatalPatch())
									{ // not the matrix or natal patch
										if (pPatch->getK() > 0.0)
										{ // suitable
											if (sett.findMate) {
												if (matePresent(pCell, othersex)) nbrlist.push_back(pCell);
											}
											else
												nbrlist.push_back(pCell);
										}
									}
								}
							}
						}
					}
				}
			}
			int listsize = (int)nbrlist.size();
			if (listsize > 0) { // there is at least one suitable neighbouring cell
				if (listsize == 1) {
					inds[i]->moveto(nbrlist[0]);
				}
				else { // select at random from the list
					int rrr = pRandom->IRandom(0, listsize - 1);
					inds[i]->moveto(nbrlist[rrr]);
				}
			}
			// else list empty - do nothing - individual retains its current location and status
		}
	}
#if RSDEBUG
	//DEBUGLOG << "Population::transfer(): 9999: ninds = " << ninds
	//	<< " ndispersers = " << ndispersers << endl;
#endif


	return ndispersers;
}

// Determine whether there is a potential mate present in a patch which a potential
// settler has reached
bool Population::matePresent(Cell* pCell, short othersex)
{
	int patch;
	Patch* pPatch;
	Population* pNewPopn;
	int popsize = 0;
	bool matefound = false;

	patch = (int)pCell->getPatch();
	if (patch != 0) {
		pPatch = (Patch*)pCell->getPatch();
		if (pPatch->getPatchNum() > 0) { // not the matrix patch
			if (pPatch->getK() > 0.0)
			{ // suitable
				pNewPopn = (Population*)pPatch->getPopn((intptr)pSpecies);
				if (pNewPopn != 0) {
					// count members of other sex already resident in the patch
					for (int stg = 0; stg < nStages; stg++) {
						popsize += pNewPopn->nInds[stg][othersex];
					}
				}
				if (popsize < 1) {
					// add any potential settlers of the other sex
					popsize += pPatch->getPossSettlers(pSpecies, othersex);
				}
			}
		}
	}
	if (popsize > 0) matefound = true;
	return matefound;
}

//---------------------------------------------------------------------------
// Determine survival and development and record in individual's status code
// Changes are NOT applied to the Population at this stage

// FOR MULTIPLE SPECIES, MAY NEED TO SEPARATE OUT THIS IDENTIFICATION STAGE,
// SO THAT IT CAN BE PERFORMED FOR ALL SPECIES BEFORE ANY UPDATING OF POPULATIONS

void Population::survival0(float localK, short option0, short option1)
{
	// option0:	0 - stage 0 (juveniles) only
	//					1 - all stages
	//					2 - stage 1 and above (all non-juveniles)
	// option1:	0 - development only (when survival is annual)
	//	  	 		1 - development and survival
	//	  	 		2 - survival only (when survival is annual)
	densDepParams ddparams = pSpecies->getDensDep();
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();

	// get surrent population size
	int ninds = (int)inds.size();
	if (ninds == 0) return;
	// set up local copies of species development and survival tables
	int nsexes;
	if (dem.repType == 0) nsexes = 1; else nsexes = 2;
	float dev[maxNbStages][maxNbSexes];
	float surv[maxNbStages][maxNbSexes];
	short minAge[maxNbStages][maxNbSexes];
	for (int stg = 0; stg < sstruct.nStages; stg++) {
		for (int sex = 0; sex < nsexes; sex++) {
			if (dem.stageStruct) {
				if (dem.repType == 1) { // simple sexual model
					// both sexes use development and survival recorded for females
					dev[stg][sex] = pSpecies->getDev(stg, 0);
					surv[stg][sex] = pSpecies->getSurv(stg, 0);
					minAge[stg][sex] = pSpecies->getMinAge(stg, 0);
				}
				else {
					dev[stg][sex] = pSpecies->getDev(stg, sex);
					surv[stg][sex] = pSpecies->getSurv(stg, sex);
					minAge[stg][sex] = pSpecies->getMinAge(stg, sex);
				}
				if (option1 == 0) surv[stg][sex] = 1.0; // development only - all survive
				if (option1 == 2) dev[stg][sex] = 0.0;  // survival only - none develops
			}
			else { // non-structured population
				if (stg == 1) { // adults
					dev[stg][sex] = 0.0; surv[stg][sex] = 0.0; minAge[stg][sex] = 0;
				}
				else { // juveniles
					dev[stg][sex] = 1.0; surv[stg][sex] = 1.0; minAge[stg][sex] = 0;
				}
			}
#if RSDEBUG
			//DEBUGLOG << "Population::survival0(): 1111 "
			//	<< " dev[" << stg << "][" << sex << "] = " << dev[stg][sex]
			//	<< " surv[" << stg << "][" << sex << "] = " << surv[stg][sex]
			//	<< endl;
#endif
		}
	}
	if (dem.stageStruct) {
		for (int stg = 0; stg < nStages; stg++) {
			for (int sex = 0; sex < nsexes; sex++) {
				if (option1 != 2 && sstruct.devDens && stg > 0) {
				// NB DD in development does NOT apply to juveniles,
					float effect = 0.0;
					if (sstruct.devStageDens) { // stage-specific density dependence
						// NOTE: matrix entries represent effect of ROW on COLUMN 
						// AND males precede females
						float weight = 0.0;
						for (int effstg = 0; effstg < nStages; effstg++) {
							for (int effsex = 0; effsex < nSexes; effsex++) {
								if (dem.repType == 2) {
									int rowincr, colincr;
									if (effsex == 0) rowincr = 1; else rowincr = 0;
									if (sex == 0) colincr = 1; else colincr = 0;
									weight = pSpecies->getDDwtDev(2 * stg + colincr, 2 * effstg + rowincr);
								}
								else {
									weight = pSpecies->getDDwtDev(stg, effstg);
								}
								effect += (float)nInds[effstg][effsex] * weight;
#if RSDEBUG
								//	DEBUGLOG << " effstg=" << effstg << " effsex=" << effsex;
								//	DEBUGLOG << " weight=" << weight << " effect=" << effect
								//		<< endl;
#endif
							}
						}
					}
					else // not stage-specific
						effect = (float)totalPop();
					if (localK > 0.0)
						dev[stg][sex] *= exp(-(ddparams.devCoeff * effect) / localK);
#if RSDEBUG
					//DEBUGLOG << "Population::survival0(): 2288 " << " effect=" << effect;
					//if (localK > 0.0)
					//	DEBUGLOG << " exp=" << exp(-(ddparams.devCoeff*effect)/localK);
					//DEBUGLOG << " dev[" << stg << "][" << sex << "] = " << dev[stg][sex]
					//	<< endl;
#endif
				} // end of if (sstruct.devDens && stg > 0)
				if (option1 != 0 && sstruct.survDens) {
#if RSDEBUG
					//	DEBUGLOG << "DD in SURVIVAL for stg=" << stg << " sex=" << sex << endl;
#endif
					float effect = 0.0;
					if (sstruct.survStageDens) { // stage-specific density dependence
						// NOTE: matrix entries represent effect of ROW on COLUMN 
						// AND males precede females
						float weight = 0.0;
						for (int effstg = 0; effstg < nStages; effstg++) {
							for (int effsex = 0; effsex < nSexes; effsex++) {
								if (dem.repType == 2) {
									int rowincr, colincr;
									if (effsex == 0) rowincr = 1; else rowincr = 0;
									if (sex == 0) colincr = 1; else colincr = 0;
									weight = pSpecies->getDDwtSurv(2 * stg + colincr, 2 * effstg + rowincr);
								}
								else {
									weight = pSpecies->getDDwtSurv(stg, effstg);
								}
								effect += (float)nInds[effstg][effsex] * weight;
#if RSDEBUG
								//	DEBUGLOG << " effstg=" << effstg << " effsex=" << effsex;
								//	DEBUGLOG << " weight=" << weight << " effect=" << effect
								//		<< endl;
#endif
							}
						}
					}
					else // not stage-specific
						effect = (float)totalPop();
					if (localK > 0.0)
						surv[stg][sex] *= exp(-(ddparams.survCoeff * effect) / localK);
#if RSDEBUG
					//DEBUGLOG << "Population::survival0(): 3333 " << " effect=" << effect;
					//if (localK > 0.0)
					//	DEBUGLOG << " exp = " << exp(-(ddparams.survCoeff*effect)/localK);
					//DEBUGLOG << " surv[" << stg << "][" << sex << "] = " << surv[stg][sex]
					//	<< endl;
#endif
				} // end of if (sstruct.survDens)
			}
		}
	}
	// identify which individuals die or develop
#if RSDEBUG
//DEBUGLOG << "Population::survival0():"  << " ninds " << ninds
//	<< endl;
#endif
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();
#if RSDEBUG
		//DEBUGLOG << "Population::survival0():"
		//	<< " i=" << i << " indId=" << inds[i]->getId()
		//	<< " stage=" << ind.stage << " status=" << ind.status << " sex=" << ind.sex
		//#if PARTMIGRN
		//	<< " migrnstatus=" << inds[i]->getMigrnStatus()
		//#endif  // PARTMIGRN 
		//	<< endl;
#endif
		if ((ind.stage == 0 && option0 < 2) || (ind.stage > 0 && option0 > 0)) {
			// condition for processing the stage is met...
			if (ind.status < 6) { // not already doomed
				if (ind.sex < sex_t::FEM || ind.sex > sex_t::MAL)
					// ?? MSVC believes it's important to bound check ind.sex
					throw runtime_error("Individual sex is out of bounds");
				double probsurv = surv[ind.stage][ind.sex];
				// does the individual survive?
				if (pRandom->Bernoulli(probsurv)) { // survives
					// does the individual develop?
					double probdev = dev[ind.stage][ind.sex];
					if (ind.stage < nStages - 1) { // not final stage
#if RSDEBUG
						//DEBUGLOG << "Population::survival0():"
						//	<< " i=" << i << " indId=" << inds[i]->getId()
						//	<< " age=" << ind.age << " minAge[stage+1]=" << minAge[ind.stage+1][ind.sex]
						//	<< " probdev=" << probdev 
						//	<< endl;
#endif
						if (ind.age >= minAge[ind.stage + 1][ind.sex]) { // old enough to enter next stage
#if RSDEBUG
							//DEBUGLOG << "Population::survival0():"
							//	<< " i=" << i << " indId=" << inds[i]->getId() << " OLD ENOUGH"
							//	<< endl;
#endif
							if (pRandom->Bernoulli(probdev)) {
								inds[i]->developing();
							}
						}
					}
				}
				else { // doomed to die
#if RSDEBUG
					//DEBUGLOG << "Population::survival0():"
					//	<< " i=" << i << " indId=" << inds[i]->getId() << " DIES"
					//	<< endl;
#endif
					inds[i]->setStatus(8);
				}
			}
		}
#if RSDEBUG
		//ind = inds[i]->getStats();
		//DEBUGLOG << "Population::survival0():"
		//	<< " i = " << i << " ID = " << inds[i]->getId()
		//	<< " stage = " << ind.stage << " status = " << ind.status
		//	<< endl;
#endif
	}
}

// Apply survival changes to the population
void Population::survival1(void)
{
	int ninds = (int)inds.size();
#if RSDEBUG
	//DEBUGLOG << "Population::survival1(): this=" << this
	//	<< " patchNum=" << pPatch->getPatchNum() << " ninds=" << ninds
	//	<< endl;
#endif
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();
#if RSDEBUG
		//DEBUGLOG << "Population::survival1(): i=" << i
		//	<< " indId=" << inds[i]->getId() << " stage=" << ind.stage << " sex=" << ind.sex
		//	<< " isDeveloping=" << ind.isDeveloping << " status=" << ind.status
		//	<< endl;
#endif
		if (ind.status > 5) { // doomed to die
			if (ind.status != 10) //not going into cold storage
				delete inds[i];
			inds[i] = NULL;
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
#if RSDEBUG
	//DEBUGLOG << "Population::survival1(): this=" << this
	//	<< " patchNum=" << pPatch->getPatchNum() << " completed individuals loop"
	//	<< endl;
#endif


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
void Population::clean(void)
{
	int ninds = (int)inds.size();
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
#if RS_RCPP
		shuffle(inds.begin(), inds.end(), pRandom->getRNG());
#else

#if !RSDEBUG
		// do not randomise individuals in RSDEBUG mode, as the function uses rand()
		// and therefore the randomisation will differ between identical runs of RS
		shuffle(inds.begin(), inds.end(), pRandom->getRNG());
#endif // !RSDEBUG

#endif // RS_RCPP
	}
}

//---------------------------------------------------------------------------
// Open population file and write header record
bool Population::outPopHeaders(int landNr, bool patchModel) {

	if (landNr == -999) { // close file
		if (outPop.is_open()) outPop.close();
		outPop.clear();
		return true;
	}
	string name;
	simParams sim = paramsSim->getSim();
	envGradParams grad = paramsGrad->getGradient();

	// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
	// ATTRIBUTES OF *ALL* SPECIES AS DETECTED AT MODEL LEVEL
	demogrParams dem = pSpecies->getDemogrParams();
	stageParams sstruct = pSpecies->getStageParams();
	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + Int2Str(sim.batchNum) + "_"
			+ "Sim" + Int2Str(sim.simulation) + "_Land" + Int2Str(landNr) + "_Pop.txt";
	}
	else {
		name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation) + "_Pop.txt";
	}
	outPop.open(name.c_str());
	outPop << "Rep\tYear\tRepSeason";
	if (patchModel) outPop << "\tPatchID\tNcells";
	else outPop << "\tx\ty";
	// determine whether environmental data need be written for populations
	bool writeEnv = false;
	if (grad.gradient) writeEnv = true;
	if (paramsStoch->envStoch()) writeEnv = true;
	if (writeEnv) outPop << "\tEpsilon\tGradient\tLocal_K";
	outPop << "\tSpecies\tNInd";
	if (dem.stageStruct) {
		if (dem.repType == 0)
		{
			for (int i = 1; i < sstruct.nStages; i++) outPop << "\tNInd_stage" << i;
			outPop << "\tNJuvs";
		}
		else {
			for (int i = 1; i < sstruct.nStages; i++)
				outPop << "\tNfemales_stage" << i << "\tNmales_stage" << i;
			outPop << "\tNJuvFemales\tNJuvMales";
		}
	}
	else {
		if (dem.repType != 0) outPop << "\tNfemales\tNmales";
	}
	outPop << endl;
	return outPop.is_open();
}

//---------------------------------------------------------------------------
// Write record to population file
void Population::outPopulation(int rep, int yr, int gen, float eps,
	bool patchModel, bool writeEnv, bool gradK)
{
	Cell* pCell;

// NEED TO REPLACE CONDITIONAL COLUMNS BASED ON ATTRIBUTES OF ONE SPECIES TO COVER
	demogrParams dem = pSpecies->getDemogrParams();

	popStats p;
	outPop << rep << "\t" << yr << "\t" << gen;
	if (patchModel) {
		outPop << "\t" << pPatch->getPatchNum();
		outPop << "\t" << pPatch->getNCells();
	}
	else {
		locn loc = pPatch->getCellLocn(0);
		outPop << "\t" << loc.x << "\t" << loc.y;
	}
	if (writeEnv) {
		if (pPatch->getPatchNum() == 0) { // matrix
			outPop << "\t0\t0\t0";
		}
		else {
			float k = pPatch->getK();
			float envval = 0.0;
			pCell = pPatch->getRandomCell();
			if (pCell != 0) envval = pCell->getEnvVal();
			outPop << "\t" << eps << "\t" << envval << "\t" << k;
		}
	}
#if RSDEBUG
	//DEBUGLOG << "Population::outPopulation(): this=" << this
	//	<< " patchNum=" << pPatch->getPatchNum()
	//	<< " totalPop()=" << totalPop()
	//	<< " nStages=" << nStages << " nSexes=" << nSexes
	//	<< endl;
#endif
	outPop << "\t" << pSpecies->getSpNum();
	if (dem.stageStruct) {
		p = getStats();
		outPop << "\t" << p.nNonJuvs;
		// non-juvenile stage totals from permanent array
		for (int stg = 1; stg < nStages; stg++) {
			for (int sex = 0; sex < nSexes; sex++) {
				outPop << "\t" << nInds[stg][sex];
			}
		}
		// juveniles from permanent array
		for (int sex = 0; sex < nSexes; sex++) {
			outPop << "\t" << nInds[0][sex];
		}
	}
	else { // non-structured population
		outPop << "\t" << totalPop();
		if (dem.repType != 0)
		{ // sexual model
			outPop << "\t" << nInds[1][0] << "\t" << nInds[1][1];
		}
	}
	outPop << endl;

	/*
	#if RS_ABC
	obsdata obs;
	if (abcYear) {
		int nobs = (int)pABCmaster->NObs();
		for (int i = 0; i < nobs; i++) {
			obs = pABCmaster->getObsData(i);
	#if RSDEBUG
	//DEBUGLOG << "Population::outPopulation(): this=" << this << " i=" << i << " yr=" << yr
	//	<< " obs.year=" << obs.year << " obs.type=" << obs.type << " obs.name=" << obs.name
	//	<< " obs.x=" << obs.x << " obs.y=" << obs.y
	//	<< endl;
	#endif
			if (obs.year == yr && obs.type == 2) {
				if (obs.name == "NInds" || obs.name == "Occupied") {
					bool match = false;
					if (patchModel) {
						if (obs.x == pPatch->getPatchNum()) {
							match = true;
	#if RSDEBUG
	//DEBUGLOG << "Population::outPopulation(): i=" << i << " PROCESS Population NInds"
	//	<< " obs.id=" << obs.id << " obs.value=" << obs.value << " obs.x=" << obs.x
	//	<< " pPatch->PatchNum()=" << pPatch->getPatchNum()
	//	<< " totalPop()=" << totalPop() << " p.nNonJuvs=" << p.nNonJuvs
	//	<< endl;
	#endif
						}
					}
					else {
						locn loc = pPatch->getCentroid();
						if (obs.x == loc.x && obs.y == loc.y) {
							match = true;
	#if RSDEBUG
	DEBUGLOG << "Population::outPopulation(): i=" << i << " PROCESS Population NInds"
		<< " obs.id=" << obs.id << " obs.value=" << obs.value << " obs.x="
		<< obs.x << " obs.y=" << obs.y << " loc.x=" << loc.x << " loc.y=" << loc.y
		<< " totalPop()=" << totalPop() << " p.nNonJuvs=" << p.nNonJuvs
		<< endl;
	#endif
						}
					}
					if (match) {
						if (obs.name == "NInds") {
							if (dem.stageStruct)
								pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,p.nNonJuvs,obs.weight);
							else
								pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,totalPop(),obs.weight);
						}
						else { // obs.name == "Occupied"
							pABCmaster->AddNewPred(sim.simulation,obs.id,rep,obs.value,p.breeding,obs.weight);
						}
					}
				}
			}
		}
	}
	#endif // ABC
	*/
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Open individuals file and write header record
void Population::outIndsHeaders(int rep, int landNr, bool patchModel)
{

	if (landNr == -999) { // close file
		if (outInds.is_open()) {
			outInds.close(); outInds.clear();
		}
		return;
	}

	string name;
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	simParams sim = paramsSim->getSim();

	if (sim.batchMode) {
		name = paramsSim->getDir(2)
			+ "Batch" + Int2Str(sim.batchNum) + "_"
			+ "Sim" + Int2Str(sim.simulation)
			+ "_Land" + Int2Str(landNr) + "_Rep" + Int2Str(rep) + "_Inds.txt";
	}
	else {
		name = paramsSim->getDir(2) + "Sim" + Int2Str(sim.simulation)
			+ "_Rep" + Int2Str(rep) + "_Inds.txt";
	}
	outInds.open(name.c_str());
	outInds << "Rep\tYear\tRepSeason\tSpecies\tIndID\tStatus";
	if (patchModel) outInds << "\tNatal_patch\tPatchID";
	else outInds << "\tNatal_X\tNatal_Y\tX\tY";
	if (dem.repType != 0) outInds << "\tSex";
	if (dem.stageStruct) outInds << "\tAge\tStage";
	if (pSpecies->getNumberOfAdaptiveTraits() > 0) outInds << "\tFitness";
	if (emig.indVar) {
		if (emig.densDep) outInds << "\tD0\tAlpha\tBeta";
		else outInds << "\tEP";
	}
	if (trfr.indVar) {
		if (trfr.moveModel) {
			if (trfr.moveType == 1) { // SMS
				outInds << "\tDP\tGB\tAlphaDB\tBetaDB";
			}
			if (trfr.moveType == 2) { // CRW
				outInds << "\tStepLength\tRho";
			}
		}
		else { // kernel
			outInds << "\tMeanDistI";
			if (trfr.twinKern) outInds << "\tMeanDistII\tPKernelI";
		}
	}
	if (sett.indVar) {
		outInds << "\tS0\tAlphaS\tBetaS";
	}
	outInds << "\tDistMoved";
#if RSDEBUG
	// ALWAYS WRITE NO. OF STEPS
	outInds << "\tNsteps";
#else
	if (trfr.moveModel) outInds << "\tNsteps";
#endif
	outInds << endl;
}

//---------------------------------------------------------------------------
// Write records to individuals file
void Population::outIndividual(Landscape* pLandscape, int rep, int yr, int gen,
	int patchNum)
{
	//int x, y, p_id;
	bool writeInd;
	pathSteps steps;
	Cell* pCell;
	landParams ppLand = pLandscape->getLandParams();
	demogrParams dem = pSpecies->getDemogrParams();
	emigRules emig = pSpecies->getEmigRules();
	transferRules trfr = pSpecies->getTransferRules();
	settleType sett = pSpecies->getSettle();
	short spNum = pSpecies->getSpNum();

	int ninds = (int)inds.size();
	for (int i = 0; i < ninds; i++) {
		indStats ind = inds[i]->getStats();
		if (yr == -1) { // write all initialised individuals
			writeInd = true;
			outInds << rep << "\t" << yr << "\t" << dem.repSeasons - 1;
		}
		else {
			if (dem.stageStruct && gen < 0) { // write status 9 individuals only
				if (ind.status == 9) {
					writeInd = true;
					outInds << rep << "\t" << yr << "\t" << dem.repSeasons - 1;
				}
				else writeInd = false;
			}
			else {
				writeInd = true;
				outInds << rep << "\t" << yr << "\t" << gen;
			}
		}
		if (writeInd) {
			outInds << "\t" << spNum << "\t" << inds[i]->getId();
			if (dem.stageStruct) outInds << "\t" << ind.status;
			else { // non-structured population
				outInds << "\t" << ind.status;
			}
			pCell = inds[i]->getLocn(1);
			locn loc = locn();
			if (pCell == 0) loc.x = loc.y = -1; // beyond boundary or in no-data cell
			else loc = pCell->getLocn();
			pCell = inds[i]->getLocn(0);
			locn natalloc = pCell->getLocn();
			//#if SEASONAL
			//		pCell = inds[i]->getLocn(2);
			//		locn prevloc = pCell->getLocn();
			//#endif
			if (ppLand.patchModel) {
				outInds << "\t" << inds[i]->getNatalPatch()->getPatchNum();
				if (loc.x == -1) outInds << "\t-1";
				else outInds << "\t" << patchNum;
			}
			else { // cell-based model
				// EITHER write co-ordinates in cell units ...
				outInds << "\t" << (float)natalloc.x << "\t" << natalloc.y;
				outInds << "\t" << (float)loc.x << "\t" << (float)loc.y;
				// ... OR write co-ordinates in real-world units
	//		outInds << "\t" << (float)natalloc.x * (float)ppLand.resol + (float)lim.minEast
	//			 << "\t" << natalloc.y * (float)ppLand.resol + (float)lim.minNorth;
	//		outInds  << "\t" << (float)loc.x * (float)ppLand.resol + (float)lim.minEast
	//			<< "\t" << (float)loc.y * (float)ppLand.resol + (float)lim.minNorth;
			}
			if (dem.repType != 0) outInds << "\t" << ind.sex;
			if (dem.stageStruct) outInds << "\t" << ind.age << "\t" << ind.stage;

			if (pSpecies->getNumberOfAdaptiveTraits() > 0) outInds << "\t" << inds[i]->getFitness();
		
			if (emig.indVar) {
				emigTraits e = inds[i]->getEmigTraits();
				if (emig.densDep) {
					outInds << "\t" << e.d0 << "\t" << e.alpha << "\t" << e.beta;
				}
				else {
					outInds << "\t" << e.d0;
				}
			} // end of if (emig.indVar)
			if (trfr.indVar) {
				if (trfr.moveModel) {
					if (trfr.moveType == 1) { // SMS
						trfrSMSTraits s = inds[i]->getSMSTraits();
						outInds << "\t" << s.dp << "\t" << s.gb;
						outInds << "\t" << s.alphaDB << "\t" << s.betaDB;
					} // end of SMS
					if (trfr.moveType == 2) { // CRW
						trfrCRWTraits c = inds[i]->getCRWTraits();
						outInds << "\t" << c.stepLength << "\t" << c.rho;
#if RSDEBUG
						//DEBUGLOG << "Population::outIndividual():"
						//	<< " patchNum=" << patchNum << " i=" << i << " ID=" << inds[i]->getId()
						//	<< " nTrfrGenes=" << nTrfrGenes << " loc[0][0].allele[0]=" << loc[0][0].allele[0]
						//	<< endl;
#endif
					} // end of CRW
				}
				else { // kernel
					trfrKernelParams k = inds[i]->getKernTraits();
					if (trfr.twinKern)
					{
						outInds << "\t" << k.meanDist1 << "\t" << k.meanDist2 << "\t" << k.probKern1;
					}
					else {
						outInds << "\t" << k.meanDist1;
					}
				}
			}

			if (sett.indVar) {
				settleTraits s = inds[i]->getSettTraits();
				outInds << "\t" << s.s0 << "\t" << s.alpha << "\t" << s.beta;
			}

			// distance moved (metres)
			if (loc.x == -1) outInds << "\t-1";
			else {
				float d = ppLand.resol * sqrt((float)((natalloc.x - loc.x) * (natalloc.x - loc.x)
					+ (natalloc.y - loc.y) * (natalloc.y - loc.y)));
				outInds << "\t" << d;
			}
#if RSDEBUG
			// ALWAYS WRITE NO. OF STEPS
			steps = inds[i]->getSteps();
			outInds << "\t" << steps.year;
#else
			if (trfr.moveModel) {
				steps = inds[i]->getSteps();
				outInds << "\t" << steps.year;
			}
#endif
			outInds << endl;
		} // end of writeInd condition

	}
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


