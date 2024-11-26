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

#include "Individual.h"
//---------------------------------------------------------------------------

int Individual::indCounter = 0;
TraitFactory Individual::traitFactory = TraitFactory();

//---------------------------------------------------------------------------

bool isAlive(indStatus s) {
	return (s != diedInTransfer
		&& s != diedInTrfrMort
		&& s != diedDemogrMort
		&& s != diedOldAge);
}

string to_string(indStatus s) {
	switch (s) {
	case initial: return "initial";
	case dispersing: return "dispersing";
	case waitSettlement: return "waitSettlement";
	case waitNextDispersal: return "waitNextDispersal";
	case settled: return "settled";
	case settledNeighbour: return "settledNeighbour";
	case diedInTransfer: return "diedInTransfer";
	case diedInTrfrMort: return "diedInTrfrMort";
	case diedDemogrMort: return "diedDemogrMort";
	case diedOldAge: return "diedOldAge";
	default: return "";
	}
}

// Individual constructor
Individual::Individual(Species* pSp, Cell* pCell, 
	Patch* pPatch, short stg, short a, short repInt,
	float probmale, bool movt, short moveType)
{
	pSpecies = pSp;
	indId = indCounter; 
	indCounter++; // unique identifier for each individual
	geneticFitness = 1.0;
	stage = stg;
	if (probmale <= 0.0) sex = FEM;
	else sex = pRandom->Bernoulli(probmale) ? MAL : FEM;
	age = a;
	status = indStatus::initial;

	if (sex == 0 && repInt > 0) { // set no. of fallow seasons for female
		fallow = pRandom->IRandom(0, repInt);
	}
	else fallow = 9999;
	isDeveloping = false;
	pPrevCell = pCurrCell = pCell;
	pNatalPatch = pPatch;
	pTrfrData = nullptr; //set to null as default
	if (movt) {

		locn loc = pCell->getLocn();
		path = new pathData;
		path->year = 0; 
		path->total = 0; 
		path->out = 0;
		path->pSettPatch = nullptr; 
		path->settleStatus = 0;

		if (moveType == 1) { // SMS
			// set up location data for SMS
			pTrfrData = make_unique<smsData>(loc, loc);
		}
		if (moveType == 2) { // CRW
			// set up continuous co-ordinates etc. for CRW movement
			float xc = ((float)pRandom->Random() * 0.999f) + (float)loc.x;
			float yc = (float)(pRandom->Random() * 0.999f) + (float)loc.y;
			float prevdrn = (float)(pRandom->Random() * 2.0 * PI);
			pTrfrData = make_unique<crwData>(prevdrn, xc, yc);
		}
	}
	else {
		path = nullptr;
		pTrfrData = make_unique<kernelData>(0.0, 0.0, 0.0);
	}
}

Individual::~Individual() {
	if (path != nullptr) delete path;
}

QuantitativeTrait* Individual::getTrait(TraitType trait) const {
	auto p = this->spTraitTable.find(trait);
	if (p == spTraitTable.end())
		throw runtime_error("Trait does not exist in trait table.");
	else return p->second.get();
}

set<TraitType> Individual::getTraitTypes() {
	auto kv = std::views::keys(this->spTraitTable);
	set< TraitType > keys{ kv.begin(), kv.end() };
	return keys;
}

//---------------------------------------------------------------------------
// Inheritance for diploid, sexual species
//---------------------------------------------------------------------------
void Individual::inheritGenes(const Individual* mother, const Individual* father) {

	int events = 0;
	const set<int> chromosomeEnds = pSpecies->getChromosomeEnds();
	const int genomeSize = pSpecies->getGenomeSize();

	int maternalStartingChromosome = pRandom->Bernoulli(0.5);
	int paternalStartingChromosome = pRandom->Bernoulli(0.5);

	set<unsigned int> maternalRecomPositions;
	set<unsigned int> paternalRecomPositions;

	// Determine which parental chromosomes are inherited
	for (int pos : chromosomeEnds) {
		if (pRandom->Bernoulli(0.5)) // switch strand for next chromosome
			maternalRecomPositions.insert(pos);
		if (pRandom->Bernoulli(0.5))
			paternalRecomPositions.insert(pos);
	}

	// Draw recombination events for maternal genome
	if (pSpecies->getRecombinationRate() > 0.0)
		events = pRandom->Binomial(genomeSize, pSpecies->getRecombinationRate());
	// if poisson exceeds genomeSize, bound to genomeSize
	int nbrCrossOvers = events + maternalRecomPositions.size();
	if (nbrCrossOvers > genomeSize) {
		nbrCrossOvers = genomeSize;
	}
	while (maternalRecomPositions.size() < nbrCrossOvers) {
		// Sample recombination sites
		maternalRecomPositions.insert(pRandom->IRandom(0, genomeSize));
	}

	// Draw recombination events for paternal genome
	if (pSpecies->getRecombinationRate() > 0.0)
		events = pRandom->Binomial(genomeSize, pSpecies->getRecombinationRate());
	nbrCrossOvers = events + paternalRecomPositions.size();
	if (nbrCrossOvers > genomeSize) {
		nbrCrossOvers = genomeSize;
	}
	while (paternalRecomPositions.size() < nbrCrossOvers) {
		paternalRecomPositions.insert(pRandom->IRandom(0, genomeSize));
	}

	// Inherit genes for each trait
	const auto& spTraits = pSpecies->getTraitTypes();
	for (auto const& trait : spTraits)
	{
		const auto motherTrait = mother->getTrait(trait);
		const auto fatherTrait = father->getTrait(trait);
		auto newTrait = motherTrait->clone(); // shallow copy pointer to species-level attributes

		// Inherit from mother first
		newTrait->inheritGenes(true, motherTrait, maternalRecomPositions, maternalStartingChromosome);
		if (newTrait->isInherited()) {
			// Inherit father trait values
			newTrait->inheritGenes(false, fatherTrait, paternalRecomPositions, paternalStartingChromosome);
			if (newTrait->getMutationRate() > 0 && pSpecies->areMutationsOn())
				newTrait->mutate();
		}
		if (trait == GENETIC_LOAD1 || trait == GENETIC_LOAD2 || trait == GENETIC_LOAD3 || trait == GENETIC_LOAD4 || trait == GENETIC_LOAD5)
			geneticFitness *= newTrait->express();

		// Add the inherited trait and genes to the newborn's list
		spTraitTable.insert(make_pair(trait, move(newTrait)));
	}
}

//---------------------------------------------------------------------------
// Inheritance for haploid, asexual species
//---------------------------------------------------------------------------
void Individual::inheritGenes(const Individual* mother) {
	set<unsigned int> recomPositions; //not used here cos haploid but need it for inherit function, not ideal 
	int startingChromosome = 0;

	const auto& spTraits = pSpecies->getTraitTypes();

	for (auto const& trait : spTraits)
	{
		const auto motherTrait = mother->getTrait(trait);
		auto newTrait = motherTrait->clone(); // shallow copy, pointer to species trait initialised and empty sequence

		newTrait->inheritGenes(true, motherTrait, recomPositions, startingChromosome);
		if (newTrait->isInherited()) {
			if (newTrait->getMutationRate() > 0 && pSpecies->areMutationsOn())
				newTrait->mutate();
		}
		if (trait == GENETIC_LOAD1 || trait == GENETIC_LOAD2 || trait == GENETIC_LOAD3 || trait == GENETIC_LOAD4 || trait == GENETIC_LOAD5)
			geneticFitness *= newTrait->express();

		// Add the inherited trait and genes to the newborn's list
		spTraitTable.insert(make_pair(trait, move(newTrait)));
	}
}

// Initialise individual trait genes from species-level traits
void Individual::setUpGenes(int resol) {

	// this way to keep spp trait table immutable i.e. not able to call getTraitTable, 
	// could pass it back by value (copy) instead but could be heavy if large map
	const auto& traitTypes = pSpecies->getTraitTypes();
	for (auto const& traitType : traitTypes)
	{
		const auto spTrait = pSpecies->getSpTrait(traitType);
		this->spTraitTable.emplace(traitType, traitFactory.Create(traitType, spTrait));
	}
	expressDispersalPhenotypes(resol);
}

void Individual::expressDispersalPhenotypes(int resol) {

	const emigRules emig = pSpecies->getEmigRules();
	const transferRules trfr = pSpecies->getTransferRules();
	const settleType sett = pSpecies->getSettle();
	const settleRules settRules = pSpecies->getSettRules(stage, sex);

	// record phenotypic traits
	if (emig.indVar)
		this->expressEmigTraits(emig.sexDep, emig.densDep);
	if (trfr.indVar)
		this->expressTransferTraits(trfr, resol);
	if (sett.indVar)
		this->expressSettlementTraits(sett.sexDep, settRules.densDep);
}

void Individual::expressTransferTraits(transferRules trfr, int resol) {
	if (trfr.usesMovtProc) {
		if (trfr.moveType == 1) {
			expressSMSTraits();
		}
		else
			expressCRWTraits();
	}
	else
		expressKernelTraits(trfr.sexDep, trfr.twinKern, resol);
}

void Individual::expressSettlementTraits(bool sexDep, bool densDep) {

	settleTraits s; s.s0 = s.alpha = s.beta = 0.0;
	if (sexDep) {
		if (this->getSex() == MAL) {
			s.s0 = getTrait(S_S0_M)->express();
			if (densDep) {
				s.alpha = getTrait(S_ALPHA_M)->express();
				s.beta = getTrait(S_BETA_M)->express();
			}
		}
		else if (this->getSex() == FEM) {
			s.s0 = getTrait(S_S0_F)->express();
			if (densDep) {
				s.alpha = getTrait(S_ALPHA_F)->express();
				s.beta = getTrait(S_BETA_F)->express();
			}
		}
		else {
			throw runtime_error("Attempt to express invalid emigration trait sex.");
		}
	}
	else {
		s.s0 = getTrait(S_S0)->express();
		if (densDep) {
			s.alpha = getTrait(S_ALPHA)->express();
			s.beta = getTrait(S_BETA)->express();
		}
	}

	pSettleTraits = make_unique<settleTraits>();
	pSettleTraits->s0 = (float)(s.s0);
	pSettleTraits->alpha = (float)(s.alpha);
	pSettleTraits->beta = (float)(s.beta);
	if (pSettleTraits->s0 < 0.0) pSettleTraits->s0 = 0.0;
	if (pSettleTraits->s0 > 1.0) pSettleTraits->s0 = 1.0;
	return;
}


// Inherit genome from parent(s), diploid
void Individual::inheritTraits(Individual* mother, Individual* father, int resol)
{
	inheritGenes(mother, father);
	expressDispersalPhenotypes(resol);
}

// Inherit genome from mother, haploid
void Individual::inheritTraits(Individual* mother, int resol)
{
	inheritGenes(mother);
	expressDispersalPhenotypes(resol);
}

//---------------------------------------------------------------------------

// Identify whether an individual is a potentially breeding female -
// if so, return her stage, otherwise return 0
int Individual::breedingFem() {
	if (sex == FEM) {
		if (status == initial 
			|| status == settled 
			|| status == settledNeighbour) 
			return stage;
		else return 0;
	}
	else return 0;
}

int Individual::getId() { return indId; }

sex_t Individual::getSex() { return sex; }

indStatus Individual::getStatus() { return status; }

float Individual::getGeneticFitness() { return geneticFitness; }

bool Individual::isViable() const {
	float probViability = geneticFitness > 1.0 ? 1.0 : geneticFitness;
	return probViability >= pRandom->Random();
}

indStats Individual::getStats() {
	indStats s;
	s.stage = stage; 
	s.sex = sex; 
	s.age = age; 
	s.status = status; 
	s.fallow = fallow;
	s.isDeveloping = isDeveloping;
	return s;
}

Cell* Individual::getLocn(const short option) {
	if (option == 0) { // return previous location
		return pPrevCell;
	}
	else { // return current location
		return pCurrCell;
	}
}

Patch* Individual::getNatalPatch() { return pNatalPatch; }

void Individual::setYearSteps(int t) {
	if (path != nullptr && t >= 0) {
		if (t >= 0) path->year = t;
		else path->year = 666;
	}
}

pathSteps Individual::getSteps() {
	pathSteps s;
	if (path == nullptr) {
		s.year = 0; 
		s.total = 0; 
		s.out = 0;
	}
	else {
		s.year = path->year; 
		s.total = path->total; 
		s.out = path->out;
	}
	return s;
}

settlePatch Individual::getSettPatch() {
	settlePatch s;
	if (path == nullptr) {
		s.pSettPatch = nullptr; 
		s.settleStatus = 0;
	}
	else {
		s.pSettPatch = path->pSettPatch; 
		s.settleStatus = path->settleStatus;
	}
	return s;
}

void Individual::setSettPatch(const settlePatch s) {
	if (path == nullptr) {
		path = new pathData;
		path->year = 0; 
		path->total = 0; 
		path->out = 0; 
		path->settleStatus = 0;
#if RS_RCPP
		path->pathoutput = 1;
#endif
	}
	if (s.settleStatus >= 0 && s.settleStatus <= 2) path->settleStatus = s.settleStatus;
	path->pSettPatch = s.pSettPatch;
}

void Individual::expressEmigTraits(bool sexDep, bool densityDep) {
	emigTraits e; 
	e.d0 = e.alpha = e.beta = 0.0;
	if (sexDep) {
		if (this->getSex() == MAL) {
			e.d0 = this->getTrait(E_D0_M)->express();
			if (densityDep) {
				e.alpha = getTrait(E_ALPHA_M)->express();
				e.beta = getTrait(E_BETA_M)->express();
			}
		}
		else if (this->getSex() == FEM) {
			e.d0 = this->getTrait(E_D0_F)->express();
			if (densityDep) {
				e.alpha = getTrait(E_ALPHA_F)->express();
				e.beta = getTrait(E_BETA_F)->express();
			}
		}
		else {
			throw runtime_error("Attempt to express invalid emigration trait sex.");
		}
	}	
	else {
		e.d0 = this->getTrait(E_D0)->express();
		if (densityDep) {
			e.alpha = getTrait(E_ALPHA)->express();
			e.beta = getTrait(E_BETA)->express();
		}
	}

	pEmigTraits = make_unique<emigTraits>();
	pEmigTraits->d0 = e.d0;
	pEmigTraits->alpha = e.alpha;
	pEmigTraits->beta = e.beta;

	// Below must never trigger, phenotype is bounded in express()
	if (pEmigTraits->d0 < 0.0) throw runtime_error("d0 value has become negative.");
	if (pEmigTraits->d0 > 1.0) throw runtime_error("d0 value has exceeded 1.");
	return;
}

// Get phenotypic emigration traits
emigTraits Individual::getIndEmigTraits() {
	emigTraits e; 
	e.d0 = e.alpha = e.beta = 0.0;
	if (pEmigTraits != nullptr) {
		e.d0 = pEmigTraits->d0;
		e.alpha = pEmigTraits->alpha;
		e.beta = pEmigTraits->beta;
	}
	return e;
}
// Set phenotypic transfer by kernel traits
void Individual::expressKernelTraits(bool sexDep, bool twinKernel, int resol) {

	trfrKernelParams k; 
	k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;

	if (sexDep) {
		if (this->sex == MAL) {
			k.meanDist1 = getTrait(KERNEL_MEANDIST_1_M)->express();

			if (twinKernel) { // twin kernel
				k.meanDist2 = getTrait(KERNEL_MEANDIST_2_M)->express();
				k.probKern1 = getTrait(KERNEL_PROBABILITY_M)->express();
			}
		}
		else if (this->sex == FEM) {
			k.meanDist1 = getTrait(KERNEL_MEANDIST_1_F)->express();

			if (twinKernel) { // twin kernel
				k.meanDist2 = getTrait(KERNEL_MEANDIST_2_F)->express();
				k.probKern1 = getTrait(KERNEL_PROBABILITY_F)->express();
			}
		}
		else {
			throw runtime_error("Attempt to express invalid kernel transfer trait sex.");
		}
	}
	else {
		k.meanDist1 = getTrait(KERNEL_MEANDIST_1)->express();

		if (twinKernel) { // twin kernel
			k.meanDist2 = getTrait(KERNEL_MEANDIST_2)->express();
			k.probKern1 = getTrait(KERNEL_PROBABILITY)->express();
		}
	}
	
	float meanDist1 = (float)(k.meanDist1);
	float meanDist2 = (float)(k.meanDist2);
	float probKern1 = (float)(k.probKern1);

	if (!pSpecies->useFullKernel()) {
		// kernel mean(s) may not be less than landscape resolution
		if (meanDist1 < resol) meanDist1 = (float)resol;
		if (meanDist2 < resol) meanDist2 = (float)resol;
	}
	if (probKern1 < 0.0) probKern1 = 0.0;
	if (probKern1 > 1.0) probKern1 = 1.0;
	auto& pKernel = dynamic_cast<kernelData&>(*pTrfrData);
	pKernel.meanDist1 = meanDist1;
	pKernel.meanDist2 = meanDist2;
	pKernel.probKern1 = probKern1;

	return;
}



// Get phenotypic emigration traits
trfrKernelParams Individual::getIndKernTraits() {
	trfrKernelParams k; 
	k.meanDist1 = k.meanDist2 = k.probKern1 = 0.0;

	if (pTrfrData != nullptr) {
		auto& pKernel = dynamic_cast<const kernelData&>(*pTrfrData);
		k.meanDist1 = pKernel.meanDist1;
		k.meanDist2 = pKernel.meanDist2;
		k.probKern1 = pKernel.probKern1;
	}

	return k;
}

void Individual::expressSMSTraits() {

	trfrSMSTraits s = pSpecies->getSpSMSTraits();

	double dp, gb, alphaDB, betaDB;
	dp = gb = alphaDB = betaDB = 0.0;
	dp = getTrait(SMS_DP)->express();
	if (s.goalType == 2) {
		gb = getTrait(SMS_GB)->express();
		alphaDB = getTrait(SMS_ALPHADB)->express();
		betaDB = getTrait(SMS_BETADB)->express();
	}

	auto& pSMS = dynamic_cast<smsData&>(*pTrfrData);
	pSMS.dp = (float)(dp);
	pSMS.gb = (float)(gb);
	if (s.goalType == 2) {
		pSMS.alphaDB = (float)(alphaDB);
		pSMS.betaDB = (int)(betaDB);
	}
	else {
		pSMS.alphaDB = s.alphaDB;
		pSMS.betaDB = s.betaDB;
	}
	if (pSMS.dp < 1.0) pSMS.dp = 1.0;
	if (pSMS.gb < 1.0) pSMS.gb = 1.0;
	if (pSMS.alphaDB <= 0.0) pSMS.alphaDB = 0.000001f;
	if (pSMS.betaDB < 1) pSMS.betaDB = 1;
	return;
}

trfrData* Individual::getTrfrData() {
	return pTrfrData.get();
}

// Get phenotypic transfer by SMS traits
trfrSMSTraits Individual::getIndSMSTraits() {

	trfrSMSTraits s; 
	s.dp = s.gb = s.alphaDB = 1.0; 
	s.betaDB = 1;
	if (pTrfrData != nullptr) {
		auto& pSMS = dynamic_cast<const smsData&>(*pTrfrData);
		s.dp = pSMS.dp; 
		s.gb = pSMS.gb;
		s.alphaDB = pSMS.alphaDB; 
		s.betaDB = pSMS.betaDB;
	}

	return s;
}


// Set phenotypic transfer by CRW traits
void Individual::expressCRWTraits() {
	trfrCRWTraits c; 
	c.stepLength = c.rho = 0.0;

	c.stepLength = getTrait(CRW_STEPLENGTH)->express();
	c.rho = getTrait(CRW_STEPCORRELATION)->express();

	auto& pCRW = dynamic_cast<crwData&>(*pTrfrData);
	pCRW.stepLength = (float)(c.stepLength);
	pCRW.rho = (float)(c.rho);
	if (pCRW.stepLength < 1.0) pCRW.stepLength = 1.0;
	if (pCRW.rho < 0.0) pCRW.rho = 0.0;
	if (pCRW.rho > 0.999) pCRW.rho = 0.999f;
	return;
}

// Get phenotypic transfer by CRW traits
trfrCRWTraits Individual::getIndCRWTraits() {

	trfrCRWTraits c; 
	c.stepLength = c.rho = 0.0;
	if (pTrfrData != 0) {
		auto& pCRW = dynamic_cast<const crwData&>(*pTrfrData);
		c.stepLength = pCRW.stepLength;
		c.rho = pCRW.rho;
	}
	return c;

}

// Get phenotypic settlement traits
settleTraits Individual::getIndSettTraits() {
	settleTraits s;
	s.s0 = s.alpha = s.beta = 0.0;
	if (pSettleTraits != nullptr) {
		s.s0 = pSettleTraits->s0;
		s.alpha = pSettleTraits->alpha;
		s.beta = pSettleTraits->beta;
	}

	return s;
}


void Individual::setStatus(indStatus s) {
	status = s;
}

void Individual::setToDevelop() {
	isDeveloping = true;
}

void Individual::develop() {
	stage++; 
	isDeveloping = false;
}

void Individual::ageIncrement(short maxage) {
	if (isAlive(status)) {
		age++;
		if (age > maxage) status = diedOldAge; // exceeds max. age - dies
		else {
			if (path != nullptr) path->year = 0;	// reset annual step count for movement models
			if (status == waitNextDispersal)
				status = dispersing;
		}
	}
}

void Individual::incFallow() { fallow++; }

void Individual::resetFallow() { fallow = 0; }

//---------------------------------------------------------------------------
// Move to a specified neighbouring cell
void Individual::moveTo(Cell* newCell) {

	// Check that location is indeed a neighbour of the current cell
	locn currloc = pCurrCell->getLocn();
	locn newloc = newCell->getLocn();
	double distance = sqrt((currloc.x - newloc.x) * (currloc.x - newloc.x)
		+ (currloc.y - newloc.y) * (currloc.y - newloc.y));
	if (distance >= 1.0 && distance < 1.5) { // ok
		pCurrCell = newCell; 
		status = settledNeighbour;
	}
}

//---------------------------------------------------------------------------
// Move to a new cell by sampling a dispersal distance from a single or double
// negative exponential kernel
// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
bool Individual::moveKernel(Landscape* pLandscape, const bool absorbing)
{
	Patch* patch;
	int patchNum = 0;
	int newX = 0, newY = 0;
	bool isDispersing = true;
	double xrand, yrand, meandist, dist, r1, rndangle, nx, ny;
	float localK;
	trfrKernelParams kern;
	Cell* pCell;
	Patch* pPatch;
	locn loc = pCurrCell->getLocn();

	landData land = pLandscape->getLandData();

	bool usefullkernel = pSpecies->useFullKernel();
	transferRules trfr = pSpecies->getTransferRules();
	settleRules sett = pSpecies->getSettRules(stage, sex);

	pCell = NULL;
	pPatch = NULL;

	if (trfr.indVar) { // get individual's kernel parameters
		kern.meanDist1 = kern.meanDist2 = kern.probKern1 = 0.0;

		auto& pKernel = dynamic_cast<const kernelData&>(*pTrfrData);

		kern.meanDist1 = pKernel.meanDist1;
		if (trfr.twinKern)
		{
			kern.meanDist2 = pKernel.meanDist2;
			kern.probKern1 = pKernel.probKern1;
		}
	}
	else { // get kernel parameters for the species
		if (trfr.sexDep) {
			if (trfr.stgDep) {
				kern = pSpecies->getSpKernTraits(stage, sex);
			}
			else {
				kern = pSpecies->getSpKernTraits(0, sex);
			}
		}
		else {
			if (trfr.stgDep) {
				kern = pSpecies->getSpKernTraits(stage, 0);
			}
			else {
				kern = pSpecies->getSpKernTraits(0, 0);
			}
		}
	}

	// scale the appropriate kernel mean to the cell size
	if (trfr.twinKern)
	{
		if (pRandom->Bernoulli(kern.probKern1))
			meandist = kern.meanDist1 / (float)land.resol;
		else
			meandist = kern.meanDist2 / (float)land.resol;
	}
	else
		meandist = kern.meanDist1 / (float)land.resol;
	// scaled mean may not be less than 1 unless emigration derives from the kernel
	// (i.e. the 'use full kernel' option is applied)
	if (!usefullkernel && meandist < 1.0) meandist = 1.0;

	int loopsteps = 0; // new counter to prevent infinite loop added 14/8/15
	do {
		do {
			do {
				// randomise the cell within the patch, provided that the individual is still in
				// its natal cell (i.e. not waiting in the matrix)
				// this is because, if the patch is very large, the individual is near the centre
				// and the (single) kernel mean is (not much more than) the cell size, an infinite
				// loop could otherwise result, as the individual never reaches the patch edge
				// (in a cell-based model, this has no effect, other than as a processing overhead)
				if (status == dispersing) {
					pCell = pNatalPatch->getRandomCell();
					if (pCell != 0) {
						loc = pCell->getLocn();
					}
				}
				// randomise the position of the individual inside the cell
				// so x and y are a corner of the cell?
				xrand = (double)loc.x + pRandom->Random() * 0.999;
				yrand = (double)loc.y + pRandom->Random() * 0.999;

				// draw factor r1 0 < r1 <= 1
				r1 = 0.0000001 + pRandom->Random() * (1.0 - 0.0000001);
				dist = (-1.0 * meandist) * log(r1);

				rndangle = pRandom->Random() * 2.0 * PI;
				nx = xrand + dist * sin(rndangle);
				ny = yrand + dist * cos(rndangle);
				if (nx < 0.0) newX = -1; else newX = (int)nx;
				if (ny < 0.0) newY = -1; else newY = (int)ny;
#ifndef NDEBUG
				if (path != 0) (path->year)++;
#endif
				loopsteps++;
			} while (loopsteps < 1000 &&
				// keep drawing if out of bounds of landscape or same cell
				((!absorbing && (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY))
					|| (!usefullkernel && newX == loc.x && newY == loc.y))
				);
			if (loopsteps < 1000) {
				if (newX < land.minX || newX > land.maxX
					|| newY < land.minY || newY > land.maxY) { // beyond absorbing boundary
					// this cannot be reached if not absorbing?
					pCell = 0;
					patch = 0;
					patchNum = -1;
				}
				else {
					pCell = pLandscape->findCell(newX, newY);
					if (pCell == 0) { // no-data cell
						patch = 0;
						patchNum = -1;
					}
					else {
						patch = pCell->getPatch();
						if (patch == 0) { // matrix
							pPatch = 0;
							patchNum = 0;
						}
						else {
							pPatch = patch;
							patchNum = pPatch->getPatchNum();
						}
					}
				}
			}
			else { // exceeded 1000 attempts
				patch = 0;
				patchNum = -1;
			}
		} while (!absorbing && patchNum < 0 && loopsteps < 1000); 			 // in a no-data region
	} while (!usefullkernel && pPatch == pNatalPatch && loopsteps < 1000); 	// still in the original (natal) patch

	if (loopsteps < 1000) {
		if (pCell == 0) { // beyond absorbing boundary or in no-data cell
			// only if absorbing=true and out of bounddaries
			pCurrCell = 0;
			status = diedInTransfer;
			isDispersing = false;
		}
		else {
			pCurrCell = pCell;
			if (pPatch == 0) localK = 0.0; // matrix
			else localK = pPatch->getK();
			if (patchNum > 0 && localK > 0.0) { // found a new patch
				status = waitSettlement; // record as potential settler
			}
			else {
				// unsuitable patch
				isDispersing = false;
				// can wait in matrix if population is stage structured ...
				if (pSpecies->stageStructured()) {
					// ... and wait option is applied ...
					if (sett.wait) { // ... it is
						status = waitNextDispersal; // waiting
					}
					else // ... it is not
						status = diedInTransfer; // dies (unless there is a suitable neighbouring cell)
				}
				else status = diedInTransfer; // dies (unless there is a suitable neighbouring cell)
			}
		}
	}
	else { // exceeded 1000 attempts
		status = diedInTransfer;
		isDispersing = false;
	}

	// apply dispersal-related mortality, which may be distance-dependent
	dist *= (float)land.resol; // re-scale distance moved to landscape scale
	if (isAlive(status) || status == diedInTransfer) {
		double dispmort;
		trfrMortParams mort = pSpecies->getMortParams();
		if (trfr.distMort) {
			dispmort = 1.0 / (1.0 + exp(-(dist - mort.mortBeta) * mort.mortAlpha));
		}
		else {
			dispmort = mort.fixedMort;
		}
		if (pRandom->Bernoulli(dispmort)) {
			status = diedInTrfrMort;
			isDispersing = false;
		}
	}

	return isDispersing;
}

//---------------------------------------------------------------------------
// Make a single movement step according to a mechanistic movement model
// Returns 1 if still dispersing (including having found a potential patch), otherwise 0
bool Individual::moveStep(Landscape* pLandscape,
	const short landIx, const bool absorbing)
{
	if (status != dispersing) return false; // not currently dispersing

	int patchNum;
	int newX, newY;
	locn loc;
	bool isDispersing = true;
	double xcnew, ycnew;
	double moveDirection;
	double probMort, rho, stepLength;
	movedata move;
	Patch* pPatch = nullptr;
	bool absorbed = false;

	landData land = pLandscape->getLandData();
	simParams sim = paramsSim->getSim();
	transferRules trfr = pSpecies->getTransferRules();
	trfrCRWTraits movt = pSpecies->getSpCRWTraits();
	settleSteps settsteps = pSpecies->getSteps(stage, sex);

	pPatch = pCurrCell->getPatch();
	if (pPatch == nullptr) { // no data
		patchNum = 0;
	} else {
		patchNum = pPatch->getPatchNum();
	}

	// Apply step-dependent mortality risk
	if (pPatch == pNatalPatch
		&& path->out == 0
		&& path->year == path->total) {
		// no mortality if ind. has not yet left natal patch
		probMort = 0.0;
	}
	else if (trfr.habMort) {
		int habIndex = pCurrCell->getHabIndex(landIx);
		probMort = habIndex < 0 ? 1.0 : pSpecies->getHabMort(habIndex);
	}
	else probMort = movt.stepMort;

	if (pRandom->Bernoulli(probMort)) {
		status = diedInTrfrMort;
		return false;
	}
	else { 
		// Take a step
		(path->year)++;
		(path->total)++;

		if (pPatch == nullptr || patchNum == 0) { // not in a patch
			// Reset path settlement status
			if (path != nullptr) path->settleStatus = 0; 
			(path->out)++;
		}
		loc = pCurrCell->getLocn();
		newX = loc.x; 
		newY = loc.y;

		switch (trfr.moveType) {

		case 1: // SMS
			move = smsMove(pLandscape, landIx, pPatch == pNatalPatch, trfr.indVar, absorbing);
			if (move.dist < 0.0) {
				// either INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
				// or individual has crossed absorbing boundary: individual dies
				status = diedInTransfer;
				isDispersing = false;
			}
			else {
				pPatch = pCurrCell->getPatch();
				if (sim.saveVisits && pPatch != pNatalPatch) {
					pCurrCell->incrVisits();
				}
			}
			break;

		case 2: // CRW

			auto& pCRW = dynamic_cast<crwData&>(*pTrfrData);

			if (trfr.indVar) {
				movt.stepLength = pCRW.stepLength;
				movt.rho = pCRW.rho;
			}
			stepLength = movt.stepLength; 
			rho = movt.rho;

			// Move in a straight line if...
			if (pPatch == pNatalPatch 
				// ... still in natal patch or
				|| (movt.straightenPath && path->settleStatus > 0) 
				// ... must straighten path to (previously determined) settlement patch 
				){
				rho = 0.99;
				path->out = 0;
			}
			int loopSteps = 0; // no infinite loop
			constexpr int maxLoopSteps = 1000;
			do {
				do {
					// Sample direction
					if (!isInLandscape(newX, newY, land)
						|| pCurrCell == nullptr) {
						// Random direction to avoid invalid area again
						moveDirection = drawDirection(pCRW.prevdrn, 0.0);
					}
					else moveDirection = drawDirection(pCRW.prevdrn, rho);

					// Get new coordinates
					xcnew = pCRW.xc + sin(moveDirection) * stepLength / land.resol;
					ycnew = pCRW.yc + cos(moveDirection) * stepLength / land.resol;
					newX = xcnew < 0.0 ? -1 : trunc(xcnew);
					newY = ycnew < 0.0 ? -1 : trunc(ycnew);

					loopSteps++;
				} while (!absorbing && 
					loopSteps < maxLoopSteps &&
					!isInLandscape(newX, newY, land));
				
				// Get cell and patch for new coordinates
				if (!isInLandscape(newX, newY, land)) pCurrCell = nullptr;
				else pCurrCell = pLandscape->findCell(newX, newY);
				if (pCurrCell == nullptr) { // no-data or beyond absorbing boundary
					pPatch = nullptr;
					if (absorbing) absorbed = true;
				}
				else pPatch = pCurrCell->getPatch();

			} while (!absorbing 
				&& pCurrCell == nullptr 
				&& loopSteps < maxLoopSteps);

			// Update individual through ref
			pCRW.prevdrn = moveDirection;
			pCRW.xc = xcnew; 
			pCRW.yc = ycnew;

			if (absorbed) { // beyond absorbing boundary or in no-data square
				status = diedInTransfer;
				isDispersing = false;
				pCurrCell = nullptr;
			}
			else if (loopSteps >= maxLoopSteps) { // unable to make a move
				// INTERNAL ERROR CONDITION - INDIVIDUAL IS IN NO-DATA SQUARE
				// NEED TO TAKE SOME FORM OF INFORMATIVE ACTION ...
				// ... individual dies as it cannot move
				status = diedInTransfer;
				isDispersing = false;
				// current cell is invalid, get back to previous cell
				pCurrCell = pPrevCell;
			}
			break;

		} // end of switch (trfr.moveType)

		// Update individual status
		if (pPatch != nullptr  // not no-data area or matrix
			&& path->total >= settsteps.minSteps) {
			if (pPatch != pNatalPatch
				&& pPatch->getK() > 0.0) {
				status = waitSettlement; // new patch is suitable
			}
		}
		if (status != waitSettlement 
			&& status != diedInTransfer) { // no suitable patch but not dead yet
			if (path->year >= settsteps.maxStepsYr) {
				status = waitNextDispersal; // try again next year
			}
			if (path->total >= settsteps.maxSteps) {
				status = diedInTransfer;
				isDispersing = false;
			}
		}
	} // end of single movement step

	return isDispersing;
}

//---------------------------------------------------------------------------

// Functions to implement the SMS algorithm

// Move to a neighbouring cell according to the SMS algorithm
movedata Individual::smsMove(Landscape* pLand, const short landIx, 
	const bool isInNatalPatch, const bool indvar, const bool absorbing)
{
	array3x3d neighbourWeights; // to hold weights/costs/probs of moving to neighbouring cells
	array3x3d goalBiasWeights;	// to hold weights for moving towards a goal location
	array3x3f habDepWeights;	// to hold weights for habitat (includes percep range)
	int newX = 0, newY = 0;
	Cell* pCell;
	Cell* pNewCell = NULL;
	double sum_nbrs = 0.0;
	movedata move;
	int cellcost, newcellcost;
	locn currLoc;

	auto& pSMS = dynamic_cast<smsData&>(*pTrfrData);
	if (pCurrCell == nullptr) {
		throw runtime_error("Individual found in a no-data cell at beginning of SMS.");
	}

	landData land = pLand->getLandData();
	trfrSMSTraits movt = pSpecies->getSpSMSTraits();
	currLoc = pCurrCell->getLocn();

	// Get directional persistence weights
	float directionalPersistence = indvar ? pSMS.dp : movt.dp;
	if ((path->out > 0 && path->out <= (movt.pr + 1))
		|| isInNatalPatch
		|| (movt.straightenPath && path->settleStatus > 0)) {
		// inflate directional persistence to help leaving the patch
		directionalPersistence *= 10.0;
	}
	neighbourWeights = getSimDirection(currLoc.x, currLoc.y, directionalPersistence);

	if (isInNatalPatch || path->settleStatus > 0) path->out = 0;

	// Get goal bias weights
	double gb;
	if (movt.goalType == 2) { // dispersal bias
		int nsteps = path->year == path->total ? 
			path->out : // first year of dispersal - use no. of steps outside natal patch
			path->total; // use total no. of steps

		float goalBias = indvar ? pSMS.gb : movt.gb;
		float alphaDB = indvar ? pSMS.alphaDB : movt.alphaDB;
		int betaDB = indvar ? pSMS.alphaDB : movt.betaDB;
		double expArg = -(nsteps - betaDB) * -alphaDB;
		if (expArg > 100.0) expArg = 100.0;
		gb = 1.0 + (goalBias - 1.0) / (1.0 + exp(expArg));
	}
	else gb = movt.gb;
	goalBiasWeights = getGoalBias(currLoc.x, currLoc.y, movt.goalType, gb);

	// Get habitat-dependent weights (mean effective costs, given perceptual range)
	habDepWeights = pCurrCell->getEffCosts();
	if (habDepWeights.cell[0][0] >= 0.0) { 
		// already calculated in previous step, skip
	} else { 
		habDepWeights = getHabMatrix(pLand, currLoc.x, currLoc.y, movt.pr, movt.prMethod,
			landIx, absorbing);
		pCurrCell->setEffCosts(habDepWeights);
	}

	// Determine effective costs for the 8 neighbours
	for (int y = 2; y > -1; y--) { // N to S
		for (int x = 0; x < 3; x++) { // W to E
			if (x == 1 && y == 1) // current cell
				neighbourWeights.cell[x][y] = 0.0;
			else {
				float stepDist = (x == 1 || y == 1) ? 1.0 : SQRT2; // adjacent or diagonal?
				neighbourWeights.cell[x][y] = stepDist 
					* neighbourWeights.cell[x][y] 
					* goalBiasWeights.cell[x][y] 
					* habDepWeights.cell[x][y];
			}
		}
	}

	// determine reciprocal of effective cost for the 8 neighbours
	for (int y = 2; y > -1; y--) {
		for (int x = 0; x < 3; x++) {
			if (neighbourWeights.cell[x][y] > 0.0)
				neighbourWeights.cell[x][y] = 1.0 / neighbourWeights.cell[x][y];
		}
	}

	// Dismiss cells outside landscape or no-data cells
	for (int y = 2; y > -1; y--) {
		for (int x = 0; x < 3; x++) {
			if (!absorbing) {
				int neighbourX = currLoc.x + x - 1;
				int neighbourY = currLoc.y + y - 1;
				if (!isInLandscape(neighbourX, neighbourY, land))
					// cell is beyond current landscape limits
					neighbourWeights.cell[x][y] = 0.0;
				else if (pLand->findCell(neighbourX, neighbourY) == nullptr) {
						neighbourWeights.cell[x][y] = 0.0; // no-data cell
				}
			}
			// Increment total for re-scaling to sum to unity
			sum_nbrs += neighbourWeights.cell[x][y];
		}
	}

	// Scale effective costs as probabilities summing to 1
	if (sum_nbrs <= 0.0) 
		throw runtime_error("SMS probabilities have summed to zero or less.");
	else {
		for (int y = 2; y > -1; y--) {
			for (int x = 0; x < 3; x++) {
				neighbourWeights.cell[x][y] = neighbourWeights.cell[x][y] / sum_nbrs;
			}
		}
	}

	// Set up cell selection probabilities
	double cumulative[9];
	int j = 0;
	cumulative[0] = neighbourWeights.cell[0][0];
	for (int y = 0; y < 3; y++) {
		for (int x = 0; x < 3; x++) {
			if (j != 0) cumulative[j] = cumulative[j - 1] + neighbourWeights.cell[x][y];
			j++;
		}
	}
	// to prevent very rare bug that random draw is greater than 0.999999999
	if (cumulative[8] != 1) cumulative[8] = 1;

	// Draw direction from selection probabilities
	// landscape boundaries and no-data cells may be reflective or absorbing
	cellcost = pCurrCell->getCost();
	int loopsteps = 0;
	constexpr int maxLoopSteps = 1000;
	do {
		do {
			double rnd = pRandom->Random();
			j = 0;

			for (int y = 0; y < 3; y++) {
				for (int x = 0; x < 3; x++) {

					if (rnd < cumulative[j]) {
						newX = currLoc.x + x - 1;
						newY = currLoc.y + y - 1;
						move.dist = (x == 1 || y == 1) ? land.resol : SQRT2;
						y = x = 999; // break from x and y loops
					}
					j++;

				}
			}
			loopsteps++;

		} while (loopsteps < maxLoopSteps
			&& (!absorbing && !isInLandscape(newX, newY, land)));

		if (loopsteps >= maxLoopSteps) pNewCell = nullptr;
		else {
			if (!isInLandscape(newX, newY, land)) {
				pNewCell = nullptr;
			}
			pNewCell = pLand->findCell(newX, newY);
		}

	} while (!absorbing 
		&& pNewCell == nullptr 
		&& loopsteps < maxLoopSteps);

	if (loopsteps >= maxLoopSteps || pNewCell == nullptr) {
		// unable to make a move or crossed absorbing boundary
		// flag individual to die
		move.dist = -123.0;
		if (pNewCell == nullptr) pCurrCell = pNewCell;
	}
	else {
		newcellcost = pNewCell->getCost();
		move.cost = move.dist * 0.5 * (cellcost + newcellcost);
		// make the selected move
		if (memory.size() == movt.memSize) {
			memory.pop(); // remove oldest memory element
		}
		memory.push(currLoc); // record previous location in memory
		pCurrCell = pNewCell;
	}
	return move;
}

// Weight neighbouring cells on basis of current movement direction
array3x3d Individual::getSimDirection(const int x, const int y, const float dirPersistence)
{
	array3x3d neighbourWeights;

	if (memory.empty()) { // no previous movement, set matrix to unity
		for (int xx = 0; xx < 3; xx++) {
			for (int yy = 0; yy < 3; yy++) {
				neighbourWeights.cell[xx][yy] = 1;
			}
		}
	}
	else { 
		// set up the matrix dependent on relationship of previous location to current
		neighbourWeights.cell[1][1] = 0;
		locn prevLoc = memory.front();
		if (x - prevLoc.x == 0 && y - prevLoc.y == 0) {
			// back to 'square 1' (first memory location) - use previous step direction only
			prevLoc = memory.back();
			if ((x - prevLoc.x) == 0 && (y - prevLoc.y) == 0) { // STILL HAVE A PROBLEM!
				for (int xx = 0; xx < 3; xx++) {
					for (int yy = 0; yy < 3; yy++) {
						neighbourWeights.cell[xx][yy] = 1.0;
					}
				}
				return neighbourWeights;
			}
		}
		double theta = atan2(x - prevLoc.x, y - prevLoc.y);
		neighbourWeights = calcWeightings(dirPersistence, theta);
	}
	return neighbourWeights;
}

// Weight neighbouring cells on basis of goal bias
//array3x3d Individual::getGoalBias(const int x,const int y,
//	const int goaltype,const float gb)
array3x3d Individual::getGoalBias(const int x, const int y,
	const int goaltype, const float gb)
{
	array3x3d neighbourGB;
	double theta;
	int xx, yy;
	auto& pSMS = dynamic_cast<const smsData&>(*pTrfrData);

	if (goaltype == 0) { // no goal set
		for (xx = 0; xx < 3; xx++) {
			for (yy = 0; yy < 3; yy++) {
				neighbourGB.cell[xx][yy] = 1.0;
			}
		}
	}
	else {
		neighbourGB.cell[1][1] = 0;
		if ((x - pSMS.goal.x) == 0 && (y - pSMS.goal.y) == 0) {
			// at goal, set matrix to unity
			for (xx = 0; xx < 3; xx++) {
				for (yy = 0; yy < 3; yy++) {
					neighbourGB.cell[xx][yy] = 1.0;
				}
			}
			return neighbourGB;
		}
		theta = atan2((x - pSMS.goal.x), (y - pSMS.goal.y));
		neighbourGB = calcWeightings(gb, theta);
	}
	return neighbourGB;
}

// Calculate weightings for neighbouring cells
array3x3d Individual::calcWeightings(const double base, const double theta) {

	array3x3d d; // 3x3 array indexed from SW corner by xx and yy
	int dx, dy, xx, yy;

	double i0 = 1.0;	// direction of theta - lowest cost bias
	double i1 = base;
	double i2 = base * base;
	double i3 = i2 * base;
	double i4 = i3 * base;		// opposite to theta - highest cost bias

	if (fabs(theta) > 7.0 * PI / 8.0) { 
		dx = 0; 
		dy = -1; 
	}
	else {
		if (fabs(theta) > 5.0 * PI / 8.0) { 
			dy = -1; 
			if (theta > 0) dx = 1; 
			else dx = -1; 
		}
		else {
			if (fabs(theta) > 3.0 * PI / 8.0) { 
				dy = 0; 
				if (theta > 0) dx = 1; 
				else dx = -1; 
			}
			else {
				if (fabs(theta) > PI / 8.0) { 
					dy = 1; 
					if (theta > 0) dx = 1; 
					else dx = -1; 
				}
				else { 
					dy = 1; 
					dx = 0; 
				}
			}
		}
	}

	d.cell[1][1] = 0; // central cell has zero weighting
	d.cell[dx + 1][dy + 1] = i0;
	d.cell[-dx + 1][-dy + 1] = i4;
	if (dx == 0 || dy == 0) { // theta points to a cardinal direction
		d.cell[dy + 1][dx + 1] = i2; 
		d.cell[-dy + 1][-dx + 1] = i2;
		if (dx == 0) { // theta points N or S
			xx = dx + 1; 
			if (xx > 1) dx -= 2; yy = dy;
			d.cell[xx + 1][yy + 1] = i1; 
			d.cell[-xx + 1][yy + 1] = i1;
			d.cell[xx + 1][-yy + 1] = i3; 
			d.cell[-xx + 1][-yy + 1] = i3;
		}
		else { // theta points W or E
			yy = dy + 1; 
			if (yy > 1) dy -= 2; xx = dx;
			d.cell[xx + 1][yy + 1] = i1; 
			d.cell[xx + 1][-yy + 1] = i1;
			d.cell[-xx + 1][yy + 1] = i3; 
			d.cell[-xx + 1][-yy + 1] = i3;
		}
	}
	else { // theta points to an ordinal direction
		d.cell[dx + 1][-dy + 1] = i2; 
		d.cell[-dx + 1][dy + 1] = i2;
		xx = dx + 1; 
		if (xx > 1) xx -= 2; 
		d.cell[xx + 1][dy + 1] = i1;
		yy = dy + 1; 
		if (yy > 1) yy -= 2;
		d.cell[dx + 1][yy + 1] = i1;
		d.cell[-xx + 1][-dy + 1] = i3; 
		d.cell[-dx + 1][-yy + 1] = i3;
	}

	return d;
}

// Weight neighbouring cells on basis of (habitat) costs
array3x3f Individual::getHabMatrix(Landscape* pLand, 
	const int currCellX, const int currCellY, const short percRange, 
	const short prMethod, const short landIx, const bool absorbing)
{
	array3x3f neighbourHabWeights; // array of effective costs to be returned
	int ncells, x4, y4;
	double weight, sumweights;
	// NW and SE corners of effective cost array relative to the current cell (x,y):
	int xmin = 0, ymin = 0, xmax = 0, ymax = 0;
	int cost, nodatacost, h;
	Cell* pCell;

	landData land = pLand->getLandData();
	if (absorbing) nodatacost = gAbsorbingNoDataCost;
	else nodatacost = gNoDataCost;

	for (int x = -1; x < 2; x++) {   // index of relative move in x direction
		for (int y = -1; y < 2; y++) { // index of relative move in x direction

			neighbourHabWeights.cell[x + 1][y + 1] = 0.0; // initialise costs array to zeroes

			// set up corners of perceptual range relative to current cell
			if (x == 0 && y == 0) { // current cell - do nothing
				xmin = 0; ymin = 0; 
				xmax = 0; ymax = 0;
			}
			else {
				if (x == 0 || y == 0) { // not diagonal (rook move)
					if (x == 0) { // vertical (N-S) move
						if (percRange % 2 == 0) { // PR even
							xmin = -percRange / 2; 
							xmax = percRange / 2; 
							ymin = y; 
							ymax = y * percRange; 
						} else { // PR odd
							xmin = -(percRange - 1) / 2; 
							xmax = (percRange - 1) / 2; 
							ymin = y;
							ymax = y * percRange; 
						} 
					}
					if (y == 0) { // horizontal (E-W) move
						if (percRange % 2 == 0) { // PR even
							xmin = x; 
							xmax = x * percRange; 
							ymin = -percRange / 2; 
							ymax = percRange / 2; 
						} else { // PR odd
							xmin = x; 
							xmax = x * percRange; 
							ymin = -(percRange - 1) / 2; 
							ymax = (percRange - 1) / 2; 
						}
					}
				}
				else { // diagonal (bishop move)
					xmin = x; 
					xmax = x * percRange; 
					ymin = y; 
					ymax = y * percRange;
				}
			}
			if (xmin > xmax) { 
				int z = xmax; 
				xmax = xmin; 
				xmin = z; 
			} // swap xmin and xmax
			if (ymin > ymax) { 
				int z = ymax; 
				ymax = ymin; 
				ymin = z; 
			} // swap ymin and ymax

			// calculate effective mean cost of cells in perceptual range
			ncells = 0; 
			weight = 0.0; 
			sumweights = 0.0;

			if (x != 0 || y != 0) { // not central cell (i.e. current cell)
				for (int x3 = xmin; x3 <= xmax; x3++) {
					for (int y3 = ymin; y3 <= ymax; y3++) {
						// if cell is out of bounds, treat landscape as a torus
						// for purpose of obtaining a cost,
						if ((currCellX + x3) < 0) 
							x4 = currCellX + x3 + land.maxX + 1;
						else { 
							if ((currCellX + x3) > land.maxX)
								x4 = currCellX + x3 - land.maxX - 1;
							else 
								x4 = currCellX + x3;
						}
						if ((currCellY + y3) < 0)
							y4 = currCellY + y3 + land.maxY + 1;
						else { 
							if ((currCellY + y3) > land.maxY)
								y4 = currCellY + y3 - land.maxY - 1;
							else y4 = currCellY + y3;
						}
						if (x4 < 0 || x4 > land.maxX || y4 < 0 || y4 > land.maxY) {
							// unexpected problem - e.g. due to ridiculously large PR
							// treat as a no-data cell
							cost = nodatacost;
						}
						else {
							// add cost of cell to total PR cost
							pCell = pLand->findCell(x4, y4);
							if (pCell == nullptr) { // no-data cell
								cost = nodatacost;
							}
							else {
								cost = pCell->getCost();
								if (cost < 0) cost = nodatacost;
								else if (cost == 0) { // cost not yet set for the cell
									h = pCell->getHabIndex(landIx);
									cost = pSpecies->getHabCost(h);
									pCell->setCost(cost);
								}
							}
						}
						if (prMethod == 1) { // arithmetic mean
							neighbourHabWeights.cell[x + 1][y + 1] += cost;
							ncells++;
						}
						if (prMethod == 2) { // harmonic mean
							if (cost > 0) {
								neighbourHabWeights.cell[x + 1][y + 1] += (1.0 / cost);
								ncells++;
							}
						}
						if (prMethod == 3) { // arithmetic mean weighted by inverse distance
							if (cost > 0) {
								// NB distance is still given by (x3,y3)
								weight = 1.0f / sqrt((pow(x3, 2) + pow(y3, 2)));
								neighbourHabWeights.cell[x + 1][y + 1] += weight * cost;
								ncells++; 
								sumweights += weight;
							}
						}

					} // end of y3 loop
				}  // end of x3 loop

				if (ncells > 0) {
					if (prMethod == 1) neighbourHabWeights.cell[x + 1][y + 1] /= ncells; // arithmetic mean
					if (prMethod == 2) neighbourHabWeights.cell[x + 1][y + 1] = ncells / neighbourHabWeights.cell[x + 1][y + 1]; // hyperbolic mean
					if (prMethod == 3 && sumweights > 0)
						neighbourHabWeights.cell[x + 1][y + 1] /= sumweights; // weighted arithmetic mean
				}
			}
			else { // central cell
				// record cost if not already recorded
				// has effect of preparing for storing effective costs for the cell
				pCell = pLand->findCell(currCellX, currCellY);
				cost = pCell->getCost();
				if (cost < 0) cost = nodatacost;
				else if (cost == 0) { // cost not yet set for the cell
					h = pCell->getHabIndex(landIx);
					cost = pSpecies->getHabCost(h);
					pCell->setCost(cost);
				}
			}

		} //end of y loop
	} //end of x loop

	return neighbourHabWeights;

}

#if RS_RCPP
//---------------------------------------------------------------------------
// Write records to movement paths file
void Individual::outMovePath(const int year)
{
	locn loc, prev_loc;

	//if (pPatch != pNatalPatch) {
	loc = pCurrCell->getLocn();
	// if still dispersing...
	if (status == dispersing) {
		// at first step, record start cell first
		if (path->total == 1) {
			prev_loc = pPrevCell->getLocn();
			outMovePaths << year << "\t" << indId << "\t"
				<< "0\t" << prev_loc.x << "\t" << prev_loc.y << "\t"
				<< "0\t"	// status at start cell is 0
				<< endl;
		}
		// then record current step
		outMovePaths << year << "\t" << indId << "\t"
			<< path->total << "\t" << loc.x << "\t" << loc.y << "\t"
			<< to_string(status) << "\t"
			<< endl;
	}
	// if not anymore dispersing...
	if (status != initial && status != dispersing) {
		prev_loc = pPrevCell->getLocn();
		// record only if this is the first step as non-disperser
		if (path->pathoutput) {
			// if this is also the first step taken at all, record the start cell first
			if (path->total == 1) {
				outMovePaths << year << "\t" << indId << "\t"
					<< "0\t" << prev_loc.x << "\t" << prev_loc.y << "\t"
					<< "0\t"	// status at start cell is 0
					<< endl;
			}
			outMovePaths << year << "\t" << indId << "\t"
				<< path->total << "\t" << loc.x << "\t" << loc.y << "\t"
				<< to_string(status) << "\t"
				<< endl;
			// current cell will be invalid (zero), so set back to previous cell
			//pPrevCell = pCurrCell;
			path->pathoutput = 0;
		}
	}
}
#endif

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

double drawDirection(double location, double rho) {
	if (rho < 0.0 || rho > 1.0) return location;
	else if (rho == 0) return pRandom->Random() * M_2PI;
	else if (rho == 1) return location;
	else return fmod(cauchy(location, -log(rho)), M_2PI);
}

double cauchy(double location, double scale) {
	if (scale < 0) return location;
	else return location + scale * tan(PI * pRandom->Random());
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

#ifndef NDEBUG
// Testing utilities

Cell* Individual::getCurrCell() const {
	return pCurrCell;
}

void Individual::setInitAngle(const float angle) {
	auto pCRW = dynamic_cast<crwData*>(pTrfrData.get());
	pCRW->prevdrn = angle;
}

// Force mutations to trigger for all traits
void Individual::triggerMutations(Species* pSp) {
	for (auto const& [trType, indTrait] : spTraitTable) {
		indTrait->mutate();
		if (trType == GENETIC_LOAD1 
			|| trType == GENETIC_LOAD2
			|| trType == GENETIC_LOAD3
			|| trType == GENETIC_LOAD4
			|| trType == GENETIC_LOAD5)
			geneticFitness *= indTrait->express();
	}
	this->expressDispersalPhenotypes(1);
}

// Shorthand function to edit a genotype with custom values
void Individual::overrideGenotype(TraitType whichTrait, const map<int, vector<shared_ptr<Allele>>>& newGenotype) {

	GeneticFitnessTrait* pGenFitTrait;
	DispersalTrait* pDispTrait;

	switch (whichTrait)
	{
	case GENETIC_LOAD1: case GENETIC_LOAD2: case GENETIC_LOAD3: case GENETIC_LOAD4: case GENETIC_LOAD5: 
		pGenFitTrait = dynamic_cast<GeneticFitnessTrait*>(this->getTrait(whichTrait));
		pGenFitTrait->getGenes() = newGenotype;
		break;
	case E_D0: case E_ALPHA: case E_BETA:
	case S_S0: case S_ALPHA: case S_BETA:
	case E_D0_F: case E_ALPHA_F: case E_BETA_F:
	case S_S0_F: case S_ALPHA_F: case S_BETA_F: 
	case E_D0_M: case E_ALPHA_M: case E_BETA_M: 
	case S_S0_M: case S_ALPHA_M: case S_BETA_M: 
	case CRW_STEPLENGTH: case CRW_STEPCORRELATION: 
	case KERNEL_MEANDIST_1: case KERNEL_MEANDIST_2: case KERNEL_PROBABILITY: 
	case KERNEL_MEANDIST_1_F: case KERNEL_MEANDIST_2_F: case KERNEL_PROBABILITY_F: 
	case KERNEL_MEANDIST_1_M: case KERNEL_MEANDIST_2_M: case KERNEL_PROBABILITY_M: 
	case SMS_DP: case SMS_GB: case SMS_ALPHADB: case SMS_BETADB:
		pDispTrait = dynamic_cast<DispersalTrait*>(this->getTrait(whichTrait));
		pDispTrait->getGenes() = newGenotype;
		break;
	default:
		throw logic_error("Wrong trait type: please choose a valid dispersal or genetic fitness trait.");
		break;
	}
};

void Individual::overrideGenotype(TraitType whichTrait, const map<int, vector<unsigned char>>& newGenotype) {

	if (!whichTrait == NEUTRAL) {
		throw logic_error("Attempt to override non-neutral trait with neutral trait genotype.\n");
	}
	NeutralTrait* pNeutralTrait;
	pNeutralTrait = dynamic_cast<NeutralTrait*>(this->getTrait(NEUTRAL));
	pNeutralTrait->getGenes() = newGenotype;
};

#endif // NDEBUG

