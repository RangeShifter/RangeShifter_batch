#include "GeneticLoad.h"

// ----------------------------------------------------------------------------------------
// for initialising population
// ----------------------------------------------------------------------------------------

GeneticLoad::GeneticLoad(SpeciesTrait* P)
{
	pSpeciesTrait = P;

	if (wildType.get() == nullptr)
		wildType = make_shared<Allele>(0.0, 0.0);

	ExpressionType expressionType = pSpeciesTrait->getExpressionType();

	_inherit_func_ptr = (pSpeciesTrait->getPloidy() == 1) ? &GeneticLoad::inheritHaploid : &GeneticLoad::inheritDiploid; //this could be changed if we wanted some alternative form of inheritance

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	map<parameter_t, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	switch (mutationDistribution) {
	case UNIFORM:
	{
		if (!mutationParameters.count(MAX))
			cout << endl << ("Error:: adaptive mutation uniform distribution parameter must contain max value (e.g. max= ) \n");

		if (!mutationParameters.count(MIN))
			cout << endl << ("Error:: adaptive mutation uniform distribution parameter must contain min value (e.g. min= ) \n");

		break;
	}
	case NORMAL:
	{

		if (!mutationParameters.count(MEAN))
			cout << endl << ("Error:: adaptive mutation distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");

		if (!mutationParameters.count(SDEV))
			cout << endl << ("Error:: adaptive mutation distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");

		break;
	}
	case GAMMA:
	{
		if (!mutationParameters.count(SHAPE))
			cout << endl << ("Error:: adaptive mutation distribution set to gamma so parameters must contain shape value (e.g. shape= ) \n");

		if (!mutationParameters.count(SCALE))
			cout << endl << ("Error:: adaptive mutation distribution set to gamma so parameters must contain scale value (e.g. scale= ) \n");

		break;
	}
	case NEGEXP:
	{
		if (!mutationParameters.count(MEAN))
			cout << endl << ("Error:: adaptive mutation distribution set to negative exponential (negative decay) so parameters must contain mean value (e.g. mean= ) \n");

		break;
	}

	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive mutation model, must be uniform/normal/gamma/negExp \n");
	}
	}

	DistributionType dominanceDistribution = pSpeciesTrait->getDominanceDistribution();
	map<parameter_t, float> dominanceParameters = pSpeciesTrait->getDominanceParameters();

	switch (dominanceDistribution) {
	case UNIFORM:
	{
		if (!dominanceParameters.count(MAX))
			cout << endl << ("Error:: adaptive dominance uniform distribution parameter must contain max value (e.g. max= ) \n");

		if (!dominanceParameters.count(MIN))
			cout << endl << ("Error:: adaptive dominance uniform distribution parameter must contain min value (e.g. min= ) \n");

		break;
	}
	case NORMAL:
	{

		if (!dominanceParameters.count(MEAN))
			cout << endl << ("Error:: adaptive dominance distribution set to normal so parameters must contain mean value (e.g. mean= ) \n");

		if (!dominanceParameters.count(SDEV))
			cout << endl << ("Error:: adaptive dominance distribution set to normal so parameters must contain sdev value (e.g. sdev= ) \n");

		break;
	}
	case GAMMA:
	{
		if (!dominanceParameters.count(SHAPE))
			cout << endl << ("Error:: adaptive dominance distribution set to gamma so parameters must contain shape value (e.g. shape= ) \n");

		if (!dominanceParameters.count(SCALE))
			cout << endl << ("Error:: adaptive dominance distribution set to gamma so parameters must contain scale value (e.g. scale= ) \n");

		break;
	}
	case NEGEXP:
	{
		if (!dominanceParameters.count(MEAN))
			cout << endl << ("Error:: adaptive dominance distribution set to negative exponential (negative decay) so parameters must contain mean value (e.g. mean= ) \n");

		break;
	}
	case SCALED:
	{
		break;
	}

	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive dominance model, must be uniform/normal/gamma/negExp/scaled \n");
		break; //should return false
	}
	}

	DistributionType initialDistribution = pSpeciesTrait->getInitialDistribution();
	map<parameter_t, float> initialParameters = pSpeciesTrait->getInitialParameters();

	//switch (expressionType) {
	//case MULTIPLICATIVE:
	//{
	//    _express_func_ptr = &AdaptiveTrait::expressMulti;
	//    break;
	//}
	//default:
	//{
	//    cout << endl << ("wrong parameter value for parameter \"expression of adaptive mutations\", must be multiplicative  \n");
	//    break; //should return false
	//}

	//}

}


// ----------------------------------------------------------------------------------------
// for creating new individuals
// ----------------------------------------------------------------------------------------

GeneticLoad::GeneticLoad(const GeneticLoad& T) : pSpeciesTrait(T.pSpeciesTrait), _inherit_func_ptr(T._inherit_func_ptr)
{}

// ----------------------------------------------------------------------------------------
// mutate uniform
// ----------------------------------------------------------------------------------------
void GeneticLoad::mutate()
{
	const int positionsSize = pSpeciesTrait->getPositionsSize();
	const auto& positions = pSpeciesTrait->getPositions();
	const short ploidy = pSpeciesTrait->getPloidy();
	const float mutationRate = pSpeciesTrait->getMutationRate();

	auto rng = pRandom->getRNG();

	for (int p = 0; p < ploidy; p++) {

		// Determine nb of mutations
		unsigned int NbMut = pRandom->Poisson(positionsSize * mutationRate);

		if (NbMut > 0) {
			vector<int> mutationPositions;
			// Draw which positions mutate
			sample(positions.begin(), positions.end(), std::back_inserter(mutationPositions),
				NbMut, rng);

			for (int m : mutationPositions) {

				float newSelectionCoef = drawSelectionCoef();
				float newDominanceCoef = drawDominance(newSelectionCoef);

				auto it = genes.find(m); //find if position in map already has mutations there

				if (it == genes.end()) {   // not found so create new entry in map with wildtype as char default
					vector<shared_ptr<Allele>> newAllelePair(2);
					newAllelePair[p] = make_shared<Allele>(newSelectionCoef, newDominanceCoef); //put new mutation value in 
					genes.insert(make_pair(m, newAllelePair));
				}
				else { //position found, already mutations there
					it->second[p] = make_shared<Allele>(newSelectionCoef, newDominanceCoef);
				}
			}
		}
	}
}
// ----------------------------------------------------------------------------------------
// get dominance value for new mutation
// ----------------------------------------------------------------------------------------


float GeneticLoad::drawDominance(float selCoef) {

	DistributionType dominanceDistribution = pSpeciesTrait->getDominanceDistribution();
	map<parameter_t, float> dominanceParameters = pSpeciesTrait->getDominanceParameters();

	float h = 1.0; //default dominance is  1

	switch (dominanceDistribution) {
	case UNIFORM:
	{
		float maxD = dominanceParameters.find(MAX)->second;
		float minD = dominanceParameters.find(MIN)->second;
		h = pRandom->FRandom(minD, maxD);
		break;
	}
	case NORMAL:
	{
		const float mean = dominanceParameters.find(MEAN)->second;
		const float sd = dominanceParameters.find(SDEV)->second;
		h = pRandom->Normal(mean, sd);
		break;
	}
	case GAMMA:
	{
		const float shape = dominanceParameters.find(SHAPE)->second;
		const float scale = dominanceParameters.find(SCALE)->second;
		h = pRandom->Gamma(shape, scale);
		break;
	}
	case NEGEXP:
	{
		const float mean = dominanceParameters.find(MEAN)->second;
		h = pRandom->NegExp(mean);
		break;
	}
	case SCALED:
	{
		const float min = 0;
		const float max = exp((-log(2 * 0.36) / 0.05) * selCoef);
		h = pRandom->FRandom(min, max);
		break;
	}

	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive dominance model, must be uniform/normal/gamma/negExp/scaled \n");
		break; //should return false
	}
	}

	return h;
}

// ----------------------------------------------------------------------------------------
// get selection coefficient for new mutation
// ----------------------------------------------------------------------------------------


float GeneticLoad::drawSelectionCoef() {

	DistributionType mutationDistribution = pSpeciesTrait->getMutationDistribution();
	map<parameter_t, float> mutationParameters = pSpeciesTrait->getMutationParameters();

	float s = 0.0; //default selection coefficient is 0

	switch (mutationDistribution) {
	case UNIFORM:
	{
		float maxD = mutationParameters.find(MAX)->second;
		float minD = mutationParameters.find(MIN)->second;
		s = pRandom->FRandom(minD, maxD);

		break;
	}
	case NORMAL:
	{
		const float mean = mutationParameters.find(MEAN)->second;
		const float sd = mutationParameters.find(SDEV)->second;
		s = pRandom->Normal(mean, sd);

		break;
	}
	case GAMMA:
	{
		const float shape = mutationParameters.find(SHAPE)->second;
		const float scale = mutationParameters.find(SCALE)->second;
		s = pRandom->Gamma(shape, scale);
		break;
	}
	case NEGEXP:
	{
		const float mean = mutationParameters.find(MEAN)->second;
		s = pRandom->NegExp(mean);
		break;
	}
	default:
	{
		cout << endl << ("Error:: wrong parameter value for adaptive mutation model, must be uniform/normal/gamma/negExp/scaled \n");
		break; //should return false
	}
	}
	return s;
}


// ----------------------------------------------------------------------------------------
// inheritance options
// ----------------------------------------------------------------------------------------


void GeneticLoad::inherit(TTrait* parentTrait, set<unsigned int> const& recomPositions, sex_t whichChromosome, int startingChromosome)
{
	auto parentCast = dynamic_cast<GeneticLoad*> (parentTrait); //horrible

	const auto& parent_seq = parentCast->get_mutations();
	if (parent_seq.size() > 0) //else nothing to inherit
		(this->*_inherit_func_ptr) (whichChromosome, parent_seq, recomPositions, startingChromosome);
}

void GeneticLoad::inheritDiploid(sex_t whichChromosome, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome) {

	auto it = recomPositions.lower_bound(parentGenes.begin()->first);

	unsigned int nextBreakpoint = *it;

	auto distance = std::distance(recomPositions.begin(), it);
	if (distance % 2 != 0)
		parentChromosome = !parentChromosome; // switch to the other one
		// use 1-parentChromosome, or switch to a sex_t ?

	for (auto const& [locus, allelePair] : parentGenes) {

		while (locus > nextBreakpoint) {
			std::advance(it, 1);
			nextBreakpoint = *it;
			parentChromosome = !parentChromosome; // switch to the other one
		}

		if (locus <= nextBreakpoint) {
			auto& allele = allelePair[parentChromosome];
			auto it = genes.find(locus);
			if (it == genes.end()) {
				// locus does not exist yet, initiate it
				vector<shared_ptr<Allele>> newAllelePair(2);
				newAllelePair[whichChromosome] = allele;
				genes.insert(make_pair(locus, newAllelePair));
			}
			else {
				// locus already exists
				// set corresponding allele
				it->second[whichChromosome] = allele;
			}
		}
	}
}

void GeneticLoad::inheritHaploid(sex_t chromosome, map<int, vector<shared_ptr<Allele>>> const& parentGenes, set<unsigned int> const& recomPositions, int parentChromosome)
{
	genes = parentGenes;
}

// ----------------------------------------------------------------------------------------
// expression options
// ----------------------------------------------------------------------------------------

float GeneticLoad::express() {

	float phenotype = 1.0;

	for (auto const& [locus, pAllelePair] : genes)
	{
		auto pAlleleLeft  = (!pAllelePair[0]) ? wildType : pAllelePair[0];
		auto pAlleleRight = (!pAllelePair[1]) ? wildType : pAllelePair[1];

		if (pAlleleLeft.get()->getId() != pAlleleRight.get()->getId()) // heterozygote
		{
			phenotype *= 1 + pAlleleLeft->getAlleleValue()  * pAlleleLeft->getDominanceCoef();
			phenotype *= 1 + pAlleleRight->getAlleleValue() * pAlleleRight->getDominanceCoef();
		}
		else { // homozygote
			phenotype *= 1 + pAlleleLeft->getAlleleValue();
			phenotype *= 1 + pAlleleRight->getAlleleValue();
		}
	}
	return phenotype;
}

// ----------------------------------------------------------------------------------------
// check if particular locus is heterozygote
// ----------------------------------------------------------------------------------------


bool GeneticLoad::isHeterozygoteAtLocus(int locus) const {

	auto it = genes.find(locus);

	if (it == genes.end()) //not found so must be wildtype homozygous
		return false;
	else {
		auto a = (!it->second[0]) ? wildType : it->second[0];
		auto b = (!it->second[1]) ? wildType : it->second[1];
		return a != b;
	}
}

// ----------------------------------------------------------------------------------------
// count heterozygote loci in genome 
// ----------------------------------------------------------------------------------------


int GeneticLoad::countHeterozygoteLoci() const {

	int count = 0;

	for (auto const& [locus, allelePair] : genes) {
		auto alleleLeft = (!allelePair[0]) ? wildType : allelePair[0];
		auto alleleRight = (!allelePair[1]) ? wildType : allelePair[1];
		count += alleleLeft != alleleRight;
	}
	return count;
}

// ----------------------------------------------------------------------------------------
// get allele value at loci 
// ----------------------------------------------------------------------------------------


float GeneticLoad::getSelectionCoefAtLoci(short whichChromosome, int position) const {

	auto it = genes.find(position);

	if (it == genes.end()) {
		return wildType->getAlleleValue(); //must still be wildtype at loci
	}
	else
		return (!it->second[whichChromosome]) ? wildType->getAlleleValue() : it->second[whichChromosome]->getAlleleValue();
}