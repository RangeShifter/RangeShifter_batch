#ifndef SNPTRAITH
#define SNPTRAITH



#include <vector>
#include <string>

#include "TTrait.h"

using namespace std;

class SNPTrait : public TTrait {

private:

	inline static int wildType = -999;
	constexpr int maxSNPAlleles = UCHAR_MAX + 1; // allele is char, can take value 0-255

	SpeciesTrait* pProtoTrait;

	//vector<vector<char>> mutations; 

	map<int, vector<unsigned char>> genes; //position <strand A , strand B>>

	void (SNPTrait::* _mutate_func_ptr) (void);
	void (SNPTrait::* _inherit_func_ptr) (sex_t chromosome, map<int, vector<unsigned char>> const& parent, set<unsigned int> const& recomPositions, int parentChromosome);

	void inheritDiploid(sex_t chromosome, map<int, vector<unsigned char>> const&, set<unsigned int> const& recomPositions, int parentChromosome);

	void inheritHaploid(sex_t chromosome, map<int, vector<unsigned char>> const& parentMutations, set<unsigned int> const& recomPositions, int parentChromosome);

	void initialiseFull(int max); //other option is that mutations map is empty until a mutation happens, default when empty is to return a 0 value for wildtype

	void mutate_KAM();
	void mutate_SSM();




public:

	//this one for species held trait table, e.g. prototype table, sets static members

	SNPTrait(SpeciesTrait* P);

	//this one for individuals, static members are not reset
	SNPTrait(const SNPTrait& T);


	virtual unique_ptr<TTrait> clone() const override { return std::make_unique<SNPTrait>(*this); }


	virtual void  inherit(TTrait* parent, set<unsigned int> const& recomPositions, sex_t chromosome, int startingChromosome) override;


	virtual void    mutate() override { (this->*_mutate_func_ptr) (); }

	virtual int getNLoci()  const override { return pProtoTrait->getPositionsSize(); }


	float getMutationRate() const override { return pProtoTrait->getMutationRate(); }

	bool isInherited() const override { return pProtoTrait->isInherited(); }

	map<int, vector<unsigned char>>& get_mutations() { return genes; } //returning reference, reciever must be const

	virtual float getSelectionCoefAtLoci(short chromosome, int position) const override;

	virtual int countHeterozygoteLoci() const;

	virtual bool isHeterozygoteAtLocus(int locus) const override;

	virtual float express() { return -9999; }

	virtual ~SNPTrait() { }
};
#endif