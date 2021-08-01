//#include <boost/serialization/export.hpp>
#include "CancerCell.h"

//BOOST_CLASS_EXPORT_IMPLEMENT(CancerCell)
#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"

#include "../compartment/Tumor.h"
#include "../compartment/LymphCentral.h"
#include "../../core/GlobalUtilities.h"


namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

using std::string;
using std::stringstream;

static int CancerCellSize = 1;

CancerCell::CancerCell()
	:_count_neighbor_Teff(0)
{
}

CancerCell::CancerCell(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _stemID(getID())
	,_divideCD(int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE)) + .5)
	, _divideFlag(true)
	, _divideCountRemaining(0)
	, _count_neighbor_Teff(0)
	, _source_CCL2(NULL)
    , _antigen_vec( {} ) //Make a function that returns the bool vector indicates the presence of antigens.
{
	_state = AgentStateEnum::CANCER_STEM;
}


CancerCell::CancerCell(const CancerCell& c)
	:Cell_Tumor(c)
	, _stemID(getID())
	,_divideCD(c._divideCD)
	, _divideFlag(c._divideFlag)
	, _divideCountRemaining(c._divideCountRemaining)
	, _count_neighbor_Teff(0)
	, _source_CCL2(NULL)
    , _antigen_vec (c._antigen_vec)
{
	setup_chem_source(_source_CCL2, CHEM_CCL2, params.getVal(PARAM_CCL2_RELEASE));
}

CancerCell::~CancerCell()
{
}

string CancerCell::toString()const{
	stringstream ss;
	ss << CellAgent::toString();
	ss << "division flag: " << _divideFlag << ", division cool down: " << _divideCD << std::endl;
	return ss.str();
}

bool CancerCell::agent_movement_step(double t, double dt, Coord& c){

	bool move = false;
	double pMove = getState() == CANCER_STEM ?
		params.getVal(PARAM_CANCER_STEM_MOVE_PROB) :
		params.getVal(PARAM_CANCER_CELL_MOVE_PROB);
	if (rng.get_unif_01() < pMove)
	{
		// move
		int idx;
		const auto shape = getCellShape();
		if (_compartment->getOneOpenVoxel(shape->getMoveDestinationVoxels(),
			shape->getMoveDirectionAnchor(), _coord, getType(), idx, rng))
		{
			move = true;
			c = getCellShape()->getMoveDirectionAnchor()[idx] + _coord;
		}
	}
	return move;
}
bool CancerCell::agent_state_step(double t, double dt, Coord& c){

	bool divide = false;
	if (!isDead() && _state == AgentStateEnum::CANCER_SENESCENT)
	{
		_life -= 1;
		if (_life <= 0)
		{
			setDead();
		}
	}

	const auto shape = getCellShape();

	// check IFNg
	Cell_Tumor::agent_state_step(t, dt, c);

	auto tumor = dynamic_cast<Tumor*>(_compartment);
	double nivo = tumor->get_Nivo();
	double ent = tumor->get_Ent();
	// decide if get killed by Tcyt
	//if (_state != AgentStateEnum::CANCER_STEM)
	{
		if (_count_neighbor_Teff > 0){
			// count neighborhood cancer cells
			int neighbor_cc = 0;
			_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
				_coord, [&](BaseAgent* ag){
				if (ag->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
					neighbor_cc += 1;
				}
				return true;
			});
			double bond = TCell::get_PD1_PDL1(_PDL1_syn, nivo);
			double supp = TCell::get_PD1_supp(bond, params.getVal(PARAM_N_PD1_PDL1));
			double NO = get_tumor().get_chem(c, CHEM_NO);
			double ArgI = get_tumor().get_chem(c, CHEM_ARGI);
			double H_mdsc_c1 = 1 - ((1 - (NO/(NO+params.getVal(PARAM_IC50_NO_CTL))))*(1-(ArgI/(ArgI+params.getVal(PARAM_IC50_ARGI_CTL))))*(1-(ent/(ent+params.getVal(PARAM_IC50_ENT_ARGI)))));
			double q = (1 - H_mdsc_c1) * double(_count_neighbor_Teff) / (_count_neighbor_Teff + neighbor_cc + params.getVal(PARAM_CELL));
			double p_kill = TCell::get_kill_prob(supp, q);
			/*
			std::cout << "T cell killing:\n"
				<< "neighbors: " << _count_neighbor_Teff << "," << neighbor_cc<< "," << q << std::endl;
			std::cout 
				<< "PDL1: " << _PDL1_syn << ", bond: " << bond << ", supp: " << supp << "\n"
				<< "pkill: " << p_kill 
				<< std::endl;
				*/
			if (rng.get_unif_01() < p_kill){
				setDead();
				tumor->inc_abm_var_exchange(Tumor::TUMEX_CC_T_KILL);
				//total_kill += 1;
				/*
				if (_state = AgentStateEnum::CANCER_PROGENITOR)
				{
					std::cout << "total kill: " << total_kill << std::endl;
				}*/
				//std::cout << "Killed: " << _state 
				// << ": " << p_kill << std::endl;
			}
		}
	}

	if (isDead())
	{
		return divide;
	}

	//printCellInfo2(t, this ,"before divide/move");
	// divide
	if (_divideCD > 0)
	{
		_divideCD--;
	}

	double H_ent_c1 = ent / (ent + params.getVal(PARAM_IC50_ENT_C));

	if (_divideFlag && _divideCD == 0)
	{
		// find location to divide
		int idx;
		if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(),
			shape->getProlifDestinationAnchor(), _coord, getType(), idx, rng))
		{
			divide = true;
			//std::cout << "Cancer cell: " << getID() << " Divided" << std::endl;
			//cout << "idx: " << idx << ", " << getCellShape()->getProlif()[idx] << endl;
			c = getCellShape()->getProlifDestinationAnchor()[idx] + _coord;

			//_divideFlag = true;
			if (_state == AgentStateEnum::CANCER_STEM)
			{
				_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) / (1-H_ent_c1) + .5);
			}
			else if (_state == AgentStateEnum::CANCER_PROGENITOR) {
				_divideCountRemaining -= 1;
				_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) / (1-H_ent_c1) + .5);
			}
		}
	}
	//printCellInfo2(t, this ,"after divide/move");

	// do other stuff
	// senescence 
	//std::cout << "Divide count remaining: " << _divideCountRemaining << ", " << _state  << "?=" << AgentStateEnum::CANCER_PROGENITOR << std::endl;
	if (_state == AgentStateEnum::CANCER_PROGENITOR && _divideCountRemaining == 0){
		//std::cout << "Cancer cell: " << getID() << ", turning senescent" << std::endl;
		setSenescent();
	}

	_count_neighbor_Teff = 0;

	return divide;
}

/*! change cell state of a CSC to progenitor.
	-# daughter cell from asymmetric division 
	-# during initialization (default state is stem)
*/
void CancerCell::setProgenitor() {
	_state = AgentStateEnum::CANCER_PROGENITOR;
	_divideCountRemaining = params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX);
	_divideCD = int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) + .5);
}

/*! set cancer cell to senescent
*/
void CancerCell::setSenescent(){
	_state = AgentStateEnum::CANCER_SENESCENT;
	_divideCD = -1;
	_divideFlag = false;
	_life = getSenescentLife();
	return;
}
/*! randomize division cooldown 
	Set to a number between [1, mean], uniform
*/
void CancerCell::randomize_div_cd(int mean){
	_divideCD = int(rng.get_unif_01()*mean) + 1;
	return;
}

/*! remaining slices for senescent cancer cells
	randomly drawn from exponential distribution.
*/
int CancerCell::getSenescentLife(void){
	double mean = params.getVal(PARAM_CANCER_SENESCENT_MEAN_LIFE);
	return int(rng.get_exponential(mean) + .5);
	//cout << "random cancer cell life: " << cLife << endl;
}

//! remarks to string
std::string CancerCell::getRemark() const{
	stringstream ss;
	ss << Cell_Tumor::getRemark()
		<< "|" << _stemID
		<< "|" << _divideCD
		<< "|" << _divideCountRemaining;
	return ss.str();
}

//! move sources (CCL2)
void CancerCell::move_all_source_sink(void)const{
	//std::cout << "moving sources: " << _coord << std::endl;
	move_source_sink(_source_CCL2);
	return;
}

//! remove sources (CCL2)
void CancerCell::remove_all_source_sink(void){
	remove_source_sink(_source_CCL2);
	return;
}
};
};
