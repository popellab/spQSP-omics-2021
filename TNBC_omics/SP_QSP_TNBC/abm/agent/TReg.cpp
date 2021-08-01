//#include <boost/serialization/export.hpp>
#include "TReg.h" 

//BOOST_CLASS_EXPORT_IMPLEMENT(Treg)

#include <iostream>
#include <sstream>

#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"
#include "TCell.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

Treg::Treg(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _divide_cd_Treg_exp(0)
	, _divide_limit_Treg_exp(params.getVal(PARAM_TREG_EXP_DIV_LIMIT))				
	//, _source_IL_10(NULL)
{
	_life = getTregLife();
}

Treg::Treg(const Treg& c)
	:Cell_Tumor(c)
	, _divide_cd_Treg_exp(c._divide_cd_Treg_exp)
	, _divide_limit_Treg_exp(c._divide_limit_Treg_exp)			
	//, _source_IL_10(NULL)
{
	_life = getTregLife();
	//setup_chem_source(_source_IL_10, CHEM_IL_10, params.getVal(PARAM_IL_10_RELEASE));
}

Treg::~Treg()
{
}

std::string Treg::toString()const{
	std::stringstream ss;
	ss << Cell_Tumor::toString();
	return ss.str();
}

/*void Treg::setDead(void)
{
	//remove_source_sink(_source_IL_10);
	CellAgent::setDead();
}*/

bool Treg::agent_movement_step(double t, double dt, Coord& c){
	bool move = false;
	/**/
	if (rng.get_unif_01() < params.getVal(PARAM_TREG_MOVE_PROB))
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

bool Treg::agent_state_step(double t, double dt, Coord& c){
	bool divide = false;
	if (!isDead())
	{
		_life--;
		if (_life == 0)
		{
			setDead();
			//std::cout << "Treg Life: " << getID() << " remaining is: " << _life << std::endl;
			//std::cout << "Treg Life: " << getID() << " remaining is dead" << std::endl;
			// remove source when cell die
			return divide;
		}
	}

	const auto shape = getCellShape();
	Cell_Tumor::agent_state_step(t, dt, c);

    auto tumor = dynamic_cast<Tumor*>(_compartment);
    double tumvol = tumor->get_Tum_Vol();
    double tregul = tumor->get_Treg();
		
	double ArgI = get_tumor().get_chem(c, CHEM_ARGI);
	if (ArgI > 0)
	{
		if (_divide_cd_Treg_exp > 0)
		{
			_divide_cd_Treg_exp--;
		}

	auto tumor = dynamic_cast<Tumor*>(_compartment);
	double ent = tumor->get_Ent();

		if (_divide_limit_Treg_exp > 0 && _divide_cd_Treg_exp == 0)
		{
			int idx;
			if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(), 
				shape->getProlifDestinationAnchor(), _coord, getType(), idx, rng))
			{
				divide = true;
				//cout << "idx: " << idx << ", " << getCellShape()->getProlif()[idx] << endl;
				c = getCellShape()->getProlifDestinationAnchor()[idx] + _coord;

				_divide_limit_Treg_exp -= 1;
				_divide_cd_Treg_exp = int(params.getVal(PARAM_TREG_EXP_INTERVAL_SLICE) / (ArgI / (ArgI + params.getVal(PARAM_EC50_ARGI_TREG)) * (1 - (tregul/ (tumvol * params.getVal(PARAM_TREGMAX)))) * (1 - ent / (ent + params.getVal(PARAM_IC50_ENT_ARGI)))) + .5);
			}

		}

	}
	return divide;
}

void Treg::move_all_source_sink(void) const
{
	//move_source_sink(_source_IL_10);
}

int Treg::getTregLife(){

	double lifeMean = params.getVal(PARAM_TREG_LIFE_MEAN_SLICE);
	double lifeSd = params.getVal(PARAM_TREG_LIFE_SD_SLICE);
	double tLifeD = lifeMean + rng.get_norm_std() * lifeSd;
	/*double lifeMean = params.getVal(PARAM_TREG_LIFE_MEAN);

	double tLifeD = rng.get_exponential(lifeMean);*/

	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	//std::cout << "random Treg life: " << tLife << std::endl;
	return tLife;
}

};
};