//#include <boost/serialization/export.hpp>
#include "TCell.h"

//BOOST_CLASS_EXPORT_IMPLEMENT(TCell)
#include "../../core/GlobalUtilities.h"
#include "../compartment/Tumor.h"


namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

//#define IFNG_ID 0

using std::string;
using std::stringstream;

static int TCellSize = 1;

TCell::TCell()
	: _count_neighbor_cancer(0)
	, _count_neighbor_Treg(0)
	, _count_neighbor_all(0)
	, _max_neighbor_PDL1(0)
{
}

TCell::TCell(SpatialCompartment* c)
	:Cell_Tumor(c)
	, _divide_flag(false)
	, _divide_cd(0)
	, _divide_limit(params.getVal(PARAM_T_DIV_LIMIT))
	, _IL2_exposure(0)
	, _IL2_release_remain(params.getVal(PARAM_IL_2_RELEASE_TIME))
	, _source_IFNg(NULL)
	, _source_IL_2(NULL)
	, _count_neighbor_cancer(0)
	, _count_neighbor_Treg(0)
	, _count_neighbor_all(0)
	, _max_neighbor_PDL1(0)
{
	_state = AgentStateEnum::T_CELL_EFF;
	_life = TCell::getTCellLife();
	_TCR_clone = 0; //This is a dummy value that will be reset later
	std::cerr << "TCell constructed: " << getID() << ", " << _TCR_clone << std::endl;
	//cout << getID() << ", " << getType() << endl;
}


TCell::TCell( const TCell & c)
	:Cell_Tumor(c)
	, _divide_flag(c._divide_flag)
	, _divide_cd(c._divide_cd)
	, _divide_limit(c._divide_limit)
	, _IL2_exposure(c._IL2_exposure)
	, _IL2_release_remain(c._IL2_release_remain)
	, _TCR_clone(c._TCR_clone)
	, _source_IFNg(NULL)
	, _source_IL_2(NULL)
	, _count_neighbor_cancer(0)
	, _count_neighbor_Treg(0)
	, _count_neighbor_all(0)
	, _max_neighbor_PDL1(0)
{
	_state = c._state;
	//cout << getID() << ", " << getType() << endl;
	_life = TCell::getTCellLife();
	if (_state == AgentStateEnum::T_CELL_CYT)
	{
		setup_chem_source(_source_IFNg, CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
		setup_chem_source(_source_IL_2, CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
	}
}


TCell::~TCell()
{
}
string TCell::toString()const{
	stringstream ss;
	ss << CellAgent::toString();
	ss << "division flag: " << _divide_flag << ", division cool down: " 
		<< _divide_cd  << std::endl;
	return ss.str();
}

bool TCell::agent_movement_step(double t, double dt, Coord& c){
	bool move = false;
	if (rng.get_unif_01() < params.getVal(PARAM_T_CELL_MOVE_PROB))
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

/*! Scan neighborhood and count neighbor cell type
	number of Teff neighbor for Cancer cell; 
	number of Cancer cell neighbor for Teff; 
	number of Treg neighbor for Teff; 
*/
void TCell::agent_state_scan(void){

	const auto shape = getCellShape();
	if (_state == AgentStateEnum::T_CELL_CYT){
		/**/
		// scan cancer cells
		//std::cout << "T cell at: "<< _coord << std::endl;
		//int nr_cancer_neighbor = 0;
		_compartment->for_each_neighbor_ag(shape->getEnvironmentLocations(),
			_coord, [&](BaseAgent* ag){
			_count_neighbor_all++;
			auto pCell = dynamic_cast<Cell_Tumor*>(ag);
			update_neighbor_PDL1(pCell->get_PDL1());
			if (ag->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
				auto pCancer = dynamic_cast<CancerCell*>(ag);
				// =============================================================================================================
				/*std::cout << "rng check: id= " << getID() << " cytotoxic t-cell with tcr: " << this->getTCR() << " at position " << this->getCoord() << std::endl;
				auto temp = pCancer->getAntigenVec();
				std::cout << "Cancer cell at position "<< pCancer->getCoord()<<" with cancer cell antigen vector: " << "[ ";
				for (size_t j = 0; j < temp.size(); j++)
				{
					std::cout << temp.at(j) << " ";
				}
				std::cout << "]" << std::endl;
				std::cout << "\n";*/
				// =============================================================================================================
				if (pCancer->getAntigenVec().at(_TCR_clone))
				{
					// ============================================================================================================
					/*std::cout << "Cytotoxic Tcell recognized Cancer Cell" << std::endl;
					std::cout << "\n";*/
					// ============================================================================================================
					inc_neighbor_cancer();						
					pCancer->inc_neighbor_Teff();
					//nr_cancer_neighbor += 1;
				}
				else {
					/*std::cout << "Cytotoxic Tcell failed to recognized Cancer Cell" << std::endl;
					std::cout << "\n";*/
				}
			}
			else if (ag->getType() == AgentTypeEnum::CELL_TYPE_TREG){
				inc_neighbor_Treg();
			}
			return false;
		});
		//std::cout << "T cell neighbor PDL1: "<< _max_neighbor_PDL1 << std::endl;
		//std::cout << "T cell neighbor Cancer: "<< nr_cancer_neighbor << std::endl;
	}
	return;
}

bool TCell::agent_state_step(double t, double dt, Coord& c){
	bool divide = false;
	if (!isDead())
	{
		_life--;
		if (_life == 0)
		{
			setDead();
			// remove source when cell die
			return divide;
		}
	}

	const auto shape = getCellShape();

	Cell_Tumor::agent_state_step(t, dt, c);

	auto tumor = dynamic_cast<Tumor*>(_compartment);

	double IL2 = get_tumor().get_chem(_coord, CHEM_IL_2);
	_IL2_exposure += params.getVal(PARAM_SEC_PER_TIME_SLICE) * IL2;

	// effector cells to proliferate on IL2 exposure
	if (_IL2_exposure > params.getVal(PARAM_IL_2_PROLIF_TH))
	{
		_divide_flag = true;
	}

	// Look for tumor Ag
	if (_state == AgentStateEnum::T_CELL_EFF)
	{
		bool cancer_found = _compartment->check_neighbor_agents(shape->getEnvironmentLocations(), // cancer_found is defined as the adjacent cancer cell has an matched TCR-antigen pair
			_coord, [&](const BaseAgent* agent) {
				if (agent->getType() == AgentTypeEnum::CELL_TYPE_CANCER && agent->getState() == AgentStateEnum::CANCER_PROGENITOR) {
					auto pCancer = dynamic_cast<const CancerCell*>(agent);
					//======================================================================================================
					/*std::cout << "RNG check: ID=" << getID() << " Effector T-cell with TCR: " << this->getTCR() << std::endl;
					auto temp = pCancer->getAntigenVec();
					std::cout << "Cancer cell Antigen Vector: " << "[ ";
					for (size_t j = 0; j < temp.size(); j++)
					{
						std::cout << temp.at(j) << " ";
					}
					std::cout << "]" << std::endl;*/
					//======================================================================================================
					if (pCancer->getAntigenVec().at(_TCR_clone))
					{
						return true;
					}
					/*std::cout << "Tcell Failed to transforming into T cyc" << std::endl;
					std::cout << "\n";*/
				}
				return false;
			});
		
		if (cancer_found)
		{
			//=========================================================
			//std::cout << "Tcell become cytotoxic from eff" << std::endl;
			//=========================================================
			_state = AgentStateEnum::T_CELL_CYT;
			_divide_flag = true;
			setup_chem_source(_source_IFNg, CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
			setup_chem_source(_source_IL_2, CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
		}
	}

	if (_state == AgentStateEnum::T_CELL_CYT){
		/**/
		double nivo = tumor->get_Nivo();

		// kill one cancer cell
		// now handled from cancer cells

		// IL-2 release time limit
		if (_IL2_release_remain > 0)
		{
			_IL2_release_remain -= params.getVal(PARAM_SEC_PER_TIME_SLICE);
		}
		else{
			// set IL-2 source to 0
			update_chem_source(_source_IL_2, 0.0);
		}
		
		// exhaustion
		if (_count_neighbor_Treg > 0)// suppresion by Treg
		{
			double bond = get_PD1_PDL1(_max_neighbor_PDL1, nivo);
			double supp = get_PD1_supp(bond, params.getVal(PARAM_N_PD1_PDL1));
			double q = double(_count_neighbor_Treg) / (_count_neighbor_all + params.getVal(PARAM_CELL));
			double p_exhaust = get_exhaust_prob_Treg(supp, q);
			if (rng.get_unif_01() < p_exhaust)
			{
				/*
				std::cout << "RNG check: ID=" << getID() <<
				"(Teff exhaust by Treg): " << rng.get_unif_01() << std::endl;
				*/
				set_suppressed();
			}
		}
		else if (_count_neighbor_all > 0)// suppresion by PDL1 
		{
			double bond = get_PD1_PDL1(_max_neighbor_PDL1, nivo);
			double supp = get_PD1_supp(bond, 1);
			double p_exhaust = get_exhaust_prob_PDL1(supp, 1.0);
			/*
			std::cout << "suppression prob: " << p_exhaust << std::endl;
			*/
			if (rng.get_unif_01() < p_exhaust)
			{
				/*
				std::cout << "RNG check: ID=" << getID() <<
				" (Teff exhaust by PDL1): " << rng.get_unif_01() << std::endl;
				*/
				set_suppressed();
			}
		}
	}

	if (_divide_cd > 0)
	{
		_divide_cd--;
	}
	
	if (_divide_limit > 0 && _divide_flag && _divide_cd == 0)
	{
		int idx;
		if (_compartment->getOneOpenVoxel(shape->getProlifDestinationVoxels(), 
			shape->getProlifDestinationAnchor(), _coord, getType(), idx, rng))
		{
			divide = true;
			//cout << "idx: " << idx << ", " << getCellShape()->getProlif()[idx] << endl;
			c = getCellShape()->getProlifDestinationAnchor()[idx] + _coord;

			_divide_flag = true;
			_divide_limit -= 1;
			_divide_cd = params.getVal(PARAM_T_DIV_INTERVAL);
		}

	}

	_count_neighbor_cancer = 0;
	_count_neighbor_Treg = 0;
	_count_neighbor_all = 0;
	_max_neighbor_PDL1 = 0;
	return divide;
}

/*! set to suppressed state
*/
void TCell::set_suppressed(void){
	_divide_flag = false;
	_divide_limit = 0;
	_state = AgentStateEnum::T_CELL_SUPP;
	remove_source_sink(_source_IFNg);
	remove_source_sink(_source_IL_2);
	return;
}

//! move sources (IFN and IL2)
void TCell::move_all_source_sink(void)const{
	//std::cout << "moving sources: " << _coord << std::endl;
	move_source_sink(_source_IFNg);
	move_source_sink(_source_IL_2);
	return;
}

//! remove sources (IFN and IL2)
void TCell::remove_all_source_sink(void){
	remove_source_sink(_source_IFNg);
	remove_source_sink(_source_IL_2);
	return;
}

void TCell::update_neighbor_PDL1(double PDL1){
	if (PDL1 > _max_neighbor_PDL1)
	{
		_max_neighbor_PDL1 = PDL1;
	}
	return;
}

//! PD1_PDL1 bond in synapse
double TCell::get_PD1_PDL1(double PDL1, double Nivo){

	static double T1 = params.getVal(PARAM_PD1_SYN);
	double T2 = PDL1;
	static double k1 = params.getVal(PARAM_PDL1_K1);
	static double k2 = params.getVal(PARAM_PDL1_K2);
	static double k3 = params.getVal(PARAM_PDL1_K3);
	/*
	double c = T1 + PDL1 + (1 + k2 * Nivo) / k1;
	double discriminant = c * c - 4 * T1 * PDL1;
	double bond = .5*(c - std::sqrt(discriminant));
	*/

	double a = 1;
	double b = (Nivo*k2/k1*(2*k3/k1-1) - 2*T2 - T1 - 1/k1)/T2;
	double c = (Nivo*k2/k1 + 1/k1  + T2 + 2*T1 )/T2;
	double d = -T1/T2;

	/*
	std::cout << k1 << "," << k2 << "," << k3 << ","
		<< T1 << "," << T2 << std::endl;
	std::cout << a << "," << b << "," << c << "," << d << std::endl;
	*/

	//Newton_Raphson_root
	int max_iter = 20;
	double tol_rel = 1E-5;
	double root = 0;
	double res, root_new, f, f1;
	int i = 0;
	while (i < max_iter){
		f = a*std::pow(root, 3) + b*std::pow(root, 2)+ c*root + d;
		f1 = 3.0*a*std::pow(root, 2) + 2.0*b*root + c;
		root_new = root - f/f1;
		res = std::abs(root_new - root) / root_new;
		if (res > tol_rel){
			i++;
			root = root_new;
		}
		else{
			break;
		}
	}
	
	double bond = T2*root;

	/*
	std::cout << "PDL1: " << PDL1 << ", PD1_syn: "<< T1<< ", nivo: " << Nivo
		<< ", root: " << root << ", bond: "<< bond << std::endl;
	*/		

	return bond;
}


//! get suppression from PD1_PDL1 bond
double TCell::get_PD1_supp(double bond, double n) {
	double k50 = params.getVal(PARAM_PD1_PDL1_HALF);
	return get_Hill_equation(bond, k50, n);
}



int TCell::getTCellLife(){
	double lifeMean = params.getVal(PARAM_T_CELL_LIFE_MEAN_SLICE);
	double lifeSd = params.getVal(PARAM_T_CELL_LIFE_SD_SLICE);
	double tLifeD = lifeMean + rng.get_norm_std() * lifeSd;
	/*
	double tLifeD = rng.get_exponential(lifeMean);
	*/

	int tLife = int(tLifeD + 0.5);
	tLife = tLife > 0 ? tLife : 0;
	//std::cout << "random T cell life: " << tLife << std::endl;
	return tLife;
}

/*! probability of Cancer cell getting killed by Teff
	\param[in] supp: suppression
	p_kill = 1 - B_esc^((1-supp)*q),
	where B_esc = exp(-tstep * k_C_death_by_T)
*/
double TCell::get_kill_prob(double supp, double q){
	return 1 - std::pow(params.getVal(PARAM_ESCAPE_BASE), q*(1-supp));
}

/*! probability of turning exhausted from PDL1-PD1 interaction
	p_exhaust_PDL1 = 1 - B^(supp*q)
*/
double TCell::get_exhaust_prob_PDL1(double supp, double q){
	return 1 - std::pow(params.getVal(PARAM_EXHUAST_BASE_PDL1), q*supp);
}

/*! probability of turning exhausted from Treg
	p_exhaust_Treg = 1 - B^((1 + supp)*q)
*/
double TCell::get_exhaust_prob_Treg(double supp, double q){
	return 1 - std::pow(params.getVal(PARAM_EXHUAST_BASE_TREG), q);
}

void TCell::set_TCR_clone(const std::vector<double> Tcell_Distribution, const double total_concentration) {

	long double avogadro_number = params.getVal(PARAM_AVOGADROS);
	std::vector<double> teff_dist = Tcell_Distribution;
	for (double teff_prob : teff_dist)
	{
		teff_prob * 1000 / total_concentration; // Scale the T Cell distribution so that the sum is 1000;
	}
	std::random_device rand_dev;
	std::default_random_engine generator(rand_dev());
	std::discrete_distribution<int> distribution(teff_dist.begin(), teff_dist.end());
	_TCR_clone = distribution(generator);
	//return true;
}
	/*if (TCR_recruit_left > 0)
	{
		return true;
	}
	else
	{
		return false;
	}

	//std::cerr << "TCell TCR reset: " << getID() << ", " << _TCR_clone << std::endl;
}
/*
bool TCell::set_TCR_clone(const std::vector<double> Tcell_Distribution, const double total_concentration, const double abm_scaler, const std::vector<double> var_abm_to_qsp_teff, double num_rec) {

	long double avogadro_number = params.getVal(PARAM_AVOGADROS);
	std::vector<double> teff_dist = Tcell_Distribution;
	for (double teff_prob : teff_dist)
	{
		teff_prob * 1000 / total_concentration; // Scale the T Cell distribution so that the sum is 1000;
	}
	std::random_device rand_dev;
	std::default_random_engine generator(rand_dev());
	std::discrete_distribution<int> distribution(teff_dist.begin(), teff_dist.end());
	_TCR_clone = distribution(generator);
	
	//! Number of tcell with this particular TCR being recruited
	double TCR_recruited = var_abm_to_qsp_teff[_TCR_clone];
	
	//! MAX Number of tcell with this particular TCR can be recruited
	double TCR_recruit_left = Tcell_Distribution[_TCR_clone] - (TCR_recruited / avogadro_number * abm_scaler) * num_rec;
	std::cout << "TCR_recruited: " << TCR_recruited << ", abmscaler / Avogadro: " << abm_scaler / avogadro_number << "ABM Scalor : " << abm_scaler <<std::endl;

	return true;
	/*if (TCR_recruit_left > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
	
	//std::cerr << "TCell TCR reset: " << getID() << ", " << _TCR_clone << std::endl;
}
*/
};
};
