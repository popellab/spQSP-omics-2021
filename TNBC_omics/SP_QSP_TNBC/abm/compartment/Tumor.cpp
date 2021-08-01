//#include <boost/serialization/export.hpp>
//#include <boost/archive/text_iarchive.hpp>
//#include <boost/archive/text_oarchive.hpp>
#include "Tumor.h"
#include "LymphCentral.h"
//BOOST_CLASS_EXPORT_IMPLEMENT(Tumor);

//#include "../agent/CancerCell.h"`
#include <algorithm>
#include <exception>
#include <numeric>
#include <queue>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

#include "../../core/GlobalUtilities.h"

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

using std::string;
using std::cout;
using std::cerr;
using std::endl;

//#include "../../MemoryUsage.h"



Tumor::Tumor(int x, int y, int z)
	: SpatialCompartment(x, y, z)
	, _stats()
	, _voxel_ic()
	, _chem()
	, _t_source()
	, _mdsc_source()	
	, _tInitDummy(new TCell(this))
	, _cInitDummy(new CancerCell(this))
	, _macInitDummy(new Mac(this))
	, _fibInitDummy(new Fib(this))
	, _TregInitDummy(new Treg(this))
	, _MDSCInitDummy(new MDSC(this))
	//, _voxel_size(params.getVal(PARAM_VOXEL_SIZE))
	, _allow_shift_grid(false)
	, _center_target()
	, _agGrid_temp(AgentGrid(x,y,z, NULL))
	, _var_abm_to_qsp()
	, _var_abm_to_qsp_teff()
	, _var_abm_to_qsp_antigen()
	, _concentration_cc(0)
	, _concentration_t_cyt(0)
	, _concentration_t_cyt_tum(0)
	, _concentration_t_reg(0)
	, _concentration_t_reg_tum(0)
	, _concentration_mdsc(0)
	, _concentration_nivo(0)
	, _concentration_ent(0)
	, _concentration_cx(0)
	, _concentration_t_exh(0)
	, _concentration_ccl2(0)
	, _tumor_volume(0)
	, _distribution_t_cyt()
	, _distribution_t_cyt_tum()
{
	//CC total, CC death total, CC death Teff, Teff recruit, Treg recruit
	_var_abm_to_qsp = std::vector<double>(TUMEX_VAR_NUM, 0);
	//update number of antigens from dead cancer cell
	_var_abm_to_qsp_antigen = std::vector<double>(17, 0); // 17 is the max number of antigen
	//update number of teff being recruited
	_var_abm_to_qsp_teff = std::vector<double>(17, 0); // 17 is the max number of TCR clones
	// dummy cells do not act as source
	//_tInitDummy->get_source_IFNg()->set_remove();
	initAgentGrid();
}

/*return the first line of the antigen data, and every row have the same 
	length, inaddition, return size() - 1 because the first element is
	cell name, so substract*/
int Tumor::getNumAntigens()const
{
	//std::cout << "The size of antigen vector is " << antigenData[0].size() - 1 << std::endl;
	return antigenData[0].size() - 1;
}


/*!
  _agGrid uses pointer to AgentGridVoxel (or its derived classes)
  as its element. This has a higher computational cost than directly
  hosting the instantications of the derived class (base won't work).
  The benefit is that this way the voxel objects can be accesses from
  the SpatialCompartment abstract class.
*/
bool Tumor::initAgentGrid() {
	//std::cout << "init agGrid:" << std::endl;
	for (int i = 0; i < _sizeX; i++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int k = 0; k < _sizeZ; k++)
			{
				_agGrid(i, j, k) = new TumorGridVoxel(i, j, k);
				//std::cout << _agGrid(i, j, k) << std::endl;
			}
		}
	}

	return true;
}

Tumor::~Tumor()
{
	//std::cout << "destructor:" << std::endl;
	for (int i = 0; i < _sizeX; i++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int k = 0; k < _sizeZ; k++)
			{
				//std::cout << _agGrid(i, j, k) << std::endl;
				delete _agGrid(i, j, k);
			}
		}
	}

	delete _tInitDummy;
	delete _cInitDummy;
	delete _macInitDummy;
	delete _fibInitDummy;
	delete _TregInitDummy;
	delete _MDSCInitDummy;
}

// add one source of entry
void Tumor::add_lymphocyte_source(const Coord& c){
	_t_source.push_back(c);
	return;
}

// add one source of entry
void Tumor::add_mdsc_source(const Coord& c){
	_mdsc_source.push_back(c);
	return;
}

static bool sortNonTCell(CellAgent * ptrCell) { return ptrCell->getType() != AgentTypeEnum::CELL_TYPE_T; }

/*! simulate generaic tumor compartment for one time slice.
	Steps:
	-# process ABM step (cellular scale events)
		-# recruit cells
		-# for each cell,
			-# progress one step and determine consequences.
			-# process consequence of the step.
		-# go through cell list, remove dead cells, add live cells to ABM snapshot stats
		-# refresh MHC
	-# process molecular scale substeps
	-# shuffle cell list.
*/
void Tumor::timeSlice(unsigned long slice) {

	static double dt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
	const double t = slice * dt;
	//static unsigned int nrMolSlice = params.getVal(PARAM_MOL_STEP_PER_SLICE);
	//static double dtMol = dt / nrMolSlice;


	// reset time step stats
	_stats.resetStats();

	// Reset to zeros
	std::fill(_var_abm_to_qsp.begin(), _var_abm_to_qsp.end(), 0);
	std::fill(_var_abm_to_qsp_teff.begin(), _var_abm_to_qsp_teff.end(), 0);
	std::fill(_var_abm_to_qsp_antigen.begin(), _var_abm_to_qsp_antigen.end(), 0);
	//------------------------   Cellular scale events   -----------------------//
	//BOOST_TEST_MESSAGE("slice:" << slice);
	//cout << "slice:" << slice << endl;

	// Recruitment
	// model comparison: find invasive front to recruit T cells and MDSC
	time_slice_recruitment();
	//cout << "after recruitment:" << endl;

	//std::cout << "RNG check (" << slice << ") rec: " << rng.get_unif_01() << std::endl;

	time_slice_movement(t, dt);
	//cout << "after movement:" << endl;

	//std::cout << "RNG check (" << slice << ") move: " << rng.get_unif_01() << std::endl;
	time_slice_state_change(t, dt);
	//cout << "after state:" << endl;
	//std::cout << "RNG check (" << slice << ") state: " << rng.get_unif_01() << std::endl;

	
	if (_allow_shift_grid && slice % params.getVal(PARAM_SHIFTGRID_INTERVAL) == 0)
	{
		shift_adjust_center();
	}
	
	unsigned int nr_cancer = 0;
	Coord c(0, 0, 0);

	for (auto p : _vecAgent) {
		if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER) {
			c = c + p->getCoord();
			nr_cancer += 1;
		}
	}

	//cout << "after shift:" << endl;

	//std::cout << "RNG check (" << slice << ") shift: " << rng.get_unif_01() << std::endl;
	time_slice_final_scan();
	//cout << "after clear dead:" << endl;

	// ------------------------  molecular scale events     ---------------------//
	/*
	for (auto p : _vecAgent)
	{
		if (p->getState() == T_CELL_CYT)
		{
			auto pt = dynamic_cast<TCell*>(p);
			if (pt->get_source_IFNg()->get_index() != _chem.get_voxel_idx(pt->getCoord()))
				if (!pt->get_source_IFNg()->is_to_remove())
					std::cout << "location mismatch (IFNg):" << pt->getCoord() << ", "
						<< _chem.idx_to_coords(pt->get_source_IFNg()->get_index())
						<< std::endl;
		}
	}*/
	time_slice_molecular(t);
	//std::cout << "RNG check (" << slice << ") mol: " << rng.get_unif_01() << std::endl;

	//cout << "after molecular:" << endl;
	// shuffle cell vector
	if (slice % params.getVal(PARAM_SHUFFLE_CELL_VEC_INTERVAL) == 0)
	{
		std::shuffle(_vecAgent.begin(), _vecAgent.end(), rng.getRNG());
	}

	/*
	std::cout << "Num sources: " 
		<< _chem.get_num_source_sink() << std::endl;
	std::cout << "CELL INFO" << std::endl;
	}*/
}

/*! Adjust the center of the tumor compartment
*/
void Tumor::shift_adjust_center(void){
	Coord3D c = get_cam_shift();
	// only shift when distance is larger than c_tol
	int c_tol = 0;
	//std::cout << c << "," << c.length() << std::endl;
	//if (c.length()>c_tol)
	if (c.z > 0)
	{
		//shift_grid(c);
		//std::cout << "Shifting: " << c << "," << c.length() << std::endl;
		Coord3D c_zshift(0, 0, c.z);
		shift_grid(c_zshift);
		std::cout << "Shifting: " << c_zshift << std::endl;
		std::cout << "Tumor CoM" << _center_target << std::endl;
	}
	return;
}

/*!
*/
Coord3D Tumor::get_cam_shift(void) {

//! get center of mass
	unsigned int nr_cancer = 0;
	Coord c(0, 0, 0);

	for (auto p : _vecAgent){
		if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
			c = c + p->getCoord();
			nr_cancer += 1;
		}
	}
	if (nr_cancer) {
		std::cout << "target: " << _center_target << "CoM: " << c / nr_cancer << "," << c.length() << std::endl;
		return c / nr_cancer - _center_target;
	}
	else {
		std::cout << "target: " << _center_target << "CoM: No cancer at all" << std::endl;
		return c;// (0,0,0)
	}
}
/*! Shift content of grid by crd
	In current version, grid size change is not allowed.
	When we need to shift center of grid by (x, y, z), 
	we move the content of the grid instead.
	The following contents need to be relocated:
	# Diffusion grid
		# concentration
		# source/sink
	# Agent grid
		# cells
		# structures
			# recruitment sources
	\param [in] crd: shifting distance (x, y, z) (of "camera")
*/
void Tumor::shift_grid(Coord3D& c_shift) {
	// k_step = 1 if k >0 else -1 (direction)
	bool x_pos = c_shift.x >= 0;
	bool y_pos = c_shift.y >= 0;
	bool z_pos = c_shift.z >= 0;

				
	int nr_shift = _sizeX*_sizeY*_sizeZ
		- (_sizeX - abs(c_shift.x))*(_sizeY - abs(c_shift.y))*(_sizeZ - abs(c_shift.z));
	int i_shift = 0;
	//cout << "shift: " << c_shift << ", nr: " << nr_shift << endl;
	CoordVec drop_out(nr_shift, Coord(0, 0, 0));
	/* shift agent grid:
		move camera by (x, y, z): 
		agents: (x0, y0, z0) -> (x0-x, y0-y, z-z0)
		*/
	// first round: save voxels to be moved out of sight to temp grid
	//	and shift remaining voxels to new location 

	for_each_grid_coord(x_pos, y_pos, z_pos, [&](Coord3D& c){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		//BOOST_TEST_MESSAGE(c);
		if (!_agGrid.inGrid(c_new)){
			//add to temp grid
			_agGrid_temp(c_target) = _agGrid(c);
			drop_out[i_shift] = c;
			i_shift++;
		}
		else{
			_agGrid(c_target) = _agGrid(c);
			for (auto pAg : _agGrid(c_target)->get_agents()){
				//BOOST_TEST_MESSAGE(c_target);
				auto pCell = dynamic_cast<Cell_Tumor*>(pAg);
				pCell->setCoord(c_target);
			}
		}
		return;
	});

	//cout << "added to drop_out: " << i_shift << endl;

	// populate new voxels (recycling out-of-site ones)
	/*
	for_each_voxel(x_pos, y_pos, z_pos, [&](Coord3D& c){
		return;
	});*/

	for (auto& c : drop_out){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		if (!_agGrid.inGrid(c_new)){
			// overwrite info and put back to grid
			_agGrid(c_target) = _agGrid_temp(c_target);
			// apply changes to voxels entering the grid (recycled)
			// remove old, 
			for (auto pAg : _agGrid(c_target)->get_agents()){
				auto pCell = dynamic_cast<Cell_Tumor*>(pAg);
				pCell->set_drop_out();
			}
			_agGrid(c_target)->remove_all_agents();
			//populate new
			populate_voxel_random(c_target);
		}
	}

	// remove dropped cells
	// now done together with dead cell

	// shift recruitment sources? 

	/* shift diffusion grid
		# shift values
		# shift point source/sink
		*/
	unsigned int num_sub = _chem.get_num_substrates();
	BioFVMGrid::chem_grid _chem_temp = _chem.get_concentrations();
	for_each_grid_coord(x_pos, y_pos, z_pos, [&](Coord3D& c){
		auto c_new = c - c_shift;
		auto c_target = _agGrid.get_toroidal(c_new);
		unsigned int idx = _chem.get_voxel_idx(c_target);
		for (size_t i = 0; i < num_sub; i++)
		{
			if (_agGrid.inGrid(c_new)){
				_chem_temp[idx][i] = _chem(c, i);
			}
			else{
				_chem_temp[idx][i] = 0;
			}
		}
		return;
	});
	_chem.setup_concentrations(_chem_temp);

	// shift source/sink: already taken care of.
	// Sources moving out of grid: set removed with agents
	// Other sources: moved together with agents.
	/* old code
	for (auto pS : _chem.get_sink_source()){
		auto c = _chem.idx_to_coords(pS->get_index());
		auto c_target = c - c_shift;
		_chem.move_source_sink(pS, c_target);
	}*/

	return;
}

void Tumor::time_slice_recruitment(){
	// CD8
	const auto dummy = _tInitDummy; //When the tumor compartment is created, Only one T-cell is created, and the T-cell serves as the prototype for the rest of T-cell
	double pRec = get_T_recruitment_prob(_concentration_t_cyt, params.getVal(PARAM_TEFF_RECRUIT_K));
	std::cout << "Central T-cell Total: " << _concentration_t_cyt << std::endl;
	std::cout << "T-cell recruit Prob: " << pRec << std::endl;
	for (auto& crd : _t_source){
		if (rng.get_unif_01() < pRec) {
			bool rec = recruitOneCellInMooreNeighborhood(dummy, crd, rng);
			if (rec)
			{
				auto pT = dynamic_cast<TCell*>(_vecAgent.back());
				pT->set_TCR_clone(_distribution_t_cyt, _concentration_t_cyt);  //for new recruited T-cell a new TCR clone is assigned pt -> get_TCR()
				_stats.incRecruit(dummy->getType(), dummy->getState());
				inc_abm_var_exchange_teff(pT->getTCR()); //collect the number of T-cell being recruited for each TCR
				/*
				bool isSufficient = pT->set_TCR_clone(_distribution_t_cyt, _concentration_t_cyt);  //for new recruited T-cell a new TCR clone is assigned pt -> get_TCR()
				if (isSufficient)
				{
					_stats.incRecruit(dummy->getType(), dummy->getState());
					inc_abm_var_exchange_teff(pT->getTCR()); //collect the number of T-cell being recruited for each TCR
				}
				else
				{
					_vecAgent.pop_back();
					removeAgentFromGrid(pT->getCoord(), pT);
					delete pT;
				}
				*/
			}
		}
	}
	// Treg
	/**/
	const auto Treg_dummy = _TregInitDummy;
	double pRec_Treg = get_T_recruitment_prob(_concentration_t_reg , params.getVal(PARAM_TREG_RECRUIT_K));
	std::cout << "Treg in central: " << _concentration_t_reg << std::endl;
	std::cout << "rec prob (Treg): " << pRec_Treg << std::endl;
	for (auto& crd : _t_source){
		if (rng.get_unif_01() < pRec_Treg) {
			bool rec = recruitOneCellInMooreNeighborhood(Treg_dummy, crd, rng);
			if (rec)
			{
				_stats.incRecruit(Treg_dummy->getType(), Treg_dummy->getState());
				auto pT = dynamic_cast<Treg*>(_vecAgent.back());
				inc_abm_var_exchange(TUMEX_TREG_REC);
			}
		}
	}

	// MDSC
	/**/
	const auto MDSC_dummy = _MDSCInitDummy;
	//std::cout << "rec prob (MDSC): " << pRec_MDSC << std::endl;
		double pBaseRec_MDSC = params.getVal(PARAM_MDSC_RECRUIT_K) / (params.getVal(PARAM_MDSC_RECRUIT_K) + params.getVal(PARAM_MDSC_RECRUIT_BY_CCL2_K) * _concentration_ccl2/(_concentration_ccl2 + params.getVal(PARAM_EC50_CCL2_REC)));
		double pRec_MDSC;	
		if (rng.get_unif_01() < pBaseRec_MDSC) {
			pRec_MDSC = get_MDSC_recruitment_prob(params.getVal(PARAM_MDSCMAX) * _tumor_volume - _concentration_mdsc, params.getVal(PARAM_MDSC_RECRUIT_K));		
			std::cout << "rec prob (MDSCbase): " << pRec_MDSC << std::endl;
		}
		else{
			pRec_MDSC = get_MDSC_recruitment_prob(params.getVal(PARAM_MDSCMAX) * _tumor_volume - _concentration_mdsc, params.getVal(PARAM_MDSC_RECRUIT_BY_CCL2_K) * _concentration_ccl2/(_concentration_ccl2 + params.getVal(PARAM_EC50_CCL2_REC)));	
			std::cout << "rec prob (MDSC): " << pRec_MDSC << std::endl;
		}	
	for (auto& crd : _mdsc_source){	
		if (rng.get_unif_01() < pRec_MDSC) {
			bool rec = recruitOneCellInMooreNeighborhood(MDSC_dummy, crd, rng);		
			if (rec)
			{
				_stats.incRecruit(MDSC_dummy->getType(), MDSC_dummy->getState());
				auto pT = dynamic_cast<MDSC*>(_vecAgent.back());
			}
		}
	}	
	
	return;
}

void Tumor::time_slice_movement(double t, double dt){
	//BOOST_TEST_MESSAGE("total cells: " << _vecAgent.size());
	for (CellVec::size_type i = 0; i < _vecAgent.size(); i++)
	{
		auto ptP = dynamic_cast<Cell_Tumor*>(_vecAgent[i]);
		Coord c(0, 0, 0);
		//cout << "move: "  <<  i << ", type: " << ptP->getType() << endl;
		//cout.flush();

		//BOOST_TEST_MESSAGE(i);
		bool move = ptP->agent_movement_step(t, dt, c);
		if (move)
		{
			// need to do this first or problems will occur when relocating other cells to current location
			try{
				removeAgentFromGrid(ptP->getCoord(), ptP);
			}
			catch (std::exception & e){
				cout << ptP->toString() << endl;
				cout << ptP->getCoord() << "->" << c << endl;
				throw(e);
			}

			//move current cell
			ptP->setCoord(c);
			addAgentToGrid(c, ptP);

			//printCellInfo(slice, ptP, "after moving cell");
			//ptP->move_all_source_sink();

			_stats.incMove(ptP->getType(), ptP->getState());
		}
	}
	return;

}
void Tumor::time_slice_state_change(double t, double dt){

	/*Scan neighbor*/
	for (auto const p: _vecAgent)
	{
		auto pCell = dynamic_cast<Cell_Tumor*>(p);
		pCell->agent_state_scan();
	}
	
	//std::cout << "RNG check "<< " scan: " << rng.get_unif_01() << std::endl;
	/*process state change*/
	CellVec::size_type originalCount = _vecAgent.size();
	for (CellVec::size_type i = 0; i < originalCount; i++)
	{
		//std::cout << "Size of _vecAgent is: " << _vecAgent.size() << std::endl;
		CellAgent *ptP = _vecAgent[i];
		Coord c(0, 0, 0);

		//std::cout << "processing: " << ptP->getID() << ", " << ptP->getType() << ", " << ptP->getState() << std::endl;
		bool divide = ptP->agent_state_step(t, dt, c);
		//std::cout << "divide: " << divide << std::endl;
		if (divide)
		{
			//std::cout << "adding daughter cell" << std::endl;
			auto ptD = ptP->createCellCopy();
			if (ptD->getType() == AgentTypeEnum::CELL_TYPE_T) // Check if the daughter cell have the same TCR as the parent cell
			{
				auto ptPTcell = dynamic_cast<TCell*>(ptP);
				auto ptDTcell = dynamic_cast<TCell*>(ptD);
				//std::cout <<"Parent Cell "<< ptPTcell->getID() << " with TCR " <<  ptPTcell->getTCR()<< " divide to daughter " << ptDTcell->getID() << " with TCR Clone "<< ptDTcell->getTCR() << std::endl;
				//std::cout << "\n";
			}
			if (ptP->getType() == AgentTypeEnum::CELL_TYPE_CANCER)
			{
				//std::cout << "cancer cell" << std::endl;

				auto ptPCancer = dynamic_cast<CancerCell*>(ptP);
				auto ptDCancer = dynamic_cast<CancerCell*>(ptD);
				ptDCancer->_stemID = ptPCancer->_stemID;
				if (ptP->getState() == AgentStateEnum::CANCER_STEM)
				{
					bool asymmetric = rng.get_unif_01() < params.getVal(PARAM_CANCER_STEM_ASYMMETRIC_DIV_PROB);
					if (asymmetric)
					{
						// asymmetric division
						ptDCancer->setProgenitor();
					}
					else{
						ptDCancer->_stemID = ptPCancer->getID();
					}
				}
			}
			
			//add daugther cell
			ptD->setCoord(c);
			addAgentToGrid(c, ptD);
			_vecAgent.push_back(ptD);

			//auto pD_tumor = dynamic_cast<Cell_Tumor*>(ptD);
			//pD_tumor->move_all_source_sink();

			_stats.incProlif(ptD->getType(), ptD->getState());
			//std::cout << "added" << std::endl;
		}
	}
	return;
}
/*! Final scan during a slice
	Process cell counts and slice statistics 
	Remove dead cells from cell vector
*/
void Tumor::time_slice_final_scan(void){

	CellVec::iterator last = stable_partition(_vecAgent.begin(), _vecAgent.end(), [&](CellAgent* p){
		auto pC = dynamic_cast<Cell_Tumor*>(p);
		return !(pC->isDead() || pC->is_drop_out());
	});// sortLiveCellAddToStats);

	// count live cells (and cancer cells)
	int nrPDL1 = 0;
	for (auto it = _vecAgent.begin(); it < last; it++)
	{
		Cell_Tumor* ptrCell = dynamic_cast<Cell_Tumor*>(*it);
		//statistics
		_stats.incCellState(ptrCell->getType(), ptrCell->getState());
		if (ptrCell->getType() == CELL_TYPE_CANCER)
		{
			inc_abm_var_exchange(TUMEX_CC);
		}
		if (ptrCell->is_PDL1_pos())
		{
			nrPDL1++;
		}
	}
	_stats.set_stats_misc(STATS_MISC_PDL1_POS, nrPDL1);
	int nr_death = 0;
	int nr_dropout = 0;
	// delete removed cells
	for (auto it = _vecAgent.begin(); it < last; it++)
	{
		auto p = dynamic_cast<Cell_Tumor*>(*it);
		if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER)
		{
			auto pCancer = dynamic_cast<CancerCell*>(p);
			//std::cout << "Cancer Cell: " << pCancer->getID() << " , state : "<< pCancer->getType() <<", Division Remain: " << pCancer->getDivCounter() << std::endl;
		}
	}
	for (CellVec::iterator it = last; it < _vecAgent.end(); it++)
	{
		auto p = dynamic_cast<Cell_Tumor*>(*it);
		if (p->isDead())
		{
			// Cancer cell dies -> antigen presentation
			if (p->getType() == AgentTypeEnum::CELL_TYPE_CANCER){
				inc_abm_var_exchange(TUMEX_CC_DEATH);
				auto pCancer = dynamic_cast<CancerCell*> (p);
				inc_abm_var_exchange_antigen(pCancer->getAntigenVec());
			}
			_stats.incDeath(p->getType(), p->getState());
			if (p->getType() == AgentTypeEnum::CELL_TYPE_TREG)
			{
				//std::cout << "Treg : " << p->getID() << " Dead" << std::endl;
			}
			if (!p->is_drop_out())
			{
				// drop out cells are already removed from voxel
				removeAgentFromGrid(p->getCoord(), p);
			}
			nr_death++;
		}
		else if (p->is_drop_out() && !p->isDead()){
			// dead and dropout are only counted as dead
			_stats.incDropOut(p->getType(), p->getState());
			nr_dropout++;
		}
		else{
			std::cout << "Error: cells neither dead nor drop out are removed" << std::endl;
		}
		// Remove from grid when step is finished so they can still function during the state_change step
		delete p;
	}
	_vecAgent.erase(last, _vecAgent.end());
	//std::cout << "Death: " << nr_death << ", dropout: " << nr_dropout << std::endl;
	return;
}

void Tumor::time_slice_molecular(double t){

	static double dt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
	_chem.remove_dead_source_sink();
	if (params.getVal(PARAM_DIFFUSION_ON) && !_chem.grid_skippable())
	{
		for (auto const pS : _chem.get_sink_source()){
			auto c = _chem.idx_to_coords(pS->get_index());
		}
		static double dt_mol = dt / params.getVal(PARAM_MOLECULAR_STEP_PER_SLICE);
		double t_mol = t;
		double t_mol_end = t + dt;
		while (t_mol < t_mol_end)
		{
			// diffusion time step, including decay/point source/sink
			_chem.timestep(dt_mol);
			t_mol += dt_mol;
		}
	}
}

//! default header for extra remark column when writing cell grid to file
std::string Tumor::getExtraRemarkHeader() const {
	std::stringstream ss;
	ss << "extra";
	return ss.str();
}

/*! Print grid snapshot to a snapshot file.
	Info includes: Element type, IL2 concentration, and MHC concentration.
*/
std::string Tumor::printGridToFile() const {
//void Tumor::printGridToFile(unsigned long slice) const {
	std::stringstream ss;
	// header
	ss << _chem.get_substrate_names() << std::endl;
	// content
	ss << _chem;
	return ss.str();
}


/*! Print state of ODE system of all T cells to ODE stats file,
	includeing a time stamp and cell ID.
	\param [in] slice: time slice
*/
void Tumor::printCellOdeToFile(unsigned long slice) const {
}

/*! print grid snapshot to screen
	\param [in] slice: time slice
*/
void Tumor::printGridToScreen(unsigned long slice)const {

	// header
	cout << "Time slice: "<< slice << endl;
	cout << "Grid size: " << getGridSize() << ", "<< _sizeX << " by "
		<< _sizeY << " by " << _sizeY <<  endl;
	//cout << getGridContent() << endl;

	cout << "x, y, z, nr_ag," << endl;
	for (int k = 0; k < _sizeZ; k++)
	{
		for (int j = 0; j < _sizeY; j++)
		{
			for (int i = 0; i < _sizeX; i++)
			{
				Coord c(i, j, k);
				cout << c.x << "," << c.y << "," << c.z << ",";
				AgentGridVoxel* voxel = _agGrid.get(i, j, k);
				TumorGridVoxel * v2 = dynamic_cast<TumorGridVoxel*>(voxel);
				cout << voxel->getNumAgents();
				//cout << v2->getDistToOrigin();
				//cout << v2->getLoc();
				cout << endl;
			}
		}
	}
}

double Tumor::get_chem(const Coord3D&c, chem_ID i)const{
	return _chem(c, i);
}
//-----------------------------------------  Protected --------------------------------------//

//-----------------------------------------  Private  --------------------------------------//

/*! Setup compartment environment. This includes: sources of T cells and MDSCs.
	create vasculature by mapping graph file to the grid.
	all locations mapped onto the graph are designated as T cell sources.
*/
void Tumor::initEnvironment() {

	int nrSource = 0;

}

template<typename numeric> // This is a function that determine if coming string can be converted to a numeric type
bool is_number(const std::string& s)
{
	numeric n;
	return((std::istringstream(s) >> n >> std::ws).eof());
}

void setClusterCancerAntigen(std::vector<std::vector<std::string>> datalist, // Use this function to fill the antigen vector for each cancer cell through a setter function.
	 std::vector<CancerCell*>& cancerVec)
{
	/*if (datalist.size() != cancerVec.size())
	{
		throw antigenCancerMismatch();
	}*/
	for (size_t i = 0; i < cancerVec.size(); i++)
	{
		std::vector<bool> antigen_vector;
		for (std::string data : datalist.at(i))
		{
			if (is_number<double>(data))
			{
				std::istringstream iss(data);
				double expression{};
				iss >> expression;
				if (expression > 0.0)
				{
					antigen_vector.push_back(true);
				}
				else {
					antigen_vector.push_back(false);
				}
			}
			//std::cout << data << " ";
		}
		//std::cout << "size: " << antigen_vector.size() << std::endl;
		if (cancerVec.size() >= 0)
		{
			cancerVec.at(i)->setAntigenVec(antigen_vector); 
		}
		//delete antigen_vector;
	}
}

//Load the antigen data from csv file
void Tumor::loadAntigenData(std::string antigenDataFile) {
	csvreader reader(antigenDataFile); 
	antigenData = reader.getdata();
	std::cout << "Antigen Data Loaded" << std::endl;
}
/*! Setup initial cells. Called when simulation starts.
	Create initial cells and set their coordinates;
	put initial cells on the grid;
	put initial cells into cell vector.
	After all intial cells are generated, go through the vector and record initial stats
	\param [in] initialCellFileName: initial arrangement of cells

	read initial cell configuration from a separate file
	structure of initial condition file:
	<initialCondition>
	  <cell>
		<x></x>
		<y></y>
		<z></z>
		<celltype></celltype>
		<cellstate></cellstate>
	  </cell>
	  <cell>
		...
	  </cell>
	</initialCondition>

	<cellCluster/> elements will have an additional subelement <count/>,
	indicating the number of cells in that cluster.
*/
void Tumor::initCell(string initialCellFileName){

	namespace pt = boost::property_tree;
	const std::string cellPath = "initialCondition";
	const std::string cellTag = "cell";
	const std::string cellClusterTag = "cellCluster";
	using std::pair;
	pt::ptree tree;

	// in case no initial condition file provided
	if (initialCellFileName.empty())
	{
		init_cell_fill_grid();
		//init_cell_single_center();
	}
	else {
		try {
			std::vector<CancerCell*> cancerInitList;

			pt::read_xml(initialCellFileName, tree, pt::xml_parser::trim_whitespace);
			// get nodes
			BOOST_FOREACH(pt::ptree::value_type const& cell, tree.get_child(cellPath)) {
				if (cell.first == cellTag) {
					icProperty ic = cell.second;
					unsigned int e = ic.get<unsigned int>("celltype");
					unsigned int s = ic.get<unsigned int>("cellstate", 0);
					unsigned int x = ic.get<unsigned int>("x");
					unsigned int y = ic.get<unsigned int>("y");
					unsigned int z = ic.get<unsigned int>("z");
					AgentType type = static_cast<AgentType>(e);
					AgentState state = static_cast<AgentState>(s);
					CellAgent* ptrCell = createOneInitCell(type, state, Coord3D(x, y, z)); //Here is a CellAngent Point and we need to use dynamic cast to safely convert that to CancerCell Pointer.
					if (type == AgentTypeEnum::CELL_TYPE_CANCER)
					{
						cancerInitList.push_back(dynamic_cast<CancerCell*>(ptrCell));
					}
				}
				else if (cell.first == cellClusterTag)
				{
					icProperty ic = cell.second;
					createClusterInitCell(ic);
				}
			}
			setClusterCancerAntigen(antigenData, cancerInitList); // Only changed the cells for createOneInitCell function, Not sure for the cluster.
		}
		catch (std::exception & e) {
			std::cerr << "Error creating initial cells" << std::endl;
			//std::cerr << e.what() << std::endl;
			throw e;
		}
	}

	// check and update stats
	for (auto ptrCell : _vecAgent)
	{
		_stats.incCellState(ptrCell->getType(), ptrCell->getState());
		if (ptrCell->isDead())
		{
			cerr << "dead cell initiated" << endl;
			exit(1);
		}
	}
}

/*! 
*/
void Tumor::init_cell_single_center(void){
	auto crd = Coord3D(_sizeX, _sizeY, _sizeZ) / 2;
	AgentType type = AgentTypeEnum::CELL_TYPE_CANCER;
	AgentState state = AgentStateEnum::CANCER_STEM;
	createOneInitCell(type, state, crd);
	//_stats.incRecruit(type, state);
	return;
}

/*! initial cell: fill grid with random population
*/
void Tumor::init_cell_fill_grid(void){
	int nr_voxel = _agGrid.getSize();
	std::vector<CancerCell*> cancerInitList; // Create a place holder for all cancer cell being created.
	/*
	// For testing
	unsigned int k = 100;
	std::vector<int> voxel_id(nr_voxel);
	std::iota(std::begin(voxel_id), std::end(voxel_id), 0); 
	//std::shuffle(voxel_id.begin(), voxel_id.begin()+k, rng.getRNG());
	rng.shuffle_first_k<int>(voxel_id, k);
	for (size_t i = 0; i < k; i++)
	{
		Coord3D crd = _agGrid.get_coord(voxel_id[i]);
		//std::cout << "Voxel_id: " << voxel_id[i] << ", " << crd << endl;
		createOneInitCell(AgentTypeEnum::CELL_TYPE_T, AgentStateEnum::T_CELL_CYT, crd);
	}*/ 

	//std::cout << "start population" << std::endl;
	Coord c_sum(0, 0, 0);
	unsigned int n = 0;
	for_each_grid_coord(true, true, true, [&](Coord3D& c){ // Start of the Lambda
		CancerCell* pCancerCell = populate_voxel_random(c);
		if (pCancerCell)
		{
			//std::cout << "CancerCell created at " << pCancerCell->getCoord() << std::endl;
			cancerInitList.push_back(pCancerCell);
			c_sum = c_sum + c;
			n += 1;
		}
		return;
	});
	std::cout << "Number of Cancer cell being Created " << cancerInitList.size() << std::endl;
	setClusterCancerAntigen(antigenData, cancerInitList); // When the entire cancer cells are created, we assign the antigen vector to all cancer cells.
	//_center_target = n ? c_sum / n : c_sum;
	_center_target = Coord3D(_sizeX / 2, _sizeY / 2, _sizeZ / 4);
	/* Testing weather antigens are successfully assigned compared with original csv file
	for(size_t i = 0; i < cancerInitList.size(); i++)
	{
		auto temp = cancerInitList.at(i)->getAntigenVec();
		std::cout << "Cancer number " << i << "[ ";
		for (size_t j = 0; j < temp.size(); j++) 
		{
			std::cout << temp.at(j); 
		}
		std::cout << " ]" << std::endl;
	} */
	//std::cout << _center_target << std::endl;
	return;
}

/*before it return the boolean
bool tumor::populate_voxel_random(const coord3d& crd){
	agenttype type = agenttypeenum::agent_dummy;
	agentstate state = agentstateenum::default_state;
	bool create_cancer_cell = false;
	int div = 0;
	if (_voxel_ic.get_type_state(crd, rng, type, state, div)){
		cellagent * ptrcell = createoneinitcell(type, state, crd);

		if (ptrcell && type == agenttypeenum::cell_type_cancer)
		{
			create_cancer_cell = true;
			pcancercell = dynamic_cast<cancercell*>(ptrcell);
			if (state == agentstateenum::cancer_progenitor)
			{
				pcancercell->setdivcounter(div);
				pcancercell->randomize_div_cd(params.getval(param_int_cancer_cell_progenitor_div_interval_slice));
			}
			else if (state == agentstateenum::cancer_stem)
			{
				pcancercell->randomize_div_cd(params.getval(param_int_cancer_cell_stem_div_interval_slice));
			}
		}
	}
	return create_cancer_cell;
}
*/

/*Current version: returning the pointer instead of the boolean*/
CancerCell* Tumor::populate_voxel_random(const Coord3D& crd) {
	AgentType type = AgentTypeEnum::AGENT_DUMMY;
	AgentState state = AgentStateEnum::DEFAULT_STATE;
	//bool create_cancer_cell = false;
	CancerCell* pCancerCell = NULL;
	int div = 0;
	if (_voxel_ic.get_type_state(crd, rng, type, state, div)) {
		CellAgent* ptrCell = createOneInitCell(type, state, crd);

		if (ptrCell && type == AgentTypeEnum::CELL_TYPE_CANCER)
		{
			//create_cancer_cell = true;
			pCancerCell = dynamic_cast<CancerCell*>(ptrCell);
			if (state == AgentStateEnum::CANCER_PROGENITOR)
			{
				pCancerCell->setDivCounter(div);
				pCancerCell->randomize_div_cd(int(params.getVal(PARAM_FLOAT_CANCER_CELL_PROGENITOR_DIV_INTERVAL_SLICE) + .5));
			}
			else if (state == AgentStateEnum::CANCER_STEM)
			{
				pCancerCell->randomize_div_cd(int(params.getVal(PARAM_FLOAT_CANCER_CELL_STEM_DIV_INTERVAL_SLICE) + .5));
			}
		}
	}
	return pCancerCell;
}


/*! recruit a cluster of cells to the lattice
	\param [in] ic: initial condition

	-# create an empty queue; count = 0
	-# push indicated location to the queue
	-# recruit first cell to lattice, count++
	-# iterate until count matched or queue is empty
	  -# deque, get coordinate
	  -# iterate until no space found
		-# push found location to queue
		-# recruit, count++
*/
void Tumor::createClusterInitCell(icProperty &ic) {

	unsigned int nrCellToCreate = ic.get<unsigned int>("count");
	//ElementType e; not working, cannot directly map string in property_tree to enum
	unsigned int e = ic.get<unsigned int>("celltype");
	unsigned int s = ic.get<unsigned int>("cellstate", 0);
	unsigned int x = ic.get<unsigned int>("x");
	unsigned int y = ic.get<unsigned int>("y");
	unsigned int z = ic.get<unsigned int>("z");


	AgentType type = static_cast<AgentType>(e);
	AgentState state = static_cast<AgentState>(s);
	const ShapeBase *shape;
	switch (type)
	{
	case AgentTypeEnum::AGENT_DUMMY:
		break;
	case AgentTypeEnum::CELL_TYPE_CANCER:
		shape = _cInitDummy->getCellShape();
		break;
	case AgentTypeEnum::CELL_TYPE_T:
		shape = _tInitDummy->getCellShape();
		break;	
	case AgentTypeEnum::CELL_TYPE_MDSC:
		shape = _MDSCInitDummy->getCellShape();
		break;	
	default:
		throw std::invalid_argument("unknown cell type in initial cells");
	}

	unsigned int count = 0;
	std::queue<Coord> nextCenter;
	auto crd = Coord(x, y, z);
	nextCenter.push(crd);
	createOneInitCell(type, state, crd);
	count++;

	while (count < nrCellToCreate || nextCenter.empty())
	{
		Coord c = nextCenter.front();
		nextCenter.pop();
		bool spaceFound = true;
		Coord cNew;
		while (count < nrCellToCreate)
		{
			// create one cell
			int idx;
			//spaceFound = getOneOpenDestinationByScan(shape->getProlifNewOccupy(),
			//	shape->getProlifSearchSeq(), shape->getProlifOccupyMap(),
			//	shape->getProlifRelocateMap(), c, ElementType(e), idx);
			spaceFound = getOneOpenVoxel(shape->getProlifDestinationVoxels(),
				shape->getProlifDestinationAnchor(), c, type, idx, rng);
			if (spaceFound)
			{
				cNew = shape->getProlifDestinationAnchor()[idx] + c;
				createOneInitCell(type, state, cNew);
				count++;
				nextCenter.push(cNew);
			}
			else {
				break;
			}
		}
	}
}

/*! recruit one cell to grid
	\param [in] e: cell type to create
	\param [in] state: cell state
	\param [in] crd: 3D coordinate, const
*/
CellAgent* Tumor::createOneInitCell(AgentType type, AgentState state, const Coord3D& crd) {

	CellAgent * ptrCell = NULL;

	//std::cout << "type: " << type << "; state: " << state << std::endl;

	// only add cell if the target voxel can take it
	if (_agGrid(crd)->isOpenToType(type))
	{
		switch (type)
		{
		case AgentTypeEnum::AGENT_DUMMY:
			throw std::invalid_argument("dummy type in initial cells");
			break;
		case AgentTypeEnum::CELL_TYPE_CANCER:
			ptrCell = _cInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_T:
			ptrCell = _tInitDummy->createCellCopy();
			break;	
		case AgentTypeEnum::CELL_TYPE_MAC:
			ptrCell = _macInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_FIB:
			ptrCell = _fibInitDummy->createCellCopy();
			break;
		case AgentTypeEnum::CELL_TYPE_MDSC:
			ptrCell = _MDSCInitDummy->createCellCopy();
			break;			
		default:
			throw std::invalid_argument("unknown cell type in initial cells");
		}

		ptrCell->setCoord(crd);
		ptrCell->setAgentState(state);
		addAgentToGrid(crd, ptrCell);
		// reset state
		// now all starts from base state; otherwise report error
		/*
		if (state != AgentStateEnum::DEFAULT_STATE)
		{
		throw std::invalid_argument("invalide initial cell state");
		}
		*/
		// randomize life
		switch (type)
		{
		case AgentTypeEnum::AGENT_DUMMY:
			break;
		case AgentTypeEnum::CELL_TYPE_CANCER:{
			auto pCancerCell = dynamic_cast<CancerCell*>(ptrCell);
			if (state == AgentStateEnum::CANCER_PROGENITOR)
			{
				pCancerCell->setProgenitor();
				pCancerCell->setDivCounter(rng.getInt(1 + params.getVal(PARAM_CANCER_CELL_PROGENITOR_DIV_MAX)));
			}
			else if (state == AgentStateEnum::CANCER_SENESCENT)
			{
				pCancerCell->setSenescent();
			}
			break;
		}
		case AgentTypeEnum::CELL_TYPE_T:{
			auto pTCell = dynamic_cast<TCell*>(ptrCell);
			int life = (int)(rng.get_unif_01() * pTCell->getCellLife() + 0.5);
			pTCell->setCellLife(life);
			if (state == AgentStateEnum::T_CELL_CYT)
			{
				//pTCell->setup_source_IFNg();
				pTCell->setup_chem_source(pTCell->get_source_IFNg(), 
					CHEM_IFN, params.getVal(PARAM_IFN_G_RELEASE));
				pTCell->setup_chem_source(pTCell->get_source_IL_2(), 
					CHEM_IL_2, params.getVal(PARAM_IL_2_RELEASE));
			}
			break;
		}
		default:
			break;
		}
		// add to cell vector
		_vecAgent.push_back(ptrCell);
		_stats.incDropIn(type, state);
	}

	return ptrCell;
}


/*! return variables needed for QSP module
	already updated during simulation.
*/
const std::vector<double>& Tumor::get_var_exchange(void){
	return _var_abm_to_qsp;
}

/*! return collected antigen variables needed for QSP module
	already updated during simulation.
*/
const std::vector<double>& Tumor::get_var_antigen_exchange(void) {
	return _var_abm_to_qsp_antigen;
}

const std::vector<double>& Tumor::get_var_teff_exchange(void) {
	return _var_abm_to_qsp_teff;
}

/*! update ABM module with variables from QSP 
*/
void Tumor::update_abm_with_qsp(const std::vector<double>& qsp_var, const std::vector<double>& qsp_var_teff, const std::vector<double>& qsp_var_eff_tum){
	_concentration_cc = qsp_var[LymphCentral::QSPEX_TUM_C] / params.getVal(PARAM_AVOGADROS);
	_distribution_t_cyt = qsp_var_teff;
	_distribution_t_cyt_tum = qsp_var_eff_tum;
	std::cout << "T-cell distribution in central comaprtment: [ ";
	for (size_t i = 0; i < qsp_var_teff.size(); i++)
	{
		std::cout << qsp_var_teff.at(i) << " ";
	}
	std::cout << " ]"<< std::endl;

	std::cout << "\n";

	std::cout << "T-cell distribution in tumor comaprtment: [ ";
	for (size_t i = 0; i < qsp_var_eff_tum.size(); i++)
	{
		std::cout << qsp_var_eff_tum.at(i) << " ";
	}
	std::cout << " ]" << std::endl;

	//_concentration_t_cyt = *std::max_element(_distribution_t_cyt.begin(), _distribution_t_cyt.end());
	_concentration_t_cyt = std::accumulate(_distribution_t_cyt.begin(), _distribution_t_cyt.end(), 0.0);
	_concentration_t_cyt_tum = std::accumulate(_distribution_t_cyt_tum.begin(), _distribution_t_cyt_tum.end(), 0.0);
	_concentration_t_reg = qsp_var[LymphCentral::QSPEX_CENT_TREG];
	_concentration_t_reg_tum = qsp_var[LymphCentral::QSPEX_TUM_TREG];
	_concentration_mdsc = qsp_var[LymphCentral::QSPEX_TUM_MDSC];
	_concentration_nivo = qsp_var[LymphCentral::QSPEX_TUM_NIVO];
	_concentration_ent = qsp_var[LymphCentral::QSPEX_TUM_ENT];
	_concentration_cx = qsp_var[LymphCentral::QSPEX_TUM_CX];
	_concentration_t_exh = qsp_var[LymphCentral::QSPEX_TUM_TEXH];
	_concentration_ccl2 = qsp_var[LymphCentral::QSPEX_TUM_CCL2];
	_tumor_volume = params.getVal(PARAM_VT_MIN) + params.getVal(PARAM_V_CELL) * (_concentration_cc + _concentration_cx) + params.getVal(PARAM_V_TCELL) * (_concentration_t_cyt_tum + _concentration_t_reg_tum + _concentration_t_exh);

	std::cout << "Tumor_Volume:" << _tumor_volume << std::endl;
	return;
}

/*! T recruitment probability
*/
double Tumor::get_T_recruitment_prob(double c, double k_rec) const{
	/*
	p = k (1/(mol * m^3)) * Vol (m^3) * Tum.T (mol)
	*/
	double num_rec = c * k_rec;
	std::cout << "num_rec: " << num_rec << std::endl;
	double p = (num_rec < 1 ? num_rec : 1);

	return p;
}

double Tumor::get_MDSC_recruitment_prob(double c, double k_rec) const {
	/*
	p = k (m^3/mol) * Tum.Vol (m^3) * (Tum.MDSCmax * Tum.Vol - Tum.MDSC) (mol)
	*/
	double num_rec = c * k_rec / _tumor_volume;
	std::cout << "num_rec_MDSC: " << num_rec << std::endl;
	double p;
	if (num_rec > 0) {
		p = (num_rec < 1 ? num_rec : 1);
	}
	else {
		p = 0;
	}

	return p;
}


void Tumor::inc_abm_var_exchange_antigen(std::vector<bool> antigen) { //Increment antigen counts to the counter whenever a cancer cell dies.
	for(size_t i = 0; i < antigen.size(); i++)
	{
		if (antigen.at(i)) 
		{
			_var_abm_to_qsp_antigen.at(i) += 1.0; //when the dead cancer cell have that antigen, _var_abm_to_qsp_antigen increment
		}
	}
}

std::string Tumor::compartment_cells_to_string(void) const {
	std::stringstream ss;
	// header
	ss << "x,y,z,Type,State," << getExtraRemarkHeader() << ",TCR"<<",CancerAntigen"<<endl;
	for (CellVec::size_type i = 0; i < _vecAgent.size(); i++)
	{
		CellAgent* ptr = _vecAgent[i];
		if (ptr->getType() == AgentTypeEnum::CELL_TYPE_T)
		{
			auto pT = dynamic_cast<TCell*>(ptr);
			ss << pT->getAdjustedX() << "," << pT->getAdjustedY()
				<< "," << pT->getAdjustedZ() << "," << pT->getType()
				<< "," << pT->getState() << "," << pT->getRemark() << "," << pT->getTCR() << "," << -1 << endl;
		}
		else if (ptr->getType() == AgentTypeEnum::CELL_TYPE_CANCER)
		{
			auto pT = dynamic_cast<CancerCell*>(ptr);
			ss << pT->getAdjustedX() << "," << pT->getAdjustedY()
				<< "," << pT->getAdjustedZ() << "," << pT->getType()
				<< "," << pT->getState() << "," << pT->getRemark() << "," << -1 <<",";

			for (bool a : pT->getAntigenVec()) {
				ss << a;
			}
			ss << endl;
		}
		else
		{
			ss << ptr->getAdjustedX() << "," << ptr->getAdjustedY()
				<< "," << ptr->getAdjustedZ() << "," << ptr->getType()
				<< "," << ptr->getState() << "," << ptr->getRemark() << "," << -1 <<"," << -1 << endl;
		}
	}
	return ss.str();

}


void Tumor::write_TumorVolume_string(unsigned long slice , std::ostream& stream)const {
	if (slice == 0)
	{
		stream << "Time" << "," << "TumorVolume" << endl;
		stream << slice << "," << get_Tum_Vol() << endl;
	}
	stream << slice << "," << get_Tum_Vol() << endl;
}

};
};
