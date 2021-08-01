#pragma once

#include "SP_QSP_shared/ABM_Base/SpatialCompartment.h"
#include "TumorGridVoxel.h"
#include "VoxelContentGen.h"
#include "../agent/TCell.h"
#include "../agent/CancerCell.h"
#include "../agent/Mac.h"
#include "../agent/Fib.h"
#include "../agent/TReg.h"
#include "../agent/MDSC.h"
#include "../../pde/DiffuseGrid.h"
#include "../../core/Stats.h"

#include "../util/csvreader.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/export.hpp>

#include <algorithm> 

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

typedef std::vector<Coord3D> CellSource;
typedef boost::property_tree::ptree icProperty;

//! Generic tumor compartment
class Tumor : public SpatialCompartment
{
public:
enum TumExVar{
	TUMEX_CC,
	TUMEX_CC_DEATH,
	TUMEX_CC_T_KILL,
	TUMEX_TEFF_REC,
	TUMEX_TREG_REC,
	TUMEX_VAR_NUM
};

enum TumExVarTeff {
	TUMEX_TEFF1_REC,
	TUMEX_TEFF2_REC,
	TUMEX_TEFF3_REC,
	TUMEX_TEFF4_REC,
	TUMEX_TEFF5_REC,
	TUMEX_TEFF6_REC,
	TUMEX_TEFF7_REC,
	TUMEX_TEFF8_REC,
	TUMEX_TEFF9_REC,
	TUMEX_TEFF10_REC,
	TUMEX_TEFF11_REC,
	TUMEX_TEFF12_REC,
	TUMEX_TEFF13_REC,
	TUMEX_TEFF14_REC,
	TUMEX_TEFF15_REC,
	TUMEX_TEFF16_REC,
	TUMEX_TEFF17_REC,
	TUMEX_VAR_NUM_TEFF
};

public:
	//! Need default constructor for serialization to work
	Tumor() {};
	Tumor(int x, int y, int z);
	virtual ~Tumor();

	// add one source of entry
	void add_lymphocyte_source(const Coord& c);

	// add one source of entry
	void add_mdsc_source(const Coord& c);	

	//! simulate for one slice
	virtual void timeSlice(unsigned long slice);

	//! recruit T Cells from invasive front 
	void recruitTCellsFront(double min, double max);

	//! get stats
	const Stats& get_stats(void) const{ return _stats; };

	//! default header for extra remark column when writing cell grid to file
	virtual std::string getExtraRemarkHeader() const;
	//! print grid snapshot to file
	std::string printGridToFile() const;
	//! print cell ODE stats to file
	void printCellOdeToFile(unsigned long slice) const;
	//! print grid snapshot to screen
	void printGridToScreen(unsigned long slice) const;
	//! get chemical grid
	DiffuseGrid & get_chem_grid(void) { return _chem; };
	//! get concentration of chemokine i.
	double get_chem(const Coord3D&c, chem_ID i)const; 

	//! add to abm var exchange counter
	void inc_abm_var_exchange(TumExVar v){ _var_abm_to_qsp[v] += 1.0;};

	//! add to abm var exchange counter
	void inc_abm_var_exchange_teff(int n) { _var_abm_to_qsp_teff[n] += 1.0; };

	//! add to abm var antigen exchange counter
	void inc_abm_var_exchange_antigen(std::vector<bool> antigen);

	//! return variables needed for QSP module 
	const std::vector<double>& get_var_exchange(void);

	//! return collected antigen variables needed for QSP module 
	const std::vector<double>& get_var_antigen_exchange(void);

	//! return collected Teff variables needed for QSP module 
	const std::vector<double>& get_var_teff_exchange(void);

	//! update ABM module with variables from QSP 
	void update_abm_with_qsp(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

	//! tumor volume
	double get_Tum_Vol(void)const{ return _tumor_volume; };
	//! treg concentration in the tumor
	double get_Treg(void)const{ return _concentration_t_reg_tum; };
	//! nivo concentration
	double get_Nivo(void)const{ return _concentration_nivo; };
    
    //! ent concentration
	double get_Ent(void)const{ return _concentration_ent; };

	//! create one initial cell and configure to explicity initial condicitons
	void TEST_AddOneCell(AgentType type, AgentState state, const Coord3D& c) {
		createOneInitCell(type, state, c);
	};

	//! set allow shift
	void set_allow_shift(bool allowed){ _allow_shift_grid = allowed; };
	
	//!Set abm scaler function
	//void setABM_scaler(double abm_scaler) { _abm_scaler = abm_scaler;  };

	VoxelContentGen _voxel_ic;
	
	// Antigen related Functions
	int getNumAntigens()const;
	void loadAntigenData(std::string antigenDataFile); // This function loads all antigen datafile provided in TNBC_sim.cpp

	//write tumor volume to the TCR file
	void write_TumorVolume_string(unsigned long slice ,std::ostream& stream)const;

	//override compartment cell to string function
	std::string compartment_cells_to_string(void) const;
protected:
	//! initialize agent grid 
	virtual bool initAgentGrid();
private:

	friend class boost::serialization::access;
	template<class Archive>
	//! boost serialization
	void serialize(Archive & ar, const unsigned int /*version*/);

	void time_slice_recruitment(void);
	void time_slice_movement(double t, double dt);
	void time_slice_state_change(double t, double dt);
	//! scan all agents for the last round in a slice
	void time_slice_final_scan(void);
	void time_slice_molecular(double t);

	//! initialize compartment: setup environmen��t
	void initEnvironment();
	//! initialize compartment: setup initial cells
	void initCell(std::string filename);
	//! initial cell: center cancer cell
	void init_cell_single_center(void);
	//! initial cell: fill grid 
	void init_cell_fill_grid(void);
	//! create one random cell
	CancerCell* populate_voxel_random(const Coord3D&);
	//! create a cluster of initial cells centered around provided coordinate
	void createClusterInitCell(icProperty &ic);
	//! create one initial cell and configure to explicity initial condicitons
	CellAgent* createOneInitCell(AgentType type, AgentState state, const Coord3D& crd);
	//! Teff/Treg recruitment probability at each entry point
	double get_T_recruitment_prob(double c, double base) const;
	//! MDSc recruitment probability at each entry point
	double get_MDSC_recruitment_prob(double c, double base) const;	
  
	//! adjust camera center by shifting grid
	void shift_adjust_center(void);
	//! get the vector of camera shift
	Coord3D get_cam_shift(void);
	//! shift contents of the grid by crd
	void shift_grid(Coord3D& crd);

	// stats
	Stats _stats;
	// Diffusion grid
	DiffuseGrid _chem;
	//! list of T cell sources
	CellSource _t_source;
	//! list of T cell sources
	CellSource _mdsc_source;
	//! Dummy T cell for initialization
	TCell * _tInitDummy;
	//! Dummy MDSC for initialization
	MDSC * _MDSCInitDummy;	
	//! Dummy Cancer cell for initialization
	CancerCell * _cInitDummy;
	Mac * _macInitDummy;
	Fib * _fibInitDummy;
	Treg * _TregInitDummy;
	//! voxel dimension, no need to serialize
	//double _voxel_size;

	/* variables used in grid shifting */

	//! allow shifting
	bool _allow_shift_grid;
	//! mass center anchor
	Coord3D _center_target;
	//! temporary AgentGrid
	AgentGrid _agGrid_temp;

	/* following are variables communicated between tumor and blood */

	//! variable values passed from QSP to ABM (no need to serialize)
	std::vector<double> _var_abm_to_qsp;
	//! collected antigen variable values passed from QSP to ABM (no need to serialize)
	std::vector<double> _var_abm_to_qsp_antigen;
	//! collected Teff values passed from QSP to ABM (no need to serialize)
	std::vector<double> _var_abm_to_qsp_teff;

	//! concentration of cancer cells(no need to serialize).
	double _concentration_cc;
	//! Distribution of All Teff in blood by their TCR.
	std::vector<double> _distribution_t_cyt;
	//! Distribution of All Teff in tumor by their TCR.
	std::vector<double> _distribution_t_cyt_tum;

	//! concentration of cytotoxic t cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_cyt;
	//! concentration of cytotoxic t cells in tumor(Unit: SI; no need to serialize).
	double _concentration_t_cyt_tum;
	//! concentration of regulatory t cells in blood(Unit: SI; no need to serialize).
	double _concentration_t_reg;
	//! concentration of regulatory t cells in the tumor(Unit: SI; no need to serialize).
	double _concentration_t_reg_tum;	
	//! concentration of mdsc(Unit: SI; no need to serialize).
	double _concentration_mdsc;
	//! concentration of dead cancer cells(Unit: SI; no need to serialize).
	double _concentration_cx;
	//! concentration of exhausted t cells(Unit: SI; no need to serialize).
	double _concentration_t_exh;			
	//! concentration of nivo in Tumor(Unit: SI; no need to serialize).
	double _concentration_nivo;
	//! concentration of ent in Tumor(Unit: SI; no need to serialize).
	double _concentration_ent;
	//! concentration of ccl2 in Tumor(Unit: SI; no need to serialize).	
	double _concentration_ccl2;	
	//! tumor volume(Unit: SI; no need to serialize).
	double _tumor_volume;
	

    double Da;	
    // Static Variable just to determine the number of antigens for tumor
	//static int num_antigens();
	std::vector<std::vector<std::string>> antigenData;

};

//BOOST_CLASS_EXPORT_KEY(Tumor);

template<class Archive>
inline void Tumor::serialize(Archive & ar, const unsigned int /* version */) {
	ar.template register_type<CancerCell>();
	ar.template register_type<TCell>();
	ar.template register_type<MDSC>();
	ar.template register_type<Mac>();
	ar.template register_type<Fib>();
	ar.template register_type<Treg>();
	ar.template register_type<TumorGridVoxel>();
	//ar.template register_type<BioFVMSinkSource>();
	// this is needed because cell._compartment is of type SpatialCompartment*
	ar.template register_type<Tumor>();
	ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SpatialCompartment);
	ar & BOOST_SERIALIZATION_NVP(_voxel_ic);
	ar & BOOST_SERIALIZATION_NVP(_stats);
	ar & BOOST_SERIALIZATION_NVP(_t_source);
	ar & BOOST_SERIALIZATION_NVP(_mdsc_source);
	ar & BOOST_SERIALIZATION_NVP(_chem);
	ar & BOOST_SERIALIZATION_NVP(_allow_shift_grid);
	ar & BOOST_SERIALIZATION_NVP(_center_target);
	/*
	ar & BOOST_SERIALIZATION_NVP(_agGrid_temp);
	ar & BOOST_SERIALIZATION_NVP(_tInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_cInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_macInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_fibInitDummy);
	ar & BOOST_SERIALIZATION_NVP(_TregInitDummy);
	*/
}

};
};

