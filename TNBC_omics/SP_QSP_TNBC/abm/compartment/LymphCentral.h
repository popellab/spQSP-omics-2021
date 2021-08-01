#pragma once

#include <boost/serialization/nvp.hpp>

#include <string>
#include <vector>
#include "SP_QSP_shared/Numerical_Adaptor/CVODE/MolecularModelCVode.h"
#include "TNBC_omics/SP_QSP_TNBC/ode/ODE_system.h"


namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

class LymphCentral
{
typedef CancerVCT::ODE_system LymphBloodQSP;
typedef CancerVCT::Param LymphBloodParam;
typedef MolecularModelCVode<LymphBloodQSP> QSP;

public:
enum QSPExVar{
	QSPEX_TUM_C,
	QSPEX_CENT_TREG,
	QSPEX_TUM_TREG,
	QSPEX_TUM_MDSC,
	QSPEX_TUM_NIVO,
	QSPEX_TUM_ENT,
	QSPEX_TUM_CX,
	QSPEX_TUM_TEXH,
	QSPEX_TUM_CCL2,
	QSPEX_VAR_NUM,
};

enum QSPExVarTeff {
	//Teff cells in the central compartment
	QSPEX_CENT_TEFF1,
	QSPEX_CENT_TEFF2,
	QSPEX_CENT_TEFF3,
	QSPEX_CENT_TEFF4,
	QSPEX_CENT_TEFF5,
	QSPEX_CENT_TEFF6,
	QSPEX_CENT_TEFF7,
	QSPEX_CENT_TEFF8,
	QSPEX_CENT_TEFF9,
	QSPEX_CENT_TEFF10,
	QSPEX_CENT_TEFF11,
	QSPEX_CENT_TEFF12,
	QSPEX_CENT_TEFF13,
	QSPEX_CENT_TEFF14,
	QSPEX_CENT_TEFF15,
	QSPEX_CENT_TEFF16,
	QSPEX_CENT_TEFF17,
	QSPEX_VAR_NUM_TEFF,
};

enum QSPExVarTeffTum {
	//Teff cells in the tumor compartment
	QSPEX_TUM_TEFF1,
	QSPEX_TUM_TEFF2,
	QSPEX_TUM_TEFF3,
	QSPEX_TUM_TEFF4,
	QSPEX_TUM_TEFF5,
	QSPEX_TUM_TEFF6,
	QSPEX_TUM_TEFF7,
	QSPEX_TUM_TEFF8,
	QSPEX_TUM_TEFF9,
	QSPEX_TUM_TEFF10,
	QSPEX_TUM_TEFF11,
	QSPEX_TUM_TEFF12,
	QSPEX_TUM_TEFF13,
	QSPEX_TUM_TEFF14,
	QSPEX_TUM_TEFF15,
	QSPEX_TUM_TEFF16,
	QSPEX_TUM_TEFF17,
	QSPEX_VAR_NUM_TEFF_TUM,
};
public:
	LymphCentral();
	~LymphCentral();

	//! setup parameters of ODE (one time)
	void setup_param(LymphBloodParam& p);

	//! time step 
	void time_step(double t, double dt);
	//! return variable values from QSP module that are needed for ABM
	const std::vector<double>& get_var_exchange(void);

	//! return variable values from QSP module that are needed for ABM for only teffs in blood
	const std::vector<double>& get_var_exchange_teff(void);

	//! return variable values from QSP module that are needed for ABM for only teffs in tumor
	const std::vector<double>& get_var_exchange_teff_tum(void);

	//! update QSP module with output from ABM
	void update_qsp_var(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);
	//! QSP headers
	std::string getQSPHeaders(void)const { return LymphBloodQSP::getHeader();};
	//! write QSP variables 
	friend std::ostream & operator<<(std::ostream &os, const LymphCentral& l);

private:

	//! boost serialization
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);


	//! QSP model excluding tumor dynamics.
	QSP _QSP_model;

	//! variable values passed from QSP to ABM
	std::vector<double> _var_qsp_to_abm;

	//! variable values passed from QSP to ABM for Teffs (central)
	std::vector<double> _var_qsp_to_abm_teff;

	//! variable values passed from QSP to ABM for Teffs (tumor)
	std::vector<double> _var_qsp_to_abm_teff_tum;
};

template<class Archive>
inline void LymphCentral::serialize(Archive & ar, const unsigned int  version) {
	//ar & BOOST_SERIALIZATION_NVP(_cancer_debris);// no need to serialize
	ar & BOOST_SERIALIZATION_NVP(_QSP_model);
	LymphBloodQSP::classSerialize(ar, version);
}

inline std::ostream & operator<<(std::ostream &os, const LymphCentral& l){
	os << l._QSP_model;
	return os;
}

};
};

