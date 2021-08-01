#include "LymphCentral.h"
#include "Tumor.h"
#include "../../core/GlobalUtilities.h"

// shorthands
// get raw value (original units)
#define GET_PARAM_RAW(x) _QSP_model.getSystem()->getParameterVal(x, true)
#define SET_PARAM_RAW(x, y) _QSP_model.getSystem()->setParameterVal(x, y, true)
#define GET_VAR_RAW(x) _QSP_model.getSystem()->getSpeciesVar(x, true)
#define SET_VAR_RAW(x, y) _QSP_model.getSystem()->setSpeciesVar(x, y, true)
// get value (SI units)
#define GET_PARAM(x) _QSP_model.getSystem()->getParameterVal(x, false)
#define SET_PARAM(x, y) _QSP_model.getSystem()->setParameterVal(x, y, false)
#define GET_VAR(x) _QSP_model.getSystem()->getSpeciesVar(x, false)
#define SET_VAR(x, y) _QSP_model.getSystem()->setSpeciesVar(x, y, false)
// parameter (SI units)
#define QSP_CONST(x) LymphBloodQSP::get_class_param(x)

// indices of parameter/variables in their vectors
// y
#define QSP_ID_TUM_C1 52

#define QSP_ID_CENT_TREG 0
#define QSP_ID_CENT_TEFF1 1
#define QSP_ID_CENT_TEFF2 2
#define QSP_ID_CENT_TEFF3 3
#define QSP_ID_CENT_TEFF4 4
#define QSP_ID_CENT_TEFF5 5
#define QSP_ID_CENT_TEFF6 6
#define QSP_ID_CENT_TEFF7 7
#define QSP_ID_CENT_TEFF8 8
#define QSP_ID_CENT_TEFF9 9
#define QSP_ID_CENT_TEFF10 10
#define QSP_ID_CENT_TEFF11 11
#define QSP_ID_CENT_TEFF12 12
#define QSP_ID_CENT_TEFF13 13
#define QSP_ID_CENT_TEFF14 14
#define QSP_ID_CENT_TEFF15 15
#define QSP_ID_CENT_TEFF16 16
#define QSP_ID_CENT_TEFF17 17

#define QSP_ID_TUM_TREG 53

#define QSP_ID_TUM_TEFF1 54
#define QSP_ID_TUM_TEFF2 55
#define QSP_ID_TUM_TEFF3 56
#define QSP_ID_TUM_TEFF4 57
#define QSP_ID_TUM_TEFF5 58
#define QSP_ID_TUM_TEFF6 59
#define QSP_ID_TUM_TEFF7 60
#define QSP_ID_TUM_TEFF8 61
#define QSP_ID_TUM_TEFF9 62
#define QSP_ID_TUM_TEFF10 63
#define QSP_ID_TUM_TEFF11 64
#define QSP_ID_TUM_TEFF12 65
#define QSP_ID_TUM_TEFF13 66
#define QSP_ID_TUM_TEFF14 67
#define QSP_ID_TUM_TEFF15 68
#define QSP_ID_TUM_TEFF16 69
#define QSP_ID_TUM_TEFF17 70

#define QSP_ID_CENT_NIVO 18
#define QSP_ID_TUM_MDSC 80
#define QSP_ID_TUM_NIVO 74
#define QSP_ID_TUM_ENT 84
#define QSP_ID_C  73

#define QSP_ID_P0 142

#define QSP_ID_P1 143
#define QSP_ID_P2 144
#define QSP_ID_P3 145
#define QSP_ID_P4 146
#define QSP_ID_P5 147
#define QSP_ID_P6 148
#define QSP_ID_P7 149
#define QSP_ID_P8 150
#define QSP_ID_P9 151
#define QSP_ID_P10 152
#define QSP_ID_P11 153
#define QSP_ID_P12 154
#define QSP_ID_P13 155
#define QSP_ID_P14 156
#define QSP_ID_P15 157
#define QSP_ID_P16 158
#define QSP_ID_P17 159


#define QSP_ID_TUM_CX 50
#define QSP_ID_TUM_TEXH 51
#define QSP_ID_TUM_CCL2 81

// class_param
#define QSP_C_MAX 15
#define QSP_CELL 9

#define QSP_P0_C1 262
#define QSP_P1_C1 277
#define QSP_P2_C1 291
#define QSP_P3_C1 305
#define QSP_P4_C1 319
#define QSP_P5_C1 333
#define QSP_P6_C1 347
#define QSP_P7_C1 361
#define QSP_P8_C1 375
#define QSP_P9_C1 389
#define QSP_P10_C1 403
#define QSP_P11_C1 417
#define QSP_P12_C1 431
#define QSP_P13_C1 445
#define QSP_P14_C1 459
#define QSP_P15_C1 473
#define QSP_P16_C1 487
#define QSP_P17_C1 501

#define QSP_DAMPS 251
#define QSP_N_T0_CLONES 19

#define QSP_N_T1_CLONES 38
#define QSP_N_T2_CLONES 50
#define QSP_N_T3_CLONES 62
#define QSP_N_T4_CLONES 74
#define QSP_N_T5_CLONES 86
#define QSP_N_T6_CLONES 98
#define QSP_N_T7_CLONES 110
#define QSP_N_T8_CLONES 122
#define QSP_N_T9_CLONES 134
#define QSP_N_T10_CLONES 146
#define QSP_N_T11_CLONES 158
#define QSP_N_T12_CLONES 170
#define QSP_N_T13_CLONES 182
#define QSP_N_T14_CLONES 194
#define QSP_N_T15_CLONES 206
#define QSP_N_T16_CLONES 218
#define QSP_N_T17_CLONES 230

#define QSP_VT_MIN 13
#define QSP_INIT_TUM_DIAM 18
#define QSP_VOL_CELL 11
#define QSP_VOL_TCELL 12

// constants
#define AVOGADROS 6.022140857E23 
#define SEC_PER_DAY 86400
#define PI 3.1416

namespace SP_QSP_IO{
namespace SP_QSP_TNBC{

LymphCentral::LymphCentral()
	: _QSP_model()
	, _var_qsp_to_abm()
{
	// Cent.Teff, Cent.Treg, Cent.Nivo
	_var_qsp_to_abm = std::vector<double>(QSPEX_VAR_NUM, 0);
	_var_qsp_to_abm_teff = std::vector<double>(QSPEX_VAR_NUM_TEFF, 0);
	_var_qsp_to_abm_teff_tum = std::vector<double>(QSPEX_VAR_NUM_TEFF_TUM, 0);
}

LymphCentral::~LymphCentral()
{
}

void LymphCentral::setup_param(LymphBloodParam& p){

	bool steadystate = true;
	unsigned int n = _QSP_model.getSystem()->get_num_variables();
	std::vector<double> ss_val(n, 0);

	if (steadystate)
	{
		QSP ss;
		LymphBloodQSP::_QSP_weight = 1;
		LymphBloodQSP::use_steady_state = true;
		LymphBloodQSP::use_resection = false;
		LymphBloodQSP::setup_class_parameters(p);
		ss.getSystem()->setup_instance_tolerance(p);
		ss.getSystem()->setup_instance_varaibles(p);
		ss.getSystem()->eval_init_assignment();

		// run to steady state or until volume condition is met
		double tss = params.getVal(PARAM_QSP_STEADYSTATE) * SEC_PER_DAY;
		double tt = 0;
		double deltatt = params.getVal(PARAM_SEC_PER_TIME_SLICE);
		double tumor_volume = QSP_CONST(QSP_VT_MIN) + (QSP_CONST(QSP_VOL_CELL) * (ss_val[52] + ss_val[50]) + QSP_CONST(QSP_VOL_TCELL) * ss_val[51]) / AVOGADROS;
		//Now Second: calculate all alive tcell (From TCR 1 to TCR 17 and Treg, 18 types of Tcell total)
		for (size_t j = 53; j < 71; j++)
		{
			tumor_volume += (QSP_CONST(QSP_VOL_TCELL) * ss_val[j]) / AVOGADROS;
		}
		double tumor_volume_ref = (PI * std::pow(QSP_CONST(QSP_INIT_TUM_DIAM),3)) / 6;
		while (tt < tss && tumor_volume < tumor_volume_ref)
		{
			ss.solve(tt, deltatt);
			for (size_t i = 0; i < n; i++)
			{
				ss_val[i] = ss.getSystem()->getSpeciesVar(i);
			}

			//Before:
			//tumor_volume = (QSP_CONST(QSP_VOL_CELL) * (ss_va[18] + ss_val[20]) + QSP_CONST(QSP_VOL_TCELL) * (ss_val[19] + ss_val[21] + ss_val[22])) / AVOGADROS; 

			//Now: First: calculate dead cancer cell, alive cancer cell and exhausted Tcell 
			tumor_volume = QSP_CONST(QSP_VT_MIN)+(QSP_CONST(QSP_VOL_CELL) * (ss_val[52] + ss_val[50]) + QSP_CONST(QSP_VOL_TCELL) * ss_val[51]) / AVOGADROS; 
			//Now Second: calculate all alive tcell (From TCR 1 to TCR 17 and Treg, 18 types of Tcell total)
			for (size_t j = 53; j < 71; j++)
			{
				tumor_volume += (QSP_CONST(QSP_VOL_TCELL) * ss_val[j]) / AVOGADROS;
			}
			tt += deltatt;
		}
		std::cout << tumor_volume << " , " << tumor_volume_ref << std::endl;
		if (tumor_volume < tumor_volume_ref)
		{		
			std::cout<<"tumor volume condition is not met"<<std::endl;
			exit(0);
		}	
	}

	// setup

	LymphBloodQSP::_QSP_weight = params.getVal(PARAM_WEIGHT_QSP);
	LymphBloodQSP::use_steady_state = false;
	LymphBloodQSP::setup_class_parameters(p);
	_QSP_model.getSystem()->setup_instance_tolerance(p);
	_QSP_model.getSystem()->setup_instance_varaibles(p);

	// load steady state
	if (steadystate)
	{
		for (size_t i = 0; i < n; i++)
		{
			_QSP_model.getSystem()->setSpeciesVar(i, ss_val[i]);
		}
	}

	_QSP_model.getSystem()->eval_init_assignment();
	_QSP_model.getSystem()->updateVar();
}

/*! solve QSP from t to t + dt
*/
void LymphCentral::time_step(double t, double dt){

	// Pharmacokinetics
	if (params.getVal(PARAM_NIVO_ON) != 0)
	{
		double week = t / (SEC_PER_DAY * params.getVal(PARAM_NIVO_DOSE_INTERVAL_TIME));
		int week_int = floor(week);
		double nivo_dose = params.getVal(PARAM_NIVO_DOSE);
		double cent_nivo = GET_VAR(QSP_ID_CENT_NIVO);

		if (week == week_int)
		{
			cent_nivo += nivo_dose;
			SET_VAR(QSP_ID_CENT_NIVO, cent_nivo);
		}
	}

	
	// solve QSP for dt
	_QSP_model.solve(t, dt);
	
	
	for (int i = 1; i <= 17; ++i) {
		assert(_QSP_model.getSystem()->getSpeciesVar(i) >= 0);
	}

	return;
}

/*! Get QSP variables for ABM.
	
	# Tum.C1 (unit: cell)
    # Cent.Teff (unit: convert from cell to mole)
	# Cent.Treg (unit: convert from cell to mole)
	# Tum.MDSC (unit: convert from cell to mole)
	# Tum.Nivo (unit: convert from 1e-6 mole/m^3 to mole/m^3)
	# Tum.ENT (unit: convert from 1e-6 mole/m^3 to mole/m^3)
    # Tum.Cx (unit: convert from cell to mole)
	# Tum.Texh (unit: convert from cell to mole)
	# Tum.CCL2 (unit: convert from 1e-6 mole/m^3 to mole/m^3)	
*/
const std::vector<double>& LymphCentral::get_var_exchange(void){

	// need to be cell count for calculating ABM scalor
	_var_qsp_to_abm[QSPEX_TUM_C] = GET_VAR_RAW(QSP_ID_TUM_C1);
	// internal SI unit for calculating rates and probabilities.
	_var_qsp_to_abm[QSPEX_CENT_TREG] = GET_VAR(QSP_ID_CENT_TREG);
	_var_qsp_to_abm[QSPEX_TUM_TREG] = GET_VAR(QSP_ID_TUM_TREG);
	_var_qsp_to_abm[QSPEX_TUM_MDSC] = GET_VAR(QSP_ID_TUM_MDSC);
	_var_qsp_to_abm[QSPEX_TUM_NIVO] = GET_VAR(QSP_ID_TUM_NIVO);
	_var_qsp_to_abm[QSPEX_TUM_ENT] = GET_VAR(QSP_ID_TUM_ENT);
	_var_qsp_to_abm[QSPEX_TUM_CX] = GET_VAR(QSP_ID_TUM_CX);
	_var_qsp_to_abm[QSPEX_TUM_TEXH] = GET_VAR(QSP_ID_TUM_TEXH);	
	_var_qsp_to_abm[QSPEX_TUM_CCL2] = GET_VAR(QSP_ID_TUM_CCL2);	

	
	return _var_qsp_to_abm;
}

const std::vector<double>& LymphCentral::get_var_exchange_teff(void) {
	// Teff from central compartment
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF1] = GET_VAR(QSP_ID_CENT_TEFF1);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF2] = GET_VAR(QSP_ID_CENT_TEFF2);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF3] = GET_VAR(QSP_ID_CENT_TEFF3);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF4] = GET_VAR(QSP_ID_CENT_TEFF4);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF5] = GET_VAR(QSP_ID_CENT_TEFF5);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF6] = GET_VAR(QSP_ID_CENT_TEFF6);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF7] = GET_VAR(QSP_ID_CENT_TEFF7);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF8] = GET_VAR(QSP_ID_CENT_TEFF8);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF9] = GET_VAR(QSP_ID_CENT_TEFF9);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF10] = GET_VAR(QSP_ID_CENT_TEFF10);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF11] = GET_VAR(QSP_ID_CENT_TEFF11);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF12] = GET_VAR(QSP_ID_CENT_TEFF12);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF13] = GET_VAR(QSP_ID_CENT_TEFF13);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF14] = GET_VAR(QSP_ID_CENT_TEFF14);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF15] = GET_VAR(QSP_ID_CENT_TEFF15);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF16] = GET_VAR(QSP_ID_CENT_TEFF16);
	_var_qsp_to_abm_teff[QSPEX_CENT_TEFF17] = GET_VAR(QSP_ID_CENT_TEFF17);


	return _var_qsp_to_abm_teff;
}

const std::vector<double>& LymphCentral::get_var_exchange_teff_tum(void) {
	// Teff from central compartment
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF1] = GET_VAR(QSP_ID_TUM_TEFF1);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF2] = GET_VAR(QSP_ID_TUM_TEFF2);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF3] = GET_VAR(QSP_ID_TUM_TEFF3);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF4] = GET_VAR(QSP_ID_TUM_TEFF4);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF5] = GET_VAR(QSP_ID_TUM_TEFF5);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF6] = GET_VAR(QSP_ID_TUM_TEFF6);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF7] = GET_VAR(QSP_ID_TUM_TEFF7);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF8] = GET_VAR(QSP_ID_TUM_TEFF8);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF9] = GET_VAR(QSP_ID_TUM_TEFF9);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF10] = GET_VAR(QSP_ID_TUM_TEFF10);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF11] = GET_VAR(QSP_ID_TUM_TEFF11);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF12] = GET_VAR(QSP_ID_TUM_TEFF12);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF13] = GET_VAR(QSP_ID_TUM_TEFF13);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF14] = GET_VAR(QSP_ID_TUM_TEFF14);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF15] = GET_VAR(QSP_ID_TUM_TEFF15);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF16] = GET_VAR(QSP_ID_TUM_TEFF16);
	_var_qsp_to_abm_teff_tum[QSPEX_TUM_TEFF17] = GET_VAR(QSP_ID_TUM_TEFF17);


	return _var_qsp_to_abm_teff_tum;
}

/*! update QSP module with output from ABM
	unit convert: item to mole
    # cancer cell death (total)
	# cancer cell death (Teff kill)
	# Teff recruitment
	# Treg recruitment
*/
void LymphCentral::update_qsp_var(const std::vector<double>& var_abm, const std::vector<double>& var_abm_antigen, const std::vector<double>& var_abm_teff){

	// Indicies for all central effector T-cells
	std::vector<int> QSP_CENT_TEFF_Indices{ QSP_ID_CENT_TEFF1,QSP_ID_CENT_TEFF2,QSP_ID_CENT_TEFF3,QSP_ID_CENT_TEFF4,QSP_ID_CENT_TEFF5,QSP_ID_CENT_TEFF6 ,QSP_ID_CENT_TEFF7, QSP_ID_CENT_TEFF8 ,QSP_ID_CENT_TEFF9 ,QSP_ID_CENT_TEFF10 ,QSP_ID_CENT_TEFF11, QSP_ID_CENT_TEFF12,QSP_ID_CENT_TEFF13, QSP_ID_CENT_TEFF14, QSP_ID_CENT_TEFF15, QSP_ID_CENT_TEFF16, QSP_ID_CENT_TEFF17 };
	// Indicies for all antigens will be transported to the lymphnode
	std::vector<int> QSP_ANTIGEN_Indices{ QSP_ID_P1 ,QSP_ID_P2 ,QSP_ID_P3 ,QSP_ID_P4 ,QSP_ID_P5,QSP_ID_P6 ,QSP_ID_P7 ,QSP_ID_P8 ,QSP_ID_P9 ,QSP_ID_P10 ,QSP_ID_P11,QSP_ID_P12, QSP_ID_P13 ,QSP_ID_P14 ,QSP_ID_P15 ,QSP_ID_P16 ,QSP_ID_P17 };
	// Indicies for all antigens will be transported to the lymphnode
	std::vector<int> QSP_TCLONES_Indices{ QSP_N_T1_CLONES, QSP_N_T2_CLONES ,QSP_N_T3_CLONES ,QSP_N_T4_CLONES ,QSP_N_T5_CLONES, QSP_N_T6_CLONES, QSP_N_T7_CLONES ,QSP_N_T8_CLONES ,QSP_N_T9_CLONES ,QSP_N_T10_CLONES ,QSP_N_T11_CLONES, QSP_N_T12_CLONES,QSP_N_T13_CLONES ,QSP_N_T14_CLONES ,QSP_N_T15_CLONES ,QSP_N_T16_CLONES ,QSP_N_T17_CLONES };
	// Indecies for antigen concentration in cancer cell clone 1
	std::vector<int> QSP_Px_C1_Indices{ QSP_P1_C1 ,QSP_P2_C1 ,QSP_P3_C1 ,QSP_P4_C1 ,QSP_P5_C1,QSP_P6_C1 ,QSP_P7_C1 ,QSP_P8_C1 ,QSP_P9_C1 ,QSP_P10_C1 ,QSP_P11_C1,QSP_P12_C1, QSP_P13_C1 ,QSP_P14_C1 ,QSP_P15_C1 ,QSP_P16_C1 ,QSP_P17_C1 };
	// convert item to internal units
	double scalar = 1 / AVOGADROS;
	std::cout << "DEAD CANCER CELL: " << var_abm[Tumor::TUMEX_CC_DEATH] << std::endl;
	std::cout << "\n";
	// CC death total, CC death Teff, Teff recruit, Treg recruit
	double cc_death_total = var_abm[Tumor::TUMEX_CC_DEATH] * scalar;
	double cc_death_Teff = var_abm[Tumor::TUMEX_CC_T_KILL] * scalar;
	double Treg_recruit = var_abm[Tumor::TUMEX_TREG_REC] * scalar;

	std::vector<double> Teff_recruit = var_abm_teff;
	for (size_t i = 0; i < Teff_recruit.size(); i++)
	{
		Teff_recruit[i] *= scalar;
	}



	// update system
	double tum_c1 = GET_VAR(QSP_ID_TUM_C1);
	double p0 = GET_VAR(QSP_ID_P0);
	double ckine = GET_VAR(QSP_ID_C);

	double factor_DAMP = QSP_CONST(QSP_DAMPS);
	ckine += cc_death_Teff * factor_DAMP;
	SET_VAR(QSP_ID_C, ckine);

	// Update for Treg epitope
	double factor_p0 = QSP_CONST(QSP_N_T0_CLONES) * QSP_CONST(QSP_P0_C1);
	p0 += cc_death_total * factor_p0;
    SET_VAR(QSP_ID_P0, p0);
	

	// FACTOR for Teff epitopes (antigens)
	std::vector<double> factor_p{};
	for (size_t i = 0; i < QSP_ANTIGEN_Indices.size(); i++)
	{
		
		factor_p.push_back(QSP_CONST(QSP_TCLONES_Indices.at(i)) * QSP_CONST(QSP_Px_C1_Indices.at(i)));
	}
	

	// EXPRESSION for epitopes
	std::vector<double> expression_p (QSP_ANTIGEN_Indices.size(), 0 );

	for (size_t i = 0; i < QSP_ANTIGEN_Indices.size(); i++)
	{
		expression_p.at(i) = GET_VAR(QSP_ANTIGEN_Indices.at(i));
		expression_p.at(i) += var_abm_antigen.at(i) * factor_p.at(i) * scalar;  //Doing extra step just to be consistent with previous code.
		
		SET_VAR(QSP_ANTIGEN_Indices.at(i), expression_p.at(i));
	}
	
	

	

	double cent_t_reg = GET_VAR(QSP_ID_CENT_TREG);
	cent_t_reg -= Treg_recruit;
	
	if (cent_t_reg < 0) 
	{
		cent_t_reg = 0;
	}
	
	SET_VAR(QSP_ID_CENT_TREG, cent_t_reg); //Reset for Treg


	std::vector <double> cent_t_eff{};
	for (int Teff_index : QSP_CENT_TEFF_Indices) 
	{
		cent_t_eff.push_back(GET_VAR(Teff_index)); // store the value of each TCR indices in the vector
	}
	
	for (size_t i = 0; i < cent_t_eff.size(); i++)
	{
		cent_t_eff.at(i) -= Teff_recruit.at(i); // subtract central compartment T-cell by T-cell being recruited for each TCR;
		
		if (cent_t_eff.at(i) < 0) 
		{
			cent_t_eff.at(i) = 0;
		}
		SET_VAR(QSP_CENT_TEFF_Indices.at(i), cent_t_eff.at(i)); //  reset the number of T-cell in central compartment
	}
	
	return;
}

};
};

