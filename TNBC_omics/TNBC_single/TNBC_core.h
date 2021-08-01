#pragma once

#include "TNBC_omics/SP_QSP_TNBC/abm/compartment/Tumor.h"
#include "TNBC_omics/SP_QSP_TNBC/abm/compartment/LymphCentral.h"
#include "FileOutputHub.h"

#include <boost/serialization/nvp.hpp>
//! Simulation central control
/*! Host all compartments (currently only one: tumor).
	Handle interaction between compartments; process global statistics output.
*/
class TNBC_Core
{
public:
	TNBC_Core();
	TNBC_Core(std::string antigenfile);
	virtual ~TNBC_Core();
	//! setup qsp module
	void setup_qsp(CancerVCT::Param& p);
	//! Initialize all compartments: from file
	void initializeSimulation(std::string);
	//! Initialize all compartments: random population 
	void initializeSimulation(void);
	//! Proceed on slice in time for all compartments
	virtual void timeSlice(const long slice);
	//! Write ABM stats header to stats file
	void write_stats_header(void) const;
	//! Write ABM stats of current time slice to stats file
	void write_stats_slice(unsigned long slice) const;
	//! write QSP content to file
	void write_QSP(unsigned long slice, bool header)const;
	//! write Tumor volume to file
	void write_TumorVolume(unsigned long slice)const;
	//! Write ODE stats of current time slice to ode stats file
	virtual void writeOde(unsigned long slice);
	//! Write grid snapshots of current time slice to file 
	virtual void writeGrids(unsigned long slice, unsigned int option);
	//! Print brief summary of current slice
	virtual void briefStats(unsigned long slice);

private:
	friend class boost::serialization::access;
	//! boost serialization 
	template<class Archive>
	void serialize(Archive & ar, const unsigned int /*version*/);

	//! tumor compartment instance

	SP_QSP_IO::SP_QSP_TNBC::Tumor _tumor;
	SP_QSP_IO::SP_QSP_TNBC::LymphCentral _lymph;
};

template<class Archive>
inline void TNBC_Core::serialize(Archive & ar, const unsigned int version){
	ar & BOOST_SERIALIZATION_NVP(_tumor);
	ar & BOOST_SERIALIZATION_NVP(_lymph);
	SP_QSP_IO::CellAgent::classSerialize(ar, version);
}