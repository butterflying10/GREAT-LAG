/**
 * @file         gsimugns.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        simulate gnss observation
 * @version      1.0
 * @date         2025-01-12
 *
 * @copyright Copyright (c) 2024, Wuhan University. All rights reserved.
 *
 */
#ifndef GSIMUGNS_H
#define GSIMUGNS_H

#include <random>
#include "gio/grtlog.h"
#include "gio/gfile.h"

#include "gset/gsetout.h"
#include "gset/gsetrec.h"
#include "gall/gallproc.h"
#include "gall/gallprec.h"
#include "gmodels/gpppmodel.h"
#include "gmodels/gprecisebiasGPP.h"
#include "gset/gsetsimu.h"
#include "gmodels/giono.h"
#include "gexport/ExportLibGREAT.h"
#include "gall/gallobs.h"
#include "gall/galltrp.h"
#include "gcoders/rinexo_leo.h"
#include "gcoders/upd.h"
using namespace std;

namespace great
{
	typedef map<par_type, map<t_gtime, double>> t_map_par;

	/**
	*@brief	   Class for recording ambiguity states and values
	*/
	class amb_value_record
	{
	public:
		/** @brief constructor */
		amb_value_record();

		/** @brief destructor */
		~amb_value_record()
		{};

		/** @brief ambiguity status */
		int status;

		/** @brief integer ambiguity value in band 1 */
		int value_L1;

		/** @brief integer ambiguity value in band 2 */
		int value_L2;

		/** @brief integer ambiguity value in band 3 */
		int value_L3;
	};

	class hw_delay_record
	{
	public:
		hw_delay_record()
		{
			inited = false;
			rand_scale = 0.0;
			phi = 0.0;
		}

		~hw_delay_record() {}

		bool inited;        ///< whether initialized for this satellite
		double rand_scale;  ///< R ~ N(0,1), fixed for each satellite
		double phi;         ///< initial phase for periodic component
	};
	/**
	*@brief	   Class for simulating GNSS observations
	*/
	class LibGREAT_LIBRARY_EXPORT t_gsimugns
	{
	public:
		/**
		* @brief Constructor
		* @note set parameter value
		* @param[in] site            site
		* @param[in] set             setbase
		* @param[in] data         all data
		* @param[in] data         random seed
		*/
		t_gsimugns(string site, t_gsetbase* set, t_gallproc* data, t_spdlog log, int random);

		/** @brief destructor */
		virtual ~t_gsimugns();

	public:
		/**
		* @brief processBatch
		*/
		bool processBatch();
	
	protected:
		t_spdlog _log;             

		t_gsetbase* _set;            
		map< GSYS, map<FREQ_SEQ, GOBSBAND> > _band_index;		///< band index
		map< GSYS, map<GOBSBAND, FREQ_SEQ> > _freq_index;		///< freq index
		map< GSYS, map<GOBSBAND, string> >	_range_order_attr; ///< range order attr
		map< GSYS, map<GOBSBAND, string> >	_phase_order_attr; ///< phase order attr
		map<string, int> _glofrq_num;          ///< glofrq number, set GLO frequency number
		map<string, map<FREQ_SEQ, double> >    _sys_wavelen;   ///< wave length
		int _frequency;								///<  frequency
		//noise setting
		map<GSYS, double> _sigCodeSimu;				///<  code noise sigma
		map<GSYS, double> _sigPhaseSimu;			///<  phase noise sigma
		map<GSYS, double> _maxCodeSimu;				///<  max code noise
		map<GSYS, double> _maxPhaseSimu;			///<  max phase noise
		map<string, double> noise_bias;
		string _site;								///<  simulated site
		set <string> _sats;							///< all sats to simu
		set <string> _systems;						///< all system to simu

		t_gallproc* _data = nullptr;				///< all data
		t_gallobs* _gobs_simu = nullptr;		  ///< all obs data include rinexo
		t_gallnav* _gall_nav = nullptr; 		  ///< all nav data include rinexn,sp3,clk
		t_gallobj* _gallobj;					 ///< all obj data include pco
		t_gupd* _gupd = nullptr;				///< all upd data 
		t_giof* _parfile = nullptr;				///< out par file

		map<t_gtime, map<string, double>> _all_simupars;
		t_gtime _simu_beg_time;			///< begin time for simu data
		t_gtime _simu_end_time;			///< end time for simu data
		int _intv;						///< intv for simu data
		bool _clk = true;			///< simu option
		bool _zwd = true;			///< simu option
		bool _ion = true;			///< simu option
		bool _isb = true;			///< simu option
		bool _ifb = true;			///< simu option
		bool _upd = true;			///< simu option

		CONSTRPAR _crd_est = CONSTRPAR::CONSTRPAR_UNDEF;

		set<string> triple_rm_sat{ "G02", "G05", "G07" , "G12" , "G13" , "G15" , "G16" , "G17" , "G19" , "G20" , "G21" , "G22" , "G28", "G29" , "G31" };


		// code hardware delay simulation
		bool _code_hw_delay = true;       ///< enable pseudorange hardware delay simulation
		bool _phase_hw_delay = true;
		double _hw_const_ns = 40;       ///< constant component amplitude in ns
		double _hw_period_ns = 0.0;       ///< periodic component amplitude in ns
		double _hw_period_sec = 6000.0;   ///< period in seconds, typical LEO value ~100 min
		t_gtime _hw_ref_epoch;            ///< reference epoch for continuous periodic simulation


	protected:
		/**
		* @brief Initializing some pars.
		*/
		t_gallpar initpars();

		/**
		* @brief get pars value from products.
		* @param[in] par_types            
		* @param[in] beg_time            begin time
		* @param[in] end_time            end time
		* @param[in] gdata_simu         all data for simu
		* @param[in] site              site
		*/
		t_map_par _extract_parval(set<par_type> par_types, t_gtime beg_time, t_gtime end_time, t_gallproc* gdata_simu = nullptr, string site = "");

		/**
		* @brief simulate one epoch observation.
		* @param[in] epoch            epoch
		* @param[in] gmodel            precise model for PPP
		* @param[in] params         all pars for simu
		*/
		bool _simuOneEpoObs(t_gtime epoch, t_gprecisebiasGPP *gmodel, t_gallpar params);

		/**
		* @brief simulate rinexo head.
		*/
		bool _simuRnxhdr();

		/**
		* @brief write simulated obs.
		*/
		bool _writeSimuObs();

		/**
		* @brief define par file out path.
		*/
		void _setout();

		/**
		* @brief write all simu pars.
		*/
		void _write_pars();

		/**
		* @brief get upd from upd products.
		* @param[in] epoch            epoch
		* @param[in] upd1            upd in band1
		* @param[in] upd2            upd in band2
		* @param[in] satdata            one sat data
		*/
		bool _get_upd_value(t_gtime epoch, double& upd1, double& upd2, t_gsatdata satdata);

		/**
		* @brief get seconds from hardware delay reference epoch.
		* @param[in] epoch            current epoch
		*/
		double _sec_from_ref(const t_gtime& epoch) const;

		/**
		* @brief simulate pseudorange hardware delay in meters.
		* @param[in] epoch            current epoch
		* @param[in] sat              satellite id
		* @param[in] obs              observation type
		*/
		double _simu_code_hw_delay(const t_gtime& epoch, const string& sat);



		std::default_random_engine      _engine;	///< random engine
		map<string, amb_value_record> _sat_amb_map;	///< map for all sat integer amb record

		map<string, hw_delay_record> _sat_hw_map;   ///< per-satellite code hardware delay states

	};
}
#endif