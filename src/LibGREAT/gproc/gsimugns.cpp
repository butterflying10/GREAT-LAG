/**
*
* @verbatim
	History
	 -1.0 xjhan		2021-05-16  creat
  @endverbatim
* Copyright (c) 2018, Wuhan University. All rights reserved.
*
* @file		    gsimugns.cpp
* @brief		obs simu processing for GNSS satellite
* @author       xjhan, Wuhan University
* @version		1.0.0
* @date		    2021-05-16
*
*/

#include "gutils/gstring.h"
#include "gsimugns.h"
#include <memory>


#define FAKE_OBS_VALUE 100.0

namespace great
{
	double interp(t_gtime epoch, map<t_gtime, double> items, int sample);
	amb_value_record::amb_value_record()
	{
		status = 0;
		value_L1 = 0;
		value_L2 = 0;
		value_L3 = 0;
	}

	t_gsimugns::t_gsimugns(string site, t_gsetbase* set, t_gallproc* data, t_spdlog log, int random) :
		_site(site),
		_set(set),
		_data(data),
		_log(log),
		_engine(random)
	{
		_sats.clear();
		_systems = dynamic_cast<t_gsetgen*>(_set)->sys();
		_sats = dynamic_cast<t_gsetgnss*>(set)->sat();
		for (auto sys : _systems)
		{
			GSYS SYS = t_gsys::str2gsys(sys);

			_band_index[SYS] = dynamic_cast<t_gsetgnss*>(set)->band_index(SYS);
			_freq_index[SYS] = dynamic_cast<t_gsetgnss*>(set)->freq_index(SYS);
			_sigCodeSimu[SYS] = dynamic_cast<t_gsetgnss*>(set)->sigC_simu(SYS);
			_maxCodeSimu[SYS] = dynamic_cast<t_gsetgnss*>(set)->maxC_simu(SYS);
			_sigPhaseSimu[SYS] = dynamic_cast<t_gsetgnss*>(set)->sigL_simu(SYS);
			_maxPhaseSimu[SYS] = dynamic_cast<t_gsetgnss*>(set)->maxL_simu(SYS);

			SPDLOG_LOGGER_INFO(log, "t_gsimugns", "range noise  ", sys , " max " , format("%16.8f", _maxCodeSimu[SYS]) , " sig " , format("%16.8f", _sigCodeSimu[SYS]));
			SPDLOG_LOGGER_INFO(log, "t_gsimugns", "phase noise  ", sys , " max " , format("%16.8f", _maxPhaseSimu[SYS]) , " sig " , format("%16.8f", _sigPhaseSimu[SYS]));
		}
		std::set<string> sats_rm = dynamic_cast<t_gsetgen*>(_set)->sat_rm();
		for (auto sat: _sats)
		{
			string tmp = t_gsys::gsys2str(t_gsys::str2gsys(sat.substr(0,1)));
			if (_systems.find(tmp) == _systems.end())
				sats_rm.insert(sat);
		}
		
		for (auto it : sats_rm) 
		{
			if (_sats.find(it) != _sats.end())
			{
				_sats.erase(it);
			}
		}
		for (auto it : _sats)
		{
			if (_sat_amb_map.find(it) == _sat_amb_map.end())
			{
				amb_value_record onerecord;
				_sat_amb_map[it] = onerecord;
			}
		}

		_range_order_attr = gnut::range_order_attr_raw;
		_phase_order_attr = gnut::phase_order_attr_raw;

		_clk = dynamic_cast<t_gsetsimu*>(_set)->clk();
		_ion = dynamic_cast<t_gsetsimu*>(_set)->ion();
		_zwd = dynamic_cast<t_gsetsimu*>(_set)->zwd();
		_upd = dynamic_cast<t_gsetsimu*>(_set)->upd();

		// set glofrq number 
		_glofrq_num = dynamic_cast<t_gallobs*>((*data)[t_gdata::ALLOBS])->glo_freq_num();
		if (_glofrq_num.empty())
		{
			_glofrq_num = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM])->glo_freq_num();
		}

		_gobs_simu = new t_gallobs();  _gobs_simu->spdlog(log); _gobs_simu->gset(set);
		_all_simupars.clear();
		_setout();
		for (auto sat : _sats)
		{
			noise_bias[sat] = std::normal_distribution<double>(0.0, 0.0001)(_engine);
		}
		_gall_nav = dynamic_cast<t_gallnav*>((*data)[t_gdata::GRP_EPHEM]);
		_gallobj = dynamic_cast<t_gallobj*>((*data)[t_gdata::ALLOBJ]);
		_gupd = dynamic_cast<t_gupd*>((*data)[t_gdata::UPD]);
	}

	t_gsimugns::~t_gsimugns()
	{
		if (_gobs_simu) { delete _gobs_simu; _gobs_simu = nullptr; }
		
		if (_parfile)
		{
			if (_parfile->is_open()) _parfile->close();
			delete _parfile;
			_parfile = nullptr;
		}
	}

	bool t_gsimugns::processBatch()
	{
		_simu_beg_time = dynamic_cast<t_gsetgen*>(_set)->beg();			///< begin time for simu data
		_simu_end_time = dynamic_cast<t_gsetgen*>(_set)->end();			///< end time for simu data
		_intv = dynamic_cast<t_gsetgen*>(_set)->sampling();				///< intv for simu data

		/*init simu params*/
		t_gallpar params = initpars();

		/*extract required params value from CLK TRP ...product*/
		set<par_type> types;
		if (_clk) types.insert(par_type::CLK);
		if (_zwd) types.insert(par_type::TRP);
		t_map_par parinfo = _extract_parval(types, _simu_beg_time, _simu_end_time, _data, _site);
		
		/*read crd xyz*/
		t_gtriple XYZ = dynamic_cast<t_gsetrec*>(_set)->get_crd_xyz(_site);
		if (double_eq(XYZ[0], 0.0) && double_eq(XYZ[1], 0.0) && double_eq(XYZ[2], 0.0))
		{
			SPDLOG_LOGGER_INFO(_log, "simugns", "ERROR : no crd for " + _site);
			return false;
		}

		/*avoid simu windup arp error*/
		t_gprecisebiasGPP gmodel(_data, _log, _set);

		/* init currunt time */
		t_gtime epoch(_simu_beg_time); 
		double sigCLK = dynamic_cast<t_gsetsimu*>(_set)->sig_simu_clk();
		/* time loop */
		while (epoch <= _simu_end_time)
		{

			/* display progress */
			if (epoch.sow() % 7200 == 0) cerr << "[great_simuobs] " << _site << " " << epoch.str_ymdhms() << "finish " << std::endl;
			if (fabs(epoch.dsec() - 1.0) < 1e-6) { epoch.reset_dsec(); epoch.add_secs(1); }

			//set simu params value
			for (int i = 0; i < params.parNumber(); i++)
			{
				string str_type = params.getPar(i).str_type();

				//skip AMB, GLO_IFB params
				if (str_type.find("AMB") != -1 || str_type.find("GLO_IFB") != -1) continue;

				//skip AMB, GLO_IFB params
				par_type parType = params.getPar(i).parType;
				switch (parType)
				{
				case par_type::CRD_X:
					if (_crd_est != CONSTRPAR::KIN) { params[i].value(XYZ[0]);  break; }						
				case par_type::CRD_Y:
					if (_crd_est != CONSTRPAR::KIN) { params[i].value(XYZ[1]);  break; }
				case par_type::CRD_Z:
					if (_crd_est != CONSTRPAR::KIN) { params[i].value(XYZ[2]);  break; }
				case par_type::CLK:
					if (_crd_est != CONSTRPAR::KIN) 
					{ 
						if (_clk)params[i].value(interp(epoch, parinfo[parType], 30));
						else params[i].value(std::normal_distribution<double>(0.0, sigCLK)(_engine));
						break; 
					}
				default:
					if (parinfo.find(parType) != parinfo.end()) params[i].value(interp(epoch, parinfo[parType], _intv));
					else params[i].value(0.0);
					break;
				}
				params[i].apriori(1e-3);
				_all_simupars[epoch][str_type] = params[i].value();
			}



			//simu one epoch gnss obs
			if (!_simuOneEpoObs(epoch, &gmodel, params))
			{
				SPDLOG_LOGGER_INFO(_log, "simugns", "ERROR: simu obs fail " + _site + " " + epoch.str_hms());
			}
			epoch.add_dsec(_intv);
		}

		/*write simulated rinexo file*/
		bool simu_hdr = _simuRnxhdr();
		bool write_simu = _writeSimuObs();

		/*write par file*/
		_write_pars();

		return true;
	}

	t_gallpar t_gsimugns::initpars()
	{
		t_gallpar params; int ipar = 0;
		params.addParam(t_gpar(_site, par_type::CRD_X, ++ipar, ""));
		params.addParam(t_gpar(_site, par_type::CRD_Y, ++ipar, ""));
		params.addParam(t_gpar(_site, par_type::CRD_Z, ++ipar, ""));
		params.addParam(t_gpar(_site, par_type::CLK, ++ipar, ""));
		params.addParam(t_gpar(_site, par_type::TRP, ++ipar, ""));

		if (_systems.find("GPS") != _systems.end() && _band_index[GSYS::GPS].size() >= 3 && _ifb)
		{
			params.addParam(t_gpar(_site, par_type::IFB_GPS, ++ipar, ""));
		}
		if (_systems.find("GAL") != _systems.end())
		{
			if (_isb) params.addParam(t_gpar(_site, par_type::GAL_ISB, ++ipar, ""));
			if (_band_index[GSYS::GAL].size() >= 3 && _ifb) params.addParam(t_gpar(_site, par_type::IFB_GAL, ++ipar, ""));
		}

		if (_systems.find("BDS") != _systems.end())
		{
			if (_isb) params.addParam(t_gpar(_site, par_type::BDS_ISB, ++ipar, ""));
			if (_band_index[GSYS::BDS].size() >= 3 && _ifb) params.addParam(t_gpar(_site, par_type::IFB_BDS, ++ipar, ""));
		}
		if (_systems.find("LEO") != _systems.end())
		{
			if (_isb) params.addParam(t_gpar(_site, par_type::LEO_ISB, ++ipar, ""));
		}
		
		for (auto sat : _sats)
		{
			vector<double> oneupd(3, 0.0);
			t_gobsgnss obs(_log,_site, sat, t_gtime());
			auto gsys = t_gsys::sat2gsys(sat);

			if (gsys==GSYS::GLO) obs.channel(_glofrq_num[sat]);

			if (_band_index.find(gsys) == _band_index.end()) continue;
			if (_band_index[gsys].find(FREQ_1) != _band_index[gsys].end())
			{
				t_gpar par_amb(_site, par_type::AMB_L1, ++ipar, sat);
				par_amb.value(0.0);
				params.addParam(par_amb);
			}
			if (_band_index[gsys].find(FREQ_2) != _band_index[gsys].end())
			{
				t_gpar par_amb(_site, par_type::AMB_L2, ++ipar, sat);
				par_amb.value(0.0);
				params.addParam(par_amb);
			}
			if (_band_index[gsys].find(FREQ_3) != _band_index[gsys].end() && (triple_rm_sat.find(sat) == triple_rm_sat.end()))
			{
				t_gpar par_amb(_site, par_type::AMB_L3, ++ipar, sat);
				par_amb.value(0.0);
				params.addParam(par_amb);
			}

			

			if (gsys == GSYS::GLO && _ifb)
			{
				t_gpar par_ifb(_site, par_type::GLO_IFB, ++ipar, sat);
				par_ifb.value((10 - _glofrq_num[sat]) * 0.05);
				params.addParam(par_ifb);
			}
		}
		return params;
	}

	t_map_par t_gsimugns::_extract_parval(set<par_type> par_types, t_gtime beg_time, t_gtime end_time, t_gallproc* gdata_simu, string site)
	{
		t_map_par parinfo;
		for (auto type : par_types)
		{
			if (type == par_type::TRP) 
			{
				t_galltrp* _galltrps = dynamic_cast<t_galltrp*>((*gdata_simu)[t_gdata::ALLTRPZTD]);
				t_gtrp* _trps = _galltrps->trp(site);
				if (_trps == nullptr)
				{
					for (t_gtime it_time = beg_time; it_time <= end_time; it_time = it_time + _intv)
					{
						parinfo[type][it_time] = 0.0;
					}
					continue;
				}
				t_gtime beg = _trps->beg();
				t_gtime last_epoch = beg;
				map<t_gtime, t_trp_record>* _trp = _trps->trp();
				for (auto &iter : (* _trp))
				{	
					t_gtime epoch = iter.first;
					double value = iter.second.TROTOT;
					if (last_epoch == epoch) continue;
					if (double_eq(value, 0.0)) continue;
					for (t_gtime it_time = last_epoch; it_time <= epoch; it_time = it_time + _intv)
					{
						parinfo[type][it_time] = value/1000;
					}
					last_epoch = epoch;
				}
			}
			else if (type == par_type::CLK)
			{
				t_gallprec* _gall_prec = dynamic_cast<t_gallprec*>((*gdata_simu)[t_gdata::ALLPREC]);
				t_gtime beg_clk = _gall_prec->beg_clk(site);
				if (beg_clk == LAST_TIME)
				{
					cout << "WARNING: simugns: " << site << " rec clk is not found in rinexc file" << endl;
					for (t_gtime it_time = beg_time; it_time <= end_time; it_time = it_time + _intv)
					{
						parinfo[type][it_time] = 0.0;
					}
					continue;
				}
				for (t_gtime it_time = beg_time; it_time <= end_time; it_time = it_time + 30)
				{
					double clk, var, dclk;
					bool chk_mask;
					double value = _gall_prec->clk(site, it_time,&clk,&var,&dclk, chk_mask);
					if (value>0)
					{
						parinfo[type][it_time] = clk;
					}
					else 
					{
						cout << "WARNING: simugns: " << it_time.mjd() << "  " << it_time.sod() << "  rec clk is not found in rinexc file" << endl;
					}
				}
			}
		}

		return parinfo;
	}

	bool t_gsimugns::_simuOneEpoObs(t_gtime epoch, t_gprecisebiasGPP* gmodel, t_gallpar params)
	{
		//reset prec interpolation points
		dynamic_cast<t_gallprec*>((*_data)[t_gdata::GRP_EPHEM])->reset_prec();

		//sat loop
		for (auto sat : _sats)
		{
			t_gallpar partmp = params;
			auto gsys = t_gsys::sat2gsys(sat);
			t_gsatdata obsdata(_log,_site, sat, epoch);
				
			vector<GOBS> obs_list;
			t_gobsgnss obsdata_simu(_log,_site, sat, epoch);
			for (FREQ_SEQ f = FREQ_1; f <= FREQ_3; f = (FREQ_SEQ)((int)f + 1))
			{
				if (_band_index[gsys].find(f) == _band_index[gsys].end()) continue;
				if (f == FREQ_3 && triple_rm_sat.find(sat) != triple_rm_sat.end()) continue;
				auto b = _band_index[gsys][f];
				obsdata.addobs(tba2gobs(GOBSTYPE::TYPE_C, b, char2gobsattr(_range_order_attr[gsys][b][0])), FAKE_OBS_VALUE);
				obsdata.addobs(tba2gobs(GOBSTYPE::TYPE_L, b, char2gobsattr(_phase_order_attr[gsys][b][0])), FAKE_OBS_VALUE / obsdata.wavelength(b));
				obs_list.push_back(tba2gobs(GOBSTYPE::TYPE_C, b, char2gobsattr(_range_order_attr[gsys][b][0])));
				obs_list.push_back(tba2gobs(GOBSTYPE::TYPE_L, b, char2gobsattr(_phase_order_attr[gsys][b][0])));
			}

			//caculate sat_crd sat_clk
			if (!gmodel->_update_obs_info_GPP(epoch, _gall_nav, _gallobj, obsdata, partmp))
			{
				SPDLOG_LOGGER_INFO(_log, "simugns", "update_obs_info fail for " + sat + " at" + epoch.str_hms());
				continue;
			}

			if (!gmodel->_prepare_obs_GPP(epoch, _gall_nav, _gallobj, partmp))
			{
				SPDLOG_LOGGER_INFO(_log, "simugns", "prepare_obs fail for " + sat + " at" + epoch.str_hms());
				continue;
			}

			// judge whether site below ground level
			t_gtriple xyz_rh = gmodel->_trs_sat_crd - gmodel->_trs_rec_crd;
			t_gtriple ell(0, 0, 0), neu_sa(0, 0, 0);
			xyz2ell(gmodel->_trs_rec_crd, ell, false);
			xyz2neu(ell, xyz_rh, neu_sa);
			if (neu_sa[2] < 0.0) {
				SPDLOG_LOGGER_INFO(_log, "simugns", "below ground level for " + sat + " at" + epoch.str_hms());
				continue;
			}

			/* simu phase noise */
			vector<double> noiseLs(3, 0);
			t_gsatdata tmp_obs = gmodel->get_crt_obs();
			double elev = tmp_obs.ele() * R2D;
			for (FREQ_SEQ f = FREQ_1; f <= FREQ_3; f = (FREQ_SEQ)((int)f + 1))
			{
				if (_band_index[gsys].find(f) == _band_index[gsys].end()) continue;

				/*quadratic functions for multipath effect*/
				double  sigma = (1.94555e-06 * elev * elev - 0.0002135 * elev + 0.0096896)/0.005 * _sigPhaseSimu[gsys];
				/* Neglect the multipath effect when the elevation angle > 30 */
				if (elev > 30) sigma = _sigPhaseSimu[gsys];

				noiseLs[(int)f - 1] = double_eq(sigma, 0.0) ? 0.0 : std::normal_distribution<double>(0.0, sigma)(_engine);

				/* reset noise when the random noise exceeds the set maximum error */
				while (fabs(noiseLs[(int)f - 1]) > _maxPhaseSimu[gsys])
				{
					noiseLs[(int)f - 1] = double_eq(sigma, 0.0) ? 0.0 : std::normal_distribution<double>(0.0, sigma)(_engine);
				}
			}

			/* simu range noise */
			double sigmaP = (_sigCodeSimu[gsys] / sin(tmp_obs.ele())) / 2;

			/* Neglect the multipath effect when the elevation angle > 30 */
			if (elev > 30) sigmaP = _sigCodeSimu[gsys];
			double noiseP = double_eq(sigmaP, 0.0) ? 0.0 : std::normal_distribution<double>(0.0, sigmaP)(_engine);

			/* reset noise when the random noise exceeds the set maximum error */
			while (fabs(noiseP) > _maxCodeSimu[gsys])
			{
				noiseP = double_eq(sigmaP, 0.0) ? 0.0 : std::normal_distribution<double>(0.0, sigmaP)(_engine);
			}
			//set site crd
			int idx_x = partmp.getParam(_site, par_type::CRD_X, ""), idx_y = partmp.getParam(_site, par_type::CRD_Y, ""), idx_z = partmp.getParam(_site, par_type::CRD_Z, "");
			t_gtriple crd(partmp[idx_x].value(), partmp[idx_y].value(), partmp[idx_z].value());

			//caculate SION at reference freq
			double ion_L1 = 0.0, ion_rms = 0.0;
			if (_ion)
			{
				t_giono_tecgrid iono_tec;
				iono_tec._band_index = _band_index;
				iono_tec.getIonoDelay(dynamic_cast<t_gionex*>((*_data)[t_gdata::IONEX]), tmp_obs, epoch, crd, ion_L1, ion_rms);
			}
			_all_simupars[epoch]["ION_L1"] = ion_L1;

			//insert amb interger for each continue arc
			double sigmaAMB = dynamic_cast<t_gsetsimu*>(_set)->sig_simu_amb();
			double upd_band1 = 0.0;
			double upd_band2 = 0.0;
			if (_sat_amb_map.find(sat) == _sat_amb_map.end())
			{
				SPDLOG_LOGGER_INFO(_log, "simugns", "ERROR : get amb failed " + sat);
			}
			else
			{
				if (_sat_amb_map[sat].status == 0) 
				{
					_sat_amb_map[sat].status = 2;
					_sat_amb_map[sat].value_L1 = int(std::normal_distribution<double>(0.0, sigmaAMB)(_engine));
					_sat_amb_map[sat].value_L2 = int(std::normal_distribution<double>(0.0, sigmaAMB)(_engine));
					_sat_amb_map[sat].value_L3 = int(std::normal_distribution<double>(0.0, sigmaAMB)(_engine));
					SPDLOG_LOGGER_INFO(_log, "simugns_amb INFO : simu amb L1 {}  value {}", sat ,_sat_amb_map[sat].value_L1);
					SPDLOG_LOGGER_INFO(_log, "simugns_amb INFO : simu amb L2 {}  value {}", sat, _sat_amb_map[sat].value_L2);
					SPDLOG_LOGGER_INFO(_log, "simugns_amb INFO : simu amb L3 {}  value {}", sat, _sat_amb_map[sat].value_L3);
				}
				if (_sat_amb_map[sat].status == 1) 
				{
					_sat_amb_map[sat].status = 2;
				}
			}
			if (_upd && !_get_upd_value(epoch, upd_band1, upd_band2, obsdata))
			{
				SPDLOG_LOGGER_INFO(_log, "simugns", "WARNING : get upd failed " + sat);
			}
			// obs type loop
			for (auto obs : obs_list)
			{
				t_gobs this_obs = t_gobs(obs);
				double omc = 0.0, wgt = 0.0;
				vector<pair<int, double> > coef;

				// PPP inverse operation to simulate obs
				if (!gmodel->_omc_obs_ALL(epoch, tmp_obs, partmp, this_obs, omc))
				{
					SPDLOG_LOGGER_INFO(_log, "simugns", "ERROR : omc obs failed " + sat);
					continue;
				}

				// simulate ion correct
				double gamma = SQR(obsdata.wavelength(this_obs.band())) / SQR(CLIGHT / obsdata.frequency(_band_index[gsys][FREQ_1]));
				double ion = gamma * ion_L1;

				double obs_value_simu = 0.0;
				double amb_val = 0.0;
				if (this_obs.is_code())
				{
					// set pseudorange noise
					double noise = 0.0;
					noise = noiseP;
					noise += ion;
					obs_value_simu = (FAKE_OBS_VALUE - omc + noise);
					if (_freq_index[gsys][this_obs.band()] == FREQ_1) _all_simupars[epoch]["noise_P1"] = noise;
					if (_freq_index[gsys][this_obs.band()] == FREQ_2) _all_simupars[epoch]["noise_P2"] = noise;
					if (_freq_index[gsys][this_obs.band()] == FREQ_3) _all_simupars[epoch]["noise_P3"] = noise;						
				}
				else if (this_obs.is_phase())
				{
					// set phase noise
					double noise = 0.0;
					par_type amb_type;

					if (_freq_index[gsys][this_obs.band()] == FREQ_1)
					{
						noise = noiseLs[0];
						amb_type = par_type::AMB_L1;
						_all_simupars[epoch]["noise_L1"] = noise;
						noise += noise_bias[sat];
						noise -= ion;
						amb_val = _sat_amb_map[sat].value_L1 * obsdata.wavelength(this_obs.band());
						amb_val -= upd_band1;
					}
					else if (_freq_index[gsys][this_obs.band()] == FREQ_2)
					{
						noise = noiseLs[1];
						amb_type = par_type::AMB_L2;
						_all_simupars[epoch]["noise_L2"] = noise;
						noise += noise_bias[sat];
						noise -= ion;
						amb_val = _sat_amb_map[sat].value_L2 * obsdata.wavelength(this_obs.band());
						amb_val -= upd_band2;
					}
					else if (_freq_index[gsys][this_obs.band()] == FREQ_3)
					{
						noise = noiseLs[2];
						amb_type = par_type::AMB_L3;
						_all_simupars[epoch]["noise_L3"] = noise;
						noise += noise_bias[sat];
						noise -= ion;
						amb_val = _sat_amb_map[sat].value_L3 * obsdata.wavelength(this_obs.band());
					}
					obs_value_simu = (FAKE_OBS_VALUE - omc + noise + amb_val) / obsdata.wavelength(this_obs.band());
					
				}
				else continue;
				obsdata_simu.addobs(obs, obs_value_simu);
			}

			_gobs_simu->addobs(make_shared<t_gobsgnss>(obsdata_simu));
		}
		for (auto sat :_sats)
		{
			if (_sat_amb_map[sat].status == 1) 
			{
				_sat_amb_map[sat].status = 0;
				_sat_amb_map[sat].value_L1 = 0;
				_sat_amb_map[sat].value_L2 = 0;
				_sat_amb_map[sat].value_L3 = 0;
			}
			if (_sat_amb_map[sat].status == 2) _sat_amb_map[sat].status = 1;
		}
		return true;
	}

	bool t_gsimugns::_simuRnxhdr()
	{
		t_rnxhdr::t_obstypes obstypes;
		for (auto gsys : _systems)
		{
			GSYS SYS = t_gsys::str2gsys(gsys);
			if (_band_index.find(SYS) == _band_index.end()) continue;
			for (FREQ_SEQ f = FREQ_1; f <= FREQ_3; f = (FREQ_SEQ)((int)f + 1))
			{
				if (_band_index[SYS].find(f) != _band_index[SYS].end())
				{
					auto b = _band_index[SYS][f];
					char csys[2]; sprintf(csys, "%c", t_gsys::gsys2char(SYS));
					obstypes[string(csys)].push_back(make_pair(tba2gobs(GOBSTYPE::TYPE_C, b, char2gobsattr(_range_order_attr[SYS][b][0])), 1));
					obstypes[string(csys)].push_back(make_pair(tba2gobs(GOBSTYPE::TYPE_L, b, char2gobsattr(_phase_order_attr[SYS][b][0])), 1));
				}
			}
		}

		auto  gobj = dynamic_cast<t_gallobj*>((*_data)[t_gdata::ID_TYPE::ALLOBJ]);
		shared_ptr<t_grec> rec_obj = dynamic_pointer_cast<t_grec> (gobj->obj(_site));
		t_grec::t_maphdr maphdr = rec_obj->gethdr();
		t_rnxhdr rnxhdr = maphdr.begin()->second;
		rnxhdr.rnxver("3.03");
		rnxhdr.interval(dynamic_cast<t_gsetgen*>(_set)->sampling());
		rnxhdr.mapobs(obstypes);
		rnxhdr.first(*_gobs_simu->epochs(_site).begin());
		rnxhdr.last(*_gobs_simu->epochs(_site).rbegin());

		rec_obj->changehdr(rnxhdr, dynamic_cast<t_gsetgen*>(_set)->beg(), "");

		return true;
	}

	bool t_gsimugns::_writeSimuObs()
	{
		string path = dynamic_cast<t_gsetout*>(_set)->outputs("rinexo");
		transform(_site.begin(), _site.end(), _site.begin(), ::tolower);
		substitute(path, "$(rec)", _site, false);
		if (path.empty()) path = _site + ".o";
		transform(_site.begin(), _site.end(), _site.begin(), ::toupper);

		t_gio* gout = new t_gfile(_log); gout->spdlog(_log); gout->path(path);
		t_gcoder* gencode = new t_rinexo_leo(_set, "", 4096);
		gencode->clear(); gencode->spdlog(_log); gencode->path(path);

		t_gallobj* gobj = new t_gallobj(_log); gobj->spdlog(_log);
		gobj->add(dynamic_cast<t_gallobj*>((*_data)[t_gdata::ID_TYPE::ALLOBJ])->obj(_site));

		gencode->add_data("ID0", _gobs_simu);
		gencode->add_data("OBJ", gobj);
		gout->coder(gencode); 
		gout->run_write();

		delete gout; gout = nullptr;
		delete gobj; gobj = nullptr;
		delete gencode; gencode = nullptr;

		return true;
	}

	void t_gsimugns::_setout()
	{
		// set par file
		string par_tmp = dynamic_cast<t_gsetout*>(_set)->outputs("par");
		if (!par_tmp.empty()) {
			substitute(par_tmp, "$(rec)", _site, false);
			if (_parfile) {
				if (_parfile->is_open()) _parfile->close();
				delete _parfile;
				_parfile = nullptr;
			}
			_parfile = new t_giof;
			_parfile->tsys(t_gtime::GPS);
			_parfile->mask(par_tmp);
			_parfile->append(dynamic_cast<t_gsetout*>(_set)->append());
		}

	}

	void t_gsimugns::_write_pars()
	{
		if (!_parfile) return;

		ostringstream os;

		// output formaat
		os << "%EPO " << fixed << setprecision(4)
			<< setw(24) << "DATE[Y-M-D h:m:s]"
			<< setw(17) << "Sec of Week[s]"
			<< setw(16) << "Sec of Day[s]"
			// write data
			<< setw(15) << "CRD_CRDX[m]"
			<< setw(15) << "CRD_CRDY[m]"
			<< setw(15) << "CRD_CRDZ[m]"
			<< setw(15) << "REC_CLK[m]"
			<< setw(15) << "TROP[m]"
			<< setw(15) << "ION_L1[m]"
			<< setw(15) << "GAL_ISB[m]"
			<< setw(15) << "BDS_ISB[m]"
			<< setw(15) << "G15CLK[m]"
			<< setw(15) << "E15CLK[m]"
			<< setw(15) << "C15CLK[m]"
			<< setw(15) << "noise_P1[m]"
			<< setw(15) << "noise_P2[m]"
			<< setw(15) << "noise_P3[m]"
			<< setw(15) << "noise_L1[m]"
			<< setw(15) << "noise_L2[m]"
			<< setw(15) << "noise_L3[m]"
			<< endl;

		map<t_gtime, map<string, double>>::iterator itpar;
		for (itpar = _all_simupars.begin(); itpar != _all_simupars.end(); itpar++)
		{
			t_gtime epo = itpar->first;
			string str_dsec = dbl2str(epo.dsec());

			//write header
			os << "=EPO " << fixed << setprecision(4)
				<< epo.str_ymdhms("", false, true) << str_dsec.substr(2)
				<< setw(17) << epo.sow() + epo.dsec()
				<< setw(16) << epo.sod() + epo.dsec();

			// write data
			os << setprecision(4) << setw(15) << itpar->second["CRD_X"]
				<< setw(15) << itpar->second["CRD_Y"]
				<< setw(15) << itpar->second["CRD_Z"]
				<< setw(15) << (itpar->second.find("CLK") != itpar->second.end() ? itpar->second["CLK"] : 0.0)
				<< setw(15) << (itpar->second.find("TRP") != itpar->second.end() ? itpar->second["TRP"] : 0.0)
				<< setw(15) << (itpar->second.find("ION_L1") != itpar->second.end() ? itpar->second["ION_L1"] : 0.0)
				<< setw(15) << (itpar->second.find("GAL_ISB") != itpar->second.end() ? itpar->second["GAL_ISB"] : 0.0)
				<< setw(15) << (itpar->second.find("BDS_ISB") != itpar->second.end() ? itpar->second["BDS_ISB"] : 0.0)
				<< setw(15) << (itpar->second.find("G15CLK") != itpar->second.end() ? itpar->second["G15CLK"] : 0.0)
				<< setw(15) << (itpar->second.find("E15CLK") != itpar->second.end() ? itpar->second["E15CLK"] : 0.0)
				<< setw(15) << (itpar->second.find("C15CLK") != itpar->second.end() ? itpar->second["C15CLK"] : 0.0)
				<< setw(15) << (itpar->second.find("noise_P1") != itpar->second.end() ? itpar->second["noise_P1"] : 0.0)
				<< setw(15) << (itpar->second.find("noise_P2") != itpar->second.end() ? itpar->second["noise_P2"] : 0.0)
				<< setw(15) << (itpar->second.find("noise_P3") != itpar->second.end() ? itpar->second["noise_P3"] : 0.0)
				<< setw(15) << (itpar->second.find("noise_L1") != itpar->second.end() ? itpar->second["noise_L1"] : 0.0)
				<< setw(15) << (itpar->second.find("noise_L2") != itpar->second.end() ? itpar->second["noise_L2"] : 0.0)
				<< setw(15) << (itpar->second.find("noise_L3") != itpar->second.end() ? itpar->second["noise_L3"] : 0.0) << endl;
		}

		_parfile->write(os.str().c_str(), os.str().size());
		_parfile->flush();
	}

	double interp(t_gtime epoch, map<t_gtime, double> items, int sample)
	{
		if (items.find(epoch) != items.end())
		{
			return items[epoch];
		}
		else if (epoch < items.begin()->first)
		{
			return items.begin()->second;
		}
		else if (epoch > (items.rbegin())->first)
		{
			return (--items.end())->second;
		}

		int down = floor(epoch.sod() / sample);
		int up = down + 1;
		t_gtime downT, upT;
		downT.from_mjd(epoch.mjd(), down * sample);
		upT.from_mjd(epoch.mjd(), up * sample);

		if (items.find(downT) == items.end())
		{
			while (1)
			{
				downT = downT - sample;
				down -= 1;
				if (items.find(downT) != items.end()) break;
			}
		}
		if (items.find(upT) == items.end())
		{
			while (1)
			{
				upT = upT + sample;
				up += 1;
				if (items.find(upT) != items.end()) break;
			}
		}

		int t0, t1;
		double val0, val1;
		t0 = epoch.sod() - down * sample;
		t1 = up * sample - epoch.sod();
		val0 = items[downT]; val1 = items[upT];
		return (val1 * t0 + val0 * t1) / (t0 + t1);
	}

	bool t_gsimugns::_get_upd_value(t_gtime epoch, double& upd1, double& upd2, t_gsatdata satdata)
	{
		one_epoch_upd tmp_wl = _gupd->get_epo_upd(UPDTYPE::WL, t_gtime(WL_IDENTIFY));
		if (tmp_wl.find(satdata.sat())== tmp_wl.end())
			return false;
		one_epoch_upd tmp_nl = _gupd->get_epo_upd(UPDTYPE::NL, epoch);
		if (tmp_nl.find(satdata.sat()) == tmp_nl.end())
			return false;
		GSYS sys = t_gsys::str2gsys(satdata.sys());
		double f1 = satdata.frequency(_band_index[sys][FREQ_1]);
		double f2 = satdata.frequency(_band_index[sys][FREQ_2]);
		upd1 = tmp_nl[satdata.sat()]->value + (f2 / (f1 - f2)) * tmp_wl[satdata.sat()]->value;
		upd1 = upd1 * satdata.wavelength(_band_index[sys][FREQ_1]);
		upd2 = tmp_nl[satdata.sat()]->value + (f1 / (f1 - f2)) * tmp_wl[satdata.sat()]->value;
		upd2 = upd2 * satdata.wavelength(_band_index[sys][FREQ_2]);
		return true;
	}
}




