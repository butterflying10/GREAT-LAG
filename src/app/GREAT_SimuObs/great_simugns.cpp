/* ----------------------------------------------------------------------
  (c) 2011-2018 Geodetic Observatory Pecny, http://www.pecny.cz (gnss@pecny.cz)
	  Research Institute of Geodesy, Topography and Cartography
	  Ondrejov 244, 251 65, Czech Republic

  Purpose: gnss data simulation
  Version: $ Rev: $

  2021-05-16 /created by xjhan

-*/

#include <random>
#include <memory>
#include <signal.h>
#include "gcfg_simugns.h"
#include "gproc/gsimugns.h"

using namespace std;
using namespace gnut;
using namespace std::chrono;
double generateRandom(double minv, double maxv);
void catch_signal(int) { cout << "Program interrupted by Ctrl-C [SIGINT,2]\n"; }

// MAIN
// ----------
int main(int argc, char** argv)
{
	// Only to cout the Reminder here
	signal(SIGINT, catch_signal);

	// Construct the gset class and init some values in the class
	t_gcfg_simugns gset;
	gset.app("GREAT/SIMUGNS", "0.9.0", "$Rev: 2448 $", "(gnss@pecny.cz)", __DATE__, __TIME__);
	// Get the arguments from the command line
	gset.arg(argc, argv, true, false);

	// Creat and set the log file : ppp.log
	auto log_type = dynamic_cast<t_gsetout*>(&gset)->log_type();
	auto log_level = dynamic_cast<t_gsetout*>(&gset)->log_level();
	auto log_name = dynamic_cast<t_gsetout*>(&gset)->log_name();
	auto log_pattern = dynamic_cast<t_gsetout*>(&gset)->log_pattern();
	spdlog::set_level(log_level);
	spdlog::set_pattern(log_pattern);
	spdlog::flush_on(spdlog::level::debug);// change err to info
	t_grtlog great_log = t_grtlog(log_type, log_level, log_name);
	auto my_logger = great_log.spdlog();

	// Prepare site list from gset
	set<string> sites = dynamic_cast<t_gsetgen*>(&gset)->recs();
	// Prepare input files list form gset
	multimap<IFMT, string> inp = gset.inputs_all();

	// DECLARATIONS/INITIALIZATIONS
	// gobs for the obs data
	t_gallobs* gobs = new t_gallobs();  gobs->spdlog(my_logger); gobs->gset(&gset);
	// gallnav for all the navigation data, gorb = gorbit data
	t_gallprec* gorb = new t_gallprec(); gorb->spdlog(my_logger);
	t_gallprec* gorb_simu = new t_gallprec(); gorb_simu->spdlog(my_logger);

	t_gdata* gdata = nullptr;
	// gpcv for the atx data read from the atx file
	t_gallpcv* gpcv = nullptr; if (gset.input_size("atx") > 0) { gpcv = new t_gallpcv;  gpcv->spdlog(my_logger); }
	// gotl for the blq data read from the blq file, which will be used for Ocean tidal corrections
	t_gallotl* gotl = nullptr; if (gset.input_size("blq") > 0) { gotl = new t_gallotl;  gotl->spdlog(my_logger); }
	// gbia for the DCB data read from the biasinex and bianern files
	t_gallbias* gbia = nullptr; if (gset.input_size("biasinex") > 0 ||
		gset.input_size("bias") > 0) {
		gbia = new t_gallbias; gbia->spdlog(my_logger);
	}
	t_gupd* gupd = nullptr; if (gset.input_size("upd") > 0) { gupd = new t_gupd;  gupd->spdlog(my_logger); }
	// gobj for the gpcv and gotl, which means the model can be used by all satellites and stations
	t_gallobj* gobj = new t_gallobj(my_logger,gpcv, gotl); gobj->spdlog(my_logger);
	
	t_gnavde* gde = new t_gnavde;
	t_gpoleut1* gerp = new t_gpoleut1;
	t_gionex* gionex = new t_gionex;
	t_galltrp* gtrp = nullptr;
	bool isZtd = gset.zwd();
	if (isZtd)
	{
		gtrp = new t_galltrp;
	}
	
	t_gtime runepoch(t_gtime::GPS);
	t_gtime lstepoch(t_gtime::GPS);

	// SET OBJECTS
	set<string>::const_iterator itOBJ;
	set<string> obj = dynamic_cast<t_gsetrec*>(&gset)->objects();
	for (itOBJ = obj.begin(); itOBJ != obj.end(); ++itOBJ) {
		string name = *itOBJ;
		shared_ptr<t_grec> rec = dynamic_cast<t_gsetrec*>(&gset)->grec(name, my_logger);
		gobj->add(rec);
	}

	gorb->use_clksp3(true);
	gorb->use_clkrnx(true);
	gorb_simu->use_clksp3(true);
	gorb_simu->use_clkrnx(true);

	// Multi gcoder for multi-thread decoding data
	vector<t_gcoder*> gcoder;
	// Multi gior for multi-thread receiving data
	vector<t_gio*> gio;

	t_gio* tgio = 0;
	t_gcoder* tgcoder = 0;

	// DATA READING
	multimap<IFMT, string>::const_iterator itINP = inp.begin();
	for (size_t i = 0; i < inp.size() && itINP != inp.end(); ++i, ++itINP)
	{
		// Get the file format/path, which will be used in decoder
		IFMT   ifmt(itINP->first);
		string path(itINP->second);
		string id("ID" + int2str(i));

		// For different file format, we prepare different data container and decoder for them.
		if (ifmt == IFMT::RINEXO_INP) { gdata = gobs; tgcoder = new t_rinexo(&gset, "", 4096); }
		else if (ifmt == IFMT::SP3_INP) { gdata = gorb; tgcoder = new t_sp3(&gset, "", 8172); }
		else if (ifmt == IFMT::RINEXC_INP) { gdata = gorb; tgcoder = new t_rinexc(&gset, "", 4096); }
		else if (ifmt == IFMT::SP3SIMU_INP) { gdata = gorb_simu; tgcoder = new t_sp3(&gset, "", 8172); }
		else if (ifmt == IFMT::RINEXCSIMU_INP) {gdata = gorb_simu; tgcoder = new t_rinexc(&gset, "", 4096);}
		else if (ifmt == IFMT::DE_INP) { gdata = gde; tgcoder = new t_dvpteph405(&gset, "", 4096); }
		else if (ifmt == IFMT::RINEXN_INP) { gdata = gorb; tgcoder = new t_rinexn(&gset, "", 4096); }
		else if (ifmt == IFMT::BLQ_INP) { gdata = gotl; tgcoder = new t_blq(&gset, "", 4096); }
		else if (ifmt == IFMT::EOP_INP) { gdata = gerp; tgcoder = new t_poleut1(&gset, "", 4096); }
		else if (ifmt == IFMT::ATX_INP) { gdata = gpcv; tgcoder = new t_atx(&gset, "", 4096); }
		else if (ifmt == IFMT::UPD_INP) { gdata = gupd; tgcoder = new t_upd(&gset, "", 4096); }
		else if (ifmt == IFMT::BIASINEX_INP) { gdata = gbia; tgcoder = new t_biasinex(&gset, "", 20480); }
		else if (ifmt == IFMT::BIAS_INP) { gdata = gbia; tgcoder = new t_biabernese(&gset, "", 20480); }
		else if (ifmt == IFMT::TRPZTD_INP) { gdata = gtrp; tgcoder = new t_trp(&gset, "", 20480); }
		else if (ifmt == IFMT::IONEX_INP) { gdata = gionex; tgcoder = new t_ionex(&gset, "", 20480); }
		else 
		{
			SPDLOG_LOGGER_INFO(my_logger, "Error: unrecognized format " + int2str(int(ifmt)));
			gdata = 0;
		}

		// Check the file path
		if (path.substr(0, 7) == "file://") {
			SPDLOG_LOGGER_INFO(my_logger, "path is file!");
		}

		// READ DATA FROM FILE
		if (tgcoder) {
			// prepare gio for the file
			tgio = new t_gfile(my_logger);
			tgio->spdlog(my_logger);
			tgio->path(path);

			// Put the file into gcoder
			tgcoder->clear();
			tgcoder->path(path);
			tgcoder->spdlog(my_logger);
			// Put the data container into gcoder
			tgcoder->add_data(id, gdata);
			tgcoder->add_data("OBJ", gobj);
			// Put the gcoder into the gio
			// Note, gcoder contain the gdata and gio contain the gcoder
			tgio->coder(tgcoder);

			runepoch = t_gtime::current_time(t_gtime::GPS);
			// Read the data from file here
			tgio->run_read();
			lstepoch = t_gtime::current_time(t_gtime::GPS);
			// Write the information of reading process to log file
			SPDLOG_LOGGER_INFO(my_logger, "main" , "READ: " , path , " time: " , dbl2str(lstepoch.diff(runepoch)) + " sec");
			// Delete 
			delete tgio; tgio = NULL;
			delete tgcoder; tgcoder = NULL;
		}
	}

	// set antennas for satllites (must be before PCV assigning)
	t_gtime beg = dynamic_cast<t_gsetgen*>(&gset)->beg();
	t_gtime end = dynamic_cast<t_gsetgen*>(&gset)->end();
    gobj->read_satinfo(beg);

	// assigning PCV pointers to objects
	gobj->sync_pcvs();
	if (gbia == NULL) gbia = new t_gallbias();

	// add all data  
	t_gallproc* data_simu = new t_gallproc(my_logger);
	if (!data_simu->Add_Data(gorb_simu)) { std::cout << "Warning: No orb simu file." << std::endl; };
	if (!data_simu->Add_Data(gobj)) { std::cout << "Warning: No obj info." << std::endl; };
	if (!data_simu->Add_Data(gotl)) { std::cout << "Warning: No OTL file." << std::endl; };
	if (!data_simu->Add_Data(gde)) { std::cout << "Warning: No DE file." << std::endl; };
	if (!data_simu->Add_Data(gerp)) { std::cout << "Warning: No POLEUT1 file." << std::endl; };
	if (!data_simu->Add_Data(gtrp)) { std::cout << "Warning: No TRP file." << std::endl; };
	if (!data_simu->Add_Data(gobs)) { std::cout << "Warning: No obs file." << std::endl; };
	if (!data_simu->Add_Data(gbia)) { std::cout << "Warning: No bias file." << std::endl; };
	if (!data_simu->Add_Data(gionex)) { std::cout << "Warning: No ionex file." << std::endl; };
	if (!data_simu->Add_Data(gupd)) { std::cout << "Warning: No upd file." << std::endl; };

	vector<t_gsimugns*> vgsimu;
	for (auto site : sites)
	{
		int random = 0; //控制随机数,不同测站仿真时模糊度值不一样,同一测站不同卫星的模糊度前后两次仿真是一样的

		for (auto c : site) random += c - 'A';

		vgsimu.push_back(0); int idx = vgsimu.size() - 1;
		vgsimu[idx] = new t_gsimugns(site, &gset, data_simu, my_logger, random);

		SPDLOG_LOGGER_INFO(my_logger, "main", site + " simulation started [" + int2str(idx) + "]");
		SPDLOG_LOGGER_INFO(my_logger, "main", site + beg.str_ymdhms("  beg: ") + end.str_ymdhms("  end: "));

		runepoch = t_gtime::current_time(t_gtime::GPS);

		// clear _prec before next site
		dynamic_cast<t_gallprec*>(gorb_simu)->reset_prec(nullptr);

		if (!vgsimu[idx]->processBatch())
		{
			SPDLOG_LOGGER_INFO(my_logger, "main", "Error: simulation fail " + site);
			continue;
		}

		lstepoch = t_gtime::current_time(t_gtime::GPS);

		SPDLOG_LOGGER_INFO(my_logger, "main", site + " simulation finished : duration  "
			+ dbl2str(lstepoch.diff(runepoch)) + " sec");

	}

	for (unsigned int i = 0; i < gio.size(); ++i) { delete gio[i]; }; gio.clear();
	for (unsigned int i = 0; i < gcoder.size(); ++i) { delete gcoder[i]; }; gcoder.clear();
	for (unsigned int i = 0; i < vgsimu.size(); ++i) { if (vgsimu[i]) { delete vgsimu[i]; vgsimu[i] = nullptr; } }

	if (gobs) delete gobs; gobs = nullptr;
	if (gpcv) delete gpcv; gpcv = nullptr;
	if (gotl) delete gotl; gotl = nullptr;
	if (gobj) delete gobj; gobj = nullptr;
	if (gorb) delete gorb; gorb = nullptr;
	if (gbia) delete gbia; gbia = nullptr;
	if (gorb_simu) delete gorb_simu; gorb_simu = nullptr;
	if (gde) delete gde; gde = nullptr;
	if (gerp) delete gerp; gerp = nullptr;
	//if (gleap) delete gleap; gleap = nullptr;
	if (gionex) delete gionex; gionex = nullptr;
	if (data_simu) delete data_simu; data_simu = nullptr;

	return 0;
}

double generateRandom(double minv, double maxv)
{
	static std::default_random_engine engine;
	std::uniform_real_distribution<double> distribution(minv, maxv);
	return distribution(engine);
}



