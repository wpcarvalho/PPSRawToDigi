#include <string>

bool updated_ntuples = true;

double timestamp0 = 1342044000;

//----------------------------------------------------------------------------------------------------
// alignment corrections
//	a: xy coupling in rad
//	b: x shift in mm
//	c: y shift in mm

double al_a_L_F = -0.02E-3, al_b_L_F = -20E-3, al_c_L_F = 231E-3;
double al_a_L_N = -1.57E-3, al_b_L_N = -24E-3, al_c_L_N = 160E-3;
double al_a_R_N = +1.71E-3, al_b_R_N = -24E-3, al_c_R_N =  26E-3;
double al_a_R_F = -0.57E-3, al_b_R_F = -21E-3, al_c_R_F =  50E-3;

bool applyTimeDependentAlignment = true;

//----------------------------------------------------------------------------------------------------
// optics and corrections

double p = 4000.;

double v_y_N = +0.02037;		// 1
double v_y_F = -0.00005;		// 1
double L_y_N = 237.7E3;			// mm
double L_y_F = 263.2E3;			// mm
double dLds_y = 4.743;			// 1
double dvds_y = -0.003801E-3;	// mm^-1

double v_x_N = -2.164;			// 1
double v_x_F = -1.866;			// 1
double L_x_N = 2.879E3;			// mm, nominal only!
double L_x_F = -0.000E3;		// mm, nominal only!
double dLds_x = -0.5359;		// 1
double dvds_x = 0.05553E-3;		// mm^-1

double v_x_rat = v_x_N / v_x_F;

double th_x_al_R = 0.;
double th_x_be_R = 0.;
double th_x_al_L = 0.;
double th_x_be_L = 0.;


//----------------------------------------------------------------------------------------------------
// selection parameters

double cut1_si, cut1_c;
double cut2_si;
double cut34_si;
double cut5_a, cut5_b, cut5_si;
double cut6_a, cut6_b, cut6_si;
double cut7_al, cut7_c, cut7_si;

//----------------------------------------------------------------------------------------------------
// binning

double t_min = 0., t_max = 1.0;
double t_min_full = 0., t_max_full = 1.2;

unsigned int N_bins = 100, N_low = 0, N_high = 50;

double t_min_fit = 0.009;


//----------------------------------------------------------------------------------------------------
// analysis parameters

double th_y_lcut_L = 0E-6, th_y_lcut_R = 0E-6;
double th_y_hcut_L = 0E-6, th_y_hcut_R = 0E-6;

double si_th_y_os = 2.16E-6, si_th_y_os_unc = 0.0E-6; // one side smearing
double si_th_x_os = 5.99E-6, si_th_x_os_unc = 0.0E-6;

double alpha_gr = 1E0;

double full_norm_corr = 0.;
double L_int_eff = 0.;	// mb^-1

double alignment_t0 = 67300;	// beginning of the first time-slice
double alignment_ts = 10*60.;	// time-slice in s
double alignment_y_min = 9.;

double eff_th_y_min = 35E-6;

std::string rp_L_N, rp_L_F, rp_R_N, rp_R_F;



//----------------------------------------------------------------------------------------------------

void Init_45b_56t()
{
	rp_L_N = "track_rp_21"; rp_L_F = "track_rp_25";
	rp_R_N = "track_rp_120"; rp_R_F = "track_rp_124";

	cut1_si = 9E-6; cut1_c = 0E-6;
	cut2_si = 3.5E-6;

	cut34_si = 0.170;

	cut5_a = 0.10735; cut5_b = -0.0E-3; cut5_si = 19E-3;
	cut6_a = 0.10705; cut6_b = +3.5E-3; cut6_si = 19E-3;

	cut7_al = -13.9;
	cut7_c = 0E-3;
	cut7_si = 9E-3;

	th_y_lcut_L = 32.5E-6; th_y_lcut_R = 30.9E-6;	
	th_y_hcut_L = 104E-6; th_y_hcut_R = 103E-6;

	full_norm_corr = 1.136;
	L_int_eff = 851.7E3;
}

//----------------------------------------------------------------------------------------------------

void Init_45t_56b()
{
	rp_L_N = "track_rp_20"; rp_L_F = "track_rp_24";
	rp_R_N = "track_rp_121"; rp_R_F = "track_rp_125";

	cut1_si = 9E-6; cut1_c = 0E-6;
	cut2_si = 3.5E-6;

	cut34_si = 0.200;

	cut5_a = 0.10735; cut5_b = -2.0E-3; cut5_si = 19E-3;
	cut6_a = 0.10705; cut6_b = +0.0E-3; cut6_si = 19E-3;

	cut7_al = -16.0;
	cut7_c = 0E-3;
	cut7_si = 9E-3;

	th_y_lcut_L = 30.6E-6; th_y_lcut_R = 31.4E-6;
	th_y_hcut_L = 106E-6; th_y_hcut_R = 106E-6;

	full_norm_corr = 1.135;
	L_int_eff = 851.7E3;
}
