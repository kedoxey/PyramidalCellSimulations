/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format off
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
// clang-format on
#include "neuron/cache/mechanism_range.hpp"
static constexpr auto number_of_datum_variables = 2;
static constexpr auto number_of_floating_point_variables = 71;
namespace {
template <typename T>
using _nrn_mechanism_std_vector = std::vector<T>;
using _nrn_model_sorted_token = neuron::model_sorted_token;
using _nrn_mechanism_cache_range = neuron::cache::MechanismRange<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_mechanism_cache_instance = neuron::cache::MechanismInstance<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_non_owning_id_without_container = neuron::container::non_owning_identifier_without_container;
template <typename T>
using _nrn_mechanism_field = neuron::mechanism::field<T>;
template <typename... Args>
void _nrn_mechanism_register_data_fields(Args&&... args) {
  neuron::mechanism::register_data_fields(std::forward<Args>(args)...);
}
}
 
#if !NRNGPU
#undef exp
#define exp hoc_Exp
#endif
 
#define nrn_init _nrn_init__Nap_Et2
#define _nrn_initial _nrn_initial__Nap_Et2
#define nrn_cur _nrn_cur__Nap_Et2
#define _nrn_current _nrn_current__Nap_Et2
#define nrn_jacob _nrn_jacob__Nap_Et2
#define nrn_state _nrn_state__Nap_Et2
#define _net_receive _net_receive__Nap_Et2 
#define rates rates__Nap_Et2 
#define states states__Nap_Et2 
 
#define _threadargscomma_ _ml, _iml, _ppvar, _thread, _globals, _nt,
#define _threadargsprotocomma_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt,
#define _internalthreadargsprotocomma_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt,
#define _threadargs_ _ml, _iml, _ppvar, _thread, _globals, _nt
#define _threadargsproto_ Memb_list* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt
#define _internalthreadargsproto_ _nrn_mechanism_cache_range* _ml, size_t _iml, Datum* _ppvar, Datum* _thread, double* _globals, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _ml->template fpfield<0>(_iml)
#define gmax_columnindex 0
#define conductance _ml->template fpfield<1>(_iml)
#define conductance_columnindex 1
#define m_instances _ml->template fpfield<2>(_iml)
#define m_instances_columnindex 2
#define m_forwardRate_rate _ml->template fpfield<3>(_iml)
#define m_forwardRate_rate_columnindex 3
#define m_forwardRate_midpoint _ml->template fpfield<4>(_iml)
#define m_forwardRate_midpoint_columnindex 4
#define m_forwardRate_scale _ml->template fpfield<5>(_iml)
#define m_forwardRate_scale_columnindex 5
#define m_reverseRate_rate _ml->template fpfield<6>(_iml)
#define m_reverseRate_rate_columnindex 6
#define m_reverseRate_midpoint _ml->template fpfield<7>(_iml)
#define m_reverseRate_midpoint_columnindex 7
#define m_reverseRate_scale _ml->template fpfield<8>(_iml)
#define m_reverseRate_scale_columnindex 8
#define m_timeCourse_TIME_SCALE _ml->template fpfield<9>(_iml)
#define m_timeCourse_TIME_SCALE_columnindex 9
#define m_timeCourse_VOLT_SCALE _ml->template fpfield<10>(_iml)
#define m_timeCourse_VOLT_SCALE_columnindex 10
#define m_steadyState_rate _ml->template fpfield<11>(_iml)
#define m_steadyState_rate_columnindex 11
#define m_steadyState_midpoint _ml->template fpfield<12>(_iml)
#define m_steadyState_midpoint_columnindex 12
#define m_steadyState_scale _ml->template fpfield<13>(_iml)
#define m_steadyState_scale_columnindex 13
#define m_q10Settings_fixedQ10 _ml->template fpfield<14>(_iml)
#define m_q10Settings_fixedQ10_columnindex 14
#define h_instances _ml->template fpfield<15>(_iml)
#define h_instances_columnindex 15
#define h_forwardRate_rate _ml->template fpfield<16>(_iml)
#define h_forwardRate_rate_columnindex 16
#define h_forwardRate_midpoint _ml->template fpfield<17>(_iml)
#define h_forwardRate_midpoint_columnindex 17
#define h_forwardRate_scale _ml->template fpfield<18>(_iml)
#define h_forwardRate_scale_columnindex 18
#define h_reverseRate_rate _ml->template fpfield<19>(_iml)
#define h_reverseRate_rate_columnindex 19
#define h_reverseRate_midpoint _ml->template fpfield<20>(_iml)
#define h_reverseRate_midpoint_columnindex 20
#define h_reverseRate_scale _ml->template fpfield<21>(_iml)
#define h_reverseRate_scale_columnindex 21
#define h_steadyState_rate _ml->template fpfield<22>(_iml)
#define h_steadyState_rate_columnindex 22
#define h_steadyState_midpoint _ml->template fpfield<23>(_iml)
#define h_steadyState_midpoint_columnindex 23
#define h_steadyState_scale _ml->template fpfield<24>(_iml)
#define h_steadyState_scale_columnindex 24
#define h_q10Settings_fixedQ10 _ml->template fpfield<25>(_iml)
#define h_q10Settings_fixedQ10_columnindex 25
#define gion _ml->template fpfield<26>(_iml)
#define gion_columnindex 26
#define m_forwardRate_x _ml->template fpfield<27>(_iml)
#define m_forwardRate_x_columnindex 27
#define m_forwardRate_r _ml->template fpfield<28>(_iml)
#define m_forwardRate_r_columnindex 28
#define m_reverseRate_x _ml->template fpfield<29>(_iml)
#define m_reverseRate_x_columnindex 29
#define m_reverseRate_r _ml->template fpfield<30>(_iml)
#define m_reverseRate_r_columnindex 30
#define m_timeCourse_V _ml->template fpfield<31>(_iml)
#define m_timeCourse_V_columnindex 31
#define m_timeCourse_ALPHA _ml->template fpfield<32>(_iml)
#define m_timeCourse_ALPHA_columnindex 32
#define m_timeCourse_BETA _ml->template fpfield<33>(_iml)
#define m_timeCourse_BETA_columnindex 33
#define m_timeCourse_t _ml->template fpfield<34>(_iml)
#define m_timeCourse_t_columnindex 34
#define m_steadyState_x _ml->template fpfield<35>(_iml)
#define m_steadyState_x_columnindex 35
#define m_q10Settings_q10 _ml->template fpfield<36>(_iml)
#define m_q10Settings_q10_columnindex 36
#define m_rateScale _ml->template fpfield<37>(_iml)
#define m_rateScale_columnindex 37
#define m_alpha _ml->template fpfield<38>(_iml)
#define m_alpha_columnindex 38
#define m_beta _ml->template fpfield<39>(_iml)
#define m_beta_columnindex 39
#define m_inf _ml->template fpfield<40>(_iml)
#define m_inf_columnindex 40
#define m_tauUnscaled _ml->template fpfield<41>(_iml)
#define m_tauUnscaled_columnindex 41
#define m_tau _ml->template fpfield<42>(_iml)
#define m_tau_columnindex 42
#define m_fcond _ml->template fpfield<43>(_iml)
#define m_fcond_columnindex 43
#define h_forwardRate_x _ml->template fpfield<44>(_iml)
#define h_forwardRate_x_columnindex 44
#define h_forwardRate_r _ml->template fpfield<45>(_iml)
#define h_forwardRate_r_columnindex 45
#define h_reverseRate_x _ml->template fpfield<46>(_iml)
#define h_reverseRate_x_columnindex 46
#define h_reverseRate_r _ml->template fpfield<47>(_iml)
#define h_reverseRate_r_columnindex 47
#define h_steadyState_x _ml->template fpfield<48>(_iml)
#define h_steadyState_x_columnindex 48
#define h_q10Settings_q10 _ml->template fpfield<49>(_iml)
#define h_q10Settings_q10_columnindex 49
#define h_rateScale _ml->template fpfield<50>(_iml)
#define h_rateScale_columnindex 50
#define h_alpha _ml->template fpfield<51>(_iml)
#define h_alpha_columnindex 51
#define h_beta _ml->template fpfield<52>(_iml)
#define h_beta_columnindex 52
#define h_fcond _ml->template fpfield<53>(_iml)
#define h_fcond_columnindex 53
#define h_inf _ml->template fpfield<54>(_iml)
#define h_inf_columnindex 54
#define h_tau _ml->template fpfield<55>(_iml)
#define h_tau_columnindex 55
#define conductanceScale _ml->template fpfield<56>(_iml)
#define conductanceScale_columnindex 56
#define fopen0 _ml->template fpfield<57>(_iml)
#define fopen0_columnindex 57
#define fopen _ml->template fpfield<58>(_iml)
#define fopen_columnindex 58
#define g _ml->template fpfield<59>(_iml)
#define g_columnindex 59
#define m_q _ml->template fpfield<60>(_iml)
#define m_q_columnindex 60
#define h_q _ml->template fpfield<61>(_iml)
#define h_q_columnindex 61
#define temperature _ml->template fpfield<62>(_iml)
#define temperature_columnindex 62
#define ena _ml->template fpfield<63>(_iml)
#define ena_columnindex 63
#define ina _ml->template fpfield<64>(_iml)
#define ina_columnindex 64
#define rate_m_q _ml->template fpfield<65>(_iml)
#define rate_m_q_columnindex 65
#define rate_h_q _ml->template fpfield<66>(_iml)
#define rate_h_q_columnindex 66
#define Dm_q _ml->template fpfield<67>(_iml)
#define Dm_q_columnindex 67
#define Dh_q _ml->template fpfield<68>(_iml)
#define Dh_q_columnindex 68
#define v _ml->template fpfield<69>(_iml)
#define v_columnindex 69
#define _g _ml->template fpfield<70>(_iml)
#define _g_columnindex 70
#define _ion_ina *(_ml->dptr_field<0>(_iml))
#define _p_ion_ina static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_dinadv *(_ml->dptr_field<1>(_iml))
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 _prop_id = _nrn_get_prop_id(_prop);
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_Nap_Et2", _hoc_setdata},
 {"rates_Nap_Et2", _hoc_rates},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_rates(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"rates", _npy_rates},
 {0, 0}
};
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"gmax_Nap_Et2", "S/cm2"},
 {"conductance_Nap_Et2", "uS"},
 {"m_forwardRate_rate_Nap_Et2", "kHz"},
 {"m_forwardRate_midpoint_Nap_Et2", "mV"},
 {"m_forwardRate_scale_Nap_Et2", "mV"},
 {"m_reverseRate_rate_Nap_Et2", "kHz"},
 {"m_reverseRate_midpoint_Nap_Et2", "mV"},
 {"m_reverseRate_scale_Nap_Et2", "mV"},
 {"m_timeCourse_TIME_SCALE_Nap_Et2", "ms"},
 {"m_timeCourse_VOLT_SCALE_Nap_Et2", "mV"},
 {"m_steadyState_midpoint_Nap_Et2", "mV"},
 {"m_steadyState_scale_Nap_Et2", "mV"},
 {"h_forwardRate_rate_Nap_Et2", "kHz"},
 {"h_forwardRate_midpoint_Nap_Et2", "mV"},
 {"h_forwardRate_scale_Nap_Et2", "mV"},
 {"h_reverseRate_rate_Nap_Et2", "kHz"},
 {"h_reverseRate_midpoint_Nap_Et2", "mV"},
 {"h_reverseRate_scale_Nap_Et2", "mV"},
 {"h_steadyState_midpoint_Nap_Et2", "mV"},
 {"h_steadyState_scale_Nap_Et2", "mV"},
 {"gion_Nap_Et2", "S/cm2"},
 {"m_forwardRate_r_Nap_Et2", "kHz"},
 {"m_reverseRate_r_Nap_Et2", "kHz"},
 {"m_timeCourse_t_Nap_Et2", "ms"},
 {"m_alpha_Nap_Et2", "kHz"},
 {"m_beta_Nap_Et2", "kHz"},
 {"m_tauUnscaled_Nap_Et2", "ms"},
 {"m_tau_Nap_Et2", "ms"},
 {"h_forwardRate_r_Nap_Et2", "kHz"},
 {"h_reverseRate_r_Nap_Et2", "kHz"},
 {"h_alpha_Nap_Et2", "kHz"},
 {"h_beta_Nap_Et2", "kHz"},
 {"h_tau_Nap_Et2", "ms"},
 {"g_Nap_Et2", "uS"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double h_q0 = 0;
 static double m_q0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void nrn_cur(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_jacob(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[2].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Nap_Et2",
 "gmax_Nap_Et2",
 "conductance_Nap_Et2",
 "m_instances_Nap_Et2",
 "m_forwardRate_rate_Nap_Et2",
 "m_forwardRate_midpoint_Nap_Et2",
 "m_forwardRate_scale_Nap_Et2",
 "m_reverseRate_rate_Nap_Et2",
 "m_reverseRate_midpoint_Nap_Et2",
 "m_reverseRate_scale_Nap_Et2",
 "m_timeCourse_TIME_SCALE_Nap_Et2",
 "m_timeCourse_VOLT_SCALE_Nap_Et2",
 "m_steadyState_rate_Nap_Et2",
 "m_steadyState_midpoint_Nap_Et2",
 "m_steadyState_scale_Nap_Et2",
 "m_q10Settings_fixedQ10_Nap_Et2",
 "h_instances_Nap_Et2",
 "h_forwardRate_rate_Nap_Et2",
 "h_forwardRate_midpoint_Nap_Et2",
 "h_forwardRate_scale_Nap_Et2",
 "h_reverseRate_rate_Nap_Et2",
 "h_reverseRate_midpoint_Nap_Et2",
 "h_reverseRate_scale_Nap_Et2",
 "h_steadyState_rate_Nap_Et2",
 "h_steadyState_midpoint_Nap_Et2",
 "h_steadyState_scale_Nap_Et2",
 "h_q10Settings_fixedQ10_Nap_Et2",
 0,
 "gion_Nap_Et2",
 "m_forwardRate_x_Nap_Et2",
 "m_forwardRate_r_Nap_Et2",
 "m_reverseRate_x_Nap_Et2",
 "m_reverseRate_r_Nap_Et2",
 "m_timeCourse_V_Nap_Et2",
 "m_timeCourse_ALPHA_Nap_Et2",
 "m_timeCourse_BETA_Nap_Et2",
 "m_timeCourse_t_Nap_Et2",
 "m_steadyState_x_Nap_Et2",
 "m_q10Settings_q10_Nap_Et2",
 "m_rateScale_Nap_Et2",
 "m_alpha_Nap_Et2",
 "m_beta_Nap_Et2",
 "m_inf_Nap_Et2",
 "m_tauUnscaled_Nap_Et2",
 "m_tau_Nap_Et2",
 "m_fcond_Nap_Et2",
 "h_forwardRate_x_Nap_Et2",
 "h_forwardRate_r_Nap_Et2",
 "h_reverseRate_x_Nap_Et2",
 "h_reverseRate_r_Nap_Et2",
 "h_steadyState_x_Nap_Et2",
 "h_q10Settings_q10_Nap_Et2",
 "h_rateScale_Nap_Et2",
 "h_alpha_Nap_Et2",
 "h_beta_Nap_Et2",
 "h_fcond_Nap_Et2",
 "h_inf_Nap_Et2",
 "h_tau_Nap_Et2",
 "conductanceScale_Nap_Et2",
 "fopen0_Nap_Et2",
 "fopen_Nap_Et2",
 "g_Nap_Et2",
 0,
 "m_q_Nap_Et2",
 "h_q_Nap_Et2",
 0,
 0};
 static Symbol* _na_sym;
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0, /* gmax */
     1e-05, /* conductance */
     3, /* m_instances */
     1.092, /* m_forwardRate_rate */
     -38, /* m_forwardRate_midpoint */
     6, /* m_forwardRate_scale */
     0.744, /* m_reverseRate_rate */
     -38, /* m_reverseRate_midpoint */
     -6, /* m_reverseRate_scale */
     1, /* m_timeCourse_TIME_SCALE */
     1, /* m_timeCourse_VOLT_SCALE */
     1, /* m_steadyState_rate */
     -52.6, /* m_steadyState_midpoint */
     4.6, /* m_steadyState_scale */
     2.95288, /* m_q10Settings_fixedQ10 */
     1, /* h_instances */
     1.33344e-05, /* h_forwardRate_rate */
     -17, /* h_forwardRate_midpoint */
     -4.63, /* h_forwardRate_scale */
     1.82522e-05, /* h_reverseRate_rate */
     -64.4, /* h_reverseRate_midpoint */
     2.63, /* h_reverseRate_scale */
     1, /* h_steadyState_rate */
     -48.8, /* h_steadyState_midpoint */
     -10, /* h_steadyState_scale */
     2.95288, /* h_q10Settings_fixedQ10 */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 71);
 	/*initialize range parameters*/
 	gmax = _parm_default[0]; /* 0 */
 	conductance = _parm_default[1]; /* 1e-05 */
 	m_instances = _parm_default[2]; /* 3 */
 	m_forwardRate_rate = _parm_default[3]; /* 1.092 */
 	m_forwardRate_midpoint = _parm_default[4]; /* -38 */
 	m_forwardRate_scale = _parm_default[5]; /* 6 */
 	m_reverseRate_rate = _parm_default[6]; /* 0.744 */
 	m_reverseRate_midpoint = _parm_default[7]; /* -38 */
 	m_reverseRate_scale = _parm_default[8]; /* -6 */
 	m_timeCourse_TIME_SCALE = _parm_default[9]; /* 1 */
 	m_timeCourse_VOLT_SCALE = _parm_default[10]; /* 1 */
 	m_steadyState_rate = _parm_default[11]; /* 1 */
 	m_steadyState_midpoint = _parm_default[12]; /* -52.6 */
 	m_steadyState_scale = _parm_default[13]; /* 4.6 */
 	m_q10Settings_fixedQ10 = _parm_default[14]; /* 2.95288 */
 	h_instances = _parm_default[15]; /* 1 */
 	h_forwardRate_rate = _parm_default[16]; /* 1.33344e-05 */
 	h_forwardRate_midpoint = _parm_default[17]; /* -17 */
 	h_forwardRate_scale = _parm_default[18]; /* -4.63 */
 	h_reverseRate_rate = _parm_default[19]; /* 1.82522e-05 */
 	h_reverseRate_midpoint = _parm_default[20]; /* -64.4 */
 	h_reverseRate_scale = _parm_default[21]; /* 2.63 */
 	h_steadyState_rate = _parm_default[22]; /* 1 */
 	h_steadyState_midpoint = _parm_default[23]; /* -48.8 */
 	h_steadyState_scale = _parm_default[24]; /* -10 */
 	h_q10Settings_fixedQ10 = _parm_default[25]; /* 2.95288 */
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 71);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ina */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _Nap_Et2_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", 1.0);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"gmax"} /* 0 */,
                                       _nrn_mechanism_field<double>{"conductance"} /* 1 */,
                                       _nrn_mechanism_field<double>{"m_instances"} /* 2 */,
                                       _nrn_mechanism_field<double>{"m_forwardRate_rate"} /* 3 */,
                                       _nrn_mechanism_field<double>{"m_forwardRate_midpoint"} /* 4 */,
                                       _nrn_mechanism_field<double>{"m_forwardRate_scale"} /* 5 */,
                                       _nrn_mechanism_field<double>{"m_reverseRate_rate"} /* 6 */,
                                       _nrn_mechanism_field<double>{"m_reverseRate_midpoint"} /* 7 */,
                                       _nrn_mechanism_field<double>{"m_reverseRate_scale"} /* 8 */,
                                       _nrn_mechanism_field<double>{"m_timeCourse_TIME_SCALE"} /* 9 */,
                                       _nrn_mechanism_field<double>{"m_timeCourse_VOLT_SCALE"} /* 10 */,
                                       _nrn_mechanism_field<double>{"m_steadyState_rate"} /* 11 */,
                                       _nrn_mechanism_field<double>{"m_steadyState_midpoint"} /* 12 */,
                                       _nrn_mechanism_field<double>{"m_steadyState_scale"} /* 13 */,
                                       _nrn_mechanism_field<double>{"m_q10Settings_fixedQ10"} /* 14 */,
                                       _nrn_mechanism_field<double>{"h_instances"} /* 15 */,
                                       _nrn_mechanism_field<double>{"h_forwardRate_rate"} /* 16 */,
                                       _nrn_mechanism_field<double>{"h_forwardRate_midpoint"} /* 17 */,
                                       _nrn_mechanism_field<double>{"h_forwardRate_scale"} /* 18 */,
                                       _nrn_mechanism_field<double>{"h_reverseRate_rate"} /* 19 */,
                                       _nrn_mechanism_field<double>{"h_reverseRate_midpoint"} /* 20 */,
                                       _nrn_mechanism_field<double>{"h_reverseRate_scale"} /* 21 */,
                                       _nrn_mechanism_field<double>{"h_steadyState_rate"} /* 22 */,
                                       _nrn_mechanism_field<double>{"h_steadyState_midpoint"} /* 23 */,
                                       _nrn_mechanism_field<double>{"h_steadyState_scale"} /* 24 */,
                                       _nrn_mechanism_field<double>{"h_q10Settings_fixedQ10"} /* 25 */,
                                       _nrn_mechanism_field<double>{"gion"} /* 26 */,
                                       _nrn_mechanism_field<double>{"m_forwardRate_x"} /* 27 */,
                                       _nrn_mechanism_field<double>{"m_forwardRate_r"} /* 28 */,
                                       _nrn_mechanism_field<double>{"m_reverseRate_x"} /* 29 */,
                                       _nrn_mechanism_field<double>{"m_reverseRate_r"} /* 30 */,
                                       _nrn_mechanism_field<double>{"m_timeCourse_V"} /* 31 */,
                                       _nrn_mechanism_field<double>{"m_timeCourse_ALPHA"} /* 32 */,
                                       _nrn_mechanism_field<double>{"m_timeCourse_BETA"} /* 33 */,
                                       _nrn_mechanism_field<double>{"m_timeCourse_t"} /* 34 */,
                                       _nrn_mechanism_field<double>{"m_steadyState_x"} /* 35 */,
                                       _nrn_mechanism_field<double>{"m_q10Settings_q10"} /* 36 */,
                                       _nrn_mechanism_field<double>{"m_rateScale"} /* 37 */,
                                       _nrn_mechanism_field<double>{"m_alpha"} /* 38 */,
                                       _nrn_mechanism_field<double>{"m_beta"} /* 39 */,
                                       _nrn_mechanism_field<double>{"m_inf"} /* 40 */,
                                       _nrn_mechanism_field<double>{"m_tauUnscaled"} /* 41 */,
                                       _nrn_mechanism_field<double>{"m_tau"} /* 42 */,
                                       _nrn_mechanism_field<double>{"m_fcond"} /* 43 */,
                                       _nrn_mechanism_field<double>{"h_forwardRate_x"} /* 44 */,
                                       _nrn_mechanism_field<double>{"h_forwardRate_r"} /* 45 */,
                                       _nrn_mechanism_field<double>{"h_reverseRate_x"} /* 46 */,
                                       _nrn_mechanism_field<double>{"h_reverseRate_r"} /* 47 */,
                                       _nrn_mechanism_field<double>{"h_steadyState_x"} /* 48 */,
                                       _nrn_mechanism_field<double>{"h_q10Settings_q10"} /* 49 */,
                                       _nrn_mechanism_field<double>{"h_rateScale"} /* 50 */,
                                       _nrn_mechanism_field<double>{"h_alpha"} /* 51 */,
                                       _nrn_mechanism_field<double>{"h_beta"} /* 52 */,
                                       _nrn_mechanism_field<double>{"h_fcond"} /* 53 */,
                                       _nrn_mechanism_field<double>{"h_inf"} /* 54 */,
                                       _nrn_mechanism_field<double>{"h_tau"} /* 55 */,
                                       _nrn_mechanism_field<double>{"conductanceScale"} /* 56 */,
                                       _nrn_mechanism_field<double>{"fopen0"} /* 57 */,
                                       _nrn_mechanism_field<double>{"fopen"} /* 58 */,
                                       _nrn_mechanism_field<double>{"g"} /* 59 */,
                                       _nrn_mechanism_field<double>{"m_q"} /* 60 */,
                                       _nrn_mechanism_field<double>{"h_q"} /* 61 */,
                                       _nrn_mechanism_field<double>{"temperature"} /* 62 */,
                                       _nrn_mechanism_field<double>{"ena"} /* 63 */,
                                       _nrn_mechanism_field<double>{"ina"} /* 64 */,
                                       _nrn_mechanism_field<double>{"rate_m_q"} /* 65 */,
                                       _nrn_mechanism_field<double>{"rate_h_q"} /* 66 */,
                                       _nrn_mechanism_field<double>{"Dm_q"} /* 67 */,
                                       _nrn_mechanism_field<double>{"Dh_q"} /* 68 */,
                                       _nrn_mechanism_field<double>{"v"} /* 69 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 70 */,
                                       _nrn_mechanism_field<double*>{"_ion_ina", "na_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_dinadv", "na_ion"} /* 1 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 2 */);
  hoc_register_prop_size(_mechtype, 71, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Nap_Et2 /home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/Nap_Et2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "Mod file for component: Component(id=Nap_Et2 type=ionChannelHH)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_internalthreadargsproto_);
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[2], _dlist1[2];
 static int states(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dm_q = rate_m_q ;
   Dh_q = rate_h_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 rates ( _threadargs_ ) ;
 Dm_q = Dm_q  / (1. - dt*( 0.0 )) ;
 Dh_q = Dh_q  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states (_internalthreadargsproto_) { {
   rates ( _threadargs_ ) ;
    m_q = m_q - dt*(- ( rate_m_q ) ) ;
    h_q = h_q - dt*(- ( rate_h_q ) ) ;
   }
  return 0;
}
 
static int  rates ( _internalthreadargsproto_ ) {
   m_forwardRate_x = ( v - m_forwardRate_midpoint ) / m_forwardRate_scale ;
   if ( m_forwardRate_x  != 0.0 ) {
     m_forwardRate_r = m_forwardRate_rate * m_forwardRate_x / ( 1.0 - exp ( 0.0 - m_forwardRate_x ) ) ;
     }
   else if ( m_forwardRate_x  == 0.0 ) {
     m_forwardRate_r = m_forwardRate_rate ;
     }
   m_reverseRate_x = ( v - m_reverseRate_midpoint ) / m_reverseRate_scale ;
   if ( m_reverseRate_x  != 0.0 ) {
     m_reverseRate_r = m_reverseRate_rate * m_reverseRate_x / ( 1.0 - exp ( 0.0 - m_reverseRate_x ) ) ;
     }
   else if ( m_reverseRate_x  == 0.0 ) {
     m_reverseRate_r = m_reverseRate_rate ;
     }
   m_timeCourse_V = v / m_timeCourse_VOLT_SCALE ;
   m_timeCourse_ALPHA = m_alpha * m_timeCourse_TIME_SCALE ;
   m_timeCourse_BETA = m_beta * m_timeCourse_TIME_SCALE ;
   if ( ( m_timeCourse_ALPHA + m_timeCourse_BETA )  == 0.0 ) {
     m_timeCourse_t = 0.0 * m_timeCourse_TIME_SCALE ;
     }
   else if ( ( m_timeCourse_ALPHA + m_timeCourse_BETA ) > ( 0.0 ) ) {
     m_timeCourse_t = ( 6.0 / ( ( m_timeCourse_ALPHA + m_timeCourse_BETA ) ) ) * m_timeCourse_TIME_SCALE ;
     }
   else {
     m_timeCourse_t = 0.0 * m_timeCourse_TIME_SCALE ;
     }
   m_steadyState_x = m_steadyState_rate / ( 1.0 + exp ( 0.0 - ( v - m_steadyState_midpoint ) / m_steadyState_scale ) ) ;
   m_q10Settings_q10 = m_q10Settings_fixedQ10 ;
   m_rateScale = m_q10Settings_q10 ;
   m_alpha = m_forwardRate_r ;
   m_beta = m_reverseRate_r ;
   m_inf = m_steadyState_x ;
   m_tauUnscaled = m_timeCourse_t ;
   m_tau = m_tauUnscaled / m_rateScale ;
   m_fcond = pow( m_q , m_instances ) ;
   h_forwardRate_x = ( v - h_forwardRate_midpoint ) / h_forwardRate_scale ;
   if ( h_forwardRate_x  != 0.0 ) {
     h_forwardRate_r = h_forwardRate_rate * h_forwardRate_x / ( 1.0 - exp ( 0.0 - h_forwardRate_x ) ) ;
     }
   else if ( h_forwardRate_x  == 0.0 ) {
     h_forwardRate_r = h_forwardRate_rate ;
     }
   h_reverseRate_x = ( v - h_reverseRate_midpoint ) / h_reverseRate_scale ;
   if ( h_reverseRate_x  != 0.0 ) {
     h_reverseRate_r = h_reverseRate_rate * h_reverseRate_x / ( 1.0 - exp ( 0.0 - h_reverseRate_x ) ) ;
     }
   else if ( h_reverseRate_x  == 0.0 ) {
     h_reverseRate_r = h_reverseRate_rate ;
     }
   h_steadyState_x = h_steadyState_rate / ( 1.0 + exp ( 0.0 - ( v - h_steadyState_midpoint ) / h_steadyState_scale ) ) ;
   h_q10Settings_q10 = h_q10Settings_fixedQ10 ;
   h_rateScale = h_q10Settings_q10 ;
   h_alpha = h_forwardRate_r ;
   h_beta = h_reverseRate_r ;
   h_fcond = pow( h_q , h_instances ) ;
   h_inf = h_steadyState_x ;
   h_tau = 1.0 / ( ( h_alpha + h_beta ) * h_rateScale ) ;
   rate_m_q = ( m_inf - m_q ) / m_tau ;
   rate_h_q = ( h_inf - h_q ) / h_tau ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for rates_Nap_Et2. Requires prior call to setdata_Nap_Et2 and that the specified mechanism instance still be in existence.", NULL);
  }
  Prop* _local_prop = _extcall_prop;
  _nrn_mechanism_cache_instance _ml_real{_local_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _local_prop ? _nrn_mechanism_access_dparam(_local_prop) : nullptr;
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 rates ( _threadargs_ );
 hoc_retpushx(_r);
}
 
static double _npy_rates(Prop* _prop) {
    double _r{0.0};
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 _nrn_mechanism_cache_instance _ml_real{_prop};
auto* const _ml = &_ml_real;
size_t const _iml{};
_ppvar = _nrn_mechanism_access_dparam(_prop);
_thread = _extcall_thread.data();
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
_nt = nrn_threads;
 _r = 1.;
 rates ( _threadargs_ );
 return(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_threadargs_);
  }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 2; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _ode_matsol1 (_threadargs_);
 }
 
static void _ode_matsol(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
   Datum* _ppvar;
   size_t _iml;   _nrn_mechanism_cache_range* _ml;   Node* _nd{};
  double _v{};
  int _cntml;
  _nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
  _ml = &_lmr;
  _cntml = _ml_arg->_nodecount;
  Datum *_thread{_ml_arg->_thread};
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _ppvar = _ml_arg->_pdata[_iml];
    _nd = _ml_arg->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  h_q = h_q0;
  m_q = m_q0;
 {
   ena = 50.0 ;
   temperature = celsius + 273.15 ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   m_q = m_inf ;
   h_q = h_inf ;
   }
 
}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 v = _v;
 initmodel(_threadargs_);
 }
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   conductanceScale = 1.0 ;
   fopen0 = m_fcond * h_fcond ;
   fopen = conductanceScale * fopen0 ;
   g = conductance * fopen ;
   gion = gmax * fopen ;
   ina = gion * ( v - ena ) ;
   }
 _current += ina;

} return _current;
}

static void nrn_cur(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_rhs = _nt->node_rhs_storage();
auto const _vec_sav_rhs = _nt->node_sav_rhs_storage();
auto const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ina += ina ;
	 _vec_rhs[_ni[_iml]] -= _rhs;
 
}
 
}

static void nrn_jacob(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_d = _nt->node_d_storage();
auto const _vec_sav_d = _nt->node_sav_d_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (_iml = 0; _iml < _cntml; ++_iml) {
  _vec_d[_ni[_iml]] += _g;
 
}
 
}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
auto* const _ml = &_lmr;
Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni;
_ni = _ml_arg->_nodeindices;
size_t _cntml = _ml_arg->_nodecount;
_thread = _ml_arg->_thread;
double* _globals = nullptr;
if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
for (size_t _iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
 {   states(_threadargs_);
  } }}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {m_q_columnindex, 0};  _dlist1[0] = {Dm_q_columnindex, 0};
 _slist1[1] = {h_q_columnindex, 0};  _dlist1[1] = {Dh_q_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/Nap_Et2.mod";
    const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=Nap_Et2 type=ionChannelHH)\n"
  "\n"
  "COMMENT\n"
  "\n"
  "    This NEURON file has been generated by org.neuroml.export (see https://github.com/NeuroML/org.neuroml.export)\n"
  "         org.neuroml.export  v1.5.3\n"
  "         org.neuroml.model   v1.5.3\n"
  "         jLEMS               v0.9.9.0\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX Nap_Et2\n"
  "    USEION na WRITE ina VALENCE 1 ? Assuming valence = 1; TODO check this!!\n"
  "    \n"
  "    RANGE gion                           \n"
  "    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!\n"
  "    RANGE conductance                       : parameter\n"
  "    \n"
  "    RANGE g                                 : exposure\n"
  "    \n"
  "    RANGE fopen                             : exposure\n"
  "    RANGE m_instances                       : parameter\n"
  "    \n"
  "    RANGE m_alpha                           : exposure\n"
  "    \n"
  "    RANGE m_beta                            : exposure\n"
  "    \n"
  "    RANGE m_tau                             : exposure\n"
  "    \n"
  "    RANGE m_inf                             : exposure\n"
  "    \n"
  "    RANGE m_rateScale                       : exposure\n"
  "    \n"
  "    RANGE m_fcond                           : exposure\n"
  "    RANGE m_forwardRate_rate                : parameter\n"
  "    RANGE m_forwardRate_midpoint            : parameter\n"
  "    RANGE m_forwardRate_scale               : parameter\n"
  "    \n"
  "    RANGE m_forwardRate_r                   : exposure\n"
  "    RANGE m_reverseRate_rate                : parameter\n"
  "    RANGE m_reverseRate_midpoint            : parameter\n"
  "    RANGE m_reverseRate_scale               : parameter\n"
  "    \n"
  "    RANGE m_reverseRate_r                   : exposure\n"
  "    RANGE m_timeCourse_TIME_SCALE           : parameter\n"
  "    RANGE m_timeCourse_VOLT_SCALE           : parameter\n"
  "    \n"
  "    RANGE m_timeCourse_t                    : exposure\n"
  "    RANGE m_steadyState_rate                : parameter\n"
  "    RANGE m_steadyState_midpoint            : parameter\n"
  "    RANGE m_steadyState_scale               : parameter\n"
  "    \n"
  "    RANGE m_steadyState_x                   : exposure\n"
  "    RANGE m_q10Settings_fixedQ10            : parameter\n"
  "    \n"
  "    RANGE m_q10Settings_q10                 : exposure\n"
  "    RANGE h_instances                       : parameter\n"
  "    \n"
  "    RANGE h_alpha                           : exposure\n"
  "    \n"
  "    RANGE h_beta                            : exposure\n"
  "    \n"
  "    RANGE h_tau                             : exposure\n"
  "    \n"
  "    RANGE h_inf                             : exposure\n"
  "    \n"
  "    RANGE h_rateScale                       : exposure\n"
  "    \n"
  "    RANGE h_fcond                           : exposure\n"
  "    RANGE h_forwardRate_rate                : parameter\n"
  "    RANGE h_forwardRate_midpoint            : parameter\n"
  "    RANGE h_forwardRate_scale               : parameter\n"
  "    \n"
  "    RANGE h_forwardRate_r                   : exposure\n"
  "    RANGE h_reverseRate_rate                : parameter\n"
  "    RANGE h_reverseRate_midpoint            : parameter\n"
  "    RANGE h_reverseRate_scale               : parameter\n"
  "    \n"
  "    RANGE h_reverseRate_r                   : exposure\n"
  "    RANGE h_steadyState_rate                : parameter\n"
  "    RANGE h_steadyState_midpoint            : parameter\n"
  "    RANGE h_steadyState_scale               : parameter\n"
  "    \n"
  "    RANGE h_steadyState_x                   : exposure\n"
  "    RANGE h_q10Settings_fixedQ10            : parameter\n"
  "    \n"
  "    RANGE h_q10Settings_q10                 : exposure\n"
  "    RANGE m_forwardRate_x                   : derived variable\n"
  "    RANGE m_reverseRate_x                   : derived variable\n"
  "    RANGE m_timeCourse_V                    : derived variable\n"
  "    RANGE m_timeCourse_ALPHA                : derived variable\n"
  "    RANGE m_timeCourse_BETA                 : derived variable\n"
  "    RANGE m_tauUnscaled                     : derived variable\n"
  "    RANGE h_forwardRate_x                   : derived variable\n"
  "    RANGE h_reverseRate_x                   : derived variable\n"
  "    RANGE conductanceScale                  : derived variable\n"
  "    RANGE fopen0                            : derived variable\n"
  "    \n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    \n"
  "    (nA) = (nanoamp)\n"
  "    (uA) = (microamp)\n"
  "    (mA) = (milliamp)\n"
  "    (A) = (amp)\n"
  "    (mV) = (millivolt)\n"
  "    (mS) = (millisiemens)\n"
  "    (uS) = (microsiemens)\n"
  "    (molar) = (1/liter)\n"
  "    (kHz) = (kilohertz)\n"
  "    (mM) = (millimolar)\n"
  "    (um) = (micrometer)\n"
  "    (umol) = (micromole)\n"
  "    (S) = (siemens)\n"
  "    \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    \n"
  "    gmax = 0  (S/cm2)                       : Will be changed when ion channel mechanism placed on cell!\n"
  "    \n"
  "    conductance = 1.0E-5 (uS)\n"
  "    m_instances = 3 \n"
  "    m_forwardRate_rate = 1.092 (kHz)\n"
  "    m_forwardRate_midpoint = -38 (mV)\n"
  "    m_forwardRate_scale = 6 (mV)\n"
  "    m_reverseRate_rate = 0.744 (kHz)\n"
  "    m_reverseRate_midpoint = -38 (mV)\n"
  "    m_reverseRate_scale = -6 (mV)\n"
  "    m_timeCourse_TIME_SCALE = 1 (ms)\n"
  "    m_timeCourse_VOLT_SCALE = 1 (mV)\n"
  "    m_steadyState_rate = 1 \n"
  "    m_steadyState_midpoint = -52.6 (mV)\n"
  "    m_steadyState_scale = 4.6 (mV)\n"
  "    m_q10Settings_fixedQ10 = 2.9528825 \n"
  "    h_instances = 1 \n"
  "    h_forwardRate_rate = 1.3334401E-5 (kHz)\n"
  "    h_forwardRate_midpoint = -17 (mV)\n"
  "    h_forwardRate_scale = -4.63 (mV)\n"
  "    h_reverseRate_rate = 1.8252202E-5 (kHz)\n"
  "    h_reverseRate_midpoint = -64.4 (mV)\n"
  "    h_reverseRate_scale = 2.63 (mV)\n"
  "    h_steadyState_rate = 1 \n"
  "    h_steadyState_midpoint = -48.8 (mV)\n"
  "    h_steadyState_scale = -10 (mV)\n"
  "    h_q10Settings_fixedQ10 = 2.9528825 \n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    \n"
  "    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel\n"
  "    v (mV)\n"
  "    celsius (degC)\n"
  "    temperature (K)\n"
  "    ena (mV)\n"
  "    ina (mA/cm2)\n"
  "    \n"
  "    \n"
  "    m_forwardRate_x                        : derived variable\n"
  "    \n"
  "    m_forwardRate_r (kHz)                  : conditional derived var...\n"
  "    \n"
  "    m_reverseRate_x                        : derived variable\n"
  "    \n"
  "    m_reverseRate_r (kHz)                  : conditional derived var...\n"
  "    \n"
  "    m_timeCourse_V                         : derived variable\n"
  "    \n"
  "    m_timeCourse_ALPHA                     : derived variable\n"
  "    \n"
  "    m_timeCourse_BETA                      : derived variable\n"
  "    \n"
  "    m_timeCourse_t (ms)                    : conditional derived var...\n"
  "    \n"
  "    m_steadyState_x                        : derived variable\n"
  "    \n"
  "    m_q10Settings_q10                      : derived variable\n"
  "    \n"
  "    m_rateScale                            : derived variable\n"
  "    \n"
  "    m_alpha (kHz)                          : derived variable\n"
  "    \n"
  "    m_beta (kHz)                           : derived variable\n"
  "    \n"
  "    m_inf                                  : derived variable\n"
  "    \n"
  "    m_tauUnscaled (ms)                     : derived variable\n"
  "    \n"
  "    m_tau (ms)                             : derived variable\n"
  "    \n"
  "    m_fcond                                : derived variable\n"
  "    \n"
  "    h_forwardRate_x                        : derived variable\n"
  "    \n"
  "    h_forwardRate_r (kHz)                  : conditional derived var...\n"
  "    \n"
  "    h_reverseRate_x                        : derived variable\n"
  "    \n"
  "    h_reverseRate_r (kHz)                  : conditional derived var...\n"
  "    \n"
  "    h_steadyState_x                        : derived variable\n"
  "    \n"
  "    h_q10Settings_q10                      : derived variable\n"
  "    \n"
  "    h_rateScale                            : derived variable\n"
  "    \n"
  "    h_alpha (kHz)                          : derived variable\n"
  "    \n"
  "    h_beta (kHz)                           : derived variable\n"
  "    \n"
  "    h_fcond                                : derived variable\n"
  "    \n"
  "    h_inf                                  : derived variable\n"
  "    \n"
  "    h_tau (ms)                             : derived variable\n"
  "    \n"
  "    conductanceScale                       : derived variable\n"
  "    \n"
  "    fopen0                                 : derived variable\n"
  "    \n"
  "    fopen                                  : derived variable\n"
  "    \n"
  "    g (uS)                                 : derived variable\n"
  "    rate_m_q (/ms)\n"
  "    rate_h_q (/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    m_q  \n"
  "    h_q  \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    ena = 50.0\n"
  "    \n"
  "    temperature = celsius + 273.15\n"
  "    \n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    m_q = m_inf\n"
  "    \n"
  "    h_q = h_inf\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=Nap_Et2 type=ionChannelHH), from conductanceScaling; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    conductanceScale = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=Nap_Et2 type=ionChannelHH), from gates; Component(id=m type=gateHHratesTauInf)\n"
  "    ? multiply applied to all instances of fcond in: <gates> ([Component(id=m type=gateHHratesTauInf), Component(id=h type=gateHHratesInf)]))\n"
  "    fopen0 = m_fcond * h_fcond ? path based, prefix = \n"
  "    \n"
  "    fopen = conductanceScale  *  fopen0 ? evaluable\n"
  "    g = conductance  *  fopen ? evaluable\n"
  "    gion = gmax * fopen \n"
  "    \n"
  "    ina = gion * (v - ena)\n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    m_q' = rate_m_q \n"
  "    h_q' = rate_h_q \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    \n"
  "    m_forwardRate_x = (v -  m_forwardRate_midpoint ) /  m_forwardRate_scale ? evaluable\n"
  "    if (m_forwardRate_x  != 0)  { \n"
  "        m_forwardRate_r = m_forwardRate_rate  *  m_forwardRate_x  / (1 - exp(0 -  m_forwardRate_x )) ? evaluable cdv\n"
  "    } else if (m_forwardRate_x  == 0)  { \n"
  "        m_forwardRate_r = m_forwardRate_rate ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    m_reverseRate_x = (v -  m_reverseRate_midpoint ) /  m_reverseRate_scale ? evaluable\n"
  "    if (m_reverseRate_x  != 0)  { \n"
  "        m_reverseRate_r = m_reverseRate_rate  *  m_reverseRate_x  / (1 - exp(0 -  m_reverseRate_x )) ? evaluable cdv\n"
  "    } else if (m_reverseRate_x  == 0)  { \n"
  "        m_reverseRate_r = m_reverseRate_rate ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    m_timeCourse_V = v /  m_timeCourse_VOLT_SCALE ? evaluable\n"
  "    m_timeCourse_ALPHA = m_alpha  *  m_timeCourse_TIME_SCALE ? evaluable\n"
  "    m_timeCourse_BETA = m_beta  *  m_timeCourse_TIME_SCALE ? evaluable\n"
  "    if (( m_timeCourse_ALPHA  +  m_timeCourse_BETA ) == 0)  { \n"
  "        m_timeCourse_t = 0.0  *  m_timeCourse_TIME_SCALE ? evaluable cdv\n"
  "    } else if (( m_timeCourse_ALPHA  +  m_timeCourse_BETA )  > ( 0 ))  { \n"
  "        m_timeCourse_t = ( 6/( ( m_timeCourse_ALPHA  +  m_timeCourse_BETA ) ) ) *  m_timeCourse_TIME_SCALE ? evaluable cdv\n"
  "    } else  { \n"
  "        m_timeCourse_t = 0.0  *  m_timeCourse_TIME_SCALE ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    m_steadyState_x = m_steadyState_rate  / (1 + exp(0 - (v -  m_steadyState_midpoint )/ m_steadyState_scale )) ? evaluable\n"
  "    m_q10Settings_q10 = m_q10Settings_fixedQ10 ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=m type=gateHHratesTauInf), from q10Settings; Component(id=null type=q10Fixed)\n"
  "    ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10Fixed)]))\n"
  "    m_rateScale = m_q10Settings_q10 ? path based, prefix = m_\n"
  "    \n"
  "    ? DerivedVariable is based on path: forwardRate/r, on: Component(id=m type=gateHHratesTauInf), from forwardRate; Component(id=null type=HHExpLinearRate)\n"
  "    m_alpha = m_forwardRate_r ? path based, prefix = m_\n"
  "    \n"
  "    ? DerivedVariable is based on path: reverseRate/r, on: Component(id=m type=gateHHratesTauInf), from reverseRate; Component(id=null type=HHExpLinearRate)\n"
  "    m_beta = m_reverseRate_r ? path based, prefix = m_\n"
  "    \n"
  "    ? DerivedVariable is based on path: steadyState/x, on: Component(id=m type=gateHHratesTauInf), from steadyState; Component(id=null type=HHSigmoidVariable)\n"
  "    m_inf = m_steadyState_x ? path based, prefix = m_\n"
  "    \n"
  "    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=m type=gateHHratesTauInf), from timeCourse; Component(id=null type=Nap_Et2_m_tau_tau)\n"
  "    m_tauUnscaled = m_timeCourse_t ? path based, prefix = m_\n"
  "    \n"
  "    m_tau = m_tauUnscaled  /  m_rateScale ? evaluable\n"
  "    m_fcond = m_q ^ m_instances ? evaluable\n"
  "    h_forwardRate_x = (v -  h_forwardRate_midpoint ) /  h_forwardRate_scale ? evaluable\n"
  "    if (h_forwardRate_x  != 0)  { \n"
  "        h_forwardRate_r = h_forwardRate_rate  *  h_forwardRate_x  / (1 - exp(0 -  h_forwardRate_x )) ? evaluable cdv\n"
  "    } else if (h_forwardRate_x  == 0)  { \n"
  "        h_forwardRate_r = h_forwardRate_rate ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    h_reverseRate_x = (v -  h_reverseRate_midpoint ) /  h_reverseRate_scale ? evaluable\n"
  "    if (h_reverseRate_x  != 0)  { \n"
  "        h_reverseRate_r = h_reverseRate_rate  *  h_reverseRate_x  / (1 - exp(0 -  h_reverseRate_x )) ? evaluable cdv\n"
  "    } else if (h_reverseRate_x  == 0)  { \n"
  "        h_reverseRate_r = h_reverseRate_rate ? evaluable cdv\n"
  "    }\n"
  "    \n"
  "    h_steadyState_x = h_steadyState_rate  / (1 + exp(0 - (v -  h_steadyState_midpoint )/ h_steadyState_scale )) ? evaluable\n"
  "    h_q10Settings_q10 = h_q10Settings_fixedQ10 ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=h type=gateHHratesInf), from q10Settings; Component(id=null type=q10Fixed)\n"
  "    ? multiply applied to all instances of q10 in: <q10Settings> ([Component(id=null type=q10Fixed)]))\n"
  "    h_rateScale = h_q10Settings_q10 ? path based, prefix = h_\n"
  "    \n"
  "    ? DerivedVariable is based on path: forwardRate/r, on: Component(id=h type=gateHHratesInf), from forwardRate; Component(id=null type=HHExpLinearRate)\n"
  "    h_alpha = h_forwardRate_r ? path based, prefix = h_\n"
  "    \n"
  "    ? DerivedVariable is based on path: reverseRate/r, on: Component(id=h type=gateHHratesInf), from reverseRate; Component(id=null type=HHExpLinearRate)\n"
  "    h_beta = h_reverseRate_r ? path based, prefix = h_\n"
  "    \n"
  "    h_fcond = h_q ^ h_instances ? evaluable\n"
  "    ? DerivedVariable is based on path: steadyState/x, on: Component(id=h type=gateHHratesInf), from steadyState; Component(id=null type=HHSigmoidVariable)\n"
  "    h_inf = h_steadyState_x ? path based, prefix = h_\n"
  "    \n"
  "    h_tau = 1/(( h_alpha + h_beta ) *  h_rateScale ) ? evaluable\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    rate_m_q = ( m_inf  -  m_q ) /  m_tau ? Note units of all quantities used here need to be consistent!\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    rate_h_q = ( h_inf  -  h_q ) /  h_tau ? Note units of all quantities used here need to be consistent!\n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "     \n"
  "    \n"
  "}\n"
  "\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
