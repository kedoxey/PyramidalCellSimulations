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
static constexpr auto number_of_datum_variables = 4;
static constexpr auto number_of_floating_point_variables = 35;
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
 
#define nrn_init _nrn_init__SK_E2
#define _nrn_initial _nrn_initial__SK_E2
#define nrn_cur _nrn_cur__SK_E2
#define _nrn_current _nrn_current__SK_E2
#define nrn_jacob _nrn_jacob__SK_E2
#define nrn_state _nrn_state__SK_E2
#define _net_receive _net_receive__SK_E2 
#define rates rates__SK_E2 
#define states states__SK_E2 
 
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
#define z_instances _ml->template fpfield<2>(_iml)
#define z_instances_columnindex 2
#define z_timeCourse_TIME_SCALE _ml->template fpfield<3>(_iml)
#define z_timeCourse_TIME_SCALE_columnindex 3
#define z_timeCourse_VOLT_SCALE _ml->template fpfield<4>(_iml)
#define z_timeCourse_VOLT_SCALE_columnindex 4
#define z_timeCourse_CONC_SCALE _ml->template fpfield<5>(_iml)
#define z_timeCourse_CONC_SCALE_columnindex 5
#define z_steadyState_TIME_SCALE _ml->template fpfield<6>(_iml)
#define z_steadyState_TIME_SCALE_columnindex 6
#define z_steadyState_VOLT_SCALE _ml->template fpfield<7>(_iml)
#define z_steadyState_VOLT_SCALE_columnindex 7
#define z_steadyState_CONC_SCALE _ml->template fpfield<8>(_iml)
#define z_steadyState_CONC_SCALE_columnindex 8
#define gion _ml->template fpfield<9>(_iml)
#define gion_columnindex 9
#define z_timeCourse_V _ml->template fpfield<10>(_iml)
#define z_timeCourse_V_columnindex 10
#define z_timeCourse_ca_conc _ml->template fpfield<11>(_iml)
#define z_timeCourse_ca_conc_columnindex 11
#define z_timeCourse_t _ml->template fpfield<12>(_iml)
#define z_timeCourse_t_columnindex 12
#define z_steadyState_V _ml->template fpfield<13>(_iml)
#define z_steadyState_V_columnindex 13
#define z_steadyState_ca_conc _ml->template fpfield<14>(_iml)
#define z_steadyState_ca_conc_columnindex 14
#define z_steadyState_x _ml->template fpfield<15>(_iml)
#define z_steadyState_x_columnindex 15
#define z_rateScale _ml->template fpfield<16>(_iml)
#define z_rateScale_columnindex 16
#define z_fcond _ml->template fpfield<17>(_iml)
#define z_fcond_columnindex 17
#define z_inf _ml->template fpfield<18>(_iml)
#define z_inf_columnindex 18
#define z_tauUnscaled _ml->template fpfield<19>(_iml)
#define z_tauUnscaled_columnindex 19
#define z_tau _ml->template fpfield<20>(_iml)
#define z_tau_columnindex 20
#define conductanceScale _ml->template fpfield<21>(_iml)
#define conductanceScale_columnindex 21
#define fopen0 _ml->template fpfield<22>(_iml)
#define fopen0_columnindex 22
#define fopen _ml->template fpfield<23>(_iml)
#define fopen_columnindex 23
#define g _ml->template fpfield<24>(_iml)
#define g_columnindex 24
#define z_q _ml->template fpfield<25>(_iml)
#define z_q_columnindex 25
#define temperature _ml->template fpfield<26>(_iml)
#define temperature_columnindex 26
#define ek _ml->template fpfield<27>(_iml)
#define ek_columnindex 27
#define ik _ml->template fpfield<28>(_iml)
#define ik_columnindex 28
#define cai _ml->template fpfield<29>(_iml)
#define cai_columnindex 29
#define cao _ml->template fpfield<30>(_iml)
#define cao_columnindex 30
#define rate_z_q _ml->template fpfield<31>(_iml)
#define rate_z_q_columnindex 31
#define Dz_q _ml->template fpfield<32>(_iml)
#define Dz_q_columnindex 32
#define v _ml->template fpfield<33>(_iml)
#define v_columnindex 33
#define _g _ml->template fpfield<34>(_iml)
#define _g_columnindex 34
#define _ion_cai *(_ml->dptr_field<0>(_iml))
#define _p_ion_cai static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_cao *(_ml->dptr_field<1>(_iml))
#define _p_ion_cao static_cast<neuron::container::data_handle<double>>(_ppvar[1])
#define _ion_ik *(_ml->dptr_field<2>(_iml))
#define _p_ion_ik static_cast<neuron::container::data_handle<double>>(_ppvar[2])
#define _ion_dikdv *(_ml->dptr_field<3>(_iml))
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
 {"setdata_SK_E2", _hoc_setdata},
 {"rates_SK_E2", _hoc_rates},
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
 {"gmax_SK_E2", "S/cm2"},
 {"conductance_SK_E2", "uS"},
 {"z_timeCourse_TIME_SCALE_SK_E2", "ms"},
 {"z_timeCourse_VOLT_SCALE_SK_E2", "mV"},
 {"z_timeCourse_CONC_SCALE_SK_E2", "mM"},
 {"z_steadyState_TIME_SCALE_SK_E2", "ms"},
 {"z_steadyState_VOLT_SCALE_SK_E2", "mV"},
 {"z_steadyState_CONC_SCALE_SK_E2", "mM"},
 {"gion_SK_E2", "S/cm2"},
 {"z_timeCourse_t_SK_E2", "ms"},
 {"z_tauUnscaled_SK_E2", "ms"},
 {"z_tau_SK_E2", "ms"},
 {"g_SK_E2", "uS"},
 {0, 0}
};
 static double delta_t = 0.01;
 static double z_q0 = 0;
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
 
#define _cvode_ieq _ppvar[4].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"SK_E2",
 "gmax_SK_E2",
 "conductance_SK_E2",
 "z_instances_SK_E2",
 "z_timeCourse_TIME_SCALE_SK_E2",
 "z_timeCourse_VOLT_SCALE_SK_E2",
 "z_timeCourse_CONC_SCALE_SK_E2",
 "z_steadyState_TIME_SCALE_SK_E2",
 "z_steadyState_VOLT_SCALE_SK_E2",
 "z_steadyState_CONC_SCALE_SK_E2",
 0,
 "gion_SK_E2",
 "z_timeCourse_V_SK_E2",
 "z_timeCourse_ca_conc_SK_E2",
 "z_timeCourse_t_SK_E2",
 "z_steadyState_V_SK_E2",
 "z_steadyState_ca_conc_SK_E2",
 "z_steadyState_x_SK_E2",
 "z_rateScale_SK_E2",
 "z_fcond_SK_E2",
 "z_inf_SK_E2",
 "z_tauUnscaled_SK_E2",
 "z_tau_SK_E2",
 "conductanceScale_SK_E2",
 "fopen0_SK_E2",
 "fopen_SK_E2",
 "g_SK_E2",
 0,
 "z_q_SK_E2",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _k_sym;
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0, /* gmax */
     1e-05, /* conductance */
     1, /* z_instances */
     1, /* z_timeCourse_TIME_SCALE */
     1, /* z_timeCourse_VOLT_SCALE */
     1e+06, /* z_timeCourse_CONC_SCALE */
     1, /* z_steadyState_TIME_SCALE */
     1, /* z_steadyState_VOLT_SCALE */
     1e+06, /* z_steadyState_CONC_SCALE */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 35);
 	/*initialize range parameters*/
 	gmax = _parm_default[0]; /* 0 */
 	conductance = _parm_default[1]; /* 1e-05 */
 	z_instances = _parm_default[2]; /* 1 */
 	z_timeCourse_TIME_SCALE = _parm_default[3]; /* 1 */
 	z_timeCourse_VOLT_SCALE = _parm_default[4]; /* 1 */
 	z_timeCourse_CONC_SCALE = _parm_default[5]; /* 1e+06 */
 	z_steadyState_TIME_SCALE = _parm_default[6]; /* 1 */
 	z_steadyState_VOLT_SCALE = _parm_default[7]; /* 1 */
 	z_steadyState_CONC_SCALE = _parm_default[8]; /* 1e+06 */
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 35);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 1); /* cai */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 2); /* cao */
 prop_ion = need_memb(_k_sym);
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ik */
 	_ppvar[3] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dikdv */
 
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

 extern "C" void _SK_E2_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", 2.0);
 	ion_reg("k", 1.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	_k_sym = hoc_lookup("k_ion");
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
                                       _nrn_mechanism_field<double>{"z_instances"} /* 2 */,
                                       _nrn_mechanism_field<double>{"z_timeCourse_TIME_SCALE"} /* 3 */,
                                       _nrn_mechanism_field<double>{"z_timeCourse_VOLT_SCALE"} /* 4 */,
                                       _nrn_mechanism_field<double>{"z_timeCourse_CONC_SCALE"} /* 5 */,
                                       _nrn_mechanism_field<double>{"z_steadyState_TIME_SCALE"} /* 6 */,
                                       _nrn_mechanism_field<double>{"z_steadyState_VOLT_SCALE"} /* 7 */,
                                       _nrn_mechanism_field<double>{"z_steadyState_CONC_SCALE"} /* 8 */,
                                       _nrn_mechanism_field<double>{"gion"} /* 9 */,
                                       _nrn_mechanism_field<double>{"z_timeCourse_V"} /* 10 */,
                                       _nrn_mechanism_field<double>{"z_timeCourse_ca_conc"} /* 11 */,
                                       _nrn_mechanism_field<double>{"z_timeCourse_t"} /* 12 */,
                                       _nrn_mechanism_field<double>{"z_steadyState_V"} /* 13 */,
                                       _nrn_mechanism_field<double>{"z_steadyState_ca_conc"} /* 14 */,
                                       _nrn_mechanism_field<double>{"z_steadyState_x"} /* 15 */,
                                       _nrn_mechanism_field<double>{"z_rateScale"} /* 16 */,
                                       _nrn_mechanism_field<double>{"z_fcond"} /* 17 */,
                                       _nrn_mechanism_field<double>{"z_inf"} /* 18 */,
                                       _nrn_mechanism_field<double>{"z_tauUnscaled"} /* 19 */,
                                       _nrn_mechanism_field<double>{"z_tau"} /* 20 */,
                                       _nrn_mechanism_field<double>{"conductanceScale"} /* 21 */,
                                       _nrn_mechanism_field<double>{"fopen0"} /* 22 */,
                                       _nrn_mechanism_field<double>{"fopen"} /* 23 */,
                                       _nrn_mechanism_field<double>{"g"} /* 24 */,
                                       _nrn_mechanism_field<double>{"z_q"} /* 25 */,
                                       _nrn_mechanism_field<double>{"temperature"} /* 26 */,
                                       _nrn_mechanism_field<double>{"ek"} /* 27 */,
                                       _nrn_mechanism_field<double>{"ik"} /* 28 */,
                                       _nrn_mechanism_field<double>{"cai"} /* 29 */,
                                       _nrn_mechanism_field<double>{"cao"} /* 30 */,
                                       _nrn_mechanism_field<double>{"rate_z_q"} /* 31 */,
                                       _nrn_mechanism_field<double>{"Dz_q"} /* 32 */,
                                       _nrn_mechanism_field<double>{"v"} /* 33 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 34 */,
                                       _nrn_mechanism_field<double*>{"_ion_cai", "ca_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_cao", "ca_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_ik", "k_ion"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"_ion_dikdv", "k_ion"} /* 3 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 4 */);
  hoc_register_prop_size(_mechtype, 35, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 SK_E2 /home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/SK_E2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "Mod file for component: Component(id=SK_E2 type=ionChannelHH)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_internalthreadargsproto_);
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[1], _dlist1[1];
 static int states(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   rates ( _threadargs_ ) ;
   Dz_q = rate_z_q ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 rates ( _threadargs_ ) ;
 Dz_q = Dz_q  / (1. - dt*( 0.0 )) ;
  return 0;
}
 /*END CVODE*/
 static int states (_internalthreadargsproto_) { {
   rates ( _threadargs_ ) ;
    z_q = z_q - dt*(- ( rate_z_q ) ) ;
   }
  return 0;
}
 
static int  rates ( _internalthreadargsproto_ ) {
   double _lcaConc ;
 _lcaConc = cai ;
   z_timeCourse_V = v / z_timeCourse_VOLT_SCALE ;
   z_timeCourse_ca_conc = _lcaConc / z_timeCourse_CONC_SCALE ;
   z_timeCourse_t = 1.0 * z_timeCourse_TIME_SCALE ;
   z_steadyState_V = v / z_steadyState_VOLT_SCALE ;
   z_steadyState_ca_conc = _lcaConc / z_steadyState_CONC_SCALE ;
   z_steadyState_x = 1.0 / ( 1.0 + pow( ( 4.3e-10 / z_steadyState_ca_conc ) , 4.8 ) ) ;
   z_rateScale = 1.0 ;
   z_fcond = pow( z_q , z_instances ) ;
   z_inf = z_steadyState_x ;
   z_tauUnscaled = z_timeCourse_t ;
   z_tau = z_tauUnscaled / z_rateScale ;
   rate_z_q = ( z_inf - z_q ) / z_tau ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
 
  if(!_prop_id) {
    hoc_execerror("No data for rates_SK_E2. Requires prior call to setdata_SK_E2 and that the specified mechanism instance still be in existence.", NULL);
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
 
static int _ode_count(int _type){ return 1;}
 
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
  cai = _ion_cai;
  cao = _ion_cao;
     _ode_spec1 (_threadargs_);
  }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  Datum* _ppvar;
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 1; ++_i) {
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
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  z_q = z_q0;
 {
   ek = - 85.0 ;
   temperature = celsius + 273.15 ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   z_q = z_inf ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel(_threadargs_);
 }
}

static double _nrn_current(_internalthreadargsprotocomma_ double _v) {
double _current=0.; v=_v;
{ {
   conductanceScale = 1.0 ;
   fopen0 = z_fcond ;
   fopen = conductanceScale * fopen0 ;
   g = conductance * fopen ;
   gion = gmax * fopen ;
   ik = gion * ( v - ek ) ;
   }
 _current += ik;

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
  cai = _ion_cai;
  cao = _ion_cao;
 auto const _g_local = _nrn_current(_threadargscomma_ _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ik += ik ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 {   states(_threadargs_);
  } }}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {z_q_columnindex, 0};  _dlist1[0] = {Dz_q_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/SK_E2.mod";
    const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=SK_E2 type=ionChannelHH)\n"
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
  "    SUFFIX SK_E2\n"
  "    USEION ca READ cai,cao VALENCE 2\n"
  "    USEION k WRITE ik VALENCE 1 ? Assuming valence = 1; TODO check this!!\n"
  "    \n"
  "    RANGE gion                           \n"
  "    RANGE gmax                              : Will be changed when ion channel mechanism placed on cell!\n"
  "    RANGE conductance                       : parameter\n"
  "    \n"
  "    RANGE g                                 : exposure\n"
  "    \n"
  "    RANGE fopen                             : exposure\n"
  "    RANGE z_instances                       : parameter\n"
  "    \n"
  "    RANGE z_tau                             : exposure\n"
  "    \n"
  "    RANGE z_inf                             : exposure\n"
  "    \n"
  "    RANGE z_rateScale                       : exposure\n"
  "    \n"
  "    RANGE z_fcond                           : exposure\n"
  "    RANGE z_timeCourse_TIME_SCALE           : parameter\n"
  "    RANGE z_timeCourse_VOLT_SCALE           : parameter\n"
  "    RANGE z_timeCourse_CONC_SCALE           : parameter\n"
  "    \n"
  "    RANGE z_timeCourse_t                    : exposure\n"
  "    RANGE z_steadyState_TIME_SCALE          : parameter\n"
  "    RANGE z_steadyState_VOLT_SCALE          : parameter\n"
  "    RANGE z_steadyState_CONC_SCALE          : parameter\n"
  "    \n"
  "    RANGE z_steadyState_x                   : exposure\n"
  "    RANGE z_timeCourse_V                    : derived variable\n"
  "    RANGE z_timeCourse_ca_conc              : derived variable\n"
  "    RANGE z_steadyState_V                   : derived variable\n"
  "    RANGE z_steadyState_ca_conc             : derived variable\n"
  "    RANGE z_tauUnscaled                     : derived variable\n"
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
  "    z_instances = 1 \n"
  "    z_timeCourse_TIME_SCALE = 1 (ms)\n"
  "    z_timeCourse_VOLT_SCALE = 1 (mV)\n"
  "    z_timeCourse_CONC_SCALE = 1000000 (mM)\n"
  "    z_steadyState_TIME_SCALE = 1 (ms)\n"
  "    z_steadyState_VOLT_SCALE = 1 (mV)\n"
  "    z_steadyState_CONC_SCALE = 1000000 (mM)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    \n"
  "    gion   (S/cm2)                          : Transient conductance density of the channel? Standard Assigned variables with ionChannel\n"
  "    v (mV)\n"
  "    celsius (degC)\n"
  "    temperature (K)\n"
  "    ek (mV)\n"
  "    ik (mA/cm2)\n"
  "    \n"
  "    cai (mM)\n"
  "    \n"
  "    cao (mM)\n"
  "    \n"
  "    \n"
  "    z_timeCourse_V                         : derived variable\n"
  "    \n"
  "    z_timeCourse_ca_conc                   : derived variable\n"
  "    \n"
  "    z_timeCourse_t (ms)                    : derived variable\n"
  "    \n"
  "    z_steadyState_V                        : derived variable\n"
  "    \n"
  "    z_steadyState_ca_conc                  : derived variable\n"
  "    \n"
  "    z_steadyState_x                        : derived variable\n"
  "    \n"
  "    z_rateScale                            : derived variable\n"
  "    \n"
  "    z_fcond                                : derived variable\n"
  "    \n"
  "    z_inf                                  : derived variable\n"
  "    \n"
  "    z_tauUnscaled (ms)                     : derived variable\n"
  "    \n"
  "    z_tau (ms)                             : derived variable\n"
  "    \n"
  "    conductanceScale                       : derived variable\n"
  "    \n"
  "    fopen0                                 : derived variable\n"
  "    \n"
  "    fopen                                  : derived variable\n"
  "    \n"
  "    g (uS)                                 : derived variable\n"
  "    rate_z_q (/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    z_q  \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    ek = -85.0\n"
  "    \n"
  "    temperature = celsius + 273.15\n"
  "    \n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    z_q = z_inf\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    ? DerivedVariable is based on path: conductanceScaling[*]/factor, on: Component(id=SK_E2 type=ionChannelHH), from conductanceScaling; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    conductanceScale = 1 \n"
  "    \n"
  "    ? DerivedVariable is based on path: gates[*]/fcond, on: Component(id=SK_E2 type=ionChannelHH), from gates; Component(id=z type=gateHHtauInf)\n"
  "    ? multiply applied to all instances of fcond in: <gates> ([Component(id=z type=gateHHtauInf)]))\n"
  "    fopen0 = z_fcond ? path based, prefix = \n"
  "    \n"
  "    fopen = conductanceScale  *  fopen0 ? evaluable\n"
  "    g = conductance  *  fopen ? evaluable\n"
  "    gion = gmax * fopen \n"
  "    \n"
  "    ik = gion * (v - ek)\n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    z_q' = rate_z_q \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    LOCAL caConc\n"
  "    \n"
  "    caConc = cai\n"
  "    \n"
  "    z_timeCourse_V = v /  z_timeCourse_VOLT_SCALE ? evaluable\n"
  "    z_timeCourse_ca_conc = caConc /  z_timeCourse_CONC_SCALE ? evaluable\n"
  "    z_timeCourse_t = 1.0  *  z_timeCourse_TIME_SCALE ? evaluable\n"
  "    z_steadyState_V = v /  z_steadyState_VOLT_SCALE ? evaluable\n"
  "    z_steadyState_ca_conc = caConc /  z_steadyState_CONC_SCALE ? evaluable\n"
  "    z_steadyState_x = 1/(1+(4.3e-10/ z_steadyState_ca_conc )^4.8) ? evaluable\n"
  "    ? DerivedVariable is based on path: q10Settings[*]/q10, on: Component(id=z type=gateHHtauInf), from q10Settings; null\n"
  "    ? Path not present in component, using factor: 1\n"
  "    \n"
  "    z_rateScale = 1 \n"
  "    \n"
  "    z_fcond = z_q ^ z_instances ? evaluable\n"
  "    ? DerivedVariable is based on path: steadyState/x, on: Component(id=z type=gateHHtauInf), from steadyState; Component(id=null type=SK_E2_z_inf_inf)\n"
  "    z_inf = z_steadyState_x ? path based, prefix = z_\n"
  "    \n"
  "    ? DerivedVariable is based on path: timeCourse/t, on: Component(id=z type=gateHHtauInf), from timeCourse; Component(id=null type=SK_E2_z_tau_tau)\n"
  "    z_tauUnscaled = z_timeCourse_t ? path based, prefix = z_\n"
  "    \n"
  "    z_tau = z_tauUnscaled  /  z_rateScale ? evaluable\n"
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
  "    rate_z_q = ( z_inf  -  z_q ) /  z_tau ? Note units of all quantities used here need to be consistent!\n"
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
