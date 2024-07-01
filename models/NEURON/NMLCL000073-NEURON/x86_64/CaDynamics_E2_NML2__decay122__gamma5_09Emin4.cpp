/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format on
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
// clang-format off
#include "neuron/cache/mechanism_range.hpp"
#include <vector>
using std::size_t;
static auto& std_cerr_stream = std::cerr;
static constexpr auto number_of_datum_variables = 7;
static constexpr auto number_of_floating_point_variables = 15;
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
 
#define nrn_init _nrn_init__CaDynamics_E2_NML2__decay122__gamma5_09Emin4
#define _nrn_initial _nrn_initial__CaDynamics_E2_NML2__decay122__gamma5_09Emin4
#define nrn_cur _nrn_cur__CaDynamics_E2_NML2__decay122__gamma5_09Emin4
#define _nrn_current _nrn_current__CaDynamics_E2_NML2__decay122__gamma5_09Emin4
#define nrn_jacob _nrn_jacob__CaDynamics_E2_NML2__decay122__gamma5_09Emin4
#define nrn_state _nrn_state__CaDynamics_E2_NML2__decay122__gamma5_09Emin4
#define _net_receive _net_receive__CaDynamics_E2_NML2__decay122__gamma5_09Emin4 
#define rates rates__CaDynamics_E2_NML2__decay122__gamma5_09Emin4 
#define states states__CaDynamics_E2_NML2__decay122__gamma5_09Emin4 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _internalthreadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
#define _internalthreadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gamma _ml->template fpfield<0>(_iml)
#define gamma_columnindex 0
#define minCai _ml->template fpfield<1>(_iml)
#define minCai_columnindex 1
#define decay _ml->template fpfield<2>(_iml)
#define decay_columnindex 2
#define depth _ml->template fpfield<3>(_iml)
#define depth_columnindex 3
#define Faraday _ml->template fpfield<4>(_iml)
#define Faraday_columnindex 4
#define currDensCa _ml->template fpfield<5>(_iml)
#define currDensCa_columnindex 5
#define concentration _ml->template fpfield<6>(_iml)
#define concentration_columnindex 6
#define extConcentration _ml->template fpfield<7>(_iml)
#define extConcentration_columnindex 7
#define cai _ml->template fpfield<8>(_iml)
#define cai_columnindex 8
#define cao _ml->template fpfield<9>(_iml)
#define cao_columnindex 9
#define ica _ml->template fpfield<10>(_iml)
#define ica_columnindex 10
#define rate_concentration _ml->template fpfield<11>(_iml)
#define rate_concentration_columnindex 11
#define Dconcentration _ml->template fpfield<12>(_iml)
#define Dconcentration_columnindex 12
#define DextConcentration _ml->template fpfield<13>(_iml)
#define DextConcentration_columnindex 13
#define _g _ml->template fpfield<14>(_iml)
#define _g_columnindex 14
#define _ion_cai *(_ml->dptr_field<0>(_iml))
#define _p_ion_cai static_cast<neuron::container::data_handle<double>>(_ppvar[0])
#define _ion_cao *(_ml->dptr_field<1>(_iml))
#define _p_ion_cao static_cast<neuron::container::data_handle<double>>(_ppvar[1])
#define _ion_ica *(_ml->dptr_field<2>(_iml))
#define _p_ion_ica static_cast<neuron::container::data_handle<double>>(_ppvar[2])
#define _ion_ca_erev *_ml->dptr_field<3>(_iml)
#define _style_ca	*_ppvar[4].get<int*>()
#define diam	(*(_ml->dptr_field<5>(_iml)))
#define area	(*(_ml->dptr_field<6>(_iml)))
 static _nrn_mechanism_cache_instance _ml_real{nullptr};
static _nrn_mechanism_cache_range *_ml{&_ml_real};
static size_t _iml{0};
static Datum *_ppvar;
 static int hoc_nrnpointerindex =  -1;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
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
 neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
_ppvar = _nrn_mechanism_access_dparam(_prop);
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", _hoc_setdata},
 {"rates_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", _hoc_rates},
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
#define iCa iCa_CaDynamics_E2_NML2__decay122__gamma5_09Emin4
 double iCa = 0;
#define initialExtConcentration initialExtConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4
 double initialExtConcentration = 0;
#define initialConcentration initialConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4
 double initialConcentration = 0;
#define surfaceArea surfaceArea_CaDynamics_E2_NML2__decay122__gamma5_09Emin4
 double surfaceArea = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"surfaceArea_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "um2"},
 {"iCa_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "nA"},
 {"initialConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "mM"},
 {"initialExtConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "mM"},
 {"minCai_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "mM"},
 {"decay_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "ms"},
 {"depth_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "um"},
 {"Faraday_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "C"},
 {"concentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "mM"},
 {"extConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "mM"},
 {"currDensCa_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", "nA / um2"},
 {0, 0}
};
 static double concentration0 = 0;
 static double delta_t = 0.01;
 static double extConcentration0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"surfaceArea_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", &surfaceArea_CaDynamics_E2_NML2__decay122__gamma5_09Emin4},
 {"iCa_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", &iCa_CaDynamics_E2_NML2__decay122__gamma5_09Emin4},
 {"initialConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", &initialConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4},
 {"initialExtConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4", &initialExtConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4},
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
 
#define _cvode_ieq _ppvar[7].literal_value<int>()
 static void _ode_synonym(_nrn_model_sorted_token const&, NrnThread&, Memb_list&, int);
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 "gamma_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 "minCai_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 "decay_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 "depth_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 "Faraday_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 0,
 "currDensCa_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 0,
 "concentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 "extConcentration_CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
 0,
 0};
 static Symbol* _morphology_sym;
 extern Node* nrn_alloc_node_;
 static Symbol* _ca_sym;
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     0.000509, /* gamma */
     0.0001, /* minCai */
     122, /* decay */
     0.1, /* depth */
     0.0964853, /* Faraday */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
   _ppvar = nrn_prop_datum_alloc(_mechtype, 8, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 15);
 	/*initialize range parameters*/
 	gamma = _parm_default[0]; /* 0.000509 */
 	minCai = _parm_default[1]; /* 0.0001 */
 	decay = _parm_default[2]; /* 122 */
 	depth = _parm_default[3]; /* 0.1 */
 	Faraday = _parm_default[4]; /* 0.0964853 */
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 15);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_morphology_sym);
 	_ppvar[5] = _nrn_mechanism_get_param_handle(prop_ion, 0); /* diam */
 	_ppvar[6] = _nrn_mechanism_get_area_handle(nrn_alloc_node_);
 prop_ion = need_memb(_ca_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0] = _nrn_mechanism_get_param_handle(prop_ion, 1); /* cai */
 	_ppvar[1] = _nrn_mechanism_get_param_handle(prop_ion, 2); /* cao */
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ica */
 	_ppvar[3] = _nrn_mechanism_get_param_handle(prop_ion, 0); // erev ca
 	_ppvar[4] = {neuron::container::do_not_search, &(_nrn_mechanism_access_dparam(prop_ion)[0].literal_value<int>())}; /* iontype for ca */
 
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

 extern "C" void _CaDynamics_E2_NML2__decay122__gamma5_09Emin4_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", 2.0);
 	_morphology_sym = hoc_lookup("morphology");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"gamma"} /* 0 */,
                                       _nrn_mechanism_field<double>{"minCai"} /* 1 */,
                                       _nrn_mechanism_field<double>{"decay"} /* 2 */,
                                       _nrn_mechanism_field<double>{"depth"} /* 3 */,
                                       _nrn_mechanism_field<double>{"Faraday"} /* 4 */,
                                       _nrn_mechanism_field<double>{"currDensCa"} /* 5 */,
                                       _nrn_mechanism_field<double>{"concentration"} /* 6 */,
                                       _nrn_mechanism_field<double>{"extConcentration"} /* 7 */,
                                       _nrn_mechanism_field<double>{"cai"} /* 8 */,
                                       _nrn_mechanism_field<double>{"cao"} /* 9 */,
                                       _nrn_mechanism_field<double>{"ica"} /* 10 */,
                                       _nrn_mechanism_field<double>{"rate_concentration"} /* 11 */,
                                       _nrn_mechanism_field<double>{"Dconcentration"} /* 12 */,
                                       _nrn_mechanism_field<double>{"DextConcentration"} /* 13 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 14 */,
                                       _nrn_mechanism_field<double*>{"_ion_cai", "ca_ion"} /* 0 */,
                                       _nrn_mechanism_field<double*>{"_ion_cao", "ca_ion"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_ica", "ca_ion"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"_ion_ca_erev", "ca_ion"} /* 3 */,
                                       _nrn_mechanism_field<int*>{"_style_ca", "#ca_ion"} /* 4 */,
                                       _nrn_mechanism_field<double*>{"diam", "diam"} /* 5 */,
                                       _nrn_mechanism_field<double*>{"area", "area"} /* 6 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 7 */);
  hoc_register_prop_size(_mechtype, 15, 8);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "#ca_ion");
  hoc_register_dparam_semantics(_mechtype, 7, "cvodeieq");
  hoc_register_dparam_semantics(_mechtype, 5, "diam");
  hoc_register_dparam_semantics(_mechtype, 6, "area");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_synonym(_mechtype, _ode_synonym);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 CaDynamics_E2_NML2__decay122__gamma5_09Emin4 /home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/CaDynamics_E2_NML2__decay122__gamma5_09Emin4.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "Mod file for component: Component(id=CaDynamics_E2_NML2__decay122__gamma5_09Emin4 type=concentrationModelHayEtAl)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates();
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[1], _dlist1[1];
 static int states(_internalthreadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargs_ ) ;
   Dconcentration = rate_concentration ;
   cai = concentration ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargs_ ) ;
 Dconcentration = Dconcentration  / (1. - dt*( 0.0 )) ;
 cai = concentration ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargs_ ) ;
    concentration = concentration - dt*(- ( rate_concentration ) ) ;
   cai = concentration ;
   }
  return 0;
}
 
static int  rates (  ) {
   surfaceArea = area ;
   iCa = - 1.0 * ( 0.01 ) * ica * surfaceArea ;
   currDensCa = iCa / surfaceArea ;
   rate_concentration = ( currDensCa * gamma / ( 2.0 * Faraday * depth ) ) - ( ( concentration - minCai ) / decay ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
  
  if(!_prop_id) {
    hoc_execerror("No data for rates_CaDynamics_E2_NML2__decay122__gamma5_09Emin4. Requires prior call to setdata_CaDynamics_E2_NML2__decay122__gamma5_09Emin4 and that the specified mechanism instance still be in existence.", NULL);
  } else {
    _setdata(_extcall_prop);
  }
   _r = 1.;
 rates (  );
 hoc_retpushx(_r);
}
 
static double _npy_rates(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r = 1.;
 rates (  );
 return(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
      Node* _nd{};
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
  ica = _ion_ica;
  cai = _ion_cai;
     _ode_spec1 ();
  _ion_cai = cai;
 }}
 
static void _ode_map(Prop* _prop, int _ieq, neuron::container::data_handle<double>* _pv, neuron::container::data_handle<double>* _pvdot, double* _atol, int _type) { 
  _ppvar = _nrn_mechanism_access_dparam(_prop);
  _cvode_ieq = _ieq;
  for (int _i=0; _i < 1; ++_i) {
    _pv[_i] = _nrn_mechanism_get_param_handle(_prop, _slist1[_i]);
    _pvdot[_i] = _nrn_mechanism_get_param_handle(_prop, _dlist1[_i]);
    _cvode_abstol(_atollist, _atol, _i);
  }
 }
 static void _ode_synonym(_nrn_model_sorted_token const& _sorted_token, NrnThread& _nt, Memb_list& _ml_arg, int _type) {
 _nrn_mechanism_cache_range _lmr{_sorted_token, _nt, _ml_arg, _type};
auto* const _ml = &_lmr;
auto const _cnt = _ml_arg._nodecount;
for (int _iml = 0; _iml < _cnt; ++_iml) {
  Datum* _ppvar = _ml_arg._pdata[_iml];
 _ion_cai =  concentration ;
   }
}
 
static void _ode_matsol_instance1(_internalthreadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
      Node* _nd{};
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
  ica = _ion_ica;
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  concentration = concentration0;
  extConcentration = extConcentration0;
 {
   initialConcentration = cai ;
   initialExtConcentration = cao ;
   rates ( _threadargs_ ) ;
   rates ( _threadargs_ ) ;
   concentration = initialConcentration ;
   extConcentration = initialExtConcentration ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 v = _v;
  cai = _ion_cai;
  cao = _ion_cao;
  ica = _ion_ica;
  cai = _ion_cai;
 initmodel();
  _ion_cai = cai;
  nrn_wrote_conc(_ca_sym, _ion_ca_erev, _ion_cai, _ion_cao, _style_ca);
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_rhs = _nt->node_rhs_storage();
auto const _vec_sav_rhs = _nt->node_sav_rhs_storage();
auto const _vec_v = _nt->node_voltage_storage();
Node *_nd; int* _ni; double _rhs, _v; int _cntml;
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 
}}

static void nrn_jacob(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type) {
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto const _vec_d = _nt->node_d_storage();
auto const _vec_sav_d = _nt->node_sav_d_storage();
auto* const _ml = &_lmr;
Node *_nd; int* _ni; int _iml, _cntml;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
  _vec_d[_ni[_iml]] += _g;
 
}}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
  cai = _ion_cai;
  cao = _ion_cao;
  ica = _ion_ica;
  cai = _ion_cai;
 { error =  states();
 if(error){
  std_cerr_stream << "at line 90 in file CaDynamics_E2_NML2__decay122__gamma5_09Emin4.mod:\n    \n";
  std_cerr_stream << _ml << ' ' << _iml << '\n';
  abort_run(error);
}
 } {
   }
  _ion_cai = cai;
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {concentration_columnindex, 0};  _dlist1[0] = {Dconcentration_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/CaDynamics_E2_NML2__decay122__gamma5_09Emin4.mod";
    const char* nmodl_file_text = 
  "TITLE Mod file for component: Component(id=CaDynamics_E2_NML2__decay122__gamma5_09Emin4 type=concentrationModelHayEtAl)\n"
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
  "    SUFFIX CaDynamics_E2_NML2__decay122__gamma5_09Emin4\n"
  "    USEION ca READ cai, cao, ica WRITE cai VALENCE 2\n"
  "    RANGE cai\n"
  "    RANGE cao\n"
  "    GLOBAL initialConcentration\n"
  "    GLOBAL initialExtConcentration\n"
  "    RANGE gamma                             : parameter\n"
  "    RANGE minCai                            : parameter\n"
  "    RANGE decay                             : parameter\n"
  "    RANGE depth                             : parameter\n"
  "    RANGE Faraday                           : parameter\n"
  "    RANGE currDensCa                        : derived variable\n"
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
  "    surfaceArea (um2)\n"
  "    iCa (nA)\n"
  "    initialConcentration (mM)\n"
  "    initialExtConcentration (mM)\n"
  "    \n"
  "    gamma = 5.09E-4 \n"
  "    minCai = 1.0E-4 (mM)\n"
  "    decay = 122 (ms)\n"
  "    depth = 0.1 (um)\n"
  "    Faraday = 0.0964853 (C / umol)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    cai (mM)\n"
  "    cao (mM)\n"
  "    ica (mA/cm2)\n"
  "    diam (um)\n"
  "    area (um2)\n"
  "    \n"
  "    currDensCa (nA / um2)                  : derived variable\n"
  "    rate_concentration (mM/ms)\n"
  "    \n"
  "}\n"
  "\n"
  "STATE {\n"
  "    concentration (mM) \n"
  "    extConcentration (mM) \n"
  "    \n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    initialConcentration = cai\n"
  "    initialExtConcentration = cao\n"
  "    rates()\n"
  "    rates() ? To ensure correct initialisation.\n"
  "    \n"
  "    concentration = initialConcentration\n"
  "    \n"
  "    extConcentration = initialExtConcentration\n"
  "    \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    \n"
  "    SOLVE states METHOD cnexp\n"
  "    \n"
  "    \n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "    rates()\n"
  "    concentration' = rate_concentration\n"
  "    cai = concentration \n"
  "    \n"
  "}\n"
  "\n"
  "PROCEDURE rates() {\n"
  "    \n"
  "    surfaceArea = area   : surfaceArea has units (um2), area (built in to NEURON) is in um^2...\n"
  "    \n"
  "    iCa = -1 * (0.01) * ica * surfaceArea :   iCa has units (nA) ; ica (built in to NEURON) has units (mA/cm2)...\n"
  "    \n"
  "    currDensCa = iCa / surfaceArea ? evaluable\n"
  "    rate_concentration = (  currDensCa   *   gamma  /(2 *  Faraday  *   depth  )) - ((  concentration   -   minCai  ) /   decay  ) ? Note units of all quantities used here need to be consistent!\n"
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
