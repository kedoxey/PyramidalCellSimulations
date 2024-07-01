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
static constexpr auto number_of_datum_variables = 6;
static constexpr auto number_of_floating_point_variables = 22;
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
 
#define nrn_init _nrn_init__MyExp2SynNMDABB
#define _nrn_initial _nrn_initial__MyExp2SynNMDABB
#define nrn_cur _nrn_cur__MyExp2SynNMDABB
#define _nrn_current _nrn_current__MyExp2SynNMDABB
#define nrn_jacob _nrn_jacob__MyExp2SynNMDABB
#define nrn_state _nrn_state__MyExp2SynNMDABB
#define _net_receive _net_receive__MyExp2SynNMDABB 
#define state state__MyExp2SynNMDABB 
 
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
#define tau1NMDA _ml->template fpfield<0>(_iml)
#define tau1NMDA_columnindex 0
#define tau2NMDA _ml->template fpfield<1>(_iml)
#define tau2NMDA_columnindex 1
#define e _ml->template fpfield<2>(_iml)
#define e_columnindex 2
#define r _ml->template fpfield<3>(_iml)
#define r_columnindex 3
#define smax _ml->template fpfield<4>(_iml)
#define smax_columnindex 4
#define sNMDAmax _ml->template fpfield<5>(_iml)
#define sNMDAmax_columnindex 5
#define Vwt _ml->template fpfield<6>(_iml)
#define Vwt_columnindex 6
#define iNMDA _ml->template fpfield<7>(_iml)
#define iNMDA_columnindex 7
#define sNMDA _ml->template fpfield<8>(_iml)
#define sNMDA_columnindex 8
#define ica _ml->template fpfield<9>(_iml)
#define ica_columnindex 9
#define g _ml->template fpfield<10>(_iml)
#define g_columnindex 10
#define A2 _ml->template fpfield<11>(_iml)
#define A2_columnindex 11
#define B2 _ml->template fpfield<12>(_iml)
#define B2_columnindex 12
#define mgblock _ml->template fpfield<13>(_iml)
#define mgblock_columnindex 13
#define factor2 _ml->template fpfield<14>(_iml)
#define factor2_columnindex 14
#define cai _ml->template fpfield<15>(_iml)
#define cai_columnindex 15
#define cao _ml->template fpfield<16>(_iml)
#define cao_columnindex 16
#define DA2 _ml->template fpfield<17>(_iml)
#define DA2_columnindex 17
#define DB2 _ml->template fpfield<18>(_iml)
#define DB2_columnindex 18
#define v _ml->template fpfield<19>(_iml)
#define v_columnindex 19
#define _g _ml->template fpfield<20>(_iml)
#define _g_columnindex 20
#define _tsav _ml->template fpfield<21>(_iml)
#define _tsav_columnindex 21
#define _nd_area *_ml->dptr_field<0>(_iml)
#define _ion_cai *(_ml->dptr_field<2>(_iml))
#define _p_ion_cai static_cast<neuron::container::data_handle<double>>(_ppvar[2])
#define _ion_cao *(_ml->dptr_field<3>(_iml))
#define _p_ion_cao static_cast<neuron::container::data_handle<double>>(_ppvar[3])
#define _ion_ica *(_ml->dptr_field<4>(_iml))
#define _p_ion_ica static_cast<neuron::container::data_handle<double>>(_ppvar[4])
#define _ion_dicadv *(_ml->dptr_field<5>(_iml))
 /* Thread safe. No static _ml, _iml or _ppvar. */
 static int hoc_nrnpointerindex =  -1;
 static _nrn_mechanism_std_vector<Datum> _extcall_thread;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static double _hoc_ghk(void*);
 static double _hoc_ghkg(void*);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {0, 0}
};
 static Member_func _member_func[] = {
 {"loc", _hoc_loc_pnt},
 {"has_loc", _hoc_has_loc},
 {"get_loc", _hoc_get_loc_pnt},
 {"ghk", _hoc_ghk},
 {"ghkg", _hoc_ghkg},
 {0, 0}
};
#define ghk ghk_MyExp2SynNMDABB
#define ghkg ghkg_MyExp2SynNMDABB
 extern double ghk( _internalthreadargsprotocomma_ double , double , double , double );
 extern double ghkg( _internalthreadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
#define fracca fracca_MyExp2SynNMDABB
 double fracca = 0.13;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"tau1NMDA", "ms"},
 {"tau2NMDA", "ms"},
 {"e", "mV"},
 {"smax", "1"},
 {"sNMDAmax", "1"},
 {"A2", "1"},
 {"B2", "1"},
 {"iNMDA", "nA"},
 {"sNMDA", "1"},
 {"ica", "nA"},
 {"g", "umho"},
 {0, 0}
};
 static double A20 = 0;
 static double B20 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"fracca_MyExp2SynNMDABB", &fracca_MyExp2SynNMDABB},
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
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(Prop*, int, neuron::container::data_handle<double>*, neuron::container::data_handle<double>*, double*, int);
static void _ode_spec(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void _ode_matsol(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 
#define _cvode_ieq _ppvar[6].literal_value<int>()
 static void _ode_matsol_instance1(_internalthreadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"MyExp2SynNMDABB",
 "tau1NMDA",
 "tau2NMDA",
 "e",
 "r",
 "smax",
 "sNMDAmax",
 "Vwt",
 0,
 "iNMDA",
 "sNMDA",
 "ica",
 "g",
 0,
 "A2",
 "B2",
 0,
 0};
 static Symbol* _ca_sym;
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     15, /* tau1NMDA */
     150, /* tau2NMDA */
     0, /* e */
     1, /* r */
     1e+09, /* smax */
     1e+09, /* sNMDAmax */
     0, /* Vwt */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
  if (nrn_point_prop_) {
    _nrn_mechanism_access_alloc_seq(_prop) = _nrn_mechanism_access_alloc_seq(nrn_point_prop_);
    _ppvar = _nrn_mechanism_access_dparam(nrn_point_prop_);
  } else {
   _ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 22);
 	/*initialize range parameters*/
 	tau1NMDA = _parm_default[0]; /* 15 */
 	tau2NMDA = _parm_default[1]; /* 150 */
 	e = _parm_default[2]; /* 0 */
 	r = _parm_default[3]; /* 1 */
 	smax = _parm_default[4]; /* 1e+09 */
 	sNMDAmax = _parm_default[5]; /* 1e+09 */
 	Vwt = _parm_default[6]; /* 0 */
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 22);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[2] = _nrn_mechanism_get_param_handle(prop_ion, 1); /* cai */
 	_ppvar[3] = _nrn_mechanism_get_param_handle(prop_ion, 2); /* cao */
 	_ppvar[4] = _nrn_mechanism_get_param_handle(prop_ion, 3); /* ica */
 	_ppvar[5] = _nrn_mechanism_get_param_handle(prop_ion, 4); /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 {0, 0}
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _MyExp2SynNMDABB_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"tau1NMDA"} /* 0 */,
                                       _nrn_mechanism_field<double>{"tau2NMDA"} /* 1 */,
                                       _nrn_mechanism_field<double>{"e"} /* 2 */,
                                       _nrn_mechanism_field<double>{"r"} /* 3 */,
                                       _nrn_mechanism_field<double>{"smax"} /* 4 */,
                                       _nrn_mechanism_field<double>{"sNMDAmax"} /* 5 */,
                                       _nrn_mechanism_field<double>{"Vwt"} /* 6 */,
                                       _nrn_mechanism_field<double>{"iNMDA"} /* 7 */,
                                       _nrn_mechanism_field<double>{"sNMDA"} /* 8 */,
                                       _nrn_mechanism_field<double>{"ica"} /* 9 */,
                                       _nrn_mechanism_field<double>{"g"} /* 10 */,
                                       _nrn_mechanism_field<double>{"A2"} /* 11 */,
                                       _nrn_mechanism_field<double>{"B2"} /* 12 */,
                                       _nrn_mechanism_field<double>{"mgblock"} /* 13 */,
                                       _nrn_mechanism_field<double>{"factor2"} /* 14 */,
                                       _nrn_mechanism_field<double>{"cai"} /* 15 */,
                                       _nrn_mechanism_field<double>{"cao"} /* 16 */,
                                       _nrn_mechanism_field<double>{"DA2"} /* 17 */,
                                       _nrn_mechanism_field<double>{"DB2"} /* 18 */,
                                       _nrn_mechanism_field<double>{"v"} /* 19 */,
                                       _nrn_mechanism_field<double>{"_g"} /* 20 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 21 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<double*>{"_ion_cai", "ca_ion"} /* 2 */,
                                       _nrn_mechanism_field<double*>{"_ion_cao", "ca_ion"} /* 3 */,
                                       _nrn_mechanism_field<double*>{"_ion_ica", "ca_ion"} /* 4 */,
                                       _nrn_mechanism_field<double*>{"_ion_dicadv", "ca_ion"} /* 5 */,
                                       _nrn_mechanism_field<int>{"_cvode_ieq", "cvodeieq"} /* 6 */);
  hoc_register_prop_size(_mechtype, 22, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 MyExp2SynNMDABB /home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/MyExp2SynNMDABB.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 0x1.78e555060882cp+16;
 static double R = 0x1.0a1013e8990bep+3;
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_internalthreadargsproto_);
/*static int _ode_matsol1(_internalthreadargsproto_);*/
 static neuron::container::field_index _slist1[2], _dlist1[2];
 static int state(_internalthreadargsproto_);
 
double ghkg ( _internalthreadargsprotocomma_ double _lv , double _lci , double _lco , double _lz ) {
   double _lghkg;
 double _lxi , _lf , _lexi , _lfxi ;
 _lf = R * ( celsius + 273.15 ) / ( _lz * ( 1e-3 ) * FARADAY ) ;
   _lxi = _lv / _lf ;
   _lexi = exp ( _lxi ) ;
   if ( fabs ( _lxi ) < 1e-4 ) {
     _lfxi = 1.0 - _lxi / 2.0 ;
     }
   else {
     _lfxi = _lxi / ( _lexi - 1.0 ) ;
     }
   _lghkg = _lf * ( ( _lci / _lco ) * _lexi - 1.0 ) * _lfxi ;
   
return _lghkg;
 }
 
static double _hoc_ghkg(void* _vptr) {
 double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _nrn_mechanism_cache_instance _ml_real{_p};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  _ppvar = _nrn_mechanism_access_dparam(_p);
  _thread = _extcall_thread.data();
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  _nt = static_cast<NrnThread*>(_pnt->_vnt);
 _r =  ghkg ( _threadargscomma_ *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 return(_r);
}
 
double ghk ( _internalthreadargsprotocomma_ double _lv , double _lci , double _lco , double _lz ) {
   double _lghk;
 double _lxi , _lf , _lexi , _lfxi ;
 _lf = R * ( celsius + 273.15 ) / ( _lz * ( 1e-3 ) * FARADAY ) ;
   _lxi = _lv / _lf ;
   _lexi = exp ( _lxi ) ;
   if ( fabs ( _lxi ) < 1e-4 ) {
     _lfxi = 1.0 - _lxi / 2.0 ;
     }
   else {
     _lfxi = _lxi / ( _lexi - 1.0 ) ;
     }
   _lghk = ( .001 ) * _lz * FARADAY * ( _lci * _lexi - _lco ) * _lfxi ;
   
return _lghk;
 }
 
static double _hoc_ghk(void* _vptr) {
 double _r;
 Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _nrn_mechanism_cache_instance _ml_real{_p};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
  _ppvar = _nrn_mechanism_access_dparam(_p);
  _thread = _extcall_thread.data();
  double* _globals = nullptr;
  if (gind != 0 && _thread != nullptr) { _globals = _thread[_gth].get<double*>(); }
  _nt = static_cast<NrnThread*>(_pnt->_vnt);
 _r =  ghk ( _threadargscomma_ *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 return(_r);
}
 
/*CVODE*/
 static int _ode_spec1 (_internalthreadargsproto_) {int _reset = 0; {
   DA2 = - A2 / tau1NMDA ;
   DB2 = - B2 / tau2NMDA ;
   }
 return _reset;
}
 static int _ode_matsol1 (_internalthreadargsproto_) {
 DA2 = DA2  / (1. - dt*( ( - 1.0 ) / tau1NMDA )) ;
 DB2 = DB2  / (1. - dt*( ( - 1.0 ) / tau2NMDA )) ;
  return 0;
}
 /*END CVODE*/
 static int state (_internalthreadargsproto_) { {
    A2 = A2 + (1. - exp(dt*(( - 1.0 ) / tau1NMDA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1NMDA ) - A2) ;
    B2 = B2 + (1. - exp(dt*(( - 1.0 ) / tau2NMDA)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2NMDA ) - B2) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  Prop* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _nrn_mechanism_cache_instance _ml_real{_pnt->_prop};
  auto* const _ml = &_ml_real;
  size_t const _iml{};
   _thread = nullptr; double* _globals = nullptr; _nt = (NrnThread*)_pnt->_vnt;   _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   double _lww ;
 _lww = _args[0] ;
   if ( r >= 0.0 ) {
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A2;
    double __primary = (A2 + factor2 * _lww * r) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1NMDA ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1NMDA ) - __primary );
    A2 += __primary;
  } else {
 A2 = A2 + factor2 * _lww * r ;
       }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B2;
    double __primary = (B2 + factor2 * _lww * r) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2NMDA ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2NMDA ) - __primary );
    B2 += __primary;
  } else {
 B2 = B2 + factor2 * _lww * r ;
       }
 }
   } }
 
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
  cai = _ion_cai;
  cao = _ion_cao;
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
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(_internalthreadargsproto_) {
  int _i; double _save;{
  A2 = A20;
  B2 = B20;
 {
   double _ltp ;
 Vwt = 0.0 ;
   if ( tau1NMDA / tau2NMDA > .9999 ) {
     tau1NMDA = .9999 * tau2NMDA ;
     }
   A2 = 0.0 ;
   B2 = 0.0 ;
   _ltp = ( tau1NMDA * tau2NMDA ) / ( tau2NMDA - tau1NMDA ) * log ( tau2NMDA / tau1NMDA ) ;
   factor2 = - exp ( - _ltp / tau1NMDA ) + exp ( - _ltp / tau2NMDA ) ;
   factor2 = 1.0 / factor2 ;
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
 _tsav = -1e20;
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
   double _liTOT ;
 mgblock = 1.0 / ( 1.0 + 0.28 * exp ( - 0.062 * v ) ) ;
   sNMDA = B2 - A2 ;
   if ( sNMDA > sNMDAmax ) {
     sNMDA = sNMDAmax ;
     }
   iNMDA = sNMDA * ( v - e ) * mgblock * ( 1.0 - fracca ) ;
   if ( fracca > 0.0 ) {
     ica = sNMDA * ghkg ( _threadargscomma_ v , cai , cao , 2.0 ) * mgblock * fracca ;
     }
   g = sNMDA * mgblock ;
   }
 _current += iNMDA;
 _current += ica;

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
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_threadargscomma_ _v);
  _ion_dicadv += (_dica - ica)/.001 * 1.e2/ (_nd_area);
 	}
 _g = (_g_local - _rhs)/.001;
  _ion_ica += ica * 1.e2/ (_nd_area);
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
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
 {   state(_threadargs_);
  } }}

}

static void terminal(){}

static void _initlists(){
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = {A2_columnindex, 0};  _dlist1[0] = {DA2_columnindex, 0};
 _slist1[1] = {B2_columnindex, 0};  _dlist1[1] = {DB2_columnindex, 0};
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/MyExp2SynNMDABB.mod";
    const char* nmodl_file_text = 
  ": $Id: MyExp2SynNMDABB.mod,v 1.4 2010/12/13 21:28:02 samn Exp $ \n"
  "NEURON {\n"
  ":  THREADSAFE\n"
  "  POINT_PROCESS MyExp2SynNMDABB\n"
  "  RANGE e, i, iNMDA, s, sNMDA, r, tau1NMDA, tau2NMDA, Vwt, smax, sNMDAmax, g\n"
  "  NONSPECIFIC_CURRENT iNMDA\n"
  "  USEION ca READ cai,cao WRITE ica\n"
  "  GLOBAL fracca\n"
  "  RANGE ica\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "  (nA) = (nanoamp)\n"
  "  (mV) = (millivolt)\n"
  "  (uS) = (microsiemens)\n"
  "  FARADAY = (faraday) (coulomb)\n"
  "  R = (k-mole) (joule/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  tau1NMDA = 15  (ms)\n"
  "  tau2NMDA = 150 (ms)\n"
  "  e        = 0	(mV)\n"
  "  r        = 1\n"
  "  smax     = 1e9 (1)\n"
  "  sNMDAmax = 1e9 (1)  \n"
  "  Vwt   = 0 : weight for inputs coming in from vector\n"
  "  fracca = 0.13 : fraction of current that is ca ions; Srupuston &al 95\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "  v       (mV)\n"
  "  iNMDA   (nA)\n"
  "  sNMDA   (1)\n"
  "  mgblock (1)\n"
  "  factor2 (1)	\n"
  "  ica	  (nA)\n"
  "  cai     (mM)\n"
  "  cao     (mM)\n"
  "  g       (umho)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "  A2 (1)\n"
  "  B2 (1)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "  LOCAL tp\n"
  "  Vwt = 0 : testing\n"
  "  if (tau1NMDA/tau2NMDA > .9999) {\n"
  "    tau1NMDA = .9999*tau2NMDA\n"
  "  }\n"
  "  A2 = 0\n"
  "  B2 = 0	\n"
  "  tp = (tau1NMDA*tau2NMDA)/(tau2NMDA - tau1NMDA) * log(tau2NMDA/tau1NMDA)\n"
  "  factor2 = -exp(-tp/tau1NMDA) + exp(-tp/tau2NMDA)\n"
  "  factor2 = 1/factor2  \n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "  LOCAL iTOT\n"
  "  SOLVE state METHOD cnexp\n"
  "  : Jahr Stevens 1990 J. Neurosci\n"
  "  mgblock = 1.0 / (1.0 + 0.28 * exp(-0.062(/mV) * v) )\n"
  "  sNMDA = B2 - A2\n"
  "  if (sNMDA>sNMDAmax) {sNMDA=sNMDAmax}: saturation\n"
  "\n"
  "  :iTOT = sNMDA * (v - e) * mgblock  \n"
  "  :iNMDA = iTOT * (1-fracca)\n"
  "  :ica = iTOT * fracca\n"
  "  \n"
  "  iNMDA = sNMDA * (v - e) * mgblock * (1-fracca)\n"
  "  if(fracca>0.0){ica =   sNMDA * ghkg(v,cai,cao,2) * mgblock * fracca}\n"
  "  g = sNMDA * mgblock\n"
  "}\n"
  "\n"
  ":::INCLUDE \"ghk.inc\"\n"
  ":::realpath /home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/ghk.inc\n"
  "COMMENT\n"
  "    GHK function that returns effective driving force\n"
  "    Slope at low voltages is 1\n"
  "    z needs to be set as a PARAMETER\n"
  "ENDCOMMENT\n"
  "\n"
  "FUNCTION ghkg(v(mV), ci(mM), co(mM), z) (mV) {\n"
  "    LOCAL xi, f, exi, fxi\n"
  "    f = R*(celsius+273.15)/(z*(1e-3)*FARADAY)\n"
  "    xi = v/f\n"
  "    exi = exp(xi)\n"
  "    if (fabs(xi) < 1e-4) {\n"
  "        fxi = 1 - xi/2\n"
  "    }else{\n"
  "        fxi = xi/(exi - 1)\n"
  "    }\n"
  "    ghkg = f*((ci/co)*exi - 1)*fxi\n"
  "}\n"
  "\n"
  "FUNCTION ghk(v(mV), ci(mM), co(mM), z) (.001 coul/cm3) {\n"
  "    LOCAL xi, f, exi, fxi\n"
  "    f = R*(celsius+273.15)/(z*(1e-3)*FARADAY)\n"
  "    xi = v/f\n"
  "    exi = exp(xi)\n"
  "    if (fabs(xi) < 1e-4) {\n"
  "        fxi = 1 - xi/2\n"
  "    }else{\n"
  "        fxi = xi/(exi - 1)\n"
  "    }\n"
  "    ghk = (.001)*z*FARADAY*(ci*exi - co)*fxi\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ":::end INCLUDE ghk.inc\n"
  "\n"
  "DERIVATIVE state {\n"
  "  A2' = -A2/tau1NMDA\n"
  "  B2' = -B2/tau2NMDA\n"
  "}\n"
  "\n"
  "NET_RECEIVE(w (uS)) {LOCAL ww\n"
  "  ww=w\n"
  "  :printf(\"NMDA Spike: %g\\n\", t)\n"
  "  if(r>=0){ : if r>=0, g = NMDA*r\n"
  "    A2 = A2 + factor2*ww*r\n"
  "    B2 = B2 + factor2*ww*r\n"
  "  }\n"
  "}\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
