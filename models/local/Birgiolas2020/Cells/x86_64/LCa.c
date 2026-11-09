/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__LCa
#define _nrn_initial _nrn_initial__LCa
#define nrn_cur _nrn_cur__LCa
#define _nrn_current _nrn_current__LCa
#define nrn_jacob _nrn_jacob__LCa
#define nrn_state _nrn_state__LCa
#define _net_receive _net_receive__LCa 
#define _f_rates _f_rates__LCa 
#define rates rates__LCa 
#define states states__LCa 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gbar _p[0]
#define ica _p[1]
#define r _p[2]
#define s _p[3]
#define Dr _p[4]
#define Ds _p[5]
#define qt _p[6]
#define _g _p[7]
#define _ion_ica	*_ppvar[0]._pval
#define _ion_dicadv	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alp(void);
 static void _hoc_bet(void);
 static void _hoc_rates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_LCa", _hoc_setdata,
 "alp_LCa", _hoc_alp,
 "bet_LCa", _hoc_bet,
 "rates_LCa", _hoc_rates,
 0, 0
};
#define alp alp_LCa
#define bet bet_LCa
 extern double alp( double , double );
 extern double bet( double , double );
 /* declare global and static user variables */
#define q10 q10_LCa
 double q10 = 3;
#define rtau rtau_LCa
 double rtau = 0;
#define rinf rinf_LCa
 double rinf = 0;
#define stau stau_LCa
 double stau = 0;
#define sinf sinf_LCa
 double sinf = 0;
#define usetable usetable_LCa
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gbar_LCa", 0, 1e+09,
 "usetable_LCa", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "stau_LCa", "ms",
 "rtau_LCa", "ms",
 "gbar_LCa", "mho/cm2",
 "ica_LCa", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double r0 = 0;
 static double s0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "q10_LCa", &q10_LCa,
 "sinf_LCa", &sinf_LCa,
 "rinf_LCa", &rinf_LCa,
 "stau_LCa", &stau_LCa,
 "rtau_LCa", &rtau_LCa,
 "usetable_LCa", &usetable_LCa,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"LCa",
 "gbar_LCa",
 0,
 "ica_LCa",
 0,
 "r_LCa",
 "s_LCa",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 8, _prop);
 	/*initialize range parameters*/
 	gbar = 0.12;
 	_prop->param = _p;
 	_prop->param_size = 8;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _LCa_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 8, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 LCa /home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/LCa.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double eca = 70;
 static double *_t_sinf;
 static double *_t_rinf;
 static double *_t_stau;
 static double *_t_rtau;
static int _reset;
static char *modelname = "LCa calcium channel with fixed reversal potential";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(double);
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(double);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Ds = ( sinf - s ) / stau ;
   Dr = ( rinf - r ) / rtau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / stau )) ;
 Dr = Dr  / (1. - dt*( ( ( ( - 1.0 ) ) ) / rtau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / stau)))*(- ( ( ( sinf ) ) / stau ) / ( ( ( ( - 1.0 ) ) ) / stau ) - s) ;
    r = r + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / rtau)))*(- ( ( ( rinf ) ) / rtau ) / ( ( ( ( - 1.0 ) ) ) / rtau ) - r) ;
   }
  return 0;
}
 
double alp (  double _lv , double _li ) {
   double _lalp;
 if ( _li  == 0.0 ) {
     _lalp = 7.5 / ( 1.0 + exp ( ( - _lv * 1.0 + 13.0 ) / 7.0 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lalp = 0.0068 / ( 1.0 + exp ( ( _lv * 1.0 + 30.0 ) / 12.0 ) ) ;
     }
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double bet (  double _lv , double _li ) {
   double _lbet;
 if ( _li  == 0.0 ) {
     _lbet = 1.65 / ( 1.0 + exp ( ( _lv * 1.0 - 14.0 ) / 4.0 ) ) ;
     }
   else if ( _li  == 1.0 ) {
     _lbet = 0.06 / ( 1.0 + exp ( - _lv * 1.0 / 11.0 ) ) ;
     }
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
   _r =  bet (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 static double _mfac_rates, _tmin_rates;
 static void _check_rates();
 static void _check_rates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_x);
    _t_sinf[_i] = sinf;
    _t_rinf[_i] = rinf;
    _t_stau[_i] = stau;
    _t_rtau[_i] = rtau;
   }
  }
 }

 static int rates(double _lv){ _check_rates();
 _n_rates(_lv);
 return 0;
 }

 static void _n_rates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  sinf = _xi;
  rinf = _xi;
  stau = _xi;
  rtau = _xi;
  return;
 }
 if (_xi <= 0.) {
 sinf = _t_sinf[0];
 rinf = _t_rinf[0];
 stau = _t_stau[0];
 rtau = _t_rtau[0];
 return; }
 if (_xi >= 200.) {
 sinf = _t_sinf[200];
 rinf = _t_rinf[200];
 stau = _t_stau[200];
 rtau = _t_rtau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 sinf = _t_sinf[_i] + _theta*(_t_sinf[_i+1] - _t_sinf[_i]);
 rinf = _t_rinf[_i] + _theta*(_t_rinf[_i+1] - _t_rinf[_i]);
 stau = _t_stau[_i] + _theta*(_t_stau[_i+1] - _t_stau[_i]);
 rtau = _t_rtau[_i] + _theta*(_t_rtau[_i+1] - _t_rtau[_i]);
 }

 
static int  _f_rates (  double _lv ) {
   double _la , _lb ;
 _la = alp ( _threadargscomma_ _lv , 0.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 0.0 ) ;
   stau = 1.0 / ( _la + _lb ) ;
   stau = stau / qt ;
   sinf = _la / ( _la + _lb ) ;
   _la = alp ( _threadargscomma_ _lv , 1.0 ) ;
   _lb = bet ( _threadargscomma_ _lv , 1.0 ) ;
   rtau = 1.0 / ( _la + _lb ) ;
   rtau = rtau / qt ;
   rinf = _la / ( _la + _lb ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
    _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  r = r0;
  s = s0;
 {
   qt = pow( q10 , ( ( celsius - 35.0 ) / 10.0 ) ) ;
   rates ( _threadargscomma_ v ) ;
   s = sinf ;
   r = rinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ica = gbar * s * r * ( v - eca ) ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  states();
 if(error){fprintf(stderr,"at line 58 in file LCa.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(s) - _p;  _dlist1[0] = &(Ds) - _p;
 _slist1[1] = &(r) - _p;  _dlist1[1] = &(Dr) - _p;
   _t_sinf = makevector(201*sizeof(double));
   _t_rinf = makevector(201*sizeof(double));
   _t_stau = makevector(201*sizeof(double));
   _t_rtau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/LCa.mod";
static const char* nmodl_file_text = 
  "TITLE LCa calcium channel with fixed reversal potential\n"
  ": Implemented in Rubin and Cleland, J. Neurophysiol 2006\n"
  ": LCa channel with parameters from US Bhalla and JM Bower,\n"
  ": J. Neurophysiol. 69:1948-1983 (1993)\n"
  ": Adapted from /usr/local/neuron/demo/release/nachan.mod - squid\n"
  ": by Andrew Davison, The Babraham Institute.\n"
  ": 25-08-98\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX LCa\n"
  "	USEION ca WRITE ica\n"
  "	RANGE gbar, ica\n"
  "	GLOBAL sinf, rinf, stau, rtau\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "}\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "CONSTANT { eca = 70 (mV) }\n"
  "\n"
  "PARAMETER {\n"
  "	v (mV)\n"
  "	dt (ms)\n"
  "	gbar	= 0.120 (mho/cm2) <0,1e9>\n"
  ":	eca = 70 (mV)\n"
  "	celsius   (degC)\n"
  "    q10 = 3\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	r s\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	ica (mA/cm2)\n"
  "	sinf\n"
  "	rinf\n"
  "	stau (ms)\n"
  "	rtau (ms)\n"
  "    qt\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    qt=q10^((celsius-35(degC))/10(degC))\n"
  "	rates(v)\n"
  "	s = sinf\n"
  "	r = rinf\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ica = gbar*s*r*(v - eca)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	rates(v)\n"
  "	s' = (sinf - s)/stau\n"
  "	r' = (rinf - r)/rtau\n"
  "}\n"
  "\n"
  "FUNCTION alp(v(mV),i) (/ms) {\n"
  "	if (i==0) {\n"
  "		alp = 7.5(/ms)/(1 + exp((-v *1(/mV) + 13)/7))\n"
  "	}else if (i==1){\n"
  "		alp = 0.0068(/ms)/(1 + exp((v *1(/mV) + 30)/12))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION bet(v(mV),i)(/ms) {\n"
  "	if (i==0) {\n"
  "		bet = 1.65(/ms)/(1 + exp((v *1(/mV) - 14)/4))\n"
  "	}else if (i==1){\n"
  "		bet = 0.06(/ms)/(1 + exp(-v* 1(/mV)/11))\n"
  "	}\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v(mV)) {LOCAL a, b\n"
  "	TABLE sinf, rinf, stau, rtau FROM -100 TO 100 WITH 200\n"
  "	a = alp(v,0)  b=bet(v,0)\n"
  "	stau = 1/(a + b)\n"
  "	stau = stau / qt\n"
  "	sinf = a/(a + b)\n"
  "	a = alp(v,1)  b=bet(v,1)\n"
  "	rtau = 1/(a + b)\n"
  "	rtau = rtau / qt\n"
  "	rinf = a/(a + b)\n"
  "}\n"
  "\n"
  ;
#endif
