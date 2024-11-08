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
 
#define nrn_init _nrn_init__KA
#define _nrn_initial _nrn_initial__KA
#define nrn_cur _nrn_cur__KA
#define _nrn_current _nrn_current__KA
#define nrn_jacob _nrn_jacob__KA
#define nrn_state _nrn_state__KA
#define _net_receive _net_receive__KA 
#define states states__KA 
#define trates trates__KA 
 
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
#define sha _p[1]
#define shi _p[2]
#define k_tauH _p[3]
#define sh_tauH _p[4]
#define ik _p[5]
#define m _p[6]
#define h _p[7]
#define ek _p[8]
#define qt _p[9]
#define Dm _p[10]
#define Dh _p[11]
#define _g _p[12]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 static void _hoc_alph(void);
 static void _hoc_alpm(void);
 static void _hoc_beth(void);
 static void _hoc_betm(void);
 static void _hoc_trates(void);
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
 "setdata_KA", _hoc_setdata,
 "alph_KA", _hoc_alph,
 "alpm_KA", _hoc_alpm,
 "beth_KA", _hoc_beth,
 "betm_KA", _hoc_betm,
 "trates_KA", _hoc_trates,
 0, 0
};
#define alph alph_KA
#define alpm alpm_KA
#define beth beth_KA
#define betm betm_KA
 extern double alph( double );
 extern double alpm( double );
 extern double beth( double );
 extern double betm( double );
 /* declare global and static user variables */
#define a0h a0h_KA
 double a0h = 0.018;
#define a0m a0m_KA
 double a0m = 0.04;
#define gmh gmh_KA
 double gmh = 0.99;
#define gmm gmm_KA
 double gmm = 0.75;
#define htau htau_KA
 double htau = 0;
#define hinf hinf_KA
 double hinf = 0;
#define mtau mtau_KA
 double mtau = 0;
#define minf minf_KA
 double minf = 0;
#define q10 q10_KA
 double q10 = 3;
#define vhalfh vhalfh_KA
 double vhalfh = -70;
#define vhalfm vhalfm_KA
 double vhalfm = -45;
#define zetah zetah_KA
 double zetah = 0.2;
#define zetam zetam_KA
 double zetam = 0.1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau_KA", "ms",
 "htau_KA", "ms",
 "gbar_KA", "mho/cm2",
 "ik_KA", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "a0m_KA", &a0m_KA,
 "vhalfm_KA", &vhalfm_KA,
 "zetam_KA", &zetam_KA,
 "gmm_KA", &gmm_KA,
 "a0h_KA", &a0h_KA,
 "vhalfh_KA", &vhalfh_KA,
 "zetah_KA", &zetah_KA,
 "gmh_KA", &gmh_KA,
 "q10_KA", &q10_KA,
 "minf_KA", &minf_KA,
 "mtau_KA", &mtau_KA,
 "hinf_KA", &hinf_KA,
 "htau_KA", &htau_KA,
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
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"KA",
 "gbar_KA",
 "sha_KA",
 "shi_KA",
 "k_tauH_KA",
 "sh_tauH_KA",
 0,
 "ik_KA",
 0,
 "m_KA",
 "h_KA",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gbar = 0.0002;
 	sha = 9.9;
 	shi = 5.7;
 	k_tauH = 1;
 	sh_tauH = -0;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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

 void _KA_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 13, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 KA /home/kedoxey/OB_Model/OlfactoryBulb/prev_ob_models/Birgiolas2020/Mechanisms/KA.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "K-A";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int trates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 trates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   trates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
   }
  return 0;
}
 
static int  trates (  double _lv ) {
   minf = 1.0 / ( 1.0 + exp ( - ( _lv / 1.0 - sha - 7.6 ) / 14.0 ) ) ;
   mtau = betm ( _threadargscomma_ _lv ) / ( qt * a0m * ( 1.0 + alpm ( _threadargscomma_ _lv ) ) ) * 1.0 ;
   hinf = 1.0 / ( 1.0 + exp ( ( _lv / 1.0 - shi + 47.4 ) / 6.0 ) ) ;
   htau = k_tauH * beth ( _threadargscomma_ _lv ) / ( qt * a0h * ( 1.0 + alph ( _threadargscomma_ _lv ) ) ) * 1.0 ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alpm (  double _lv ) {
   double _lalpm;
 _lalpm = exp ( zetam * ( _lv / 1.0 - vhalfm ) ) ;
   
return _lalpm;
 }
 
static void _hoc_alpm(void) {
  double _r;
   _r =  alpm (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double betm (  double _lv ) {
   double _lbetm;
 _lbetm = exp ( zetam * gmm * ( _lv / 1.0 - vhalfm ) ) ;
   
return _lbetm;
 }
 
static void _hoc_betm(void) {
  double _r;
   _r =  betm (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double alph (  double _lv ) {
   double _lalph;
 _lalph = exp ( zetah * ( _lv / 1.0 - vhalfh - sh_tauH ) ) ;
   
return _lalph;
 }
 
static void _hoc_alph(void) {
  double _r;
   _r =  alph (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double beth (  double _lv ) {
   double _lbeth;
 _lbeth = exp ( zetah * gmh * ( _lv / 1.0 - vhalfh - sh_tauH ) ) ;
   
return _lbeth;
 }
 
static void _hoc_beth(void) {
  double _r;
   _r =  beth (  *getarg(1) );
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
  ek = _ion_ek;
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
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   qt = pow( q10 , ( ( celsius - 24.0 ) / 10.0 ) ) ;
   trates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
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
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ik = gbar * m * h * ( v - ek ) ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
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
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 56 in file KA.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/kedoxey/OB_Model/OlfactoryBulb/prev_ob_models/Birgiolas2020/Mechanisms/KA.mod";
static const char* nmodl_file_text = 
  "TITLE K-A\n"
  ": K-A current for Mitral Cells from Wang et al (1996)\n"
  ": M.Migliore Jan. 2002\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX KA\n"
  "	USEION k READ ek WRITE ik\n"
  "	RANGE  gbar, ik, m, h, sha,shi, k_tauH,sh_tauH\n"
  "	GLOBAL minf, mtau, hinf, htau\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gbar = 0.0002   	(mho/cm2)	\n"
  "								\n"
  "	celsius  (degC)\n"
  "	ek	= -70	(mV)            : must be explicitly def. in hoc\n"
  "	v 		(mV)\n"
  "	a0m=0.04\n"
  "	vhalfm=-45\n"
  "	zetam=0.1\n"
  "	gmm=0.75\n"
  "\n"
  "	a0h=0.018\n"
  "	vhalfh=-70\n"
  "	zetah=0.2\n"
  "	gmh=0.99\n"
  "\n"
  "	sha=9.9\n"
  "	shi=5.7\n"
  "	\n"
  "	q10=3\n"
  "	\n"
  "	k_tauH=1.0    : 2.5; added by GL\n"
  "	sh_tauH=-0    : -20; added by GL\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(pS) = (picosiemens)\n"
  "	(um) = (micron)\n"
  "} \n"
  "\n"
  "ASSIGNED {\n"
  "	ik 		(mA/cm2)\n"
  "	minf 		mtau (ms)	 	\n"
  "	hinf 		htau (ms)\n"
  "	qt\n"
  "}\n"
  " \n"
  "\n"
  "STATE { m h}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ik = gbar*m*h*(v - ek)\n"
  "} \n"
  "\n"
  "INITIAL {\n"
  "	qt=q10^((celsius-24(degC))/10(degC))\n"
  "	trates(v)\n"
  "	m=minf  \n"
  "	h=hinf  \n"
  "}\n"
  "\n"
  "DERIVATIVE states {   \n"
  "    trates(v)\n"
  "    m' = (minf-m)/mtau\n"
  "    h' = (hinf-h)/htau\n"
  "}\n"
  "\n"
  "PROCEDURE trates(v (mV)) {  \n"
  "	minf = 1/(1 + exp(-(v/1(mV)-sha-7.6)/14))\n"
  "	mtau = betm(v)/(qt*a0m*(1+alpm(v)))*1(ms)\n"
  "\n"
  "	hinf = 1/(1 + exp((v/1(mV)-shi+47.4)/6))\n"
  "	htau = k_tauH*beth(v)/(qt*a0h*(1+alph(v)))*1(ms)\n"
  "}\n"
  "\n"
  "FUNCTION alpm(v(mV)) {\n"
  "	alpm = exp(zetam*(v/1(mV)-vhalfm))\n"
  "}\n"
  "\n"
  "FUNCTION betm(v(mV)) {\n"
  "    betm = exp(zetam*gmm*(v/1(mV)-vhalfm))\n"
  "}\n"
  "\n"
  "FUNCTION alph(v(mV)) {\n"
  "    alph = exp(zetah*(v/1(mV)-vhalfh-sh_tauH))\n"
  "}\n"
  "\n"
  "FUNCTION beth(v(mV)) {\n"
  "    beth = exp(zetah*gmh*(v/1(mV)-vhalfh-sh_tauH))\n"
  "}\n"
  ;
#endif
