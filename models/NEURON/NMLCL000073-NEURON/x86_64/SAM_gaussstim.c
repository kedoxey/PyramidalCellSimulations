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
 
#define nrn_init _nrn_init__SAM_GaussStim
#define _nrn_initial _nrn_initial__SAM_GaussStim
#define nrn_cur _nrn_cur__SAM_GaussStim
#define _nrn_current _nrn_current__SAM_GaussStim
#define nrn_jacob _nrn_jacob__SAM_GaussStim
#define nrn_state _nrn_state__SAM_GaussStim
#define _net_receive _net_receive__SAM_GaussStim 
#define event_time event_time__SAM_GaussStim 
#define init_sequence init_sequence__SAM_GaussStim 
#define seed seed__SAM_GaussStim 
 
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
#define interval _p[0]
#define number _p[1]
#define start _p[2]
#define factor _p[3]
#define refrac _p[4]
#define k _p[5]
#define fm _p[6]
#define t_shift _p[7]
#define x _p[8]
#define N_forward _p[9]
#define N_backward _p[10]
#define N_normal _p[11]
#define N_total _p[12]
#define amp_mean _p[13]
#define amp_std _p[14]
#define rand _p[15]
#define event _p[16]
#define on _p[17]
#define end _p[18]
#define m _p[19]
#define diff _p[20]
#define _tsav _p[21]
#define _nd_area  *_ppvar[0]._pval
 
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
 /* declaration of user functions */
 static double _hoc_event_time();
 static double _hoc_init_sequence();
 static double _hoc_invl();
 static double _hoc_seed();
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

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "event_time", _hoc_event_time,
 "init_sequence", _hoc_init_sequence,
 "invl", _hoc_invl,
 "seed", _hoc_seed,
 0, 0
};
#define invl invl_SAM_GaussStim
 extern double invl( double );
 /* declare global and static user variables */
#define pi pi_SAM_GaussStim
 double pi = 3.14159;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "fm", 0, 1000,
 "interval", 1e-09, 1e+09,
 "k", 0, 1,
 "number", 0, 1e+09,
 "refrac", 0, 3,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "interval", "ms",
 "start", "ms",
 "t_shift", "ms",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "pi_SAM_GaussStim", &pi_SAM_GaussStim,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"SAM_GaussStim",
 "interval",
 "number",
 "start",
 "factor",
 "refrac",
 "k",
 "fm",
 "t_shift",
 0,
 "x",
 "N_forward",
 "N_backward",
 "N_normal",
 "N_total",
 "amp_mean",
 "amp_std",
 "rand",
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 22, _prop);
 	/*initialize range parameters*/
 	interval = 10;
 	number = 10000;
 	start = 0;
 	factor = 4;
 	refrac = 0.5;
 	k = 1;
 	fm = 200;
 	t_shift = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 22;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _SAM_gaussstim_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 22, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
 add_nrn_has_net_event(_mechtype);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 SAM_GaussStim /home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/SAM_gaussstim.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int event_time();
static int init_sequence(double);
static int seed(double);
 
static int  seed (  double _lx ) {
   set_seed ( _lx ) ;
    return 0; }
 
static double _hoc_seed(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 seed (  *getarg(1) );
 return(_r);
}
 
static int  init_sequence (  double _lt ) {
   if ( number > 0.0 ) {
     on = 1.0 ;
     event = _lt ;
     end = _lt + 1e-6 + interval * ( number - 1.0 ) ;
     }
    return 0; }
 
static double _hoc_init_sequence(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 init_sequence (  *getarg(1) );
 return(_r);
}
 
double invl (  double _lmean ) {
   double _linvl;
 double _lstd ;
 if ( _lmean <= 0. ) {
     _lmean = .01 ;
     }
   _lstd = _lmean / factor ;
   _linvl = normrand ( _lmean , _lstd ) ;
   if ( _linvl >= interval ) {
     _linvl = fmod ( _linvl , interval ) ;
     N_forward = N_forward + 1.0 ;
     }
   else if ( _linvl < 0.0 ) {
     _linvl = fmod ( _linvl , interval ) + interval ;
     N_backward = N_backward + 1.0 ;
     }
   else {
     N_normal = N_normal + 1.0 ;
     }
   diff = interval - _linvl ;
   
return _linvl;
 }
 
static double _hoc_invl(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r =  invl (  *getarg(1) );
 return(_r);
}
 
static int  event_time (  ) {
   double _ldiff2 , _lT , _lrnd ;
 _ldiff2 = diff ;
   if ( number > 0.0 ) {
     _lT = invl ( _threadargscomma_ m ) ;
     _lrnd = _lT ;
     if ( _lT  == 0.0  && _ldiff2  == 0.0 ) {
       _lT = _lT + dt ;
       }
     event = _lT + event + _ldiff2 ;
     N_total = N_total + 1.0 ;
     }
   if ( event > end ) {
     on = 0.0 ;
     }
    return 0; }
 
static double _hoc_event_time(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 event_time (  );
 return(_r);
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 0.0 ) {
     if ( _args[0] > 0.0  && on  == 0.0 ) {
       init_sequence ( _threadargscomma_ t ) ;
       net_send ( _tqitem, _args, _pnt, t +  0.0 , 1.0 ) ;
       }
     else if ( _args[0] < 0.0  && on  == 1.0 ) {
       on = 0.0 ;
       }
     }
   if ( _lflag  == 3.0 ) {
     if ( on  == 0.0 ) {
       init_sequence ( _threadargscomma_ t ) ;
       net_send ( _tqitem, _args, _pnt, t +  0.0 , 1.0 ) ;
       }
     }
   if ( _lflag  == 1.0  && on  == 1.0 ) {
     if ( x  == 0.0 ) {
       rand = normrand ( amp_mean , amp_std ) ;
       x = rand * ( 1.0 + k * cos ( ( 0.001 ) * 2.0 * pi * fm * ( t - t_shift ) ) ) / ( 1.0 + k ) ;
       net_event ( _pnt, t ) ;
       event_time ( _threadargs_ ) ;
       if ( event - t <= refrac + dt  && event - t >= refrac ) {
         event = event + dt ;
         }
       if ( on  == 1.0 ) {
         net_send ( _tqitem, _args, _pnt, t +  event - t , 1.0 ) ;
         }
       net_send ( _tqitem, _args, _pnt, t +  refrac , 2.0 ) ;
       }
     else if ( x  != 0.0 ) {
       net_event ( _pnt, t ) ;
       event_time ( _threadargs_ ) ;
       if ( on  == 1.0 ) {
         net_send ( _tqitem, _args, _pnt, t +  event - t , 1.0 ) ;
         }
       }
     }
   if ( _lflag  == 2.0 ) {
     x = 0.0 ;
     }
   } }

static void initmodel() {
  int _i; double _save;_ninits++;
{
 {
   on = 0.0 ;
   x = 0.0 ;
   diff = 0.0 ;
   m = interval / 2.0 ;
   N_forward = 0.0 ;
   N_backward = 0.0 ;
   N_normal = 0.0 ;
   N_total = 0.0 ;
   if ( start >= 0.0  && number > 0.0 ) {
     event = start + invl ( _threadargscomma_ m ) ;
     if ( event < 0.0 ) {
       event = 0.0 ;
       }
     net_send ( _tqitem, (double*)0, _ppvar[1]._pvoid, t +  event , 3.0 ) ;
     }
   }

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
 _tsav = -1e20;
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

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

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
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/SAM_gaussstim.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "generate sinosoidal amplitude modulated pulse trains;\n"
  "\n"
  "x=Xmax*[1+kcos(2*pi*fm(t-t_shift))], where x is time-varying envelope of a pulse train\n"
  "\n"
  "Xmax=N(mean,std)\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "NEURON	{ \n"
  "  POINT_PROCESS SAM_GaussStim\n"
  "  RANGE x\n"
  "  RANGE interval, number, start,factor,refrac\n"
  "  \n"
  "  RANGE N_backward,N_forward,N_normal,N_total\n"
  "\n"
  "  RANGE k,fm,t_shift\n"
  "  RANGE amp_mean,amp_std, rand\n"
  "  \n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	interval	= 10 (ms) <1e-9,1e9>: time between spikes (msec)\n"
  "	number	= 10000 <0,1e9>	: number of spikes\n"
  "	start		= 0 (ms)	: start of first spike\n"
  "	factor =4  : portion of std to mean\n"
  "	refrac		=0.5 <0,3>    : absolute refractory period, up limit of freq=1kHz\n"
  "						: refrac >=dt, otherwise x can't be reset to 0 \n"
  "        pi=3.1415927\n"
  "	k=1		<0,1>   : modulation depth\n"
  "	fm=200		<0, 1000>   : modulation frequency\n"
  "	t_shift=0 (ms)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	x\n"
  "	event (ms)\n"
  "	\n"
  "	on\n"
  "	end (ms)\n"
  "	m (ms)            : mean of Gaussian\n"
  "	diff (ms)\n"
  "	N_forward  : swap spike whose value exceed forwardly one interval \n"
  "	N_backward  : swap spike whose value exceed backwardly one interval \n"
  "	N_normal\n"
  "	N_total\n"
  "	\n"
  " 	amp_mean\n"
  "	amp_std\n"
  "	rand\n"
  "}\n"
  "\n"
  "PROCEDURE seed(x) {\n"
  "	set_seed(x)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	on = 0\n"
  "	x = 0\n"
  "	diff=0\n"
  "	m=interval/2  : each T has normal distribution N(interval/2, interval/2/factor)\n"
  "	N_forward=0\n"
  "	N_backward=0\n"
  "	N_normal=0\n"
  "	N_total=0\n"
  "	if (start >= 0 && number > 0) {\n"
  "		\n"
  "		event = start + invl(m) \n"
  "		\n"
  "		if (event < 0) {\n"
  "			event = 0\n"
  "		}\n"
  "		net_send(event, 3)\n"
  "	}\n"
  "}	\n"
  "\n"
  "PROCEDURE init_sequence(t(ms)) {\n"
  "	if (number > 0) {\n"
  "		on = 1\n"
  "		event = t\n"
  "		end = t + 1e-6 + interval*(number-1)\n"
  "		:printf(\"next event=%g\\n\",event)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION invl(mean (ms)) (ms) { LOCAL std\n"
  "	if (mean <= 0.) {\n"
  "		mean = .01 (ms) : I would worry if it were 0.\n"
  "	}\n"
  "	std=mean/factor  : std=T/(2*factor)\n"
  "	invl = normrand(mean, std)  : relative to current interval \n"
  "	\n"
  "	if(invl>=interval) { \n"
  "		:printf(\"original=%g\\n\",invl)\n"
  "		invl=fmod(invl,interval)\n"
  "		:printf(\"now=%g\\n\",invl)\n"
  "\n"
  "		N_forward=N_forward+1\n"
  "\n"
  "\n"
  "		}else if(invl<0) { \n"
  "\n"
  "			:printf(\"original=%g\\n\",invl)\n"
  "			invl=fmod(invl,interval)+interval\n"
  "			:printf(\"now=%g\\n\",invl)\n"
  "\n"
  "			N_backward=N_backward+1\n"
  "			}else {\n"
  "			N_normal=N_normal+1\n"
  "			}\n"
  "		\n"
  "		diff=interval-invl\n"
  "	\n"
  "	:printf(\"invl=%g\\n\",invl)\n"
  "}\n"
  "\n"
  "PROCEDURE event_time() {LOCAL diff2,T,rnd\n"
  "        diff2=diff\n"
  "	if (number > 0) {\n"
  "	   T=invl(m)\n"
  "	   rnd=T\n"
  "	   if(T==0 && diff2==0) { T=T+dt } : previous and current spikes overlapped\n"
  "	   \n"
  "	   :printf(\"event=%g  diff=%g\\n\",event,diff2)\n"
  "	   event = T+event + diff2    :compute absolute event time, relative to 0ms\n"
  "	   \n"
  " 	   N_total=N_total+1\n"
  " 	}\n"
  " 			\n"
  "	if (event > end) {\n"
  "		on = 0\n"
  "	}\n"
  "}\n"
  "\n"
  "NET_RECEIVE (w) {\n"
  "	if (flag == 0) { : external event\n"
  "		if (w > 0 && on == 0) { : turn on spike sequence\n"
  "			init_sequence(t)\n"
  "			net_send(0, 1)\n"
  "		}else if (w < 0 && on == 1) { : turn off spiking\n"
  "			on = 0\n"
  "		}\n"
  "	}\n"
  "	if (flag == 3) { : from INITIAL\n"
  "		if (on == 0) {\n"
  "			init_sequence(t)\n"
  "			net_send(0, 1)\n"
  "		}\n"
  "	}\n"
  "	if (flag == 1 && on == 1) {\n"
  "		if(x == 0){  : after refractory\n"
  "				rand=normrand(amp_mean, amp_std) : normal distribution\n"
  "				x=rand*(1+k*cos((0.001)*2*pi*fm*(t-t_shift)))/(1+k)\n"
  "				net_event(t)\n"
  "				event_time()  : after each spike call next spike time\n"
  "\n"
  "		  		 if (event-t <= refrac+dt && event-t >= refrac) {\n"
  "		      		 	event=event+dt    : happen when next event at reset edge time\n"
  "		       			:printf(\"next spike on edge %g\\n\", event)\n"
  "		  		 }\n"
  "\n"
  "		  		 if (on==1) {\n"
  "						net_send(event - t, 1) : self-feed next event\n"
  "		  				 }\n"
  "				net_send(refrac, 2)  : refractory period\n"
  "		\n"
  "		} else if (x!=0) {\n"
  "			net_event(t)\n"
  "			event_time() : although this spike will be ignored, still call next spike : independent cycles \n"
  "			:printf(\"inside refrac\\n\")\n"
  "\n"
  "			if (on == 1) {\n"
  "				net_send(event - t, 1)\n"
  "					  }\n"
  "		}\n"
  "		:printf(\"x=1 at %g\\n\",t)\n"
  "	}\n"
  "	if (flag == 2) {\n"
  "		:printf(\"x=0 at %g\\n\",t)\n"
  "		x = 0\n"
  "		\n"
  "	}\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
