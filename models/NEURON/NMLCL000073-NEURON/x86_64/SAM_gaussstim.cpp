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
static constexpr auto number_of_datum_variables = 2;
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
#define interval _ml->template fpfield<0>(_iml)
#define interval_columnindex 0
#define number _ml->template fpfield<1>(_iml)
#define number_columnindex 1
#define start _ml->template fpfield<2>(_iml)
#define start_columnindex 2
#define factor _ml->template fpfield<3>(_iml)
#define factor_columnindex 3
#define refrac _ml->template fpfield<4>(_iml)
#define refrac_columnindex 4
#define k _ml->template fpfield<5>(_iml)
#define k_columnindex 5
#define fm _ml->template fpfield<6>(_iml)
#define fm_columnindex 6
#define t_shift _ml->template fpfield<7>(_iml)
#define t_shift_columnindex 7
#define x _ml->template fpfield<8>(_iml)
#define x_columnindex 8
#define N_forward _ml->template fpfield<9>(_iml)
#define N_forward_columnindex 9
#define N_backward _ml->template fpfield<10>(_iml)
#define N_backward_columnindex 10
#define N_normal _ml->template fpfield<11>(_iml)
#define N_normal_columnindex 11
#define N_total _ml->template fpfield<12>(_iml)
#define N_total_columnindex 12
#define amp_mean _ml->template fpfield<13>(_iml)
#define amp_mean_columnindex 13
#define amp_std _ml->template fpfield<14>(_iml)
#define amp_std_columnindex 14
#define rand _ml->template fpfield<15>(_iml)
#define rand_columnindex 15
#define event _ml->template fpfield<16>(_iml)
#define event_columnindex 16
#define on _ml->template fpfield<17>(_iml)
#define on_columnindex 17
#define end _ml->template fpfield<18>(_iml)
#define end_columnindex 18
#define m _ml->template fpfield<19>(_iml)
#define m_columnindex 19
#define diff _ml->template fpfield<20>(_iml)
#define diff_columnindex 20
#define _tsav _ml->template fpfield<21>(_iml)
#define _tsav_columnindex 21
#define _nd_area *_ml->dptr_field<0>(_iml)
 static _nrn_mechanism_cache_instance _ml_real{nullptr};
static _nrn_mechanism_cache_range *_ml{&_ml_real};
static size_t _iml{0};
static Datum *_ppvar;
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_event_time(void*);
 static double _hoc_init_sequence(void*);
 static double _hoc_invl(void*);
 static double _hoc_seed(void*);
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
 neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
_ppvar = _nrn_mechanism_access_dparam(_prop);
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
 {"event_time", _hoc_event_time},
 {"init_sequence", _hoc_init_sequence},
 {"invl", _hoc_invl},
 {"seed", _hoc_seed},
 {0, 0}
};
#define invl invl_SAM_GaussStim
 extern double invl( double );
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
#define pi pi_SAM_GaussStim
 double pi = 3.14159;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {"fm", 0, 1000},
 {"interval", 1e-09, 1e+09},
 {"k", 0, 1},
 {"number", 0, 1e+09},
 {"refrac", 0, 3},
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {"interval", "ms"},
 {"start", "ms"},
 {"t_shift", "ms"},
 {0, 0}
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"pi_SAM_GaussStim", &pi_SAM_GaussStim},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
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
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
     10, /* interval */
     10000, /* number */
     0, /* start */
     4, /* factor */
     0.5, /* refrac */
     1, /* k */
     200, /* fm */
     0, /* t_shift */
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
  if (nrn_point_prop_) {
    _nrn_mechanism_access_alloc_seq(_prop) = _nrn_mechanism_access_alloc_seq(nrn_point_prop_);
    _ppvar = _nrn_mechanism_access_dparam(nrn_point_prop_);
  } else {
   _ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
    _nrn_mechanism_access_dparam(_prop) = _ppvar;
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 22);
 	/*initialize range parameters*/
 	interval = _parm_default[0]; /* 10 */
 	number = _parm_default[1]; /* 10000 */
 	start = _parm_default[2]; /* 0 */
 	factor = _parm_default[3]; /* 4 */
 	refrac = _parm_default[4]; /* 0.5 */
 	k = _parm_default[5]; /* 1 */
 	fm = _parm_default[6]; /* 200 */
 	t_shift = _parm_default[7]; /* 0 */
  }
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 22);
 	_nrn_mechanism_access_dparam(_prop) = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 
#define _tqitem &(_ppvar[2])
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _SAM_gaussstim_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nullptr, nullptr, nullptr, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype,
                                       _nrn_mechanism_field<double>{"interval"} /* 0 */,
                                       _nrn_mechanism_field<double>{"number"} /* 1 */,
                                       _nrn_mechanism_field<double>{"start"} /* 2 */,
                                       _nrn_mechanism_field<double>{"factor"} /* 3 */,
                                       _nrn_mechanism_field<double>{"refrac"} /* 4 */,
                                       _nrn_mechanism_field<double>{"k"} /* 5 */,
                                       _nrn_mechanism_field<double>{"fm"} /* 6 */,
                                       _nrn_mechanism_field<double>{"t_shift"} /* 7 */,
                                       _nrn_mechanism_field<double>{"x"} /* 8 */,
                                       _nrn_mechanism_field<double>{"N_forward"} /* 9 */,
                                       _nrn_mechanism_field<double>{"N_backward"} /* 10 */,
                                       _nrn_mechanism_field<double>{"N_normal"} /* 11 */,
                                       _nrn_mechanism_field<double>{"N_total"} /* 12 */,
                                       _nrn_mechanism_field<double>{"amp_mean"} /* 13 */,
                                       _nrn_mechanism_field<double>{"amp_std"} /* 14 */,
                                       _nrn_mechanism_field<double>{"rand"} /* 15 */,
                                       _nrn_mechanism_field<double>{"event"} /* 16 */,
                                       _nrn_mechanism_field<double>{"on"} /* 17 */,
                                       _nrn_mechanism_field<double>{"end"} /* 18 */,
                                       _nrn_mechanism_field<double>{"m"} /* 19 */,
                                       _nrn_mechanism_field<double>{"diff"} /* 20 */,
                                       _nrn_mechanism_field<double>{"_tsav"} /* 21 */,
                                       _nrn_mechanism_field<double*>{"_nd_area", "area"} /* 0 */,
                                       _nrn_mechanism_field<Point_process*>{"_pntproc", "pntproc"} /* 1 */,
                                       _nrn_mechanism_field<void*>{"_tqitem", "netsend"} /* 2 */);
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
static const char *modelname = "";

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
    auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _setdata(_p);
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
    auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _setdata(_p);
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
    auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _setdata(_p);
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
    auto* const _pnt = static_cast<Point_process*>(_vptr);
  auto* const _p = _pnt->_prop;
  if (!_p) {
    hoc_execerror("POINT_PROCESS data instance not valid", NULL);
  }
   _setdata(_p);
 _r = 1.;
 event_time (  );
 return(_r);
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{   neuron::legacy::set_globals_from_prop(_pnt->_prop, _ml_real, _ml, _iml);
    _ppvar = _nrn_mechanism_access_dparam(_pnt->_prop);
  if (_tsav > t){ hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = nullptr;}
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
     net_send ( _tqitem, nullptr, _ppvar[1].get<Point_process*>(), t +  event , 3.0 ) ;
     }
   }

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
 _tsav = -1e20;
   _v = _vec_v[_ni[_iml]];
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

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
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON/SAM_gaussstim.mod";
    const char* nmodl_file_text = 
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
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
