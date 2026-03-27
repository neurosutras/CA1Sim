/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
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
 
#define nrn_init _nrn_init__GABAb
#define _nrn_initial _nrn_initial__GABAb
#define nrn_cur _nrn_cur__GABAb
#define _nrn_current _nrn_current__GABAb
#define nrn_jacob _nrn_jacob__GABAb
#define nrn_state _nrn_state__GABAb
#define _net_receive _net_receive__GABAb 
#define bindkin bindkin__GABAb 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gmax _p[0]
#define gmax_columnindex 0
#define i _p[1]
#define i_columnindex 1
#define g _p[2]
#define g_columnindex 2
#define lastrelease _p[3]
#define lastrelease_columnindex 3
#define R _p[4]
#define R_columnindex 4
#define G _p[5]
#define G_columnindex 5
#define C _p[6]
#define C_columnindex 6
#define Gn _p[7]
#define Gn_columnindex 7
#define on _p[8]
#define on_columnindex 8
#define DR _p[9]
#define DR_columnindex 9
#define DG _p[10]
#define DG_columnindex 10
#define v _p[11]
#define v_columnindex 11
#define _g _p[12]
#define _g_columnindex 12
#define _tsav _p[13]
#define _tsav_columnindex 13
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
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
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
 _extcall_prop = _prop;
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
 0, 0
};
 /* declare global and static user variables */
#define Cdur Cdur_GABAb
 double Cdur = 1;
#define Cmax Cmax_GABAb
 double Cmax = 1;
#define Deadtime Deadtime_GABAb
 double Deadtime = 1;
#define Erev Erev_GABAb
 double Erev = -95;
#define KD KD_GABAb
 double KD = 100;
#define K4 K4_GABAb
 double K4 = 0.034;
#define K3 K3_GABAb
 double K3 = 0.18;
#define K2 K2_GABAb
 double K2 = 0.0012;
#define K1 K1_GABAb
 double K1 = 0.09;
#define n n_GABAb
 double n = 4;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Cmax_GABAb", "mM",
 "Cdur_GABAb", "ms",
 "Deadtime_GABAb", "ms",
 "K1_GABAb", "/ms",
 "K2_GABAb", "/ms",
 "K3_GABAb", "/ms",
 "K4_GABAb", "/ms",
 "Erev_GABAb", "mV",
 "gmax", "umho",
 "i", "nA",
 "g", "umho",
 "lastrelease", "ms",
 0,0
};
 static double G0 = 0;
 static double R0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Cmax_GABAb", &Cmax_GABAb,
 "Cdur_GABAb", &Cdur_GABAb,
 "Deadtime_GABAb", &Deadtime_GABAb,
 "K1_GABAb", &K1_GABAb,
 "K2_GABAb", &K2_GABAb,
 "K3_GABAb", &K3_GABAb,
 "K4_GABAb", &K4_GABAb,
 "KD_GABAb", &KD_GABAb,
 "n_GABAb", &n_GABAb,
 "Erev_GABAb", &Erev_GABAb,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABAb",
 "gmax",
 0,
 "i",
 "g",
 "lastrelease",
 0,
 "R",
 "G",
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
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 14;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _gabab_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABAb /Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/gabab.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "minimal model of GABAB receptors (NET_RECEIVE version)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int bindkin(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   if ( on  && ( t - lastrelease ) > Cdur ) {
     C = 0.0 ;
     on = 0.0 ;
     }
   DR = K1 * C * ( 1.0 - R ) - K2 * R ;
   DG = K3 * R - K4 * G ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 if ( on  && ( t - lastrelease ) > Cdur ) {
   C = 0.0 ;
   on = 0.0 ;
   }
 DR = DR  / (1. - dt*( ( K1 * C )*( ( ( - 1.0 ) ) ) - ( K2 )*( 1.0 ) )) ;
 DG = DG  / (1. - dt*( ( - ( K4 )*( 1.0 ) ) )) ;
  return 0;
}
 /*END CVODE*/
 static int bindkin (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   if ( on  && ( t - lastrelease ) > Cdur ) {
     C = 0.0 ;
     on = 0.0 ;
     }
    R = R + (1. - exp(dt*(( K1 * C )*( ( ( - 1.0 ) ) ) - ( K2 )*( 1.0 ))))*(- ( ( ( K1 )*( C ) )*( ( 1.0 ) ) ) / ( ( ( K1 )*( C ) )*( ( ( - 1.0 ) ) ) - ( K2 )*( 1.0 ) ) - R) ;
    G = G + (1. - exp(dt*(( - ( K4 )*( 1.0 ) ))))*(- ( ( K3 )*( R ) ) / ( ( - ( K4 )*( 1.0 ) ) ) - G) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   if ( ( t - lastrelease ) > ( Cdur + Deadtime ) ) {
     C = Cmax ;
     on = 1.0 ;
     lastrelease = t ;
     }
   } }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  G = G0;
  R = R0;
 {
   R = 0.0 ;
   G = 0.0 ;
   C = 0.0 ;
   on = 0.0 ;
   lastrelease = - 1e9 ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   Gn = pow( G , n ) ;
   g = gmax * Gn / ( Gn + KD ) ;
   i = g * ( v - Erev ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
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
 {   bindkin(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = R_columnindex;  _dlist1[0] = DR_columnindex;
 _slist1[1] = G_columnindex;  _dlist1[1] = DG_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/gabab.mod";
static const char* nmodl_file_text = 
  "TITLE minimal model of GABAB receptors (NET_RECEIVE version)\n"
  "\n"
  "COMMENT\n"
  "Minimal kinetic model for GABA-B receptors (Destexhe & Sejnowski 1995).\n"
  "Rewritten to use NET_RECEIVE for compatibility with NetCon/VecStim,\n"
  "replacing the old POINTER pre / voltage-threshold approach.\n"
  "\n"
  "Kinetics:\n"
  "  dR/dt = K1 * T * (1-R) - K2 * R\n"
  "  dG/dt = K3 * R - K4 * G\n"
  "  Open fraction: O = G^n / (G^n + KD)\n"
  "  i = gmax * O * (v - Erev)\n"
  "\n"
  "Written by Destexhe, Laval University, 1995.\n"
  "NET_RECEIVE adaptation: compatible with Synapse() / VecStim.\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "    POINT_PROCESS GABAb\n"
  "    RANGE R, G, g, gmax, lastrelease\n"
  "    NONSPECIFIC_CURRENT i\n"
  "    GLOBAL Cmax, Cdur, Deadtime\n"
  "    GLOBAL K1, K2, K3, K4, KD, n, Erev\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (nA) = (nanoamp)\n"
  "    (mV) = (millivolt)\n"
  "    (umho) = (micromho)\n"
  "    (mM)  = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    Cmax     = 1      (mM)   : max transmitter concentration\n"
  "    Cdur     = 1      (ms)   : transmitter pulse duration\n"
  "    Deadtime = 1      (ms)   : refractory period between events\n"
  "\n"
  "    K1  = 0.09   (/ms mM)   : receptor binding rate\n"
  "    K2  = 0.0012 (/ms)      : receptor unbinding rate\n"
  "    K3  = 0.18   (/ms)      : G-protein production rate\n"
  "    K4  = 0.034  (/ms)      : G-protein decay rate\n"
  "    KD  = 100               : K+ channel dissociation constant\n"
  "    n   = 4                 : G-protein binding cooperativity\n"
  "    Erev = -95    (mV)      : reversal potential (E_K)\n"
  "    gmax = 0      (umho)    : maximum conductance (set by user)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    v           (mV)\n"
  "    i           (nA)\n"
  "    g           (umho)\n"
  "    C           (mM)        : current transmitter concentration\n"
  "    Gn                      : G^n\n"
  "    lastrelease (ms)        : time of last NET_RECEIVE event\n"
  "    on                      : 1 = transmitter being released\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    R                       : fraction of activated receptors\n"
  "    G                       : fraction of activated G-protein\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "    R           = 0\n"
  "    G           = 0\n"
  "    C           = 0\n"
  "    on          = 0\n"
  "    lastrelease = -1e9\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    SOLVE bindkin METHOD cnexp\n"
  "    Gn = G^n\n"
  "    g  = gmax * Gn / (Gn + KD)\n"
  "    i  = g * (v - Erev)\n"
  "}\n"
  "\n"
  "DERIVATIVE bindkin {\n"
  "    : transmitter pulse: on for Cdur ms after each event\n"
  "    if (on && (t - lastrelease) > Cdur) {\n"
  "        C  = 0\n"
  "        on = 0\n"
  "    }\n"
  "    R' = K1 * C * (1 - R) - K2 * R\n"
  "    G' = K3 * R - K4 * G\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight) {\n"
  "    : weight is ignored - gmax controls amplitude\n"
  "    : only fire if outside the dead time\n"
  "    if ((t - lastrelease) > (Cdur + Deadtime)) {\n"
  "        C           = Cmax\n"
  "        on          = 1\n"
  "        lastrelease = t\n"
  "    }\n"
  "}\n"
  ;
#endif
