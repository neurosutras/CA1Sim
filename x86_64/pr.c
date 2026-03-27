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
 
#define nrn_init _nrn_init__Pr
#define _nrn_initial _nrn_initial__Pr
#define nrn_cur _nrn_cur__Pr
#define _nrn_current _nrn_current__Pr
#define nrn_jacob _nrn_jacob__Pr
#define nrn_state _nrn_state__Pr
#define _net_receive _net_receive__Pr 
#define setRandObjRef setRandObjRef__Pr 
 
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
#define P0 _p[0]
#define P0_columnindex 0
#define f _p[1]
#define f_columnindex 1
#define tau_F _p[2]
#define tau_F_columnindex 2
#define d1 _p[3]
#define d1_columnindex 3
#define tau_D1 _p[4]
#define tau_D1_columnindex 4
#define P _p[5]
#define P_columnindex 5
#define random _p[6]
#define random_columnindex 6
#define F _p[7]
#define F_columnindex 7
#define D1 _p[8]
#define D1_columnindex 8
#define tlast _p[9]
#define tlast_columnindex 9
#define v _p[10]
#define v_columnindex 10
#define _tsav _p[11]
#define _tsav_columnindex 11
#define _nd_area  *_ppvar[0]._pval
#define randObjPtr	*_ppvar[2]._pval
#define _p_randObjPtr	_ppvar[2]._pval
 
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
 static int hoc_nrnpointerindex =  2;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_randGen(void*);
 static double _hoc_setRandObjRef(void*);
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
 "randGen", _hoc_randGen,
 "setRandObjRef", _hoc_setRandObjRef,
 0, 0
};
#define randGen randGen_Pr
 extern double randGen( _threadargsproto_ );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "P0", 0, 1,
 "d1", 0, 1,
 "f", 0, 1e+09,
 "tau_D1", 1e-09, 1e+09,
 "tau_F", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "P0", "1",
 "f", "1",
 "tau_F", "ms",
 "d1", "1",
 "tau_D1", "ms",
 "tlast", "ms",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"Pr",
 "P0",
 "f",
 "tau_F",
 "d1",
 "tau_D1",
 0,
 "P",
 "random",
 "F",
 "D1",
 "tlast",
 0,
 0,
 "randObjPtr",
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
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	P0 = 0.2;
 	f = 1.769;
 	tau_F = 67.351;
 	d1 = 0.878;
 	tau_D1 = 92.918;
  }
 	_prop->param = _p;
 	_prop->param_size = 12;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _pr_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,(void*)0, (void*)0, (void*)0, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
 add_nrn_has_net_event(_mechtype);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Pr /Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/pr.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "stochastic release probability";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int setRandObjRef(_threadargsproto_);
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   F = 1.0 + ( F - 1.0 ) * exp ( - ( t - tlast ) / tau_F ) ;
   D1 = 1.0 - ( 1.0 - D1 ) * exp ( - ( t - tlast ) / tau_D1 ) ;
   if ( P0 * F * D1 > 1.0 ) {
     P = 1.0 ;
     }
   else {
     P = P0 * F * D1 ;
     }
   random = randGen ( _threadargs_ ) ;
   if ( random <= P ) {
     net_event ( _pnt, t ) ;
     }
   tlast = t ;
   F = F + f ;
   D1 = D1 * d1 ;
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       double* _p = _pnt->_prop->param;
    Datum* _ppvar = _pnt->_prop->dparam;
    Datum* _thread = (Datum*)0;
    NrnThread* _nt = (NrnThread*)_pnt->_vnt;
 tlast = t ;
   }
 
/*VERBATIM*/
#include "nrnran123.h"
 
double randGen ( _threadargsproto_ ) {
   double _lrandGen;
 
/*VERBATIM*/
   if (_p_randObjPtr) {
      /*
      :Supports separate independent but reproducible streams for
      : each instance. However, the corresponding hoc Random
      : distribution MUST be set to Random.uniform(0,1)
      */
      _lrandGen = nrn_random_pick((Rand*)_p_randObjPtr);
   }else{
      hoc_execerror("Random object ref not set correctly for randObjPtr"," only via hoc Random");
   }
 
return _lrandGen;
 }
 
static double _hoc_randGen(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  randGen ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
static int  setRandObjRef ( _threadargsproto_ ) {
   
/*VERBATIM*/
   Rand** pv4 = (Rand**)(&_p_randObjPtr);
   if (ifarg(1)) {
      *pv4 = nrn_random_arg(1);
   }else{
      *pv4 = (Rand*)0;
   }
  return 0; }
 
static double _hoc_setRandObjRef(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 setRandObjRef ( _p, _ppvar, _thread, _nt );
 return(_r);
}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
 {
   P = P0 ;
   random = 1.0 ;
   F = 1.0 ;
   D1 = 1.0 ;
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

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{
} return _current;
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
}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/pr.mod";
static const char* nmodl_file_text = 
  "TITLE stochastic release probability\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Milstein 2015. When this point process receives a spike, it requests a random number from a random number generator and\n"
  "compares it to an internal release probability variable 'P' to decide if a single vesicle should be released at this\n"
  "synapse. The release probability is then updated to follow specified dynamics of facilitation and depression. In order\n"
  "to make use of this event, additional synaptic mechanisms must be connected to this point process via a NetCon object.\n"
  "\n"
  "Dynamics based on:\n"
  "\n"
  "Implementation of the model of short-term facilitation and depression described in\n"
  "  Varela, J.A., Sen, K., Gibson, J., Fost, J., Abbott, L.R., and Nelson, S.B.\n"
  "  A quantitative description of short-term plasticity at excitatory synapses\n"
  "  in layer 2/3 of rat primary visual cortex\n"
  "  Journal of Neuroscience 17:7926-7940, 1997\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS Pr\n"
  "	RANGE P, P0, random, f, tau_F, d1, tau_D1, F, D1, tlast\n"
  "    THREADSAFE\n"
  "    POINTER randObjPtr\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    : the (1) is needed for the range limits to be effective\n"
  "    P0 = 0.200          (1)     < 0, 1 >        : basal release probability\n"
  "    f = 1.769           (1)     < 0, 1e9 >      : additive facilitation per spike\n"
  "    tau_F = 67.351      (ms)    < 1e-9, 1e9 >   : rate of decay back to baseline following facilitation\n"
  "    d1 = 0.878          (1)     < 0, 1 >        : multiplicative fast depression per spike\n"
  "    tau_D1 = 92.918     (ms)    < 1e-9, 1e9 >   : rate of decay back to baseline following fast depression\n"
  "    :d2 = 0.975         (1)     < 0, 1 >        : multiplicative slow depression per spike\n"
  "    :tau_D2 = 9200      (ms)    < 1e-9, 1e9 >   : rate of decay back to baseline following slow depression\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	P				        : instantaneous release probability\n"
  "    randObjPtr              : pointer to a hoc random number generator Random.uniform(0,1)\n"
  "    random                  : individual instance of random number\n"
  "    F                       : current level of facilitation\n"
  "    D1                      : current level of fast depression\n"
  "    :D2                      : current level of slow depression\n"
  "    tlast (ms)              : time of last spike\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	P = P0\n"
  "    random = 1\n"
  "    F = 1\n"
  "    D1 = 1\n"
  "    :D2 = 1\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight) {\n"
  "    INITIAL {\n"
  "        tlast = t\n"
  "    }\n"
  "    F = 1 + (F-1)*exp(-(t - tlast)/tau_F)\n"
  "    D1 = 1 - (1-D1)*exp(-(t - tlast)/tau_D1)\n"
  ":    D2 = 1 - (1-D2)*exp(-(t - tlast)/tau_D2)\n"
  ":    if (P0*F*D1*D2 > 1) {\n"
  "    if (P0*F*D1 > 1) {\n"
  "        P = 1\n"
  "    } else {\n"
  ":        P = P0*F*D1*D2\n"
  "        P = P0*F*D1\n"
  "    }\n"
  "    random = randGen()\n"
  "    if (random <= P) {\n"
  "        net_event(t)\n"
  "    }\n"
  "    tlast = t\n"
  "    F = F + f\n"
  "    D1 = D1 * d1\n"
  ":    D2 = D2 * d2\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "#include \"nrnran123.h\"\n"
  "ENDVERBATIM\n"
  "\n"
  "FUNCTION randGen() {\n"
  "VERBATIM\n"
  "   if (_p_randObjPtr) {\n"
  "      /*\n"
  "      :Supports separate independent but reproducible streams for\n"
  "      : each instance. However, the corresponding hoc Random\n"
  "      : distribution MUST be set to Random.uniform(0,1)\n"
  "      */\n"
  "      _lrandGen = nrn_random_pick((Rand*)_p_randObjPtr);\n"
  "   }else{\n"
  "      hoc_execerror(\"Random object ref not set correctly for randObjPtr\",\" only via hoc Random\");\n"
  "   }\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE setRandObjRef() {\n"
  "VERBATIM\n"
  "   Rand** pv4 = (Rand**)(&_p_randObjPtr);\n"
  "   if (ifarg(1)) {\n"
  "      *pv4 = nrn_random_arg(1);\n"
  "   }else{\n"
  "      *pv4 = (Rand*)0;\n"
  "   }\n"
  "ENDVERBATIM\n"
  "}\n"
  ;
#endif
