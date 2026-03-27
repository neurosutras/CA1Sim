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
 
#define nrn_init _nrn_init__NMDA_KIN5
#define _nrn_initial _nrn_initial__NMDA_KIN5
#define nrn_cur _nrn_cur__NMDA_KIN5
#define _nrn_current _nrn_current__NMDA_KIN5
#define nrn_jacob _nrn_jacob__NMDA_KIN5
#define nrn_state _nrn_state__NMDA_KIN5
#define _net_receive _net_receive__NMDA_KIN5 
#define kstates kstates__NMDA_KIN5 
#define mgblock mgblock__NMDA_KIN5 
 
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
#define Cmax _p[0]
#define Cmax_columnindex 0
#define Cdur _p[1]
#define Cdur_columnindex 1
#define kon _p[2]
#define kon_columnindex 2
#define koff _p[3]
#define koff_columnindex 3
#define CC _p[4]
#define CC_columnindex 4
#define CO _p[5]
#define CO_columnindex 5
#define Beta _p[6]
#define Beta_columnindex 6
#define Alpha _p[7]
#define Alpha_columnindex 7
#define Erev _p[8]
#define Erev_columnindex 8
#define Kd _p[9]
#define Kd_columnindex 9
#define gamma _p[10]
#define gamma_columnindex 10
#define mg _p[11]
#define mg_columnindex 11
#define gmax _p[12]
#define gmax_columnindex 12
#define kin_scale _p[13]
#define kin_scale_columnindex 13
#define i _p[14]
#define i_columnindex 14
#define g _p[15]
#define g_columnindex 15
#define B _p[16]
#define B_columnindex 16
#define Ru _p[17]
#define Ru_columnindex 17
#define Rb _p[18]
#define Rb_columnindex 18
#define Rc _p[19]
#define Rc_columnindex 19
#define Ro _p[20]
#define Ro_columnindex 20
#define C _p[21]
#define C_columnindex 21
#define kin_factor _p[22]
#define kin_factor_columnindex 22
#define scale _p[23]
#define scale_columnindex 23
#define DRu _p[24]
#define DRu_columnindex 24
#define DRb _p[25]
#define DRb_columnindex 25
#define DRc _p[26]
#define DRc_columnindex 26
#define DRo _p[27]
#define DRo_columnindex 27
#define v _p[28]
#define v_columnindex 28
#define _g _p[29]
#define _g_columnindex 29
#define _tsav _p[30]
#define _tsav_columnindex 30
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
 static double _hoc_mgblock(void*);
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
 "mgblock", _hoc_mgblock,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Cmax", "mM",
 "Cdur", "ms",
 "kon", "/mM/ms",
 "koff", "/ms",
 "CC", "/ms",
 "CO", "/ms",
 "Beta", "/ms",
 "Alpha", "/ms",
 "Erev", "mV",
 "Kd", "mM",
 "gamma", "/mV",
 "mg", "mM",
 "gmax", "umho",
 "kin_scale", "1",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double Ro0 = 0;
 static double Rc0 = 0;
 static double Rb0 = 0;
 static double Ru0 = 0;
 static double delta_t = 0.01;
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
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"NMDA_KIN5",
 "Cmax",
 "Cdur",
 "kon",
 "koff",
 "CC",
 "CO",
 "Beta",
 "Alpha",
 "Erev",
 "Kd",
 "gamma",
 "mg",
 "gmax",
 "kin_scale",
 0,
 "i",
 "g",
 "B",
 0,
 "Ru",
 "Rb",
 "Rc",
 "Ro",
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
 	_p = nrn_prop_data_alloc(_mechtype, 31, _prop);
 	/*initialize range parameters*/
 	Cmax = 1;
 	Cdur = 0.3;
 	kon = 86.89;
 	koff = 0.69;
 	CC = 9.64;
 	CO = 2.6;
 	Beta = 0.68;
 	Alpha = 0.079;
 	Erev = 0;
 	Kd = 9.98;
 	gamma = 0.101;
 	mg = 1;
 	gmax = 0.003026;
 	kin_scale = 1.83;
  }
 	_prop->param = _p;
 	_prop->param_size = 31;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
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
 
#define _tqitem &(_ppvar[2]._pvoid)
 static void _net_receive(Point_process*, double*, double);
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _nmda_kin5_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 3,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 31, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "netsend");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 NMDA_KIN5 /Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/nmda_kin5.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int mgblock(_threadargsprotocomma_ double);
 extern double *_nrn_thread_getelm(SparseObj*, int, int);
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[4], _dlist1[4]; static double *_temp1;
 static int kstates();
 
static int kstates (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<4;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 /* ~ Ru <-> Rb ( C * kon , koff )*/
 f_flux =  C * kon * Ru ;
 b_flux =  koff * Rb ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  C * kon ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  koff ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ Rb <-> Rc ( CC , CO * kin_factor )*/
 f_flux =  CC * Rb ;
 b_flux =  CO * kin_factor * Rc ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 0) += (f_flux - b_flux);
 
 _term =  CC ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  CO * kin_factor ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ Rc <-> Ro ( Beta , Alpha * kin_factor )*/
 f_flux =  Beta * Rc ;
 b_flux =  Alpha * kin_factor * Ro ;
 _RHS1( 0) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  Beta ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 3 ,0)  -= _term;
 _term =  Alpha * kin_factor ;
 _MATELM1( 0 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
    } return _reset;
 }
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t;   if (_lflag == 1. ) {*(_tqitem) = 0;}
 {
   if ( _lflag  == 0.0 ) {
     C = Cmax ;
     scale = _args[0] ;
     net_send ( _tqitem, _args, _pnt, t +  Cdur , 1.0 ) ;
     }
   else {
     C = 0. ;
     }
   } }
 
static int  mgblock ( _threadargsprotocomma_ double _lv ) {
   B = 1. / ( 1. + exp ( gamma * ( - _lv ) ) * ( mg / Kd ) ) ;
   kin_factor = ( 1.0 - kin_scale ) * B + kin_scale ;
    return 0; }
 
static double _hoc_mgblock(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r = 1.;
 mgblock ( _p, _ppvar, _thread, _nt, *getarg(1) );
 return(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<4;_i++) _p[_dlist1[_i]] = 0.0;}
 /* ~ Ru <-> Rb ( C * kon , koff )*/
 f_flux =  C * kon * Ru ;
 b_flux =  koff * Rb ;
 DRu -= (f_flux - b_flux);
 DRb += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ Rb <-> Rc ( CC , CO * kin_factor )*/
 f_flux =  CC * Rb ;
 b_flux =  CO * kin_factor * Rc ;
 DRb -= (f_flux - b_flux);
 DRc += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ Rc <-> Ro ( Beta , Alpha * kin_factor )*/
 f_flux =  Beta * Rc ;
 b_flux =  Alpha * kin_factor * Ro ;
 DRc -= (f_flux - b_flux);
 DRo += (f_flux - b_flux);
 
 /*REACTION*/
    } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<4;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 /* ~ Ru <-> Rb ( C * kon , koff )*/
 _term =  C * kon ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  koff ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ Rb <-> Rc ( CC , CO * kin_factor )*/
 _term =  CC ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  CO * kin_factor ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ Rc <-> Ro ( Beta , Alpha * kin_factor )*/
 _term =  Beta ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 3 ,0)  -= _term;
 _term =  Alpha * kin_factor ;
 _MATELM1( 0 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
    } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 4;}
 
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
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 4, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  Rc = Rc0;
  Rb = Rb0;
  Ru = Ru0;
  Ro = Ro0;
 {
   C = 0. ;
   Ru = 1. ;
   Rb = 0. ;
   Rc = 0. ;
   Ro = 0. ;
   scale = 1. ;
   mgblock ( _threadargscomma_ v ) ;
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
   mgblock ( _threadargscomma_ v ) ;
   g = scale * gmax * B * Ro ;
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
 {  sparse_thread(&_thread[_spth1]._pvoid, 4, _slist1, _dlist1, _p, &t, dt, kstates, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = Rc_columnindex;  _dlist1[0] = DRc_columnindex;
 _slist1[1] = Rb_columnindex;  _dlist1[1] = DRb_columnindex;
 _slist1[2] = Ru_columnindex;  _dlist1[2] = DRu_columnindex;
 _slist1[3] = Ro_columnindex;  _dlist1[3] = DRo_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/nmda_kin5.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "GluN-R (NMDA-type Glutamate Receptors)\n"
  "\n"
  "Postsynaptic current is additionally constrained by voltage-dependent channel\n"
  "block by extracellular Mg.\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "Synaptic mechanism based on a simplified model of transmitter binding to\n"
  "postsynaptic receptors.\n"
  "\n"
  "Modified by Aaron Milstein, 2015.\n"
  "\n"
  "Modification of original code by:\n"
  "\n"
  "A. Destexhe & Z. Mainen, The Salk Institute, March 12, 1993.\n"
  "Last modif. Sept 8, 1993.\n"
  "\n"
  "Reference:\n"
  "\n"
  "Destexhe, A., Mainen, Z. and Sejnowski, T.J.  An efficient method for\n"
  "computing synaptic conductances based on a kinetic model of receptor binding.\n"
  "Neural Computation, 6: 14-18, 1994.\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "Upon arrival of a presynaptic spike (a net_event), the concentration of the\n"
  "neurotransmitter C in the synaptic cleft is briefly stepped to concentration\n"
  "Cmax for duration Cdur.\n"
  "\n"
  "C     _____ . . . . . . Cmax\n"
  "      |   |\n"
  " _____|   |______ . . . 0\n"
  "     t0   t0 + Cdur\n"
  "\n"
  "The receptors then bind transmitter, change conformation, and open their\n"
  "channel according to the following kinetic scheme:\n"
  "\n"
  "       ---[C] * kon-->        -------CC----->        -----Beta----->\n"
  "C + Ru <-----koff-----   Rb   <------CO------   Rc   <----Alpha-----   Ro\n"
  "\n"
  "where Ru, Rb, Rc, and Ro are respectively the fraction of channels in the\n"
  "unbound, closed bound, closed cleft, and open states of the postsynaptic\n"
  "receptor. kon and koff are the binding and unbinding rate constants, CC and\n"
  "CO are the closing and opening rates of the receptor ligand-binding cleft,\n"
  "and Beta and Alpha are the opening and closing rates of the channel.\n"
  "\n"
  "The maximal conductance of the channel can be set for each instance of the\n"
  "synaptic mechanism by specifying gmax, and the relative weight of events\n"
  "can be set through the weight parameter of the associated netcon object.\n"
  "\n"
  "The postsynaptic current is given by:\n"
  "\n"
  "i = weight * gmax * Ro * (V-Erev)\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS NMDA_KIN5\n"
  "	RANGE Cmax, Cdur, kon, koff, CC, CO, Beta, Alpha, Erev, Kd, gamma, mg, gmax, g, B, kin_scale, Ro\n"
  "	NONSPECIFIC_CURRENT i\n"
  "}\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "\n"
  "	Cmax        = 1.        (mM)        : transmitter concentration during release event\n"
  "	Cdur	    = 0.3       (ms)		: transmitter duration (rising phase)\n"
  "	kon         = 86.89     (/mM/ms)    : unbound receptor ligand-binding rate\n"
  "    koff        = 0.69      (/ms)       : bound receptor ligand-unbinding rate\n"
  "    CC          = 9.64      (/ms)       : bound receptor cleft closing rate\n"
  "    CO          = 2.60      (/ms)       : bound receptor cleft opening rate\n"
  "    Beta	    = 0.68      (/ms)	    : channel opening rate\n"
  "    Alpha       = 0.079     (/ms)       : open channel closing rate\n"
  "	Erev	    = 0.        (mV)		: reversal potential\n"
  "	Kd          = 9.98      (mM)        : modulate Mg concentration dependence\n"
  "    gamma       = 0.101     (/mV)       : modulate slope of Mg sensitivity\n"
  "    mg          = 1.0       (mM)        : extracellular Mg concentration\n"
  "    gmax	    = 0.003026  (umho)	    : maximum conductance\n"
  "    kin_scale   = 1.83      (1)         : scale voltage sensitivity of decay kinetics\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)		: postsynaptic voltage\n"
  "	i 		(nA)		: current = g * (v - Erev)\n"
  "	g 		(umho)		: conductance\n"
  "	C       (mM)        : unbound transmitter concentration\n"
  "    B                   : fraction of channels not blocked by extracellular Mg\n"
  "    kin_factor          : exponential scaling of decay rates with voltage\n"
  "    scale               : allow netcon weight to scale conductance\n"
  "}\n"
  "\n"
  "STATE {\n"
  "    Ru                  : fraction of receptors not bound to transmitter\n"
  "    Rb                  : fraction of receptors bound to transmitter\n"
  "    Rc                  : fraction of receptors in closed cleft state\n"
  "    Ro                  : fraction of channels in open state\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	C = 0.\n"
  "    Ru = 1.\n"
  "    Rb = 0.\n"
  "    Rc = 0.\n"
  "    Ro = 0.\n"
  "    scale = 1.\n"
  "    mgblock(v)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    mgblock(v)\n"
  "    SOLVE kstates METHOD sparse\n"
  "	g = scale * gmax * B * Ro\n"
  "    i = g * (v - Erev)\n"
  "}\n"
  "\n"
  "KINETIC kstates {\n"
  "    ~ Ru <-> Rb     (C * kon, koff)\n"
  "    ~ Rb <-> Rc     (CC, CO * kin_factor)\n"
  "    ~ Rc <-> Ro     (Beta, Alpha * kin_factor)\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight) {\n"
  "    if (flag == 0) { : a new spike received\n"
  "        C = Cmax\n"
  "        scale = weight\n"
  "        net_send(Cdur, 1)\n"
  "    } else {    : a self event\n"
  "        C = 0.\n"
  "    }\n"
  "}\n"
  "\n"
  "PROCEDURE mgblock(v(mV)) {\n"
  "	: from Jahr & Stevens\n"
  "    B = 1. / (1. + exp(gamma * (-v)) * (mg / Kd))\n"
  "    kin_factor = (1 - kin_scale) * B + kin_scale\n"
  "}\n"
  ;
#endif
