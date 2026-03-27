/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
 
#define nrn_init _nrn_init__cat
#define _nrn_initial _nrn_initial__cat
#define nrn_cur _nrn_cur__cat
#define _nrn_current _nrn_current__cat
#define nrn_jacob _nrn_jacob__cat
#define nrn_state _nrn_state__cat
#define _net_receive _net_receive__cat 
#define rates rates__cat 
#define states states__cat 
 
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
#define gcatbar _p[0]
#define gcatbar_columnindex 0
#define iCa _p[1]
#define iCa_columnindex 1
#define minf _p[2]
#define minf_columnindex 2
#define hinf _p[3]
#define hinf_columnindex 3
#define taum _p[4]
#define taum_columnindex 4
#define tauh _p[5]
#define tauh_columnindex 5
#define m _p[6]
#define m_columnindex 6
#define h _p[7]
#define h_columnindex 7
#define cai _p[8]
#define cai_columnindex 8
#define cao _p[9]
#define cao_columnindex 9
#define Dm _p[10]
#define Dm_columnindex 10
#define Dh _p[11]
#define Dh_columnindex 11
#define gcat _p[12]
#define gcat_columnindex 12
#define _g _p[13]
#define _g_columnindex 13
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_iCa	*_ppvar[2]._pval
#define _ion_diCadv	*_ppvar[3]._pval
 
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
 static void _hoc_KTF(void);
 static void _hoc_alpm(void);
 static void _hoc_alph(void);
 static void _hoc_betm(void);
 static void _hoc_beth(void);
 static void _hoc_efun(void);
 static void _hoc_ghk(void);
 static void _hoc_h2(void);
 static void _hoc_rates(void);
 static void _hoc_states(void);
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
 "setdata_cat", _hoc_setdata,
 "KTF_cat", _hoc_KTF,
 "alpm_cat", _hoc_alpm,
 "alph_cat", _hoc_alph,
 "betm_cat", _hoc_betm,
 "beth_cat", _hoc_beth,
 "efun_cat", _hoc_efun,
 "ghk_cat", _hoc_ghk,
 "h2_cat", _hoc_h2,
 "rates_cat", _hoc_rates,
 "states_cat", _hoc_states,
 0, 0
};
#define KTF KTF_cat
#define _f_betm _f_betm_cat
#define _f_alpm _f_alpm_cat
#define _f_beth _f_beth_cat
#define _f_alph _f_alph_cat
#define alpm alpm_cat
#define alph alph_cat
#define betm betm_cat
#define beth beth_cat
#define efun efun_cat
#define ghk ghk_cat
#define h2 h2_cat
 extern double KTF( double );
 extern double _f_betm( double );
 extern double _f_alpm( double );
 extern double _f_beth( double );
 extern double _f_alph( double );
 extern double alpm( double );
 extern double alph( double );
 extern double betm( double );
 extern double beth( double );
 extern double efun( double );
 extern double ghk( double , double , double );
 extern double h2( double );
 /* declare global and static user variables */
#define eca eca_cat
 double eca = 140;
#define ki ki_cat
 double ki = 0.001;
#define tfi tfi_cat
 double tfi = 0.68;
#define tfa tfa_cat
 double tfa = 1;
#define tBase tBase_cat
 double tBase = 23.5;
#define usetable usetable_cat
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_cat", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tBase_cat", "degC",
 "ki_cat", "mM",
 "gcatbar_cat", "mho/cm2",
 "iCa_cat", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tBase_cat", &tBase_cat,
 "ki_cat", &ki_cat,
 "tfa_cat", &tfa_cat,
 "tfi_cat", &tfi_cat,
 "eca_cat", &eca_cat,
 "usetable_cat", &usetable_cat,
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
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"cat",
 "gcatbar_cat",
 0,
 "iCa_cat",
 "minf_cat",
 "hinf_cat",
 "taum_cat",
 "tauh_cat",
 0,
 "m_cat",
 "h_cat",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _Ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gcatbar = 0;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 prop_ion = need_memb(_Ca_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* iCa */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_diCadv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _cat_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	ion_reg("Ca", 2.0);
 	_ca_sym = hoc_lookup("ca_ion");
 	_Ca_sym = hoc_lookup("Ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "Ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "Ca_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cat /Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/cat.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double KTOMV = .0853;
 static double *_t_alph;
 static double *_t_beth;
 static double *_t_alpm;
 static double *_t_betm;
 static double _zfacm , _zfach ;
static int _reset;
static char *modelname = "t-type calcium channel with high threshold for activation";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
static int states();
 static double _n_betm(double);
 static double _n_alpm(double);
 static double _n_beth(double);
 static double _n_alph(double);
 
double h2 (  double _lcai ) {
   double _lh2;
 _lh2 = ki / ( ki + _lcai ) ;
   
return _lh2;
 }
 
static void _hoc_h2(void) {
  double _r;
   _r =  h2 (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghk (  double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lnu , _lf ;
 _lf = KTF ( _threadargscomma_ celsius ) / 2.0 ;
   _lnu = _lv / _lf ;
   _lghk = - _lf * ( 1. - ( _lci / _lco ) * exp ( _lnu ) ) * efun ( _threadargscomma_ _lnu ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   _r =  ghk (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double KTF (  double _lcelsius ) {
   double _lKTF;
 _lKTF = ( ( 25. / 293.15 ) * ( _lcelsius + 273.15 ) ) ;
   
return _lKTF;
 }
 
static void _hoc_KTF(void) {
  double _r;
   _r =  KTF (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double efun (  double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   _r =  efun (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_alph, _tmin_alph;
 static void _check_alph();
 static void _check_alph() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_alph =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_alph)/200.; _mfac_alph = 1./_dx;
   for (_i=0, _x=_tmin_alph; _i < 201; _x += _dx, _i++) {
    _t_alph[_i] = _f_alph(_x);
   }
  }
 }

 double alph(double _lv){ _check_alph();
 return _n_alph(_lv);
 }

 static double _n_alph(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_alph(_lv); 
}
 _xi = _mfac_alph * (_lv - _tmin_alph);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_alph[0];
 }
 if (_xi >= 200.) {
 return _t_alph[200];
 }
 _i = (int) _xi;
 return _t_alph[_i] + (_xi - (double)_i)*(_t_alph[_i+1] - _t_alph[_i]);
 }

 
double _f_alph (  double _lv ) {
   double _lalph;
 _lalph = 1.6e-4 * exp ( - ( _lv + 57.0 ) / 19.0 ) ;
   
return _lalph;
 }
 
static void _hoc_alph(void) {
  double _r;
    _r =  alph (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_beth, _tmin_beth;
 static void _check_beth();
 static void _check_beth() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_beth =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_beth)/200.; _mfac_beth = 1./_dx;
   for (_i=0, _x=_tmin_beth; _i < 201; _x += _dx, _i++) {
    _t_beth[_i] = _f_beth(_x);
   }
  }
 }

 double beth(double _lv){ _check_beth();
 return _n_beth(_lv);
 }

 static double _n_beth(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_beth(_lv); 
}
 _xi = _mfac_beth * (_lv - _tmin_beth);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_beth[0];
 }
 if (_xi >= 200.) {
 return _t_beth[200];
 }
 _i = (int) _xi;
 return _t_beth[_i] + (_xi - (double)_i)*(_t_beth[_i+1] - _t_beth[_i]);
 }

 
double _f_beth (  double _lv ) {
   double _lbeth;
 _lbeth = 1.0 / ( exp ( ( - _lv + 15.0 ) / 10.0 ) + 1.0 ) ;
   
return _lbeth;
 }
 
static void _hoc_beth(void) {
  double _r;
    _r =  beth (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_alpm, _tmin_alpm;
 static void _check_alpm();
 static void _check_alpm() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_alpm =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_alpm)/200.; _mfac_alpm = 1./_dx;
   for (_i=0, _x=_tmin_alpm; _i < 201; _x += _dx, _i++) {
    _t_alpm[_i] = _f_alpm(_x);
   }
  }
 }

 double alpm(double _lv){ _check_alpm();
 return _n_alpm(_lv);
 }

 static double _n_alpm(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_alpm(_lv); 
}
 _xi = _mfac_alpm * (_lv - _tmin_alpm);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_alpm[0];
 }
 if (_xi >= 200.) {
 return _t_alpm[200];
 }
 _i = (int) _xi;
 return _t_alpm[_i] + (_xi - (double)_i)*(_t_alpm[_i+1] - _t_alpm[_i]);
 }

 
double _f_alpm (  double _lv ) {
   double _lalpm;
 _lalpm = 0.1967 * ( - 1.0 * _lv + 19.88 ) / ( exp ( ( - 1.0 * _lv + 19.88 ) / 10.0 ) - 1.0 ) ;
   
return _lalpm;
 }
 
static void _hoc_alpm(void) {
  double _r;
    _r =  alpm (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_betm, _tmin_betm;
 static void _check_betm();
 static void _check_betm() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_betm =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_betm)/200.; _mfac_betm = 1./_dx;
   for (_i=0, _x=_tmin_betm; _i < 201; _x += _dx, _i++) {
    _t_betm[_i] = _f_betm(_x);
   }
  }
 }

 double betm(double _lv){ _check_betm();
 return _n_betm(_lv);
 }

 static double _n_betm(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_betm(_lv); 
}
 _xi = _mfac_betm * (_lv - _tmin_betm);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_betm[0];
 }
 if (_xi >= 200.) {
 return _t_betm[200];
 }
 _i = (int) _xi;
 return _t_betm[_i] + (_xi - (double)_i)*(_t_betm[_i+1] - _t_betm[_i]);
 }

 
double _f_betm (  double _lv ) {
   double _lbetm;
 _lbetm = 0.046 * exp ( - _lv / 22.73 ) ;
   
return _lbetm;
 }
 
static void _hoc_betm(void) {
  double _r;
    _r =  betm (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int  states (  ) {
   rates ( _threadargscomma_ v ) ;
   m = m + _zfacm * ( minf - m ) ;
   h = h + _zfach * ( hinf - h ) ;
   
/*VERBATIM*/
        return 0;
  return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _la ;
 _la = alpm ( _threadargscomma_ _lv ) ;
   taum = 1.0 / ( tfa * ( _la + betm ( _threadargscomma_ _lv ) ) ) ;
   minf = _la / ( _la + betm ( _threadargscomma_ _lv ) ) ;
   _zfacm = ( 1.0 - exp ( - dt / taum ) ) ;
   _la = alph ( _threadargscomma_ _lv ) ;
   tauh = 1.0 / ( tfi * ( _la + beth ( _threadargscomma_ _lv ) ) ) ;
   hinf = _la / ( _la + beth ( _threadargscomma_ _lv ) ) ;
   _zfach = ( 1.0 - exp ( - dt / tauh ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("cat", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_Ca_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   gcat = gcatbar * m * m * h * h2 ( _threadargscomma_ cai ) ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gcat = gcatbar * m * m * h * h2 ( _threadargscomma_ cai ) ;
   iCa = gcat * ghk ( _threadargscomma_ v , cai , cao ) ;
   }
 _current += iCa;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _diCa;
  _diCa = iCa;
 _rhs = _nrn_current(_v);
  _ion_diCadv += (_diCa - iCa)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_iCa += iCa ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
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

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
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
  cai = _ion_cai;
  cao = _ion_cao;
 { error =  states();
 if(error){fprintf(stderr,"at line 60 in file cat.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_alph = makevector(201*sizeof(double));
   _t_beth = makevector(201*sizeof(double));
   _t_alpm = makevector(201*sizeof(double));
   _t_betm = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/samgritz/Desktop/Rutgers/Milstein_Lab/Code/Rutgers-Neuroscience-PhD/Biophysical_model_sims/CA1Sim-minimal/cat.mod";
static const char* nmodl_file_text = 
  "TITLE t-type calcium channel with high threshold for activation\n"
  ": used in somatic and dendritic regions \n"
  ": it calculates I_Ca using channel permeability instead of conductance\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "\n"
  "	FARADAY = 96520 (coul)\n"
  "	R = 8.3134 (joule/degK)\n"
  "	KTOMV = .0853 (mV/degC)\n"
  "}\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "PARAMETER {           :parameters that can be entered when function is called in cell-setup \n"
  "        dt            (ms)\n"
  "	v             (mV)\n"
  "        tBase = 23.5  (degC)\n"
  "	celsius = 22  (degC)\n"
  "	gcatbar = 0   (mho/cm2)  : initialized conductance\n"
  "	ki = 0.001    (mM)\n"
  "	cai = 5.e-5   (mM)       : initial internal Ca++ concentration\n"
  "	cao = 2       (mM)       : initial external Ca++ concentration\n"
  "        tfa = 1                  : activation time constant scaling factor\n"
  "        tfi = 0.68               : inactivation time constant scaling factor\n"
  "        eca = 140                : Ca++ reversal potential\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX cat\n"
  "	USEION ca READ cai,cao \n"
  "        USEION Ca WRITE iCa VALENCE 2\n"
  "        : The T-current does not activate calcium-dependent currents.\n"
  "        : The construction with dummy ion Ca prevents the updating of the \n"
  "        : internal calcium concentration. \n"
  "        RANGE gcatbar, hinf, minf, taum, tauh, iCa\n"
  "}\n"
  "\n"
  "STATE {	m h }  : unknown activation and inactivation parameters to be solved in the DEs \n"
  "\n"
  "ASSIGNED {     : parameters needed to solve DE\n"
  "	iCa (mA/cm2)\n"
  "        gcat  (mho/cm2) \n"
  "        minf\n"
  "        hinf\n"
  "        taum\n"
  "        tauh\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  ":        tadj = 3^((celsius-tBase)/10)   : assume Q10 of 3\n"
  "	rates(v)\n"
  "        m = minf\n"
  "        h = hinf\n"
  "	gcat = gcatbar*m*m*h*h2(cai)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "	gcat = gcatbar*m*m*h*h2(cai) : maximum channel permeability\n"
  "	iCa = gcat*ghk(v,cai,cao)    : dummy calcium current induced by this channel\n"
  "\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "FUNCTION h2(cai(mM)) {\n"
  "	h2 = ki/(ki+cai)\n"
  "}\n"
  "\n"
  "FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) { LOCAL nu,f\n"
  "        f = KTF(celsius)/2\n"
  "        nu = v/f\n"
  "        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)\n"
  "}\n"
  "\n"
  "FUNCTION KTF(celsius (degC)) (mV) {   : temperature-dependent adjustment factor\n"
  "        KTF = ((25./293.15)*(celsius + 273.15))\n"
  "}\n"
  "\n"
  "FUNCTION efun(z) {\n"
  "	if (fabs(z) < 1e-4) {\n"
  "		efun = 1 - z/2\n"
  "	}else{\n"
  "		efun = z/(exp(z) - 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION alph(v(mV)) {\n"
  "	TABLE FROM -150 TO 150 WITH 200\n"
  "	alph = 1.6e-4*exp(-(v+57)/19)\n"
  "}\n"
  "\n"
  "FUNCTION beth(v(mV)) {\n"
  "        TABLE FROM -150 TO 150 WITH 200\n"
  "	beth = 1/(exp((-v+15)/10)+1.0)\n"
  "}\n"
  "\n"
  "FUNCTION alpm(v(mV)) {\n"
  "	TABLE FROM -150 TO 150 WITH 200\n"
  "	alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)\n"
  "}\n"
  "\n"
  "FUNCTION betm(v(mV)) {\n"
  "	TABLE FROM -150 TO 150 WITH 200\n"
  "	betm = 0.046*exp(-v/22.73)\n"
  "}\n"
  "\n"
  "UNITSON\n"
  "LOCAL facm,fach\n"
  "\n"
  ":if state_cagk is called from hoc, garbage or segmentation violation will\n"
  ":result because range variables won't have correct pointer.  This is because\n"
  ": only BREAKPOINT sets up the correct pointers to range variables.\n"
  "PROCEDURE states() {     : exact when v held constant; integrates over dt step\n"
  "        rates(v)\n"
  "        m = m + facm*(minf - m)\n"
  "        h = h + fach*(hinf - h)\n"
  "        VERBATIM\n"
  "        return 0;\n"
  "        ENDVERBATIM\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v (mV)) { :callable from hoc\n"
  "        LOCAL a\n"
  "        a = alpm(v)\n"
  "        taum = 1/(tfa*(a + betm(v))) : estimation of activation tau\n"
  "        minf =  a/(a+betm(v))        : estimation of activation steady state\n"
  "        facm = (1 - exp(-dt/taum))\n"
  "        a = alph(v)\n"
  "        tauh = 1/(tfi*(a + beth(v))) : estimation of inactivation tau\n"
  "        hinf = a/(a+beth(v))         : estimation of inactivation steady state\n"
  "        fach = (1 - exp(-dt/tauh))\n"
  "}\n"
  ;
#endif
