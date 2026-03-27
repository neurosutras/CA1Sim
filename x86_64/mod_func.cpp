#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _ampa_kin_reg(void);
extern void _cad_reg(void);
extern void _cagk_reg(void);
extern void _cal_reg(void);
extern void _calH_reg(void);
extern void _car_reg(void);
extern void _cat_reg(void);
extern void _gaba_a_kin_reg(void);
extern void _gabab_reg(void);
extern void _h_reg(void);
extern void _kad_reg(void);
extern void _kap_reg(void);
extern void _kca_reg(void);
extern void _kdr_reg(void);
extern void _km_reg(void);
extern void _km2_reg(void);
extern void _nas_reg(void);
extern void _nax_reg(void);
extern void _nmda_kin5_reg(void);
extern void _pr_reg(void);
extern void _vecevent_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"./ampa_kin.mod\"");
    fprintf(stderr, " \"./cad.mod\"");
    fprintf(stderr, " \"./cagk.mod\"");
    fprintf(stderr, " \"./cal.mod\"");
    fprintf(stderr, " \"./calH.mod\"");
    fprintf(stderr, " \"./car.mod\"");
    fprintf(stderr, " \"./cat.mod\"");
    fprintf(stderr, " \"./gaba_a_kin.mod\"");
    fprintf(stderr, " \"./gabab.mod\"");
    fprintf(stderr, " \"./h.mod\"");
    fprintf(stderr, " \"./kad.mod\"");
    fprintf(stderr, " \"./kap.mod\"");
    fprintf(stderr, " \"./kca.mod\"");
    fprintf(stderr, " \"./kdr.mod\"");
    fprintf(stderr, " \"./km.mod\"");
    fprintf(stderr, " \"./km2.mod\"");
    fprintf(stderr, " \"./nas.mod\"");
    fprintf(stderr, " \"./nax.mod\"");
    fprintf(stderr, " \"./nmda_kin5.mod\"");
    fprintf(stderr, " \"./pr.mod\"");
    fprintf(stderr, " \"./vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _ampa_kin_reg();
  _cad_reg();
  _cagk_reg();
  _cal_reg();
  _calH_reg();
  _car_reg();
  _cat_reg();
  _gaba_a_kin_reg();
  _gabab_reg();
  _h_reg();
  _kad_reg();
  _kap_reg();
  _kca_reg();
  _kdr_reg();
  _km_reg();
  _km2_reg();
  _nas_reg();
  _nax_reg();
  _nmda_kin5_reg();
  _pr_reg();
  _vecevent_reg();
}

#if defined(__cplusplus)
}
#endif
