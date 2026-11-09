#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _AmpaNmdaSyn_reg(void);
extern void _CaPool_reg(void);
extern void _CaT_reg(void);
extern void _GabaSyn_reg(void);
extern void _gapjunction_reg(void);
extern void _Ih_reg(void);
extern void _KA_reg(void);
extern void _KCa_reg(void);
extern void _Kd_reg(void);
extern void _KM_reg(void);
extern void _Kslow_reg(void);
extern void _LCa_reg(void);
extern void _Na_reg(void);
extern void _VecStim_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"./AmpaNmdaSyn.mod\"");
    fprintf(stderr," \"./CaPool.mod\"");
    fprintf(stderr," \"./CaT.mod\"");
    fprintf(stderr," \"./GabaSyn.mod\"");
    fprintf(stderr," \"./gapjunction.mod\"");
    fprintf(stderr," \"./Ih.mod\"");
    fprintf(stderr," \"./KA.mod\"");
    fprintf(stderr," \"./KCa.mod\"");
    fprintf(stderr," \"./Kd.mod\"");
    fprintf(stderr," \"./KM.mod\"");
    fprintf(stderr," \"./Kslow.mod\"");
    fprintf(stderr," \"./LCa.mod\"");
    fprintf(stderr," \"./Na.mod\"");
    fprintf(stderr," \"./VecStim.mod\"");
    fprintf(stderr, "\n");
  }
  _AmpaNmdaSyn_reg();
  _CaPool_reg();
  _CaT_reg();
  _GabaSyn_reg();
  _gapjunction_reg();
  _Ih_reg();
  _KA_reg();
  _KCa_reg();
  _Kd_reg();
  _KM_reg();
  _Kslow_reg();
  _LCa_reg();
  _Na_reg();
  _VecStim_reg();
}
