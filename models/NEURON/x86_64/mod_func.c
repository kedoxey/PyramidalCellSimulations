#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _CaDynamics_E2_NML2__decay122__gamma5_09Emin4_reg(void);
extern void _CaDynamics_E2_NML2__decay460__gamma5_01Emin4_reg(void);
extern void _Ca_HVA_reg(void);
extern void _Ca_LVAst_reg(void);
extern void _Ih_reg(void);
extern void _Im_reg(void);
extern void _K_Pst_reg(void);
extern void _K_Tst_reg(void);
extern void _MyExp2SynBB_reg(void);
extern void _MyExp2SynNMDABB_reg(void);
extern void _Nap_Et2_reg(void);
extern void _NaTa_t_reg(void);
extern void _pas_nml2_reg(void);
extern void _ProbAMPA2_reg(void);
extern void _ProbAMPANMDA2_reg(void);
extern void _ProbNMDA2_reg(void);
extern void _ProbUDFsyn2_reg(void);
extern void _SAM_gaussstim_reg(void);
extern void _SK_E2_reg(void);
extern void _SKv3_1_reg(void);
extern void _vecevent_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"./CaDynamics_E2_NML2__decay122__gamma5_09Emin4.mod\"");
    fprintf(stderr," \"./CaDynamics_E2_NML2__decay460__gamma5_01Emin4.mod\"");
    fprintf(stderr," \"./Ca_HVA.mod\"");
    fprintf(stderr," \"./Ca_LVAst.mod\"");
    fprintf(stderr," \"./Ih.mod\"");
    fprintf(stderr," \"./Im.mod\"");
    fprintf(stderr," \"./K_Pst.mod\"");
    fprintf(stderr," \"./K_Tst.mod\"");
    fprintf(stderr," \"./MyExp2SynBB.mod\"");
    fprintf(stderr," \"./MyExp2SynNMDABB.mod\"");
    fprintf(stderr," \"./Nap_Et2.mod\"");
    fprintf(stderr," \"./NaTa_t.mod\"");
    fprintf(stderr," \"./pas_nml2.mod\"");
    fprintf(stderr," \"./ProbAMPA2.mod\"");
    fprintf(stderr," \"./ProbAMPANMDA2.mod\"");
    fprintf(stderr," \"./ProbNMDA2.mod\"");
    fprintf(stderr," \"./ProbUDFsyn2.mod\"");
    fprintf(stderr," \"./SAM_gaussstim.mod\"");
    fprintf(stderr," \"./SK_E2.mod\"");
    fprintf(stderr," \"./SKv3_1.mod\"");
    fprintf(stderr," \"./vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _CaDynamics_E2_NML2__decay122__gamma5_09Emin4_reg();
  _CaDynamics_E2_NML2__decay460__gamma5_01Emin4_reg();
  _Ca_HVA_reg();
  _Ca_LVAst_reg();
  _Ih_reg();
  _Im_reg();
  _K_Pst_reg();
  _K_Tst_reg();
  _MyExp2SynBB_reg();
  _MyExp2SynNMDABB_reg();
  _Nap_Et2_reg();
  _NaTa_t_reg();
  _pas_nml2_reg();
  _ProbAMPA2_reg();
  _ProbAMPANMDA2_reg();
  _ProbNMDA2_reg();
  _ProbUDFsyn2_reg();
  _SAM_gaussstim_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
  _vecevent_reg();
}
