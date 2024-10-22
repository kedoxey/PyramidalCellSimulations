#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _CaDynamics_E2_NML2__cADpyr_231_axonal_reg(void);
extern void _CaDynamics_E2_NML2__cADpyr_231_somatic_reg(void);
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
extern void _NaTs2_t_reg(void);
extern void _pas_nml2_reg(void);
extern void _ProbAMPA2_reg(void);
extern void _ProbAMPANMDA2_reg(void);
extern void _ProbNMDA2_reg(void);
extern void _ProbUDFsyn2_reg(void);
extern void _SK_E2_reg(void);
extern void _SKv3_1_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/CaDynamics_E2_NML2__cADpyr_231_axonal.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/CaDynamics_E2_NML2__cADpyr_231_somatic.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/Ca_HVA.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/Ca_LVAst.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/Ih.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/Im.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/K_Pst.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/K_Tst.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/MyExp2SynBB.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/MyExp2SynNMDABB.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/Nap_Et2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/NaTa_t.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/NaTs2_t.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/pas_nml2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/ProbAMPA2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/ProbAMPANMDA2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/ProbNMDA2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/ProbUDFsyn2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/SK_E2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000660-NEURON/SKv3_1.mod\"");
    fprintf(stderr, "\n");
  }
  _CaDynamics_E2_NML2__cADpyr_231_axonal_reg();
  _CaDynamics_E2_NML2__cADpyr_231_somatic_reg();
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
  _NaTs2_t_reg();
  _pas_nml2_reg();
  _ProbAMPA2_reg();
  _ProbAMPANMDA2_reg();
  _ProbNMDA2_reg();
  _ProbUDFsyn2_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
}

#if defined(__cplusplus)
}
#endif
