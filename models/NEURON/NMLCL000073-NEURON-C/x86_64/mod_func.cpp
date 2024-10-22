#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _CaDynamics_E2_NML2__decay122__gamma5_09Emin4_reg(void);
extern void _CaDynamics_E2_NML2__decay460__gamma5_01Emin4_reg(void);
extern void _Ca_HVA_reg(void);
extern void _Ca_LVAst_reg(void);
extern void _Ih_reg(void);
extern void _Im_reg(void);
extern void _K_Pst_reg(void);
extern void _K_Tst_reg(void);
extern void _Nap_Et2_reg(void);
extern void _NaTa_t_reg(void);
extern void _pas_nml2_reg(void);
extern void _SK_E2_reg(void);
extern void _SKv3_1_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/CaDynamics_E2_NML2__decay122__gamma5_09Emin4.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/CaDynamics_E2_NML2__decay460__gamma5_01Emin4.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/Ca_HVA.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/Ca_LVAst.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/Ih.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/Im.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/K_Pst.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/K_Tst.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/Nap_Et2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/NaTa_t.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/pas_nml2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/SK_E2.mod\"");
    fprintf(stderr, " \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/NEURON/NMLCL000073-NEURON-C/SKv3_1.mod\"");
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
  _Nap_Et2_reg();
  _NaTa_t_reg();
  _pas_nml2_reg();
  _SK_E2_reg();
  _SKv3_1_reg();
}

#if defined(__cplusplus)
}
#endif
