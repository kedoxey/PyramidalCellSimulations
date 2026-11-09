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

    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/AmpaNmdaSyn.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/CaPool.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/CaT.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/GabaSyn.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/gapjunction.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/Ih.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/KA.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/KCa.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/Kd.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/KM.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/Kslow.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/LCa.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/Na.mod\"");
    fprintf(stderr," \"/home/kedoxey/CRCNS/PyramidalCellSimulations/models/local/Birgiolas2020/Mechanisms/VecStim.mod\"");
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
