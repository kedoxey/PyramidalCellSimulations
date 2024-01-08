from neuron.units import mV, ms
import matplotlib.pyplot as plt
import os
import json
import model_helpers as mh
from neuron import h, load_mechanisms

test_name = 'soma_pas'

code_version = 'Hay'
flag = ''  # -Converted for NML->NEURON with pnml

model_version = 'NEURON'

nmldb_id = 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
model_name = f'{nmldb_id}-{model_version}{flag[:2]}'

cell_model = 'Hay2011'  # 'Hay2011'
cell_type = 'PYR'
cell_name = 'L5PC'  # 'L5PC'

hoc_name = cell_name  # 'L5PC'

cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models')
model_dir = os.path.join(models_dir, model_version, model_name)  # 'L5bPCmodelsEH')
hocs_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir,'models')
mod_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir, 'mod')

output_dir = os.path.join(model_dir,'output')

mh.compile_mechs(cwd,hocs_dir,mod_dir)
load_mechanisms(model_dir)

hoc_file = os.path.join(hocs_dir, f'{cell_name}_{test_name}.hoc')

h.load_file('stdrun.hoc')
h.load_file(hoc_file)

h.cvode_active(0)
h.tstop = 3000 * ms
h.celsius = 34
h.steps_per_ms = 100
h.dt = 1.0 / h.steps_per_ms
h.v_init = -80 * mV
# h.finitialize(-80 * mV)

hoc_cell = getattr(h, cell_name)
hoc_cell = hoc_cell()

cell = hoc_cell
soma = cell.soma_0

t = h.Vector().record(h._ref_t)
v = h.Vector().record(soma(0.5)._ref_v)

h.run()

fig, axs = plt.subplots(1, 1, figsize=(8,8), sharex=True)
# axs = axs.ravel()

axs.plot(t,v)
axs.set_xlabel('time (ms)')
axs.set_ylabel('membrane potential (mV)')

fig.savefig(os.path.join(output_dir,f'{cell_name}-{test_name}-NEURON.jpg'), dpi=300)