import matplotlib.pyplot as plt
import os
import json
import model_helpers as mh
from neuron import h, load_mechanisms
from netpyne import specs, sim


### Parameters for simulation ###
test_name = ''
hoc_fname = 'L5PC'
vinit = -80

input_amps = [0.35477, 0.44346, 0.53215, 1.0643]
amp_idx = 0  # if -1 then no current injection
input_amp = input_amps[amp_idx]
test_label = f'input_{amp_idx}-{test_name}'

### Model information ###
code_version = 'Hay'
model_version = 'NEURON'

nmldb_id = 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
model_name = f'{nmldb_id}-{model_version}'

cell_model = 'Hay2011'  # 'Hay2011'
cell_type = 'PYR'
cell_name = 'L5PC'  # 'L5PC'
cell_label = cell_name+'_hoc'

### Define paths ###
cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models') 
model_dir = os.path.join(models_dir, model_version, model_name)  # 'L5bPCmodelsEH')
hocs_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir,'models')
mod_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir, 'mod')

hoc_file = os.path.join(hocs_dir, f'{cell_name}{test_name}.hoc')

### Download model ###
mh.download_from_nmldb(nmldb_id, model_version)

output_dir = os.path.join(model_dir,'output')
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

### Compile mechs ###
mh.compile_mechs(cwd,hocs_dir,mod_dir)
load_mechanisms(model_dir)

### Network Parameters ###
netParams = specs.NetParams()

### Import cell into NetPyNE ###
importedCellParams = netParams.importCellParams(label=cell_label,
                                                conds={'cellType': cell_type, 'cellModel': cell_model},
                                                fileName=hoc_file,
                                                cellName=cell_name
                                                )

for sec in importedCellParams['secs']:
    importedCellParams[sec]['vinit'] = vinit

### Create population ###
pop_label = cell_label+'_pop'
netParams.popParams[pop_label] = {'cellType': cell_type, 
                                  'cellModel': cell_model,
                                  'numCells': 1}

### Add input ###
if amp_idx > -1:
    netParams.stimSourceParams['Input_IC'] = {
        'type': 'IClamp',
        'del': 700,
        'dur': 2000,
        'amp': input_amp  
    }

    netParams.stimTargetParams['Input_IC->Soma'] = {
        'source': 'Input_IC',
        'sec': 'soma_0',
        'loc': 0.5,
        'conds': {'pop': pop_label}
    }

    print(f'current clamp added at {input_amp} pA')
else:
    print('no input')

### Simulation configuration ###
cfg = specs.SimConfig()					            # object of class SimConfig to store simulation configuration
cfg.duration = 3000 						            # Duration of the simulation, in ms
cfg.dt = 0.01								                # Internal integration timestep to use
cfg.verbose = True							                # Show detailed messages
cfg.recordTraces = {'V_soma':{'sec':'soma_0','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStep = 0.01
cfg.filename = os.path.join(output_dir,cell_name+'_'+test_label) 			# Set file output name
cfg.saveJson = False
cfg.analysis['plotTraces'] = {'include': ['all'], 'saveFig': True} # Plot recorded traces for this list of cells
cfg.hParams['celsius'] = 34.0 
cfg.hParams['v_init'] = vinit

### Run simulation ###
sim.createSimulateAnalyze(netParams = netParams, simConfig = cfg)

### Plot morphology ###
sim.analysis.plotShape(showSyns=True, dist=0.8, saveFig=True, axisLabels=True)
