import matplotlib.pyplot as plt
import os
import json
import model_helpers as mh
from neuron import h, load_mechanisms
from netpyne import specs, sim
import time

time_flag = False
start_t = time.time()

### Parameters for simulation ###
run_NML = False
sim_name = 'test_syn-a1'
hoc_fname = 'L5PC'
vinit = -80

### Input parameters ###
input_amps = [0, 0.35477, 0.44346, 0.53215, 1.0643]
amp_idx = 0  
input_amp = input_amps[amp_idx]
in_delay = 700
in_dur = 2000

sim_label = f'input_{round(input_amp,1)}-{sim_name}'

### Model information ###
code_version = 'Hay'
model_version = 'NeuroML' if run_NML else 'NEURON'

nmldb_id = 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
model_name = f'{nmldb_id}-{model_version}'

cell_model = 'Hay2011'  # 'Hay2011'
cell_type = 'PYR'
cell_name = 'L5PC'  # 'L5PC'
cell_label = cell_name+'_hoc'
pop_label = cell_name+'_Pop'

### Define paths ###
cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models') 
model_dir = os.path.join(models_dir, model_version, model_name)  # 'L5bPCmodelsEH')
hocs_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir,'models')
mod_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir, 'mod')

hoc_file = os.path.join(hocs_dir, f'{hoc_fname}.hoc')

### Download model ###
mh.download_from_nmldb(nmldb_id, model_version)

### Get output directories ###
output_dir, sim_dir = mh.create_output_dirs(sim_name, model_dir)

### Generate network if running NeuroML ###
if run_NML:
    net_nml_file = mh.generate_network(model_dir, cell_name, pop_label, 
                                       force=True, 
                                       input_amp=input_amp, 
                                       start=in_delay, 
                                       stop=in_delay+in_dur)

### Compile mechs ###
mh.compile_mechs(cwd,hocs_dir,mod_dir)
load_mechanisms(model_dir)

### Simulation configuration ###
cfg = specs.SimConfig()					                    # object of class SimConfig to store simulation configuration
cfg.duration = 3000 						                # Duration of the simulation, in ms
cfg.dt = 0.01								                # Internal integration timestep to use
cfg.verbose = True							                # Show detailed messages
cfg.recordTraces = {'V_soma':{'sec':'soma_0','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStep = 0.01
cfg.filename = os.path.join(sim_dir,cell_name+'_'+sim_label) 	# Set file output name
cfg.saveJson = False
cfg.analysis['plotTraces'] = {'include': [pop_label], 'saveFig': True}  # Plot recorded traces for this list of cells
cfg.hParams['celsius'] = 34.0 
cfg.hParams['v_init'] = vinit

### Import cell ###
if run_NML:
    gid, netParams = sim.importNeuroML2(net_nml_file, simConfig=cfg, simulate=False, analyze=False, return_net_params_also=True)
    # sim.importNeuroML2SimulateAnalyze(net_nml_file, simConfig=cfg)
else:
    netParams = specs.NetParams()

    importedCellParams = netParams.importCellParams(label=cell_label,
                                                    conds={'cellType': cell_type, 'cellModel': cell_model},
                                                    fileName=hoc_file,
                                                    cellName=cell_name
                                                    )

    for sec in importedCellParams['secs']:
        importedCellParams[sec]['vinit'] = vinit

    ### Create population ###
    netParams.popParams[pop_label] = {'cellType': cell_type, 
                                    'cellModel': cell_model,
                                    'numCells': 1}
    
# mh.get_components(netParams.cellParams[cell_label])
### Get all sections ###
syn_secs = mh.get_components(importedCellParams, 'all')

### Add AMPA/NMDA synapse ###
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}

exc_syns = ['AMPA', 'NMDA']

### Add synaptic input ###
syn_method = 'cell'  # 'stim'

if 'cell' in syn_method:
    netParams.popParams['vecstim'] = {
        'cellModel': 'VecStim',
        'numCells': 1,
        'spikePattern': {'type': 'poisson',
                         'start': 10,
                         'stop': -1,
                         'frequency': 50}
    }
    
    netParams.connParams[f'vecstim->{pop_label}'] = {
        'preConds': {'pop': 'vecstim'},
        'postConds': {'pop': pop_label},
        'sec': syn_secs,
        'synMech': exc_syns,
        'weight': 1,
        'delay': 5,
        'probability': 1.0
    }

    netParams.subConnParams[f'vecstim->{pop_label}'] = {
        'preConds': {'pop': 'vecstim'},
        'postConds': {'pop': pop_label},
        'sec': all,
        'groupSynMech': exc_syns,
        'density': 'uniform'
    }
else:
    f = 60  # Hz, frequency of input
    p_mean = 1000/f  # ms, mean time between spikes
    netParams.stimSourceParams['Input_syn'] = {
        'type': 'NetStim',
        'interval': f'poisson({p_mean})',
        'number': 1e9,
        'start': 0,
        'noise': 1
        # 'rate': 100,
        # 'noise': 0.5
    }
    netParams.stimTargetParams['Input_syn->soma'] = {
        'source': 'Input_syn',
        'conds': {'pop': pop_label},
        'weight': 1,
        'delay': 5,
        'synMech': 'AMPANMDA',
        'sec': all
    }

### Add input ###
netParams.stimSourceParams['Input_IC'] = {
    'type': 'IClamp',
    'del': in_delay,
    'dur': in_dur,
    'amp': input_amp  
}

netParams.stimTargetParams['Input_IC->Soma'] = {
    'source': 'Input_IC',
    'sec': 'soma_0',
    'loc': 0.5,
    'conds': {'pop': pop_label}
}

### Run simulation ###
sim.createSimulateAnalyze(netParams=netParams, simConfig=cfg, output=False)

### Plot morphology ###
# sim.analysis.plotShape(showSyns=True, dist=0.8, saveFig=True, axisLabels=True)

### Simulation time ###
if time_flag:
    end_t = time.time()
    final_t = end_t - start_t
    m, s = divmod(final_t, 60)
    time_path = os.path.join(sim_dir,f'{model_version}_{sim_label}_sim-time.txt')

    with open(time_path, "w") as text_file:
        text_file.write(f'Simulation time = {m} m {s} s')
    print(f'Simulation time = {m} m {s} s')
