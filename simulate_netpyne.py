import matplotlib.pyplot as plt
import os
import json
import model_helpers as mh
import numpy as np
from neuron import h, load_mechanisms
from netpyne import specs, sim
from netpyne import cell
from netpyne import support
import time

time_flag = False
start_t = time.time()

### Parameters for simulation ###
run_NML = False
plot_morphology = True
enable_syns = True
record_LFP = True

spk_type = 'poisson'  # poisson or gaussian

syns_lb = 0.1  # 0 if not distibuted by distance from soma
syns_ub = 1
syns_type = 'all'  # soma, basal, apical, basal_apical, basal_soma, apical_soma, basal_apical_soma, all
num_syns = 50

sim_name = f'MEETING-{syns_type}'
hoc_fname = 'L5PC'
vinit = -80

### Input parameters ###
input_amps = [0, 0.35477, 0.44346, 0.53215, 1.0643]
amp_idx = 0
input_amp = input_amps[amp_idx]
stim_delay = 700
stim_dur = 2000
sim_dur = 3000

# Poisson spike trains
spk_freq = 20

# Gaussian spike trains
gauss_mean = 1.4
gauss_std = 0.4

sim_label = f'{sim_name}-dist_{syns_lb}_{syns_ub}'

sim_message = 'Simulations for meeting with Sharon and Yi, varying parameters and syns_dist_scale and syns_type'

# if 'poisson' in spk_type:
#     sim_label = f'{spk_type}-freq_{spk_freq}-{sim_name}'
# else:
#     sim_label = f'{spk_type}-mean_{gauss_mean}-std_{gauss_std}-{sim_name}'

# sim_label = f'input_{round(input_amp,2)}-{sim_name}'

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

### Wrtie README containing simulation description
mh.create_sim_description(sim_dir, run_NML=run_NML, spk_type=spk_type, syns_lb=syns_lb, syns_ub=syns_ub, syns_type=syns_type, num_syns=num_syns, vinit=vinit, spk_freq=spk_freq, sim_message=sim_message)

### Generate network if running NeuroML ###
if run_NML:
    net_nml_file = mh.generate_network(model_dir, cell_name, pop_label, 
                                       force=True, 
                                       input_amp=input_amp, 
                                       start=stim_delay, 
                                       stop=stim_delay+stim_dur)

### Compile mechs ###
mh.compile_mechs(cwd,hocs_dir,mod_dir)
load_mechanisms(model_dir)

### Simulation configuration ###
cfg = specs.SimConfig()					                    # object of class SimConfig to store simulation configuration
cfg.duration = sim_dur 						                # Duration of the simulation, in ms
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


### Get sections ###
# basal, apical, basal_apical, basal_soma, apical_soma, basal_apical_soma, all
if syns_lb > 0:
    syn_secs = mh.get_secs_from_dist(hoc_file, cell_name, syns_lb, syns_ub)
else:
    syn_secs = mh.get_components(importedCellParams, syns_type)

### Add AMPA/NMDA synapse ###
netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0}

exc_syns = ['AMPA', 'NMDA']

### Add synaptic input ###
if enable_syns:
    syn_method = 'cell'  # 'stim' or 'cell'

    if 'cell' in syn_method:

        # Poisson spike pattern
        if 'poisson' in spk_type:
            netParams.popParams['vecstim'] = {
                'cellModel': 'VecStim',
                'numCells': num_syns,  # int(len(syn_secs)/4),
                'spikePattern': {'type': 'poisson',
                                'start': stim_delay,
                                'stop': stim_delay+stim_dur,
                                'frequency': spk_freq}
            }
        # Gaussian spike pattern
        else:
            netParams.popParams['vecstim'] = {
                'cellModel': 'VecStim',
                'numCells': 50,  # int(len(syn_secs)/4),
                'spikePattern': {'type': 'gauss',
                                'mu': gauss_mean,
                                'sigma': gauss_std}
            }

        netParams.connParams[f'vecstim->{pop_label}'] = {
            'preConds': {'pop': 'vecstim'},
            'postConds': {'pop': pop_label},
            'sec': syn_secs,
            'synsPerConn': 1,
            'synMech': exc_syns,
            'weight': 0.001,  # 0.5,
            # 'synMechWeightFactor': [0.5,0.5],
            'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
            'probability': 1.0,
        }

        netParams.subConnParams[f'vecstim->{pop_label}'] = {
            'preConds': {'pop': 'vecstim'},
            'postConds': {'pop': pop_label},
            'sec': syn_secs,
            'groupSynMech': exc_syns,
            'density': 'uniform'
        }
    else:
        f = 50  # Hz, frequency of input
        p_mean = 1000/f  # ms, mean time between spikes
        netParams.stimSourceParams['Input_syn'] = {
            'type': 'NetStim',
            'interval': f'poisson({p_mean})',
            'number': 1e9,
            'start': 10,
            'noise': 1
            # 'rate': 100,
            # 'noise': 0.5
        }
        netParams.stimTargetParams['Input_syn->soma'] = {
            'source': 'Input_syn',
            'conds': {'pop': pop_label},
            'weight': 1,
            'delay': 5,
            'synMech': exc_syns,
            'sec': syn_secs
        }

### Add input ###
netParams.stimSourceParams['Input_IC'] = {
    'type': 'IClamp',
    'del': stim_delay,
    'dur': stim_dur,
    'amp': input_amp  
}

netParams.stimTargetParams['Input_IC->Soma'] = {
    'source': 'Input_IC',
    'sec': 'soma_0',
    'loc': 0.5,
    'conds': {'pop': pop_label}
}

### Add linear probe ###
if record_LFP:
    probe_L = 1280//5
    channels = 1
    depths = 10
    elec_dist = probe_L//depths  # microns
    disp = 100  # 150

    elec_pos = [[x*elec_dist, (y*elec_dist - disp)*-1, 0] for x in range(channels) for y in range(depths)]  # 
    # -x is left and -y is above soma
    # elec_pos.reverse()

    cfg.recordLFP = elec_pos
    cfg.analysis['plotLFP'] = {'saveFig': True}

### Run simulation ###
sim.createSimulateAnalyze(netParams=netParams, simConfig=cfg, output=False)

### Plot morphology ###
if plot_morphology:
    sim.analysis.plotShape(showSyns=True, dist=0.8, includePre=[None], includePost=[pop_label], 
                           saveFig=True, axisLabels=True, returnPlotter=True, synColor='darkcyan') 

### Simulation time ###
if time_flag:
    end_t = time.time()
    final_t = end_t - start_t
    m, s = divmod(final_t, 60)
    time_path = os.path.join(sim_dir,f'{model_version}_{sim_label}_sim-time.txt')

    with open(time_path, "w") as text_file:
        text_file.write(f'Simulation time = {m} m {s} s')
    print(f'Simulation time = {m} m {s} s')
