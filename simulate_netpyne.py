import matplotlib.pyplot as plt
import os
import secrets
import json
import model_helpers as mh
import numpy as np
from neuron import h, load_mechanisms
from netpyne import specs, sim
from netpyne import cell
from netpyne import support
from argparse import Namespace
import time

time_flag = False
start_t = time.time()

### Import simulation config ###
config_name = 'layer_config'
params = Namespace(**mh.load_config(config_name))

### Model information ###
model_version = 'NeuroML' if params.run_NML else 'NEURON'

nmldb_id =  params.nmldb_id  # 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
model_name = f'{nmldb_id}-{model_version}'

### Download model ###
mh.download_from_nmldb(nmldb_id, model_version)

### Define paths ###
cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models') 
model_dir = os.path.join(models_dir, model_version, model_name)  # 'L5bPCmodelsEH')
hocs_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir,'models')
mod_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir, 'mod')

cell_model = params.publication  # 'Hay2011'
cell_type = 'PYR'
cell_name = mh.get_cell_name(model_dir)  # 'L5PC'
cell_label = cell_name+'_hoc'
pop_label = cell_name+'_Pop'

hoc_file = os.path.join(hocs_dir, f'{cell_name}.hoc')

### Get output directories ###
output_dir, sim_dir = mh.create_output_dirs(params.sim_name, model_dir)

### Wrtie README containing simulation description
mh.save_config(sim_dir,params.sim_label,config_name)
# mh.create_sim_description(sim_dir, run_NML=run_NML, spk_type=spk_type, syns_lb=syns_lb, syns_ub=syns_ub, syns_type=syns_type, num_syns=num_syns, vinit=vinit, spk_freq=spk_freq, sim_message=sim_message)

### Generate network if running NeuroML ###
if params.run_NML:
    net_nml_file = mh.generate_network(model_dir, cell_name, pop_label, 
                                       force=True, 
                                       input_amp=params.input_amp, 
                                       start=params.stim_delay, 
                                       stop=params.stim_delay+params.stim_dur)

### Compile mechs ###
mh.compile_mechs(cwd,hocs_dir,mod_dir)
load_mechanisms(model_dir)

### Simulation configuration ###
cfg = specs.SimConfig()					                    # object of class SimConfig to store simulation configuration
cfg.duration = params.sim_dur 						                # Duration of the simulation, in ms
cfg.dt = params.dt								                # Internal integration timestep to use
cfg.verbose = True							                # Show detailed messages
cfg.recordTraces = {'V_soma':{'sec':'soma_0','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStep = params.recordStep
cfg.filename = os.path.join(sim_dir,cell_name+'_'+params.sim_label) 	# Set file output name
cfg.saveJson = False
cfg.analysis['plotTraces'] = {'include': [pop_label], 'saveFig': True}  # Plot recorded traces for this list of cells
cfg.hParams['celsius'] = 34.0 
cfg.hParams['v_init'] = params.vinit

### Import cell ###
if params.run_NML:
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
        importedCellParams[sec]['vinit'] = params.vinit

    ### Create population ###
    netParams.popParams[pop_label] = {'cellType': cell_type, 
                                      'cellModel': cell_model,
                                      'numCells': 1}


### Get sections ###
# basal, apical, basal_apical, basal_soma, apical_soma, basal_apical_soma, all
if params.syns_lb > 0:
    syn_secs = mh.get_secs_from_dist(hoc_file, cell_name, params.syns_lb, params.syns_ub)
    if params.add_soma:
        syn_secs.append('soma_0')
else:
    syn_secs = mh.get_components(importedCellParams, params.syns_type)

# Layer 1 input
syn_secs_L1 = mh.get_secs_from_dist(hoc_file, cell_name, 0.9, 1)

syn_secs_E, syn_secs_I = mh.get_rand_secs(syn_secs, params.num_syns)

cfg.recordTraces['V_syn'] = {'sec':secrets.choice(syn_secs),'loc':0.5,'var':'v'}

### Add AMPA/NMDA synapse ###
if 'HS' in params.syns_source:
    # Hay & Segev 2015
    netParams.synMechParams['AMPA'] = {'mod':'ProbAMPA2', 'tau_r_AMPA': 0.2, 'tau_d_AMPA': 1.7, 'e': 0, 'gmax': 0.0004}
    netParams.synMechParams['NMDA'] = {'mod': 'ProbNMDA2', 'tau_r_NMDA': 0.29, 'tau_d_NMDA': 43, 'e': 0, 'gmax': 0.0004}
    netParams.synMechParams['GABAA'] = {'mod': 'ProbUDFsyn2', 'tau_r': 1, 'tau_d': 20, 'e': -80, 'gmax': 0.001}
    inh_syns = ['GABAA']
else:
    # Dura-Bernal et al. 2024
    netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
    netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0} 
    netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
    netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93}
    inh_syns = ['GABAA', 'GABAB']

exc_syns = ['AMPA', 'NMDA']

### Add synaptic input ###
if params.enable_syns:
    syn_method = 'cell'  # 'stim' or 'cell'

    if 'cell' in syn_method:

        # Poisson spike pattern
        if 'poisson' in params.spk_type:
            ### Excitatory synapses ###
            netParams.popParams['vecstim'] = {
                'cellModel': 'VecStim',
                'numCells': params.num_syns,  # int(len(syn_secs)/4),
                'spikePattern': {'type': 'poisson',
                                'start': params.stim_delay,
                                'stop': params.stim_delay+params.stim_dur,
                                'frequency': params.spk_freq}
            }
            ## Layer 1 ##
            netParams.popParams['vecstimL1'] = {
                'cellModel': 'VecStim',
                'numCells': params.num_syns//2,  # int(len(syn_secs)/4),
                'spikePattern': {'type': 'poisson',
                                'start': params.stim_delay,
                                'stop': params.stim_delay+params.stim_dur,
                                'frequency': params.spk_freq}
            }
            ### Inhibitory synapses ###
            if params.add_inh:
                netParams.popParams['vecstimI'] = {
                    'cellModel': 'VecStim',
                    'numCells': params.num_syns//5,  # int(len(syn_secs)/4),
                    'spikePattern': {'type': 'poisson',
                                    'start': params.stim_delay,
                                    'stop': params.stim_delay+params.stim_dur,
                                    'frequency': params.spk_freq}
                }
        # Gaussian spike pattern
        else:
            netParams.popParams['vecstim'] = {
                'cellModel': 'VecStim',
                'numCells': params.num_syns,  # int(len(syn_secs)/4),
                'spikePattern': {'type': 'gauss',
                                'mu': params.gauss_mean,
                                'sigma': params.gauss_std}
            }

        ### Excitatory synapses ###
        netParams.connParams[f'vecstim->{pop_label}'] = {
            'preConds': {'pop': 'vecstim'},
            'postConds': {'pop': pop_label},
            'sec': syn_secs_E,
            'synsPerConn': params.synsPerConn,
            'synMech': exc_syns,
            'weight': params.syns_weight,  # 
            # 'synMechWeightFactor': [0.5,0.5],
            'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
            'probability': 1.0,
        }

        netParams.subConnParams[f'vecstim->{pop_label}'] = {
            'preConds': {'pop': 'vecstim'},
            'postConds': {'pop': pop_label},
            'sec': syn_secs_E,
            'groupSynMech': exc_syns,
            'density': 'uniform'
        }
        ## Layer 1 ##
        netParams.connParams[f'vecstimL1->{pop_label}'] = {
            'preConds': {'pop': 'vecstimL1'},
            'postConds': {'pop': pop_label},
            'sec': syn_secs_L1,
            'synsPerConn': params.synsPerConn,
            'synMech': exc_syns,
            'weight': params.syns_weight,  # 
            # 'synMechWeightFactor': [0.5,0.5],
            'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
            'probability': 1.0,
        }

        netParams.subConnParams[f'vecstimL1->{pop_label}'] = {
            'preConds': {'pop': 'vecstimL1'},
            'postConds': {'pop': pop_label},
            'sec': syn_secs_L1,
            'groupSynMech': exc_syns,
            'density': 'uniform'
        }

        ### Inhibitory synapses ###
        if params.add_inh:
            netParams.connParams[f'vecstimI->{pop_label}'] = {
                'preConds': {'pop': 'vecstimI'},
                'postConds': {'pop': pop_label},
                'sec': syn_secs_I,
                'synsPerConn': params.synsPerConn,
                'synMech': inh_syns,
                'weight': params.syns_weight,  # 
                # 'synMechWeightFactor': [0.5,0.5],
                'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
                'probability': 1.0,
            }

            netParams.subConnParams[f'vecstimI->{pop_label}'] = {
                'preConds': {'pop': 'vecstimI'},
                'postConds': {'pop': pop_label},
                'sec': syn_secs_I,
                'groupSynMech': inh_syns,
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
    'del': params.stim_delay,
    'dur': params.stim_dur,
    'amp': params.input_amps[params.amp_idx] 
}

netParams.stimTargetParams['Input_IC->Soma'] = {
    'source': 'Input_IC',
    'sec': 'soma_0',
    'loc': 0.5,
    'conds': {'pop': pop_label}
}

### Add linear probe ###
if params.record_LFP:
    probe_L = 1280//5
    channels = 1
    elec_dist = probe_L//params.depths  # microns
    disp = 100  # 150

    elec_pos = [[x*elec_dist, (y*elec_dist - disp)*-1, 0] for x in range(channels) for y in range(params.depths)]  # 
    # -x is left and -y is above soma
    # elec_pos.reverse()

    cfg.recordLFP = elec_pos
    cfg.analysis['plotLFP'] = {'saveFig': True}

### Run simulation ###
(pops, cells, conns, stims, simData) = sim.createSimulateAnalyze(netParams=netParams, simConfig=cfg, output=True)

### Process LFP ###
# lfp_bp_low, lfp_bp_spikes = mh.get_filtered_signal(simData['LFP'], cfg.dt)
# spkt = list(simData['spkt'])
mh.plot_lfp(simData, params.recordStep, params.sim_label, sim_dir)

### Plot morphology ###
if params.plot_morphology:
    sim.analysis.plotShape(showSyns=True, dist=0.8, includePre=[None], includePost=[pop_label], 
                           saveFig=True, axisLabels=True, returnPlotter=True, synColor='darkcyan') 


### Simulation time ###
if time_flag:
    end_t = time.time()
    final_t = end_t - start_t
    m, s = divmod(final_t, 60)
    time_path = os.path.join(sim_dir,f'{model_version}_{params.sim_label}_sim-time.txt')

    with open(time_path, "w") as text_file:
        text_file.write(f'Simulation time = {m} m {s} s')
    print(f'Simulation time = {m} m {s} s')
