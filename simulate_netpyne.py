import matplotlib.pyplot as plt
import os
import secrets
import json
import pickle
import model_helpers as mh
import numpy as np
from neuron import h, load_mechanisms
from netpyne import specs, sim, conversion
import argparse as ap
import time


def run_sim(config_name, *batch_params):

    ### Import simulation config ###
    # config_name = args.config_name
    params = ap.Namespace(**mh.load_config(config_name))

    for batch_param, batch_value in batch_params[0].items():
        setattr(params, batch_param, batch_value)

    ### Set simulation name and label ###
    params.sim_name = f'{params.sim_name}-{params.syns_type}'
    if isinstance(params.syns_weight, list):
        params.sim_label = f'{params.sim_name}_{params.num_syns_E}Ex{params.syns_weight[0]}AMPAx{params.syns_weight[1]}NMDA_{params.num_poisson}x{params.spk_freq}Hz'
    else:
        params.sim_label = f'{params.sim_name}_{params.num_syns_E}Ex{params.syns_weight}_{params.num_poisson}x{params.spk_freq}Hz'
    
    if params.add_bkg:
        params.sim_label += '+bkg'

    params.sim_label += f'{params.sim_flag}'

    ### Model information ###
    model_version = 'NeuroML' if params.run_NML else 'NEURON'

    nmldb_id =  params.nmldb_id  # 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
    model_name = f'{nmldb_id}-{model_version}'

    ### Download model ###
    cell_model = mh.download_from_nmldb(nmldb_id, model_version)  # AuthorYear

    ### Define paths ###
    cwd = os.getcwd()
    models_dir = os.path.join(cwd, 'models') 
    model_dir = os.path.join(models_dir, model_version, model_name)  # 'L5bPCmodelsEH')
    hocs_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir,'models')
    mod_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir, 'mod')

    cell_type = 'PYR'
    cell_name = mh.get_cell_name(model_dir)  # 'L5PC'
    cell_label = cell_name+'_hoc'
    pop_label = cell_name+'_Pop'

    hoc_file = os.path.join(hocs_dir, f'{cell_name}.hoc')

    ### Copy synapses ###
    # mh.copy_synapses(model_dir)

    ### Get output directories ###
    output_dir, sim_dir = mh.create_output_dirs(params.sim_name, params.sim_label, model_dir)

    ### Wrtie README containing simulation description
    mh.write_config(params,sim_dir,params.sim_label,config_name)
    # mh.create_sim_description(sim_dir, run_NML=run_NML, spk_type=spk_type, syns_lb=syns_lb, syns_ub=syns_ub, syns_type=syns_type, num_syns=num_syns, vinit=vinit, spk_freq=spk_freq, sim_message=sim_message)

    ### Generate network if running NeuroML ###
    if params.run_NML:
        net_nml_file = mh.generate_network(model_dir, cell_name, pop_label, 
                                        force=True, 
                                        input_amp=params.input_amp, 
                                        start=params.stim_delay, 
                                        stop=params.stim_delay+params.stim_dur)

    ### Compile mechs ###
    mh.compile_mechs(cwd,hocs_dir,mod_dir)  #,force=True)
    load_mechanisms(model_dir)

    ### Simulation configuration ###
    cfg = specs.SimConfig()					                    # object of class SimConfig to store simulation configuration
    cfg.duration = params.sim_dur 						                # Duration of the simulation, in ms
    cfg.dt = params.dt								                # Internal integration timestep to use
    cfg.verbose = True							                # Show detailed messages
    cfg.recordTraces = {'V_soma': {'sec': 'soma_0', 'loc': 0.5, 'var': 'v'}}  # Dict with traces to record
    cfg.recordStep = params.recordStep
    # cfg.recordStim = True
    cfg.filename = os.path.join(sim_dir,cell_name+'_'+params.sim_label) 	# Set file output name
    cfg.savePickle = False
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
        
        channel_secs = mh.get_compartments(hoc_file, importedCellParams, cell_name, params.channel_secs)
        importedCellParams = mh.toggle_channels(importedCellParams, channel_secs, params.channel_toggles)  #,'Na',params.soma_na_toggle)

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
        syn_secs = mh.get_compartments(hoc_file, importedCellParams, cell_name, params.syns_type)

    # Layer inhibitory sections
    layer_bounds = {'L1': {'lb': 5/6, 'ub': 1},
                    'L2': {'lb': 2/3, 'ub': 5/6},
                    'L4': {'lb': 1/6, 'ub': 5/12}}
    layer_secs = {'L1': [],
                'L2': [],
                'L4': []}

    for layer, bounds in layer_bounds.items():
        layer_secs[layer] = mh.get_secs_from_dist(hoc_file, cell_name, bounds['lb'], bounds['ub'], secs_lim='apic')

    layer_secs['L5'] = ['soma_0']

    # syn_secs_L1 = mh.get_secs_from_dist(hoc_file, cell_name, 0.9, 1)

    syn_secs_E, syn_secs_I = mh.get_rand_secs(syn_secs, params.num_syns_E, params.num_syns_I, params.seed)

    ### Add AMPA/NMDA synapse ###
    if 'HS' in params.syns_source:
        # Hay & Segev 2015
        netParams.synMechParams['AMPA_NMDA'] = {'mod':'ProbAMPANMDA2', 'tau_r_AMPA': 0.3, 'tau_d_AMPA': 3, 'tau_r_NMDA': 2, 'tau_d_NMDA': 65, 'e': 0, 'gmax': 0.0004}
        exc_syns =  ['AMPA_NMDA']
        exc_syn_locs = [0.5]
        # netParams.synMechParams['AMPA'] = {'mod':'ProbAMPA2', 'tau_r_AMPA': 0.3, 'tau_d_AMPA': 3, 'e': 0, 'gmax': 0.0004}
        # netParams.synMechParams['NMDA'] = {'mod': 'ProbNMDA2', 'tau_r_NMDA': 2, 'tau_d_NMDA': 65, 'e': 0, 'gmax': 0.0004}
        netParams.synMechParams['GABAA'] = {'mod': 'ProbUDFsyn2', 'tau_r': 1, 'tau_d': 20, 'e': -80, 'gmax': 0.001}
        inh_syns = ['GABAA']
    else:
        # Dura-Bernal et al. 2024
        netParams.synMechParams['AMPA'] = {'mod':'MyExp2SynBB', 'tau1': 0.05, 'tau2': 5.3, 'e': 0}
        netParams.synMechParams['NMDA'] = {'mod': 'MyExp2SynNMDABB', 'tau1NMDA': 15, 'tau2NMDA': 150, 'e': 0} 
        exc_syns = ['AMPA', 'NMDA']
        exc_syn_locs = [0.5, 0.5]
        netParams.synMechParams['GABAA'] = {'mod':'MyExp2SynBB', 'tau1': 0.07, 'tau2': 18.2, 'e': -80}
        netParams.synMechParams['GABAB'] = {'mod':'MyExp2SynBB', 'tau1': 3.5, 'tau2': 260.9, 'e': -93}
        inh_syns = ['GABAA', 'GABAB']


    for syn_sec_E in syn_secs_E:
        cfg.recordTraces[f'V_{syn_sec_E}'] = {'sec':syn_sec_E,'loc':0.5,'var':'v'}
        cfg.recordTraces[f'I_{syn_sec_E}_ampa'] = {'sec':syn_sec_E,'loc':exc_syn_locs[0],'synMech':'AMPA_NMDA','var':'g_AMPA'}
        cfg.recordTraces[f'I_{syn_sec_E}_nmda'] = {'sec':syn_sec_E,'loc':exc_syn_locs[0],'synMech':'AMPA_NMDA','var':'g_NMDA'}

    if 'apic_32' not in syn_sec_E:
        cfg.recordTraces[f'V_apic_32'] = {'sec':'apic_32','loc':0.5,'var':'v'}
        cfg.recordTraces[f'I_apic_32_ampa'] = {'sec':'apic_32','loc':exc_syn_locs[0],'synMech':'AMPA_NMDA','var':'g_AMPA'}
        cfg.recordTraces[f'I_apic_32_nmda'] = {'sec':'apic_32','loc':exc_syn_locs[0],'synMech':'AMPA_NMDA','var':'g_NMDA'}


    ### Add synaptic input ###
    if params.num_syns_E == 0:
        params.enable_syns = False

    if params.enable_syns:
        
        # Poisson spike pattern

        ### Layer inhibitory input
        num_I_each = params.num_syns_I // len(layer_secs.keys())

        if params.num_syns_I > 0:
            for layer, layer_secs in layer_secs.items():
                netParams.popParams[f'vecstim_I{layer}'] = {
                        'cellModel': 'VecStim',
                        'numCells': num_I_each,  # int(len(syn_secs)/4),
                        'spikePattern': {'type': 'poisson',
                                        'start': params.stim_delay,
                                        'stop': params.stim_delay+params.stim_dur,
                                        'frequency': params.spk_freq }  # np.random.randint(params.spk_freq_lb, params.spk_freq_ub, 1)[0]}
                    }
                
                netParams.connParams[f'vecstim_I{layer}->{pop_label}'] = {
                        'preConds': {'pop': f'vecstim_I{layer}'},
                        'postConds': {'pop': pop_label},
                        'sec': layer_secs,
                        'synsPerConn': params.synsPerConn,
                        'synMech': inh_syns,
                        'weight': params.syns_weight,  # 
                        # 'synMechWeightFactor': [0.5,0.5],
                        'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
                        'probability': 1.0,
                    }

                netParams.subConnParams[f'vecstimI{layer}->{pop_label}'] = {
                    'preConds': {'pop': f'vecstimI{layer}'},
                    'postConds': {'pop': pop_label},
                    'sec': layer_secs,
                    'groupSynMech': inh_syns,
                    'density': 'uniform'
                }

        ### Excitatory synapses ###
        num_E_each = params.num_syns_E // params.num_poisson

        for i_poisson in range(params.num_poisson):

            netParams.popParams[f'vecstim_E{i_poisson}'] = {
                'cellModel': 'VecStim',
                'numCells': num_E_each,  # int(len(syn_secs)/4),
                'spikePattern': {'type': 'poisson',
                                'start': params.stim_delay,
                                'stop': params.stim_delay+params.stim_dur,
                                'frequency': params.spk_freq}  #  np.random.randint(params.spk_freq_lb, params.spk_freq_ub, 1)[0]}
            }

            netParams.connParams[f'vecstim_E{i_poisson}->{pop_label}'] = {
                'preConds': {'pop': f'vecstim_E{i_poisson}'},
                'postConds': {'pop': pop_label},
                'sec': syn_secs_E,
                'synsPerConn': params.synsPerConn,
                'synMech': exc_syns,
                'weight': params.syns_weight,  # 
                # 'synMechWeightFactor': [0.5,0.5],
                'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
                'probability': 1.0,
                'loc': exc_syn_locs
            }


    ### Add input ###
    netParams.stimSourceParams['Input_IC'] = {
        'type': 'IClamp',
        'del': params.stim_delay,
        'dur': params.stim_dur,
        'amp': params.input_amp 
    }

    netParams.stimTargetParams['Input_IC->Soma'] = {
        'source': 'Input_IC',
        'sec': 'soma_0',
        'loc': 0.5,
        'conds': {'pop': pop_label}
    }

     ### Background input ###
    if params.add_bkg:
        netParams.stimSourceParams['bkg'] = {'type': 'NetStim', 'rate': 100, 'noise': 1}
        # netParams.stimTargetParams['bkg->ALL'] = {'source': 'bkg', 'conds': {'cellType': [cell_label]}, 
        #                                           'weight': 0.01, 'delay': 'max(1, normal(5,2))', 'synMech': 'AMPA_NMDA'}
        netParams.stimTargetParams['bkg->ALL'] = {'source': 'bkg', 'sec': 'soma_0', 'loc': 0.5, 
                                                'conds': {'pop': pop_label}, 'weight': 8, 
                                                'delay': 'max(1, normal(5,2))', 'synMech': 'AMPA_NMDA'}


    ### Add linear probe ###
    if params.record_LFP:
        probe_L = 300
        channels = 1
        elec_dist = probe_L//params.depths  # microns
        disp = 130  # 150

        elec_pos = [[x*elec_dist, (y*elec_dist - disp)*-1, 0] for x in range(channels) for y in range(params.depths)]

        apic_pos = [[0, -930-(y*elec_dist - disp), 0] for y in range(2)]

        elec_pos.extend(apic_pos)  # 
        # -x is left and -y is above soma
        # elec_pos.reverse()

        cfg.recordLFP = elec_pos
        cfg.analysis['plotLFP'] = {'saveFig': True}


    ### Run simulation ###
    (pops, cells, conns, stims, simData) = sim.createSimulateAnalyze(netParams=netParams, simConfig=cfg, output=True)

    # mh.save_simData(simData, params.sim_label, sim_dir)

    synColors = ('firebrick','darkcyan')
    colormapE, colormapI = mh.get_colormaps(params.num_syns_E, params.num_syns_I)
    secSynColors = mh.get_syn_sec_colors(cells[0], (colormapE, colormapI), synColors)

    if params.enable_syns:
        spikeTrains = mh.plot_pre_spike_trains(cells, conns, params.sim_label, sim_dir)

        if len(syn_secs_E) < 175:
            mh.plot_secs(simData, spikeTrains, params.sim_label, sim_dir, secSynColors)

    mh.plot_soma(simData, params.sim_label, sim_dir)

    if params.record_LFP:
        mh.plot_isolated_LFP(simData, params.syns_type, params.sim_label, sim_dir)

    ### Plot morphology ###
    if params.plot_morphology:
        sim.analysis.plotShape(showSyns=True, dist=0.8, includePre=[None], includePost=[pop_label], 
                            saveFig=True, axisLabels=True, returnPlotter=True, secSynColors=secSynColors, colormaps=(colormapE,colormapI), synColors=synColors)




# run_sim('exc_config', ('syns_weight',0.5))
