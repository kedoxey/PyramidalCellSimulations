import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import secrets
import json
import pickle
import model_helpers as mh
import numpy as np
from neuron import h, load_mechanisms
from netpyne import specs, sim
import argparse as ap
import time

import utils


def run_sim(config_name, *batch_params):

    ### Import simulation config ###
    # config_name = args.config_name
    params = ap.Namespace(**utils.config.load_config(config_name))

    for batch_param, batch_value in batch_params[0].items():
        setattr(params, batch_param, batch_value)

    ### Get model parameters from neuroml-db ###
    neuron_details = utils.wrangle.get_model_details(params.nmldb_id)
    params.dt = neuron_details['model']['Optimal_DT']
    params.recordStep = 3*params.dt
    rheobase = neuron_details['model']['Rheobase_High']
    bias_curr = neuron_details['model']['Bias_Current']
    params.vinit = neuron_details['model']['Resting_Voltage']

    if not isinstance(params.input_amp, int):
        if 'rheobase' in params.input_amp:
            params.input_amp = 3*rheobase

    ### Set simulation name and label ###
    params.sim_name = f'{params.sim_name}_{params.syns_type}'
    params.sim_label = f'{params.sim_name}'

    # if params.input_amp != 0:
    params.sim_label += f'-{round(params.input_amp,3)}nA'

    if params.enable_syns:
        if isinstance(params.syns_weight, list):
            params.sim_label += f'-{params.num_syns_E}Ex{params.syns_weight[0]}AMPAx{params.syns_weight[1]}NMDA-{params.num_poisson}x{params.spk_freq}Hz'
        else:
            params.sim_label += f'-{params.num_syns_E}Ex{params.syns_weight}-{params.num_poisson}x{params.spk_freq}Hz'
    
    if params.use_probes:
        params.sim_label += f'-{params.num_probes}x{params.total_channels}elec'

    if params.add_bkg:
        params.sim_label += '+bkg'

    params.sim_label += f'-{params.sim_flag}'

    ### Model information ###
    if params.local:
        model_version = 'local'

        model_name = params.model_name

    else:
        model_version = 'NeuroML' if params.run_NML else 'NEURON'

        nmldb_id =  params.nmldb_id  # 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
        model_name = f'{nmldb_id}-{model_version}'


    ### Define paths ###
    cwd = os.getcwd()
    models_dir = os.path.join(cwd, 'models')
    model_dir = os.path.join(models_dir, model_version, model_name)

    hocs_dir = os.path.join(model_dir, params.hocs_dname) if params.hocs_dname else model_dir
    mod_dir = os.path.join(model_dir, params.mod_dname) if params.mod_dname else model_dir


    ### Download or define model ###
    cell_model = params.cell_model if params.cell_model else utils.wrangle.download_from_nmldb(nmldb_id, model_version)
    cell_name = params.cell_name if params.cell_name else utils.wrangle.get_cell_name(model_dir)
    if params.local: params.sim_label += f'-{cell_name}'

    cell_type = params.cell_type
    cell_label = cell_name+'_hoc'
    pop_label = cell_name+'_Pop'

    hoc_file = os.path.join(hocs_dir, f'{cell_name}.hoc')

    ### Copy synapses ###
    # config.copy_synapses(model_dir)

    ### Get output directories ###
    output_dir, sim_dir = utils.config.create_output_dirs(params.sim_name, params.sim_label, model_dir)

    ### Wrtie README containing simulation description
    utils.config.write_config(params,sim_dir,params.sim_label,config_name)
    # config.create_sim_description(sim_dir, run_NML=run_NML, spk_type=spk_type, syns_lb=syns_lb, syns_ub=syns_ub, syns_type=syns_type, num_syns=num_syns, vinit=vinit, spk_freq=spk_freq, sim_message=sim_message)

    ### Generate network if running NeuroML ###
    if params.run_NML:
        net_nml_file = utils.wrangle.generate_network(model_dir, cell_name, pop_label, 
                                                force=True, 
                                                input_amp=params.input_amp, 
                                                start=params.stim_delay, 
                                                stop=params.stim_delay+params.stim_dur)

    ### Compile mechs ###
    utils.wrangle.compile_mechs(cwd,hocs_dir,mod_dir)  #,force=True)
    load_mechanisms(mod_dir)

    ### Instantiate simulation configuration ###
    cfg = specs.SimConfig()					                    # object of class SimConfig to store simulation configuration

    ### Import cell ###
    netParams = specs.NetParams()

    importedCellParams = netParams.importCellParams(label=cell_label,
                                                    conds={'cellType': cell_type, 'cellModel': cell_model},
                                                    fileName=hoc_file,
                                                    cellName=cell_name
                                                    )

    # netParams.defaultThreshold = -20
    for sec in importedCellParams['secs']:
        importedCellParams[sec]['vinit'] = params.vinit
    
    channel_secs = utils.setup.get_compartments(hoc_file, importedCellParams, cell_name, params.channel_secs)
    importedCellParams = utils.setup.toggle_channels(importedCellParams, channel_secs, params.channel_toggles)  #,'Na',params.soma_na_toggle)

    if params.local:
        importedCellParams = utils.setup.update_cell_params(importedCellParams, cell_name, os.path.join(hocs_dir, f'{params.cell_type}_model_params.pkl'))

    ### Define geometry
    netpyne_geometry = True
    if netpyne_geometry:
        netParams.propVelocity = 100.0
        netParams.probLengthConst = 150.0

        neuron_morpho = utils.wrangle.get_morphophetrics(params.nmldb_id)
        for metric_dict in neuron_morpho:
            if metric_dict['Metric_ID'] == 'Height':
                model_height = metric_dict['Maximum']
            if metric_dict['Metric_ID'] == 'Width':
                model_width = metric_dict['Maximum']
            if metric_dict['Metric_ID'] == 'Depth':
                model_depth = metric_dict['Maximum']

        channel_spacing = 50.
        buffer_dim = 50.

        x_dim =  channel_spacing*np.floor((model_width+buffer_dim)/channel_spacing)+channel_spacing
        y_dim =  channel_spacing*np.floor((model_height+buffer_dim)/channel_spacing)+channel_spacing
        z_dim =  channel_spacing*np.floor((model_depth+buffer_dim)/channel_spacing)+channel_spacing

        if x_dim>z_dim:
            z_dim = x_dim
        else:
            x_dim = z_dim

        netParams.sizeX = x_dim # x-dimension (horizontal length) size in um
        netParams.sizeY = y_dim # y-dimension (vertical height or cortical depth) size in um
        netParams.sizeZ = z_dim # z-dimension (horizontal length) size in um

        cfg.pt3dRelativeToCellLocation = False  # Make grid as function from soma location (at origin) (default: True)
        cfg.invertedYCoord = False  # make y-axis coordinate negative so they represent depth when visualized (0 at the top) (default: True)

    ### Create population ###
    netParams.popParams[pop_label] = {'cellType': cell_type, 
                                    'cellModel': cell_model,
                                    'numCells': 1}

    for sec_key in importedCellParams['secs']:
        if 'soma' in sec_key:
            soma_name = sec_key

    ### Add synaptic input ###
    if params.num_syns_E == 0:
        params.enable_syns = False

    if params.enable_syns:

            ### Get sections ###
        # basal, apical, basal_apical, basal_soma, apical_soma, basal_apical_soma, all
        if params.syns_lb > 0:
            syn_secs = utils.setup.get_secs_from_dist(hoc_file, cell_name, soma_name, params.syns_lb, params.syns_ub)
            if params.add_soma:
                syn_secs.append(soma_name)
        else:
            syn_secs = utils.setup.get_compartments(hoc_file, importedCellParams, cell_name, soma_name, params.syns_type)

        # syn_secs_L1 = mh.get_secs_from_dist(hoc_file, cell_name, 0.9, 1)

        syn_secs_E, syn_secs_I = utils.setup.get_rand_secs(syn_secs, params.num_syns_E, params.num_syns_I, params.seed)

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

        if 'apic_32' not in syn_secs_E:
            cfg.recordTraces[f'V_apic_32'] = {'sec':'apic_32','loc':0.5,'var':'v'}
            cfg.recordTraces[f'I_apic_32_ampa'] = {'sec':'apic_32','loc':exc_syn_locs[0],'synMech':'AMPA_NMDA','var':'g_AMPA'}
            cfg.recordTraces[f'I_apic_32_nmda'] = {'sec':'apic_32','loc':exc_syn_locs[0],'synMech':'AMPA_NMDA','var':'g_NMDA'}
        
        # Poisson spike pattern
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



            netParams.connParams[f'vecstim_E{i_poisson}!->{pop_label}'] = {
                'preConds': {'pop': f'vecstim_E{i_poisson}'},
                'postConds': {'pop': pop_label},
                'sec': syn_secs_E,
                'loc': exc_syn_locs*num_E_each,
                'synsPerConn': params.synsPerConn,
                'synMech': exc_syns,
                'weight': params.syns_weight,  # 
                # 'synMechWeightFactor': [0.5,0.5],
                'delay': 5,  # 'defaultDelay + dist_2D/propVelocity',
                'probability': 1.0
            }

        ### Layer inhibitory input
        if params.num_syns_I > 0:
            # Layer inhibitory sections
            layer_bounds = {'L1': {'lb': 5/6, 'ub': 1},
                            'L2': {'lb': 2/3, 'ub': 5/6},
                            'L4': {'lb': 1/6, 'ub': 5/12}}
            layer_secs = {'L1': [],
                        'L2': [],
                        'L4': []}

            for layer, bounds in layer_bounds.items():
                layer_secs[layer] = utils.setup.get_secs_from_dist(hoc_file, cell_name, soma_name, bounds['lb'], bounds['ub'], secs_lim='apic')

            layer_secs['L5'] = [soma_name]

            num_I_each = params.num_syns_I // len(layer_secs.keys())

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


    ### Add input ###
    # if params.input_amp != 0:

    if 'soma' in params.input_sec:
        params.input_sec = soma_name

    netParams.stimSourceParams['SSInput'] = {
        'type': 'IClamp',
        'del': 0,
        'dur': params.stim_delay + params.stim_dur,
        'amp': bias_curr  # bias current
    }
    netParams.stimTargetParams[f'SSInput->{params.input_sec}'] = {
        'source': 'SSInput',
        'sec': params.input_sec,
        'loc': 0.5,
        'conds': {'pop': pop_label}
    }

    netParams.stimSourceParams['Input_IC'] = {
        'type': 'IClamp',
        'del': params.stim_delay,
        'dur': params.stim_dur,
        'amp': params.input_amp 
    }
    netParams.stimTargetParams[f'Input_IC->{params.input_sec}'] = {
        'source': 'Input_IC',
        'sec': params.input_sec,
        'loc': 0.5,
        'conds': {'pop': pop_label}
    }


     ### Background input ###
    if params.add_bkg:
        netParams.stimSourceParams['bkg'] = {'type': 'NetStim', 'rate': 100, 'noise': 1}
        # netParams.stimTargetParams['bkg->ALL'] = {'source': 'bkg', 'conds': {'cellType': [cell_label]}, 
        #                                           'weight': 0.01, 'delay': 'max(1, normal(5,2))', 'synMech': 'AMPA_NMDA'}
        netParams.stimTargetParams['bkg->ALL'] = {'source': 'bkg', 'sec': soma_name, 'loc': 0.5, 
                                                'conds': {'pop': pop_label}, 'weight': 8, 
                                                'delay': 'max(1, normal(5,2))', 'synMech': 'AMPA_NMDA'}

    rec_electrode = None
    ### Add linear probe ###
    if params.record_LFP:

        if params.use_probes:

            if params.detect_limit:
                rec_electrode, r_thetas = utils.setup.define_detect_lim_geom(params.sim_label, sim_dir)
            else:
                r_thetas = None
                sampling_dist = 'adjusted_uniform'
                rec_electrode = utils.setup.define_electrode_geom(params.num_probes, params.total_channels, sampling_dist, params.sim_label, sim_dir, model_dir)

            cfg.recordLFP = rec_electrode
            cfg.analysis['plotLFP'] = {'plots': ['locations'], 'electrodes': ['all'], 'saveFig': True, 'showFig': False}
        else:
            probe_L = 300
            channels = 1
            elec_dist = probe_L//params.depths  # microns
            disp = 130  # 150

            elec_pos = [[x*elec_dist, (y*elec_dist - disp)*-1, 0] for x in range(channels) for y in range(params.depths)]

            if params.apical_depths > 0:
                apic_pos = [[0, -930-(y*elec_dist - disp), 0] for y in range(params.apical_depths)]
                elec_pos.extend(apic_pos)  # 
            # -x is left and -y is above soma
            elec_pos.reverse()
            cfg.recordLFP = elec_pos

        
    ### Simulation configuration ###
    cfg.duration = params.sim_dur 						                # Duration of the simulation, in ms
    cfg.dt = params.dt	
    cfg.recordStep = params.recordStep							                # Internal integration timestep to use
    cfg.verbose = True		
    cfg.recordCells = ['all']					                # Show detailed messages
    cfg.recordTraces[f'V_{soma_name}'] = {'sec': soma_name, 'loc': 0.5, 'var': 'v'}  # Dict with traces to record
    # cfg.recordStim = True
    cfg.filename = os.path.join(sim_dir,cell_name+'_'+params.sim_label) 	# Set file output name
    cfg.savePickle = params.save_pickle
    # cfg.analysis['plotTraces'] = {'include': [pop_label], 'saveFig': False}  # Plot recorded traces for this list of cells
    cfg.hParams['celsius'] = 34.0 
    cfg.hParams['v_init'] = params.vinit

    ### Run simulation ###
    (pops, cells, conns, stims, simData) = sim.createSimulateAnalyze(netParams=netParams, simConfig=cfg, output=True)

    # mh.save_simData(simData, params.sim_label, sim_dir)

    ### Save LFP data ###
    if params.record_LFP:

        waveforms_df = utils.process.reformat_data(simData, rec_electrode, soma_name, params.stim_delay, params.sim_dur, params.nmldb_id, params.sim_label, sim_dir, r_thetas)

        if params.detect_limit:

            utils.process.save_detect_limit(waveforms_df, model_dir)

    ### Plot sections ###
    synColors = {'E': 'firebrick', 'I': 'darkcyan'}
    colormapE, colormapI = utils.plotting.get_colormaps(params.num_syns_E, params.num_syns_I)
    secSynColors = utils.plotting.get_syn_sec_colors(cells[0], params.use_colormaps, (colormapE, colormapI), synColors)

    if params.enable_syns:
        spikeTrains = utils.plotting.plot_pre_spike_trains(cells, conns, params.sim_label, sim_dir)

        if len(syn_secs_E) < 175:
            utils.plotting.plot_secs(simData, soma_name, spikeTrains, params.sim_label, sim_dir, secSynColors)

        utils.plotting.plot_syns_traces(simData, syn_secs_E, params.sim_label, sim_dir, synColors)

    ### Plot somatic spiking ###
    utils.plotting.plot_soma(simData, soma_name, params.sim_label, sim_dir)
        
    ### Plot LFP ###
    if params.record_LFP:
        if not params.use_probes:
            utils.plotting.plot_isolated_LFP(simData, soma_name, params.syns_type, params.num_syns_E, params.sim_label, sim_dir, output_dir)
            utils.plotting.plot_isolated_syn_traces(simData, soma_name, syn_secs, params.syns_type, params.num_syns_E, params.sim_label, sim_dir, output_dir, synColors)
            utils.plotting.plot_isolated_traces(simData, soma_name, syn_secs, params.syns_type, params.num_syns_E, params.sim_label, sim_dir, output_dir, synColors)
            utils.plotting.plot_isolated_soma_pot(simData, soma_name, params.syns_type, params.num_syns_E, params.sim_label, sim_dir, output_dir)
        else:
            utils.plotting.plot_eap_kernel('closest', sim_dir, params.sim_label)
            utils.plotting.plot_eap_kernel('farthest', sim_dir, params.sim_label)
            
        sim.analysis.plotLFP(plots=['locations'], saveFig=True, showFig=False)

    if params.log_firing_rate:
        utils.process.save_firing_rate(simData, soma_name, params.sim_dur, params.syns_type, params.num_syns_E, output_dir)

    ### Plot morphology ###
    if params.plot_morphology:
        sim.analysis.plotShape(showSyns=True, dist=0.8, includePre=[None], includePost=[pop_label], axisLabels=False, includeGrid=False,
                               saveFig=True, fontSize=10, returnPlotter=True, bkgColor=mpl.colors.to_rgba('w'), 
                               secSynColors=secSynColors, colormaps=(colormapE,colormapI), synColors=synColors)
        
