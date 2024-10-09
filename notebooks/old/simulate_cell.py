import matplotlib.pyplot as plt
import os
import json
import model_helpers as mh
from neuron import h, load_mechanisms
from netpyne import specs, sim

test_name = 'soma_pas'

code_version = 'Hay'

model_version = 'NEURON'

if 'Hay' in code_version:
    nmldb_id = 'NMLCL000073'  # 'NMLCL000073' (Hay et al. 2011)
    model_name = f'{nmldb_id}-{model_version}'

    cell_model = 'Hay2011'  # 'Hay2011'
    cell_type = 'PYR'
    cell_name = 'L5PC'  # 'L5PC'
else:
    nmldb_id = 'NMLCL001535'  # 'NMLCL000073' (Hay et al. 2011)
    model_name = f'{nmldb_id}-{model_version}'

    cell_model = 'Allen'  # 'Hay2011'
    cell_type = 'PYR'
    cell_name = 'Cell_473863035'  # 'L5PC'

### Define paths ###
hoc_name = cell_name  # 'L5PC'

cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models')
model_dir = os.path.join(models_dir, model_version, model_name)  # 'L5bPCmodelsEH')
hocs_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir,'models')
mod_dir = model_dir if 'biophys' not in model_name else os.path.join(model_dir, 'mod')

### Download model ###
mh.download_from_nmldb(nmldb_id, model_version)

output_dir = os.path.join(model_dir,'output')
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

### Simulation configuration ###
if 'Hay' in code_version:
    input_amps = [0.35477, 0.44346, 0.53215, 1.0643]
else:
    input_amps = [0.10574, 0.13218, 0.15862, 0.31723]

amp_idx = -1  # if -1 then no current injection
input_amp = input_amps[amp_idx]
test_label = f'input_{amp_idx}-{test_name}'

vinit = -80

## cfg
cfg = specs.SimConfig()					            # object of class SimConfig to store simulation configuration
cfg.duration = 3000 						            # Duration of the simulation, in ms
cfg.dt = 0.01								                # Internal integration timestep to use
cfg.verbose = True							                # Show detailed messages
cfg.recordTraces = {'V_soma':{'sec':'soma_0','loc':0.5,'var':'v'}}  # Dict with traces to record
cfg.recordStim = True
cfg.recordStep = 0.01
cfg.filename = os.path.join(output_dir,cell_name+'_'+test_label) 			# Set file output name
cfg.saveJson = False
cfg.analysis['plotTraces'] = {'include': ['all'], 'saveFig': True} # Plot recorded traces for this list of cells
cfg.hParams['celsius'] = 34.0 
cfg.hParams['v_init'] = vinit

netParams = specs.NetParams()

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
        'conds': {'cellType': cell_name}
    }

    print(f'current clamp added at {input_amp} pA')
else:
    print('no input')


if 'NeuroML' in model_version:
    mh.compile_mechs(cwd,hocs_dir,mod_dir)
    load_mechanisms(model_dir)
    h.load_file('stdrun.hoc')

    nml2_file_name = os.path.join(model_dir, cell_name+'.cell.nml')

    nml_net_fname = os.path.join(model_dir, cell_name) + '_Network.net.nml'

    sim.importNeuroML2SimulateAnalyze(nml_net_fname, cfg)

elif 'NEURON' in model_version:

    ### Compile mechs ###
    mh.compile_mechs(cwd,hocs_dir,mod_dir)
    load_mechanisms(model_dir)

    ### Import cell into NetPyNE ###
    hoc_file = os.path.join(hocs_dir, f'{cell_name}_{test_name}.hoc')

    ### Network parameters ###
    netParams = specs.NetParams()

    cell_label = cell_name+'_hoc'

    importedCellParams = netParams.importCellParams(label=cell_label,
                                                    conds={'cellType': cell_type, 'cellModel': cell_model},
                                                    fileName=hoc_file,
                                                    cellName=cell_name,
                                                    importSynMechs=False,
                                                    somaAtOrigin=True
                                                )

    for sec in importedCellParams['secs']:
        importedCellParams[sec]['vinit'] = vinit

    ### Create population ###
    pop_label = cell_label+'_pop'
    netParams.popParams[pop_label] = {'cellType': cell_type, 
                                      'cellModel': cell_model,
                                      'numCells': 1}

    ### Add linear probe ###

    ####################
    # Geometry
    # --------
    # (Cartesian axes)
    #          y    z
    #          ^  ^
    #          | /
    #   x <--- o --
    #         /|
    ####################

    # probe_L = 1280
    # channels = 1
    # depths = 10
    # elec_dist = probe_L//depths  # microns

    # elec_pos = [[x*elec_dist, y*elec_dist - 150, 0] for x in range(channels) for y in range(depths)]

    # cfg.recordLFP = elec_pos
    # cfg.analysis['plotLFP'] = {'saveFig': True}

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

    ### Run simulation ###
    sim.createSimulateAnalyze(netParams = netParams, simConfig = cfg)

    ### Plot morphology ###
    sim.analysis.plotShape(showSyns=True, dist=0.8, saveFig=True, axisLabels=True)
