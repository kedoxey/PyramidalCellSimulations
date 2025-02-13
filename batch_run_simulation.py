import os
import model_helpers as mh
from sklearn.model_selection import ParameterGrid
from itertools import chain
from simulate_cell import run_sim


config_name = 'default_config'

sim_name = 'INJ'

# sim_dur = 1000 if 'FR' in sim_name else 1000
# stim_delay = 0 if 'FR' in sim_name else 400
# stim_dur = sim_dur

stim_dur = 60
stim_delay = 1000
sim_dur = stim_dur + stim_delay + 40

param_sets = {#'active': {'log_firing_rate': True}
              'probes_mod': {'use_probes': [True],
                         'num_probes': [100],
                         'total_channels': [64]}}  # total_channels must be greater than 1
              # 'pas_soma': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'soma'},
              # 'pas_all': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'all'}}
              
group_num_syns = {'soma': [0]}
                  #'soma': [50*i for i in [1,4,7]], # [5*i for i in range(61, 111)]
                  # 'basal': [50]}
                  # 'apical_distal': [1800]}

paramGrids = []

for sim_flag, param_set in param_sets.items():
  for syns_type, num_syns in group_num_syns.items():

    paramGrid = {'sim_name': [sim_name],
                  'sim_flag': [sim_flag],
                  'nmldb_id': ['NMLCL000073'],
                  'enable_syns': [False],
                  'syns_type': [syns_type], #, 'basal', 'apical_distal'],  
                  'num_syns_E': num_syns,
                  'add_bkg': [False], 
                  'record_LFP': [True],
                  # 'depths': [4],
                  # 'apical_depths': [2],
                  'sim_dur': [sim_dur],  # 5000 or 1000
                  'stim_dur': [stim_dur], # 4900 or 400
                  'stim_delay': [stim_delay],  # 100 or 400
                  'input_amp': ['rheobase'],
                  'save_pickle': [True],
                  'dt': [0.05],
                  'recordStep': [0.05]}
    
    for param_name, param in param_set.items():
      paramGrid[param_name] = param

    paramGrids.append(paramGrid)
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  run_sim(config_name, batchParams)

  print(f"!!! Simulation ran for {batchParams['num_syns_E']} {batchParams['syns_type']} synpases !!!")
