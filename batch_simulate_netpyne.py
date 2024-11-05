import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid
from itertools import chain


config_name = 'default_config'

sim_name = 'FR'

param_sets = {'active': {'log_firing_rate': True}}
              # 'pas_soma': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'soma'},
              # 'pas_all': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'all'}}
              
group_num_syns = {#'soma': [5*i for i in range(91, 111)],
                  # 'basal': [5*i for i in range(83,101)]}
                  'apical_distal': [1000*i for i in range(6,11)]}

paramGrids = []

for sim_flag, param_set in param_sets.items():
  for syns_type, num_syns in group_num_syns.items():

    paramGrid = {'sim_name': [sim_name],
                  'sim_flag': [sim_flag],
                  'nmldb_id': ['NMLCL000073'],
                  'enable_syns': [True],
                  'syns_type': [syns_type], #, 'basal', 'apical_distal'],  
                  'num_syns_E': num_syns,
                  'add_bkg': [False], 
                  'record_LFP': [False],
                  'depths': [4],
                  'apical_depths': [2],
                  'sim_dur': [5000],  # 5000 or 1000
                  'stim_dur': [4900], # 4900 or 400
                  'stim_delay': [1000],  # 1000 or 400
                  'save_pickle': [False]}
    
    for param_name, param in param_set.items():
      paramGrid[param_name] = [param]

    paramGrids.append(paramGrid)
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  simulate_netpyne.run_sim(config_name, batchParams)

  print(f"!!! Simulation ran for {batchParams['num_syns_E']} {batchParams['syns_type']} synpases !!!")
