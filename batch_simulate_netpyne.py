import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid
from itertools import chain


config_name = 'default_config'

sim_name = 'FR'

sim_dur = 5000 if 'FR' in sim_name else 1000
stim_delay = 100 if 'FR' in sim_name else 400
stim_dur = sim_dur - stim_delay

param_sets = {'active': {'log_firing_rate': True}}
              # 'pas_soma': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'soma'},
              # 'pas_all': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'all'}}
              
group_num_syns = {'basal': [50*i for i in range(8,16)]} # [5*i for i in range(61, 111)]
                  #  'basal': [200, 500]}
                  # 'apical_distal': [1000*i for i in range(6,11)]} [200*i for i in range(2, 50)]

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
                  'sim_dur': [sim_dur],  # 5000 or 1000
                  'stim_dur': [stim_dur], # 4900 or 400
                  'stim_delay': [stim_delay],  # 100 or 400
                  'save_pickle': [False]}
    
    for param_name, param in param_set.items():
      paramGrid[param_name] = [param]

    paramGrids.append(paramGrid)
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  simulate_netpyne.run_sim(config_name, batchParams)

  print(f"!!! Simulation ran for {batchParams['num_syns_E']} {batchParams['syns_type']} synpases !!!")
