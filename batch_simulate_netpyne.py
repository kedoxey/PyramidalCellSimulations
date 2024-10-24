import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid
from itertools import chain


config_name = 'default_config'

sim_name = 'FR'

param_sets = {'active': {'log_firing_rate': True}}#,
              # 'pas_soma': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'soma'},
              # 'pas_all': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'all'}}
              
# group_num_syns = {'soma': [25, 50, 100, 200, 300],
#                   'basal': [25, 50, 100, 200, 300],
#                   'apical_distal': [25, 50, 100, 200, 300]}

paramGrids = []

for sim_flag, param_set in param_sets.items():

  paramGrid = {'sim_name': [sim_name],
                'sim_flag': [sim_flag],
                'nmldb_id': ['NMLCL000073'],
                'enable_syns': [True],
                'syns_type': ['soma', 'basal', 'apical_distal'],  
                'num_syns_E': [25, 50, 100, 200, 300],
                'add_bkg': [False], 
                'record_LFP': [True],
                'depths': [4],
                'apical_depths': [2],
                'sim_dur': [1000],
                'stim_dur': [400],
                'stim_delay': [400],
                'save_pickle': [False]}
  
  for param_name, param in param_set.items():
    paramGrid[param_name] = [param]

  paramGrids.append(paramGrid)
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  simulate_netpyne.run_sim(config_name, batchParams)

      # cmd_line_txt = f'python simulate_netpyne.py {config_name} {batch_param} {batch_value}'
      
      # os.system(cmd_line_txt)
