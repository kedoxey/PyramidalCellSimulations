import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid
from itertools import chain


config_name = 'default_config'

sim_name = 'LFP'
param_sets = {'syns': {'enable_syns': True}}  #,
              # '-pass_soma': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'soma'},
              # '-pass_all': {'channel_toggles': {'Na': 0, 'Ca': 0, 'K': 0},
              #                'channel_secs': 'all'}}
              # '_no-Na': {'channel_toggles': {'Na': 0}}}
group_num_syns = {#'soma': [25, 50, 100, 200, 300],
                  # 'basal': [25, 50, 100, 200, 300],
                  # 'apical_distal': [25, 50, 100, 200, 300]}
  #  'soma': [10, 50, 100, 200, 300],  # [5*i for i in range(50,61)],
  #                 'basal': list(chain([10, 50, 100, 200, 300],[5*i for i in range(61,101)])),
                  'apical_distal': [5*i for i in range(228,320)]}

paramGrids = []

for sim_flag, param_set in param_sets.items():
  for syn_group, num_syns in group_num_syns.items():
    for num_syn in num_syns:

      paramGrid = {'sim_name': [sim_name],
                    'sim_flag': [sim_flag],
                    'nmldb_id': ['NMLCL000073'],
                    'syns_type': [syn_group],  
                    'num_syns_E': [num_syn],
                    'add_bkg': [False], 
                    'record_LFP': [True],
                    'depths': [4],
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
