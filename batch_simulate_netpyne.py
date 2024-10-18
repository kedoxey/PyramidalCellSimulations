import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid
from itertools import chain


configName = 'exc_config'

param_sets = {'-syns': {'enable_syns': True}}
              # '_no-Na': {'channel_toggles': {'Na': 0}}}
group_num_syns = {'soma': [10, 50, 100, 200, 300],  # [5*i for i in range(50,61)],
                  'basal': list(chain([10, 50, 100, 200, 300],[5*i for i in range(61,101)])),
                  'apical_distal': list(chain([10, 50, 100, 200, 300],[5*i for i in range(181,241)]))}
                  # 'basal': 30}
                  # 'basal_soma': 30,
                  # 'apical_proximal': 35,
                  # 'apical_distal': 700}
paramGrids = []

for sim_flag, param_set in param_sets.items():
    for syn_group, num_syns in group_num_syns.items():
      for num_syn in num_syns:

        paramGrid = {'sim_name': ['LFP'],
                    'sim_flag': [sim_flag],
                    'nmldb_id': ['NMLCL000073'],  # NMLCL000660, NMLCL000073
                    'syns_type': [syn_group],  #, 'apical_proximal', 'apical_distal', 'basal'],  #, 'apical_distal'],  #'soma', 'basal', 'apical'],  # 'basal_apical', 'basal_soma', 'apical_soma', 'basal_apical_soma', 'all'
                    'num_syns_E': [num_syn], # [50, 70],  # 90, 110],
                    'add_bkg': [False], 
                    'record_LFP': [True],
                    'depths': [4],
                    'sim_dur': [1000],
                    'stim_dur': [400],
                    'stim_delay': [400]}  #[0, 1.0], [1.0, 0]]}
        
        for param_name, param in param_set.items():
          paramGrid[param_name] = [param]

        paramGrids.append(paramGrid)
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

    simulate_netpyne.run_sim(configName, batchParams)

        # cmd_line_txt = f'python simulate_netpyne.py {config_name} {batch_param} {batch_value}'
        
        # os.system(cmd_line_txt)
