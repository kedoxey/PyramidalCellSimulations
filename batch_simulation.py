import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid


config_name = 'default_config'

paramGrids = {'sim_name': ['EXC'],
              'sim_flag': ['syns'],
              'nmldb_id': ['NMLCL000073'],
              'syns_type': ['basal', 'soma'],  
              'num_syns_E': [50, 100],
              'enable_syns': [True],
              'record_LFP': [False],
              'sim_dur': [1000],
              'stim_dur': [600],
              'stim_delay': [200]}
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  simulate_netpyne.run_sim(config_name, batchParams)
