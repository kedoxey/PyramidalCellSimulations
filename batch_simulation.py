import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid


config_name = 'default_config'

paramGrids = {'sim_name': ['INJ'],
              'sim_flag': ['iClamp'],

              # 'local': [True],
              'nmldb_id': ['NMLCL000073'],
              # 'publication': ['Markram2015'],
              # 'model_name': ['Birgiolas2020'],
              # 'hocs_dname': ['Cells'],
              # 'mod_dname': ['Mechanisms'],
              # 'cell_type': ['GC'],
              # 'cell_model': ['Birgiolas2020'],
              # 'cell_name': ['GC1'],

              'enable_syns': [False],
              # 'log_firing_rate': [True],
              # 'syns_type': ['apical'],
              # 'num_syns_E': [0, 50],
              'input_amp': [0.15],
              'input_sec': ['soma_0'],
              'record_LFP': [False],
              'sim_dur': [1000],
              'stim_dur': [1000-100],
              'stim_delay': [50]}
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  simulate_netpyne.run_sim(config_name, batchParams)

# TODO: plot firing rates
