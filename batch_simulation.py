import os
import simulate_netpyne
import model_helpers as mh
from sklearn.model_selection import ParameterGrid


config_name = 'default_config'

paramGrids = {'sim_name': ['INJ'],
              'sim_flag': ['iclamp'],

              'local': [True],
              'model_name': ['Birgiolas2020'],
              'hocs_dname': ['Cells'],
              'mod_dname': ['Mechanisms'],
              'cell_type': ['MIT'],
              'cell_model': ['Birgiolas2020'],
              'cell_name': ['MC5'],

              'enable_syns': [False],
              'input_amp': [0.25],
              'input_sec': ['soma'],
              'record_LFP': [False],
              'sim_dur': [3000],
              'stim_dur': [2500],
              'stim_delay': [200]}
      
batchParamsList = list(ParameterGrid(paramGrids))

for batchParams in batchParamsList:

  simulate_netpyne.run_sim(config_name, batchParams)
