# Characterize single cell dynamics in NetPyNE.

In [batch_simulation.py](https://github.com/kedoxey/PyramidalCellSimulations/blob/dev/batch_simulation.py), simulations can be run one at a time or in batch. Parameter sets are handled by the [ParameterGrid](https://scikit-learn.org/dev/modules/generated/sklearn.model_selection.ParameterGrid.html) module from [scikit-learn](https://scikit-learn.org/dev/index.html). A dictionary consisting of lists of parameter values is passed to `ParamaterGrid` and a list of all possible paramater value combinations is returned. The default parameter set can be found in [config/default-config.py](https://github.com/kedoxey/PyramidalCellSimulations/blob/dev/config/default_config.yml). All parameters defined with `ParameterGrid` will overwrite the default values.


#### Current supported models:
* Hay et al. 2011 Layer 5b Pyramidal Cell from [NeuroML-DB](https://neuroml-db.org/model_info?model_id=NMLCL000073 "NeuorML-DB")
* Hay and Segev 2015 excitatory (AMPA, NMDA) and inhibitory (GABAA, GABAB) synpases from [ModelDB](https://modeldb.science/156780 "ModelDB")

## Output Location
All output from the simulation will be saved in the following path:
```
models/NEURON/<nmldb-id>-NEURON/output/<sim_name>/<sim_label>
```

## Example Parameter Sets
### Injected Current
#### Single Simulation
Run a simulation for *1000 ms* with injected current to the soma with an amplitude of *0.44 nA* for *600 ms* with a *200 ms* delay.
```
paramGrids = {'sim_name': ['INJ'],
              'sim_flag': ['iclamp'],
              'nmldb_id': ['NMLCL000073'],
              'enable_syns': [False],
              'input_amp': [0.44],
              'input_sec': ['soma_0'],
              'sim_dur': [1000],
              'stim_dur': [600],
              'stim_delay': [200]}
```
#### Batch Simulations
Run simulations for *1000 ms* each with injected current to the soma with an amplitudes of *0.44 nA*, *0.63 nA*, and *0.8 nA* for *600 ms* with a *200 ms* delay. `ParameterGrid` will return 3 parameter sets with unique values for the parameter `input_amp`. 
```
paramGrids = {'sim_name': ['INJ'],
              'sim_flag': ['iclamp'],
              'nmldb_id': ['NMLCL000073'],
              'enable_syns': [False],
              'input_amp': [0.44, 0.63, 0.8],
              'input_sec': ['soma_0'],
              'sim_dur': [1000],
              'stim_dur': [600],
              'stim_delay': [200]}
```
### Synaptic Current
#### Single Simulation
Run a simulation for *1000 ms* with *50* excitatory synapses on basal dendrites that are stimulated by Poisson spike trains with frequency *40 Hz* for *600 ms* with a *200 ms* delay.
```
paramGrids = {'sim_name': ['EXC'],
             'sim_flag': ['syns'],
             'nmldb_id': ['NMLCL000073'],
             'syns_type': ['basal'],  
             'num_syns_E': [50],
             'sim_dur': [1000],
             'stim_dur': [600],
             'stim_delay': 200]}
```
#### Batch Simulations
Run simulations for *1000 ms* each with *50* and *100* excitatory synapses on basal and apical distal dendrites that are stimulated by Poisson spike trains with frequency *40 Hz* for *600 ms* with a *200 ms* delay. `ParameterGrid` will return 4 parameter sets with unique combinations of the parameters `syns_type` and `num_syns_E`. 
```
paramGrids = {'sim_name': ['EXC'],
             'sim_flag': ['syns'],
             'nmldb_id': ['NMLCL000073'],
             'syns_type': ['basal', 'apical_distal'],  
             'num_syns_E': [50, 100],
             'enable_syns': [True],
             'record_LFP': [False],
             'sim_dur': [1000],
             'stim_dur': [600],
             'stim_delay': [200]}
```
### Combined
Run simulations with different input types.
```
# unique parameters for each parameter set
param_sets = {'syns': {'enable_syns': [True],
                       'syns_type': ['basal', 'apical_distal'],
                       'num_syns_E': [50, 100]},
              'iclamp': {'enable_syns': [False],
                         'input_amp': [0.44, 0.63, 0.8],
                         'input_sec': ['soma_0']}}

paramGrids = []

for sim_flag, param_set in param_sets.items():
  # common parameters between parameter sets
  paramGrid = {'sim_name': ['IN'],
               'sim_flag': [sim_flag],
               'nmldb_id': ['NMLCL000073'],
               'record_LFP': [False],
               'sim_dur': [1000],
               'stim_dur': [600],
               'stim_delay': [200]}
	            
  for param_name, param in param_set.items():
    paramGrid[param_name] = [param]
		
  paramGrids.append(paramGrid)
```
### Different number of synpases for each location
Run simulations with different numbers of synapses for each group of synapse locations.
```
# unique parameters for each parameter set
param_sets = {'syns': {'enable_syns': [True]}
group_num_syns = {'basal': [50, 100],
                  'apical_distal': [200, 300]}

paramGrids = []

for sim_flag, param_set in param_sets.items():
  for syn_group, num_syns in group_num_syns.items():
    # common parameters between parameter sets
    paramGrid = {'sim_name': ['IN'],
                 'sim_flag': [sim_flag],
                 'nmldb_id': ['NMLCL000073'],
                 'syns_type': [syn_group],
                 'num_syns_E': num_syns,
                 'record_LFP': [False],
                 'sim_dur': [1000],
                 'stim_dur': [600],
                 'stim_delay': [200]}
  	            
    for param_name, param in param_set.items():
      paramGrid[param_name] = [param]
  		
    paramGrids.append(paramGrid)
```
