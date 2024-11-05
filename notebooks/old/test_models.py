import matplotlib
import os
import json
import model_helpers as mh
from netpyne import specs, sim


# nmldb_id = "NMLCL000073"
model_source = 'ModelDB'
model_name = "L5PCbiophys"
model_version = '1'  # '1', '2', '3', '4'

cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models')
model_dir = os.path.join(models_dir, model_source, 'L5bPCmodelsEH')
hocs_dir = os.path.join(model_dir,'models')
mod_dir = os.path.join(model_dir, 'mod')

mh.compile_mechs(cwd, hocs_dir, mod_dir)

hoc_file = os.path.join(hocs_dir, model_name+model_version+'.hoc')

from neuron import h, load_mechanisms

load_mechanisms('/home/kedoxey/CRCNS/PyramidalCellSimulations/models/ModelDB/L5bPCmodelsEH/mod')

# h.load_file(hoc_file)

# Network parameters
netParams = specs.NetParams()

netParams.importCellParams(
    label='L5PCbiophys1_hoc',
    conds={'cellType': 'PYR', 'cellModel': 'Hay2011'},
    fileName=hoc_file,
    cellName=model_name,
    importSynMechs=False
)

temp = 5
