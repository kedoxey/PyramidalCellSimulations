import requests
import json
import os
import pyneuroml as pynml
import neuroml as nml
import neuroml.writers as writers
import numpy as np
import re
from neuroml.utils import component_factory
from zipfile import ZipFile
from urllib.request import urlopen

from neuromllite import Network, Cell, Population, InputSource, Input, NetworkGenerator

cwd = os.getcwd()


def compile_mechs(cwd, hocs_dir, mod_dir, force=False):
    if force or not os.path.exists(os.path.join(hocs_dir, 'x86_64')):
        os.chdir(hocs_dir)
        os.system(f'nrnivmodl {mod_dir}')
        os.chdir(cwd)
    else:
        print('Mechanisms already compiled!')


# download specified version of model from neuroml-db 
def download_from_nmldb(model_id, version):
    
    zip_url = f'https://neuroml-db.org/GetModelZip?modelID={model_id}&version={version}'

    unzip_path = os.path.join(cwd,'models',version,f'{model_id}-{version}')
    zip_path = os.path.join(cwd,'models','zips',model_id+'.zip')

    if not os.path.exists(unzip_path):
        # download model if it has not been already
        nmldb_response = requests.get(zip_url)
        open(zip_path,'wb').write(nmldb_response.content)

        os.makedirs(unzip_path)

        with ZipFile(zip_path,'r') as zObject:
            zObject.extractall(path=unzip_path)

        print(f'Model {model_id} successfully downloaded!')
    else:
        print(f'Model {model_id} already downloaded.')

    
def generate_network(nml_dir, cell_name, pop_label, pop_size=1, force=False, **input_args):
    
    cell_nml_path = os.path.join(nml_dir, f'{cell_name}.cell.nml')
    net_nml_path = os.path.join(nml_dir, f'{cell_name}_Net.net.nml')

    if force or not os.path.exists(net_nml_path):

        net = Network(id=f'{cell_name}_Net')
        net.notes = 'Hay et al. 2011 model of Neocortex Layer 5 Pyramidal Cell for import into NetPyNE'
        net.parameters = {'input_amp': input_args['input_amp']}

        cell = Cell(id=cell_name, neuroml2_source_file=cell_nml_path)
        net.cells.append(cell)

        pop = Population(id=pop_label,
                         size=pop_size,
                         component=cell.id)
        net.populations.append(pop)

        NetworkGenerator.generate_neuroml2_from_network(net,
                                                        nml_file_name=net_nml_path,
                                                        format='xml')
        
        print(f'Network for {cell_name} successfully generated!')
    else:
        print(f'Network for {cell_name} already generated!')

    return net_nml_path

def create_output_dirs(test_name, model_dir):

    output_dir = os.path.join(model_dir,'output')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    sim_dir = os.path.join(output_dir, test_name)
    if not os.path.exists(sim_dir):
        os.mkdir(sim_dir)

    return output_dir, sim_dir

def get_components(cell):
    secs = np.array(list(cell['secs'].keys()))

    apics = secs[np.where(np.char.find(secs,'apic')>=0)]
    dends = secs[np.where(np.char.find(secs,'dend')>=0)]

    apics = sorted(apics, key=lambda s: int(re.search(r'\d+', s).group()))
    dends = sorted(dends, key=lambda s: int(re.search(r'\d+', s).group()))

    temp = 1
