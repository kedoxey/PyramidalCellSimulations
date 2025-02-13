import json
import os
import requests

import numpy as np

from neuromllite import Network, Cell, Population, NetworkGenerator
from zipfile import ZipFile

cwd = os.getcwd()


### Functions for compiling and parsing model hoc files and downloading models from NeuroML-DB ###
def compile_mechs(cwd, hocs_dir, mod_dir, force=False):
    if force or not os.path.exists(os.path.join(hocs_dir, 'x86_64')):
        os.chdir(hocs_dir)
        os.system(f'nrnivmodl {mod_dir}')
        os.chdir(cwd)
    else:
        print('Mechanisms already compiled!')


# get the name of the cell from the hoc file
def get_cell_name(model_dir):
    
    model_files = np.array(os.listdir(model_dir))
    hoc_file = model_files[np.char.endswith(model_files,'.hoc')][0]
    cell_name = hoc_file.split('.')[0]
    # for model_file in model_files:
    #     if '.hoc' in model_file:
    #         cell_name = model_file.split('.')[0]
    #         break
    return cell_name


# download specified version of model from neuroml-db 
def download_from_nmldb(model_id, version):

    # get author and year of corresponding publication to population the cell_model parameter
    detail_url = f'https://neuroml-db.org/api/model?id={model_id}'
    detail_response = requests.get(detail_url, verify=False)
    detail_dict = json.loads(detail_response.content)

    author = detail_dict['publication']['authors'][0]['Person_Last_Name']
    year = detail_dict['publication']['record']['Year']
    author_year = f'{author}{year}'

    zip_url = f'https://neuroml-db.org/GetModelZip?modelID={model_id}&version={version}'

    unzip_path = os.path.join(cwd,'models',version,f'{model_id}-{version}')
    zip_path = os.path.join(cwd,'models','zips',model_id+'.zip')

    if not os.path.exists(unzip_path):
        # download model if it has not been already
        nmldb_response = requests.get(zip_url, verify=False)
        open(zip_path,'wb').write(nmldb_response.content)

        os.makedirs(unzip_path)

        with ZipFile(zip_path,'r') as zObject:
            zObject.extractall(path=unzip_path)

        # copy_synapses(unzip_path)

        print(f'Model {model_id} successfully downloaded!')
    else:
        print(f'Model {model_id} already downloaded.')

    return author_year


def get_model_details(nmldb_id):

    nmldb_url = f'http://neuroml-db.org/api/model?id={nmldb_id}'

    nmldb_response = requests.get(nmldb_url, verify=False)

    return nmldb_response.json()


def get_morphophetrics(nmldb_id):

    nmldb_url = f'http://neuroml-db.org/api/morphometrics?id={nmldb_id}'

    nmldb_response = requests.get(nmldb_url, verify=False)

    return nmldb_response.json()


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
