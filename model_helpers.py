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
from neuron import h

from neuromllite import Network, Cell, Population, InputSource, Input, NetworkGenerator

cwd = os.getcwd()

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

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


def create_sim_description(output_dir, **params):

    file_name = 'README.md'
    file_path = os.path.join(output_dir, file_name)

    if not os.path.exists(file_path):

        message = '--- Simulation Parameters ---\n'

        for key, value in params.items():
            if 'message' in key:
                message = f'--- Simulation Description ---\n{value}\n' + message
            else:
                message += f'{key} = {value}\n'

    f = open(file_path, 'w')
    f.write(message)
    f.close()
    
    print("Simulation description saved!")


def get_components(cell, group_name='all'):

    secs = np.array(list(cell['secs'].keys()))
    basal_group = list(secs[np.where(np.char.find(secs,'dend')>=0)])
    apic_group = list(secs[np.where(np.char.find(secs,'apic')>=0)])
    all_group = list(secs)
    soma = ['soma_0']

    if group_name == 'basal':
        return basal_group
    elif group_name == 'apical':
        return apic_group
    elif group_name == 'basal_apical':
        return basal_group + apic_group
    elif group_name == 'basal_soma':
        return basal_group + soma
    elif group_name == 'apical_soma':
        return apic_group + soma
    elif group_name == 'basal_apical_soma':
        return basal_group + apic_group + soma
    elif group_name == 'soma':
        return soma
    else:
        return all_group


def get_hoc_cell(filename, cell_name):

    h.load_file(filename)

    hoc_cell = getattr(h, cell_name)
    hoc_cell = hoc_cell()
    cell = hoc_cell

    return cell


def get_total_soma_distance(cell):
    
    soma = cell.soma_0

    dist_from_soma = 0
    furthest_sec = ''

    for sec in h.allsec():
        sec_name = sec.name().split('.')[1]
        
        dist = h.distance(soma(0.5), sec(0.5))
        if dist > dist_from_soma:
            dist_from_soma = dist
            furthest_sec = sec_name

        temp = 6

    return dist_from_soma, furthest_sec


def get_secs_from_dist(filename, cell_name, lb, ub=1):
    cell = get_hoc_cell(filename, cell_name)
    soma = cell.soma_0

    total_distance, _ = get_total_soma_distance(cell)
    distance = lb*total_distance
    dist_bound = ub*total_distance

    secs_from_dist = []

    for sec in h.allsec():
        sec_name = sec.name().split('.')[1]
        
        dist = h.distance(soma(0.5), sec(0.5))
        if (distance < dist) and (dist < dist_bound):
            secs_from_dist.append(sec_name)

    return secs_from_dist
