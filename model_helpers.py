import requests
import json
import os
import yaml
import shutil
import pickle
import pyneuroml as pynml
import neuroml as nml
import neuroml.writers as writers
import numpy as np
import matplotlib.pyplot as plt
import re
from neuroml.utils import component_factory
from zipfile import ZipFile
from urllib.request import urlopen
from scipy.signal import butter, lfilter
from neuron import h

from neuromllite import Network, Cell, Population, InputSource, Input, NetworkGenerator

cwd = os.getcwd()

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


### Model files ###
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

def copy_synapses(model_dir):

    synapse_dir = os.path.join(cwd,'models','synapses')
    synapse_files = os.listdir(synapse_dir)
    for synapse_file in synapse_files:
        src = os.path.join(synapse_dir,synapse_file)
        dst = os.path.join(model_dir,synapse_file)

        if not os.path.exists(dst):
            with open(src, 'r') as f_in:
                data_in = f_in.read()
            with open(dst, 'w') as f_out:
                f_out.write(data_in)
        # shutil.copy2(src, dst)


# download specified version of model from neuroml-db 
def download_from_nmldb(model_id, version):
    
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


### For config file ###
def join(loader,node):
    seq = loader.construct_sequence(node)
    return ''.join([str(i) for i in seq])

def load_config(config_name='default_config'):
    cwd = os.getcwd()
    config_dir = os.path.join(cwd, 'config')
    config_file = os.path.join(config_dir, config_name+'.yml')

    yaml.add_constructor('!join', join)
    with open(config_file) as f:
        config_params = yaml.full_load(f)

    config_params['input_amp'] = config_params['input_amps'][config_params['amp_idx']]
    
    return config_params

def save_config(sim_dir, sim_label, config_name='default_config'):
    cwd = os.getcwd()
    config_dir = os.path.join(cwd, 'config')
    config_file = os.path.join(config_dir, config_name+'.yml')

    shutil.copy2(config_file, os.path.join(sim_dir,'config-'+sim_label+'.yml'))
    
    print('Config saved!')


### File system ###
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


### Model analysis ###
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

def get_rand_secs(sec_list, num_syns_E, num_syns_I):

    sec_list = np.array(sec_list)

    locs_E = np.random.randint(len(sec_list),size=num_syns_E)
    locs_I = np.random.randint(len(sec_list),size=num_syns_I)

    return list(sec_list[locs_E]), list(sec_list[locs_I])


### Signal processing ###
def butter_bandpass(lowcut, highcut, fs, order=5):
 nyq = 0.5 * fs
 low = lowcut / nyq
 high = highcut / nyq
 b, a = butter(order, [low, high], btype='band')
 return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
   b, a = butter_bandpass(lowcut, highcut, fs, order=order)
   y = lfilter(b, a, data)
   return y

def get_filtered_signal(lfp, dt):

    lfp = [lfp_d[0] for lfp_d in lfp]

    lfp_bp_low = butter_bandpass_filter(lfp, 10, 300, 1/dt*1000, order=4)
    lfp_bp_spikes = butter_bandpass_filter(lfp, 301, 700, 1/dt*1000, order=4)

    return lfp_bp_low, lfp_bp_spikes

def plot_lfp(simData, dt, sim_label, output_dir):
    with open(os.path.join(output_dir,f'{sim_label}_simData.pkl'),'wb') as fp:
        pickle.dump(simData,fp)

    t = list(simData['t'])
    with open(os.path.join(output_dir,f'{sim_label}_t.pkl'),'wb') as fp:
        pickle.dump(t,fp)

    spkt = list(simData['spkt'])
    with open(os.path.join(output_dir,f'{sim_label}_spkt.pkl'),'wb') as fp:
        pickle.dump(spkt,fp)

    # lfp_bp_low, lfp_bp_spikes = get_filtered_signal(simData['LFP'], dt)

    # with open(os.path.join(output_dir,f'{sim_label}lfp_bp_low.pkl'),'wb') as fp:
    #     pickle.dump(lfp_bp_low,fp)
    
    # with open(os.path.join(output_dir,f'{sim_label}lfp_bp_spikes.pkl'),'wb') as fp:
    #     pickle.dump(lfp_bp_spikes,fp)

    # fig, ax = plt.subplots(1, 1, figsize=(20,12))

    # ax.plot(t[0:len(lfp_bp_low)],lfp_bp_low*10000-200,color='cornflowerblue',label='BP filtered: 10-300Hz')
    # ax.plot(t[0:len(lfp_bp_spikes)],lfp_bp_spikes*10000-800,color='orchid',label='BP filtered: >300Hz')
    # ax.vlines(spkt, [195], [205], 'k')
    # ax.legend(loc='upper right')
    # ax.set_xlabel('Time (ms)')
    # ax.set_ylabel('LFP')      
    # ax.set_title('LFPs')


    # fig.savefig(os.path.join(output_dir,'lfp_fig.png'),bbox_inches='tight',dpi=300)

