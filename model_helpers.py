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
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import re
from neuroml.utils import component_factory
from zipfile import ZipFile
from urllib.request import urlopen
from scipy import stats
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

    # config_params['input_amp'] = config_params['input_amp']
    
    return config_params

def write_config(config, sim_dir, sim_label, config_name='default_config'):
    cwd = os.getcwd()

    with open(os.path.join(sim_dir,'config-'+sim_label+'.yml'), 'w') as outfile:
        yaml.dump(config, outfile)

def save_config(sim_dir, sim_label, config_name='default_config'):
    cwd = os.getcwd()
    config_dir = os.path.join(cwd, 'config')
    config_file = os.path.join(config_dir, config_name+'.yml')

    shutil.copy2(config_file, os.path.join(sim_dir,'config-'+sim_label+'.yml'))
    
    print('Config saved!')


### File system ###
def create_output_dirs(sim_name, sim_label, model_dir):

    output_dir = os.path.join(model_dir,'output')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    sim_top_dir = os.path.join(output_dir, sim_name)
    if not os.path.exists(sim_top_dir):
        os.mkdir(sim_top_dir)

    sim_dir = os.path.join(sim_top_dir, sim_label)
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
def get_compartments(hoc_fname, cell, cell_name, group_name='all'):

    secs = np.array(list(cell['secs'].keys()))
    basal_group = list(secs[np.where(np.char.find(secs,'dend')>=0)])
    apic_group = list(secs[np.where(np.char.find(secs,'apic')>=0)])
    all_group = list(secs)
    soma = ['soma_0']

    if group_name == 'basal':
        return basal_group
    elif group_name == 'apical':
        return apic_group
    elif group_name == 'apical_distal':
        return get_secs_from_dist(hoc_fname, cell_name, 0.5, 1, secs_lim='apic')
    elif group_name == 'apical_proximal':
        return get_secs_from_dist(hoc_fname, cell_name, 0, 0.5, secs_lim='apic')
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


def toggle_channels(cell_params, sec_names, toggles):  #, ion_name, toggle):

    for sec_name in sec_names:
        sec = cell_params['secs'][sec_name]

        for ion_name, toggle in toggles.items():
            for mech_name, mech in list(sec['mechs'].items()):
                if (ion_name in mech_name) and ('gmax' in mech.keys()):
                    mech['gmax'] *= toggle
    
    return cell_params


def get_secs_from_dist(filename, cell_name, lb, ub=1, secs_lim='all'):
    cell = get_hoc_cell(filename, cell_name)
    soma = cell.soma_0

    total_distance, _ = get_total_soma_distance(cell)
    dist_lb = lb*total_distance
    dist_ub = ub*total_distance

    secs_from_dist = []

    for sec in h.allsec():
        sec_name = sec.name().split('.')[1]

        get_dist = True
        if 'all' not in secs_lim:
            if secs_lim not in sec_name:
                get_dist = False

        if get_dist:
            dist = h.distance(soma(0.5), sec(0.5))
            if (dist_lb < dist) and (dist < dist_ub):
                secs_from_dist.append(sec_name)


    return secs_from_dist

def get_rand_secs(sec_list, num_syns_E, num_syns_I, seed=74):

    sec_list = np.array(sec_list)

    np.random.seed(seed)

    locs_E = np.random.randint(len(sec_list),size=num_syns_E)
    locs_I = np.random.randint(len(sec_list),size=num_syns_I)

    return list(sec_list[locs_E]), list(sec_list[locs_I])

def get_rand_sec(sec_list, seed=74):
    
    sec_list = np.array(sec_list)

    np.random.seed(seed)

    return sec_list[np.random.randint(len(sec_list))]


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


### Figure Customization ###
def get_colormaps(numSynsE=None, numSynsI=None):
    
    if numSynsE:
        cmap = mpl.colormaps['YlOrRd']
        colorsE = cmap(np.linspace(0,1,numSynsE))
    else:
        colorsE = None

    if numSynsI:
        cmap = mpl.colormaps['GnBu']
        colorsI = cmap(np.linspace(0,1,numSynsI))
    else:
        colorsI = None
    
    return colorsE, colorsI

def get_syn_sec_colors(cell, colormaps, synColors):
    secSynCount = 0
    secSynColors = {}
    for secName, sec in cell.secs.items():
        for synMech in sec['synMechs']:
            secSynColors[secName] = {'E': synColors['E'],
                                     'I': synColors['I']}
            if 'GABA' in synMech['label']:
                secSynColors[secName]['I'] = colormaps[1][secSynCount]
            else:
                secSynColors[secName]['E'] = colormaps[0][secSynCount]

        if bool(sec['synMechs']):
            secSynCount += 1

    return secSynColors


### Plotting ###
def plot_isolated_LFP(simData, syns_type, num_syns, sim_label, sim_dir, output_dir):
    t = np.array(simData['t'])
    dt = t[1] - t[0]

    V_soma = np.array(simData['V_soma']['cell_0'])
    t_spikes = t[np.where(V_soma>10)]
    # t_bound = 600 if len(np.where(t_spikes>600)[0])>0 else 500
    
    # t_bound = 500 if 'distal' in syns_type else 600

    plot_flag = True

    if len(t_spikes) > 0:
        t_bound = int(100*np.floor(t_spikes[-1]/100))
        t_spike = t_spikes[np.where(t_spikes>t_bound)[0][0]]
        slice_start = int((t_spike - 2.25)/dt)
        slice_end = int((t_spike + 4.5)/dt)

        save_eap_time(syns_type, num_syns, slice_start, slice_end, t_spike, output_dir)
    else:
        time_window_path = os.path.join(output_dir,'eap_time_windows.pkl')
        with open(time_window_path,'rb') as fp:
            time_windows = pickle.load(fp)
        
        try:
            slice_start = time_windows[syns_type][num_syns][0]
            slice_end = time_windows[syns_type][num_syns][1]
            t_spike = time_windows[syns_type][num_syns][2]
        except KeyError:
            plot_flag = False

    if plot_flag:
        t_slice = t[slice_start:slice_end]
        
        lfp = simData['LFP']
        num_elecs = len(lfp[0])
        lfp_elec = {i: np.zeros(len(lfp)) for i in range(num_elecs)}
        for t_i, lfp_t in enumerate(lfp):
            for elec_i, lfp_t_i in enumerate(lfp_t):
                lfp_elec[elec_i][t_i] = lfp_t_i

        colors = {0: (91/255,154/255,204/255),
                1: (223/255,183/255,10/255),
                2: (89/255,199/255,128/255),
                3: (221/255,59/255,5/255),
                4: (72/255,156/255,155/255),
                5: (223/255,132/255,7/255)}

        fig, axs = plt.subplots(num_elecs//2, 2, figsize=(5,10))
        axs = axs.ravel()

        for ax_i, i in enumerate(reversed(range(num_elecs))):

            lfp_i = lfp_elec[i]
            lfp_slice = lfp_i[slice_start:slice_end]

            height = max(lfp_slice) - min(lfp_slice)
            lfp_slice_norm = lfp_slice/height
            lfp_slice_norm -= lfp_slice_norm[0]

            axs[ax_i].plot(t_slice,lfp_slice_norm,color=colors[i],zorder=12)
            axs[ax_i].axvline(t_spike,alpha=0.4,color='k',zorder=1,linestyle='--')
            xticks = [(int(t_spike.round(0))-2)+2*i for i in range(4)]
            axs[ax_i].set_xticks(xticks)
            axs[ax_i].set_title(i)

        axs[-1].set_xlabel('Time (ms)')
        axs[-2].set_xlabel('Time (ms)')

        fig.suptitle('Normalized Amplitude')
        fig.tight_layout()

        fig.savefig(os.path.join(sim_dir,f'{sim_label}-isolated_LFP.png'),bbox_inches='tight',dpi=300)


def save_eap_time(syns_type, num_syns, slice_start, slice_end, t_spike, output_dir):

    file_name = 'eap_time_windows.pkl'
    file_path = os.path.join(output_dir,file_name)

    if os.path.exists(file_path):
        with open(file_path,'rb') as fp:
            time_windows = pickle.load(fp)

        if syns_type not in time_windows.keys():
            time_windows[syns_type] = {num_syns: (slice_start, slice_end, t_spike)}
        else:
            time_windows[syns_type][num_syns] = (slice_start, slice_end, t_spike)
    else:
        time_windows = {syns_type: {num_syns: (slice_start, slice_end, t_spike)}}

    with open(file_path,'wb') as fp:
        pickle.dump(time_windows,fp)

def save_firing_rate(simData, stim_delay, stim_dur, syns_type, num_syns, output_dir):
    
    V_soma = np.array(simData['V_soma']['cell_0'])
    t = np.array(simData['t'])
    t_spikes = t[np.where(V_soma>0)]
    
    if len(t_spikes) > 0:
        t_spikes_stim = t_spikes[np.where((t_spikes > stim_delay) & (t_spikes < stim_delay + stim_dur))]
        firing_rate = len(t_spikes_stim)
    else:
        firing_rate = 0
    
    file_name = 'firing_rates.pkl'
    file_path = os.path.join(output_dir, file_name)
    
    if os.path.exists(file_path):
        with open(file_path,'rb') as fp:
            firing_rates = pickle.load(fp)

        try:
            firing_rates[syns_type][num_syns] = firing_rate
        except KeyError:
            firing_rates[syns_type] = {num_syns: firing_rate}

    else:
        firing_rates = {syns_type: {num_syns: firing_rate}}
        
    with open(file_path,'wb') as fp:
        pickle.dump(firing_rates,fp)
                

def plot_pre_spike_trains(cells, conns, sim_label, sim_dir):

    spike_trains = {}

    fig, axs = plt.subplots(figsize=(5,8))

    for conn in conns[0]:
        sec_name = conn['sec']
        if conn['hObj'].pre():
            pre_gid = int(re.findall(r'\d+', conn['hObj'].pre().hname())[0])
            train_name = f'{sec_name}_{pre_gid}'

            if train_name not in spike_trains.keys():
                pre_spikes = list(cells[pre_gid+1].hSpkTimes)
                spike_trains[train_name] = pre_spikes

                axs.vlines(pre_spikes,pre_gid-0.25,pre_gid+0.25)
        
    axs.set_yticks(range(len(spike_trains.keys())))
    axs.set_yticklabels(spike_trains.keys())
    axs.set_xlim([0,1500])

    fig.savefig(os.path.join(sim_dir,f'{sim_label}-presyn_spike_trains.png'),bbox_inches='tight',dpi=300)

    return spike_trains


def plot_soma(simData, sim_label, sim_dir):

    t = np.array(simData['t'])
    V_soma = np.array(simData['V_soma']['cell_0'])

    ### PLOT SOMA MEMBRANE POTENTIAL ###
    fig, axs = plt.subplots(figsize=(8,5))
    axs.plot(t, V_soma)
    axs.set_title(f'Soma Membrane Potential')
    axs.set_ylabel('Voltage (mV)')
    axs.set_xlabel('Time (ms)')

    fig.tight_layout()
    fig.savefig(os.path.join(sim_dir,f'{sim_label}-soma_pot.png'),bbox_inches='tight',dpi=300)

    ### PLOT SOMA AND APICAL DISTAL COMPARTMENT MEMBRANE POTENTIALS FOR BACKPROPAGATION ###
    fig, axs = plt.subplots(figsize=(8,5))
    axs.plot(t, V_soma, label='soma')
    axs.plot(t, np.array(simData['V_apic_32']['cell_0']), label='apic_32')
    axs.legend(loc='upper right')
    axs.set_title(f'Membrane Potentials')
    axs.set_ylabel('Voltage (mV)')
    axs.set_xlabel('Time (ms)')

    fig.tight_layout()
    fig.savefig(os.path.join(sim_dir,f'{sim_label}-backpropagation.png'),bbox_inches='tight',dpi=300)


def plot_secs(simData, spike_trains, sim_label, sim_dir, sec_syn_colors):

    sec_traces = list(simData.keys())
    [sec_traces.remove(key) for key in ['spkt', 'spkid', 't', 'V_soma', 'avgRate', '__dict__'] if key in sec_traces]
    num_secs = len(sec_traces)//3

    t = np.array(simData['t'])
    V_soma = np.array(simData['V_soma']['cell_0'])
    avg_V_soma = np.average(V_soma)
    t_spikes = t[np.where(V_soma>-10)]
    
    ### PLOT COMPARTMENT MEMBRANE POTENTIALS AND SYNAPTIC CONDUCTANCES ###
    ## combined synaptic conductances ###
    if num_secs > 0:
        fig, axs = plt.subplots(num_secs, 2, figsize=(10,num_secs*2))
        axs.ravel()

        for i in range(num_secs):

            sec_i = 3*i

            if num_secs > 1:
                idx_0 = (i, 0)
                idx_1 = (i, 1)
            else:
                plt_i = 2*i
                idx_0 = (plt_i)
                idx_1 = (plt_i+1)
        
            sec_name = sec_traces[sec_i].split('V_')[1]
            
            spk_train_names = [key for key in spike_trains.keys() if sec_name in key]
            for i, spk_train_name in enumerate(spk_train_names):
                shift = 4*i  # gid%len(spk_train_names)
                axs[idx_0].vlines(spike_trains[spk_train_name], (avg_V_soma/2-1)+shift, (avg_V_soma/2+2)+shift, color='royalblue', zorder=12)

            axs[idx_0].plot(t, V_soma, label='soma', zorder=1)
            axs[idx_0].plot(t, list(simData[sec_traces[sec_i]]['cell_0']), label=f'{sec_name}', zorder=1)
            axs[idx_0].set_title(f'{sec_name} Membrane Potential')
            axs[idx_0].set_ylabel('Voltage (mV)')
            axs[idx_0].legend(loc='lower right')

            plot_spike_times = False

            I_ampa = list(simData[sec_traces[sec_i+1]]['cell_0'])
            I_nmda = list(simData[sec_traces[sec_i+2]]['cell_0'])
            if len(I_ampa) > 100:
                axs[idx_1].plot(t, I_ampa, label='AMPA', zorder=1, color='firebrick')
                if max(I_ampa) > 0:
                    axs[idx_1].vlines(t_spikes, [(max(I_ampa)/2)-0.00005], [(max(I_ampa)/2)+0.00005], 'k', zorder=12)
                    for i, spk_train_name in enumerate(spk_train_names):
                        shift = (4*i+1)*0.0001  # gid%len(spk_train_names)
                        axs[idx_1].vlines(spike_trains[spk_train_name], [((max(I_ampa)/4)-0.00005)+shift], [((max(I_ampa)/4)+0.00005)+shift], color='royalblue', zorder=13)

                    plot_spike_times = True
            if len(I_nmda) > 100:
                axs[idx_1].plot(t, I_nmda, label='NMDA', zorder=2, color=sec_syn_colors[sec_name]['E'], path_effects=[path_effects.SimpleLineShadow(offset=(0.5,-0.5)),path_effects.Normal()])  # 'salmon'
                
                if not plot_spike_times:
                    axs[idx_1].vlines(t_spikes, [(max(I_nmda)/2)-0.00005], [(max(I_nmda)/2)+0.00005], 'k', zorder=12)
                    for i, spk_train_name in enumerate(spk_train_names):
                        shift = (4*i+1)*0.0001  # gid%len(spk_train_names)
                        axs[idx_1].vlines(spike_trains[spk_train_name], [((max(I_nmda)/2)-0.00005)+shift], [((max(I_nmda)/2)+0.00005)+shift], color='royalblue', zorder=13)

                axs[idx_1].legend(loc='upper right')
                
            axs[idx_1].set_title(f'{sec_name} Synaptic Conductance')
            axs[idx_1].set_ylabel('g (uS)')
        
        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-secs.png'),bbox_inches='tight',dpi=300)

        ## separate synaptic conductances ##
        fig = plt.figure(figsize=(15,num_secs*5))
        outer_grid = fig.add_gridspec(num_secs,2)

        for i in range(num_secs):

            axL = fig.add_subplot(outer_grid[i,0])
            inner_grid = outer_grid[i,1].subgridspec(nrows=2, ncols=1)
            (axRa, axRb) = inner_grid.subplots()

            sec_i = 3*i

            if num_secs > 1:
                idx_0 = (i, 0)
                idx_1 = (i, 1)
            else:
                plt_i = 2*i
                idx_0 = (plt_i)
                idx_1 = (plt_i+1)
            
            sec_name = sec_traces[sec_i].split('V_')[1]

            for i, spk_train_name in enumerate(spk_train_names):
                shift = 4*i  # gid%len(spk_train_names)
                axL.vlines(spike_trains[spk_train_name], (avg_V_soma/2-1)+shift, (avg_V_soma/2+2)+shift, color='royalblue', zorder=12)

            axL.plot(t, simData[f'V_soma']['cell_0'], label='soma')
            axL.plot(t, simData[f'V_{sec_name}']['cell_0'], label=f'{sec_name}')
            axL.set_title(f'{sec_name} Membrane Potential')
            axL.set_ylabel('Voltage (mV)')
            axL.set_xlabel('Time (ms)')
            axL.legend(loc='upper right')

            I_ampa = list(simData[f'I_{sec_name}_ampa']['cell_0'])
            I_nmda = list(simData[f'I_{sec_name}_nmda']['cell_0'])


            if len(I_ampa) > 100:
                axRa.plot(t, simData[f'I_{sec_name}_ampa']['cell_0'], label='AMPA', color='firebrick', zorder=1)
                if max(I_ampa) > 0:
                    axRa.vlines(t_spikes, [(max(I_ampa)/2)-0.00005], [(max(I_ampa)/2)+0.00005], 'k', zorder=12)
                    for i, spk_train_name in enumerate(spk_train_names):
                        shift = (4*i+1)*0.0001  # gid%len(spk_train_names)
                        axRa.vlines(spike_trains[spk_train_name], [((max(I_ampa)/4)-0.00005)+shift], [((max(I_ampa)/4)+0.00005)+shift], color='royalblue', zorder=13)

                axRa.set_title(f'{sec_name} AMPA Synaptic Conductance')
                # axRa.set_xlim([800,1200])
                axRa.set_xticks([])
                axRa.set_ylabel('g (uS)')

            if len(I_nmda) > 100:
                axRb.plot(t, simData[f'I_{sec_name}_nmda']['cell_0'], label='NMDA', color='salmon', zorder=1)
                if max(I_nmda) > 0:
                    axRb.vlines(t_spikes, [(max(I_nmda)/2)-0.00005], [(max(I_nmda)/2)+0.00005], 'k', zorder=12)
                    for i, spk_train_name in enumerate(spk_train_names):
                        shift = (4*i+1)*0.0001  # gid%len(spk_train_names)
                        axRb.vlines(spike_trains[spk_train_name], [((max(I_ampa)/4)-0.00005)+shift], [((max(I_ampa)/4)+0.00005)+shift], color='royalblue', zorder=13)

                axRb.set_title(f'{sec_name} NMDA Synaptic Conductance')
                # axRb.set_xlim([800,1200])
                axRb.set_ylabel('g (uS)')
                axRb.set_xlabel('Time (ms)')

            # ax2.legend(loc='upper right')

        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-secs_synsSep.png'),bbox_inches='tight')  # ,dpi=200)

    return "Section membrane potential and synaptic conductances plotted!"


def plot_syns_traces(simData, syn_secs, sim_label, sim_dir, syn_colors):

    t = np.array(simData['t'])

    ### PLOT SYNAPSE MEMBRANE POTENTIALS ###
    fig, axs = plt.subplots(figsize=(8,5))
    for syn_sec in syn_secs:
        V_sec = np.array(simData[f'V_{syn_sec}']['cell_0'])
        
        axs.plot(t, V_sec, color=syn_colors['E'])

    axs.set_title(f'Synapse Location Membrane Potentials')
    axs.set_ylabel('Voltage (mV)')
    axs.set_xlabel('Time (ms)')

    fig.tight_layout()
    fig.savefig(os.path.join(sim_dir,f'{sim_label}-syns_pot.png'),bbox_inches='tight',dpi=300)


def plot_isoalted_syn_traces(simData, syn_secs, syns_type, num_syns, sim_label, sim_dir, output_dir, syn_colors):

    t = np.array(simData['t'])
    time_window_path = os.path.join(output_dir,'eap_time_windows.pkl')
    with open(time_window_path,'rb') as fp:
        time_windows = pickle.load(fp)
    
    plot_flag = True
    try:
        slice_start = time_windows[syns_type][num_syns][0]
        slice_end = time_windows[syns_type][num_syns][1]
    except KeyError:
        plot_flag = False


    ### PLOT SYNAPSE LOCATION MEMBRANE POTENTIALS OF ISOLATED EAP ###
    if plot_flag:
        t_window = t[slice_start:slice_end]

        fig, axs = plt.subplots(figsize=(8,5))
        for syn_sec in syn_secs:
            # TODO: support for empty V_syn_sec
            if f'V_{syn_sec}' in simData.keys():
                V_sec = np.array(simData[f'V_{syn_sec}']['cell_0'])
                V_sec_window = V_sec[slice_start:slice_end]
            
                axs.plot(t_window, V_sec_window, color=syn_colors['E'])

        axs.set_title(f'Synapse Location Membrane Potentials')
        axs.set_ylabel('Voltage (mV)')
        axs.set_xlabel('Time (ms)')

        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-isolated_syns_pot.png'),bbox_inches='tight',dpi=300)

