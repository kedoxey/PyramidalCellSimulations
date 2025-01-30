import os
import pickle

import numpy as np
import pandas as pd


#### Functions for processing data post simulation ###

def reformat_data(simData, rec_electrode, soma_name, delay, nmldb_id, sim_label, sim_dir):

    columns = ['Model_ID','t','vm','ve','x_bar','y_bar','z_bar',
                   'num_spikes','did_spike','first_spkt']
    waveforms_df = pd.DataFrame(columns=columns)

    try:
        did_spike = True
        num_spikes = len(np.array(simData['spkt']))
        spkt = np.array(simData['spkt'])[0]
    except IndexError:
        did_spike = False
        num_spikes = 0
        spkt = np.nan

    t = np.array(simData['t'])
    Vm = np.array(simData[f'V_{soma_name}']['cell_0'])
    Ve = np.array(simData['LFP'])

    temp_Ve = Ve.T

    max_t = 1040

    num_channels, _ = np.shape(rec_electrode)

    for chan_i in range(num_channels):

        x_bar = rec_electrode[chan_i][0]
        y_bar = rec_electrode[chan_i][1]
        z_bar = rec_electrode[chan_i][2]

        ve_i = temp_Ve[chan_i]
        
        this_interval = (t >= delay) & (t <= max_t)

        try:
            reduced_ve_i = ve_i[this_interval]
        except IndexError:
            this_interval = this_interval[:-1]
            reduced_ve_i = ve_i[this_interval]

        try:
            reduced_t = t[this_interval]
            reduced_Vm = Vm[this_interval]
        except IndexError:
            this_interval = (t >= delay) & (t <= max_t)
            reduced_t = t[this_interval]
            reduced_Vm = Vm[this_interval]

        df = pd.DataFrame(columns=columns)

        df['Model_ID'] = [nmldb_id]
        df['t'] = [reduced_t]
        df['vm'] = [reduced_Vm]
        df['ve'] = [reduced_ve_i]
        df['x_bar'] = [x_bar]
        df['y_bar'] = [y_bar]
        df['z_bar'] = [z_bar]
        df['first_spkt'] = [spkt]
        df['did_spike'] = [did_spike]
        df['num_spikes'] = [num_spikes]

        join_frames = [waveforms_df, df]
        waveforms_df = pd.concat(join_frames, ignore_index=True)

    file_name = f'{sim_label}-simulated_eaps.pkl'
    file_path = os.path.join(sim_dir, file_name)

    waveforms_df.to_pickle(file_path, protocol=3)


def save_eap_time(syns_type, num_syns, slice_start, slice_end, t_spike, output_dir):

    file_name = 'eap_time_windows.pkl'
    file_path = os.path.join(output_dir,'time_windows', file_name)

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

def save_firing_rate(simData, soma_name, sim_dur, syns_type, num_syns, output_dir):
    
    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])
    t = np.array(simData['t'])
    # t_adapt = t[np.where(t > stim_delay*3)]
    t_spikes = t[np.where(V_soma>0)]

    spkt = np.array(simData['spkt'])
    spkid = np.array(simData['spkid'])

    # TODO: adapting threshold
    adapt_delay = 200
    
    if len(t_spikes) > 0:
        # t_spikes = t_spikes[np.where(t_spikes > stim_delay*4)]   # only include spikes after initial burst and once cell has adapted
        # firing_rate = len(t_spikes) / ((stim_dur - stim_delay*4) / 1000)  # calculate firing rate based on excluding time for adaptation and scale to Hz (spikes/s)
        # isis = [t[i+1] - t[i] for i in range(len(t_spikes))]
        # f_isi = 1/isis[0]
        # f_ss = 1/np.mean(isis[-3:])

        spike_times = spkt[np.where(spkid == 0)]
        spike_times = spike_times[np.where(spike_times > adapt_delay)]
        num_spikes = len(spike_times)
        num_isi = num_spikes - 1
        # if num_isi > 0:
        #     msf = num_isi / (spike_times[-1] - spike_times[0]) * 1000
        firing_rate = num_spikes / ((sim_dur - adapt_delay) / 1000)
    else:
        firing_rate = 0
    
    file_name = 'firing_rates.pkl'
    file_path = os.path.join(output_dir, 'firing_rates', file_name)
    
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
                

def get_isolated_time_window(syns_type, num_syns, output_dir):

    time_window_path = os.path.join(output_dir,'time_windows', 'eap_time_windows.pkl')
    with open(time_window_path,'rb') as fp:
        time_windows = pickle.load(fp)

    slice_groups = [syns_type, 'soma']
    
    for slice_group in slice_groups:
        try:
            slice_start = time_windows[slice_group][num_syns][0]
            slice_end = time_windows[slice_group][num_syns][1]
            t_spike = time_windows[slice_group][num_syns][2]
            plot_flag = True
            break
        except KeyError:
            slice_start = None
            slice_end = None
            t_spike = None
            plot_flag = False

    return plot_flag, slice_start, slice_end, t_spike

