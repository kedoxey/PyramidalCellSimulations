import os
import pickle

import numpy as np
import pandas as pd
import colorednoise as cn

from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline


#### Functions for processing data post simulation ###

def reformat_data(simData, rec_electrode, soma_name, delay, sim_dur, nmldb_id, sim_label, sim_dir, r_thetas=None):

    columns = ['Model_ID','t','vm','ve','x_bar','y_bar','z_bar',
                   'num_spikes','did_spike','first_spkt']
    if r_thetas:
        columns.extend(['r', 'theta'])
    waveforms_df = pd.DataFrame(columns=columns)

    spkt = np.array(simData['spkt'])
    spkid = np.array(simData['spkid'])

    try:
        did_spike = True
        num_spikes = len(spkt[np.where(spkid == 0)])
        first_spike = spkt[np.where(spkid == 0)][0]
    except IndexError:
        did_spike = False
        num_spikes = 0
        first_spike = np.nan

    t = np.array(simData['t'])
    Vm = np.array(simData[f'V_{soma_name}']['cell_0'])
    Ve = np.array(simData['LFP'])

    temp_Ve = Ve.T

    max_t = sim_dur-5

    num_channels, _ = np.shape(rec_electrode)

    for chan_i in range(num_channels):

        x_bar = rec_electrode[chan_i][0]
        y_bar = rec_electrode[chan_i][1]
        z_bar = rec_electrode[chan_i][2]

        if r_thetas:
            r = r_thetas[chan_i][0]
            theta = r_thetas[chan_i][1]

        ve_i = temp_Ve[chan_i]
        
        this_interval = (t >= delay) & (t <= max_t)

        try:
            reduced_ve_i = ve_i[this_interval]
        except IndexError:
            time_diff = len(this_interval) - len(ve_i)
            if time_diff < 0:
                ve_i = ve_i[:time_diff]
            else:
                this_interval = this_interval[:-time_diff]
            # this_interval = this_interval[:-1]
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
        df['first_spkt'] = [first_spike]
        df['did_spike'] = [did_spike]
        df['num_spikes'] = [num_spikes]
        if r_thetas:
            df['r'] = [r]
            df['theta'] = [theta]

        join_frames = [waveforms_df, df]
        waveforms_df = pd.concat(join_frames, ignore_index=True)

    file_name = f'{sim_label}-simulated_eaps.pkl'
    file_path = os.path.join(sim_dir, file_name)

    waveforms_df.to_pickle(file_path, protocol=3)

    return waveforms_df


def save_detect_limit(data_df, model_dir):

    rs = np.unique(data_df.r.values)

    r_amps = {}
    within_limit = {}

    for r in rs:

        r_df = data_df[data_df.r == r]

        amps = []

        for i, row in r_df.iterrows():

            ve = row.ve
            amp = np.abs(np.max(ve) - np.min(ve))*1000
            amps.append(amp)

        avg_amp = np.average(amps)
        r_amps[r] = avg_amp
        if avg_amp > 20:
            within_limit[r] = avg_amp

    detect_limit = np.array(list(within_limit.keys())[-1])
    
    np.save(os.path.join(model_dir, 'detectability_limit'), detect_limit)


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


def find_n_closest_probes(probes, n_closest, max_channel_i = 25):
    
    num_probes, _, _ = np.shape(probes)
    distances = []
    
    for i in range(num_probes):
        probe = probes[i,:,:]
        x,y,z = probe[max_channel_i]
        dist = np.sqrt(x**2+y**2+z**2)
        distances.append(dist)
        
    print('Closest probe found at %s microns...'%np.min(distances))
    return np.argsort(distances)[:n_closest]  # np.argmin(distances)


def find_n_farthest_probes(probes, n_farthest, max_channel_i = 25):
    
    num_probes, _, _ = np.shape(probes)
    distances = []
    
    for i in range(num_probes):
        probe = probes[i,:,:]
        x,y,z = probe[max_channel_i]
        dist = np.sqrt(x**2+y**2+z**2)
        distances.append(dist)
        
    print('Farthest probe found at %s microns...'%np.max(distances))
    return np.argsort(distances)[-n_farthest:]  # np.argmax(distances)


def get_probe_max_channel(probe_df):

    amps = []

    for _, row in probe_df.iterrows():

        amps.append(np.max(row.ve) - np.min(row.ve))

    max_channel = np.argmax(amps)

    return max_channel


def load_eap_probe_data(sim_dir, sim_label):

    eap_file_path = os.path.join(sim_dir, f'{sim_label}-simulated_eaps.pkl')
    probe_file_path = os.path.join(sim_dir, f'{sim_label}-recording_probe_locs.pkl')

    data_df = pd.read_pickle(eap_file_path)
    probes = np.array(pd.read_pickle(probe_file_path))

    num_probes, num_channels, _ = np.shape(probes)

    probe_nums = []
    channel_nums = []

    for i, row in data_df.iterrows():
        probe_nums.append(int(i/num_channels))
        channel_nums.append(i%num_channels)

    data_df['probe_num'] = probe_nums
    data_df['channel_num'] = channel_nums

    return data_df, probes


def get_avgerage_waves(include_channels, prep_probes_df):

    avg_waves = []

    for channel_num in include_channels:

        channel_df = prep_probes_df[prep_probes_df.channel_num == channel_num]

        waves = []

        for _, row in channel_df.iterrows():

            ve = row.ve
            amp = np.max(ve) - np.min(ve)
            norm_ve = np.divide(ve, amp)
            wave = list(norm_ve)
            waves.append(wave)


        waves = np.array(waves)
        channel_avg = np.average(waves, 0)

        avg_waves.append(channel_avg)

    avg_waves = list(reversed(avg_waves))

    return avg_waves


def preprocess_probes(max_channel, include_channels, probes, probes_df):
    
    prep_probes_df = probes_df.copy()
    
    probe_amps = []

    for probe_num in probes:

        probe_df = probes_df[probes_df.probe_num == probe_num]

        prep_probe_df = preprocess_probe_data(probe_df, max_channel_i=max_channel, resample=512, conversion_factor=1e3, num_noisy_spikes=200)
        
        prep_probes_df.update(prep_probe_df)

        probe_amp = []

        for _, row in prep_probe_df.iterrows():

            if row.channel_num in include_channels:

                ve = row.ve
                amp = np.max(ve) - np.min(ve)
                probe_amp.append(amp)

        probe_amps.append(probe_amp)

    return prep_probes_df, probe_amps


def preprocess_probe_data(probe_df,max_channel_i=15,downsample=82,resample=128,conversion_factor=1e3, scaling_factor=10,num_noisy_spikes=100):
    
    preprocessed_probe_df = probe_df.copy()
    new_ves = []
    all_amps = []
    
    duration = 5.6
    peak_shift = 1.4 # duration/4
    

    t= probe_df.t.iloc[0]
    dt = t[1]-t[0]
    Fs = 1000./dt
    
    # approximate times for aligning waveforms    
    max_chan_ve = probe_df[probe_df['channel_num']==max_channel_i]['ve'].iloc[0]
    peak_time = t[np.argmin(max_chan_ve)]
    
    start_time = peak_time-peak_shift
    end_time = start_time+duration
    
    # approximate reduced timeseries condition
    t_cond = (t>=start_time)&(t<end_time)
    
    # find actual maximum channel
    temp_amps = []
    for _, row_df in probe_df.iterrows():
        wave = row_df['ve']
        wave = wave[t_cond]
        amp = np.max(wave)-np.min(wave)
        
        temp_amps.append(amp)
        
    # actual times for aligning waveforms
    temp_amps = np.array(temp_amps)
    max_channel_i = np.argmax(temp_amps)
    # print(max_channel_i)
    
    max_chan_ve = probe_df[probe_df['channel_num']==max_channel_i]['ve'].iloc[0]
    peak_time = t[np.argmin(max_chan_ve)]
    
    # readjust after finding real max
    start_time = peak_time-peak_shift
    end_time = start_time+duration

    # reduced timeseries
    t_cond = (t>=start_time)&(t<end_time)
    reduced_t = t[t_cond]
    
    T = len(reduced_t)
    
    # interpolation grid
    interp_t = np.linspace(0,duration,T)
    
    # downsample and upsample grids
    down_t = np.linspace(0,duration,downsample)
    if resample is False:
        re_t = down_t
    else:
        re_t = np.linspace(0,duration,resample)
    
    
    # iterate over channels
    for _, row_df in probe_df.iterrows():
        
        # extract signal 
        ve = row_df['ve']
        reduced_ve = ve[t_cond]
        
        # scale signal
        reduced_ve = conversion_factor*reduced_ve 
        
        # downsample signal
        cs = CubicSpline(interp_t,reduced_ve,axis=0)
        new_ve = cs(down_t)
        
        # noise details
        noise_scale = 10 # 15, 25
        noise_generator = cn.powerlaw_psd_gaussian
        p1 = 1
        p2 = (num_noisy_spikes,len(new_ve))
        args = [p1,p2]

        
        # add "single-trial" noise
        noise_signals = scaling_factor*noise_generator(*args)
        uncorr_noise_signal = np.mean(noise_signals,axis=0)

        synth_ve = np.add(new_ve,uncorr_noise_signal)
        
        
        # resample signal
        if resample is False:
            new_ve = synth_ve
        else:
            cs = CubicSpline(down_t,synth_ve,axis=0)
            new_ve = cs(re_t)
            
        # center channel
        median_ve = np.median(new_ve,axis=0)
        new_ve = np.subtract(new_ve,median_ve) 
            
        amp = np.max(new_ve)-np.min(new_ve)
        
        
        new_ves.append(new_ve)
        all_amps.append(amp)
    
    preprocessed_probe_df['ve'] = new_ves
    preprocessed_probe_df['t'] = [re_t for _ in range(len(new_ves))]
    preprocessed_probe_df['amplitude'] = all_amps
    
    return preprocessed_probe_df
        

