import os
import pickle
import re

import matplotlib as mpl
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import utils.process
# from utils.process import save_eap_time, get_isolated_time_window, find_n_closest_probes, find_n_farthest_probes, get_probe_max_channel, preprocess_probe_data


### Functions for plotting simulation data ###

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


def get_syn_sec_colors(cell, use_colormaps, colormaps, synColors):
    secSynCount = 0
    secSynColors = {}
    for secName, sec in cell.secs.items():
        for synMech in sec['synMechs']:
            if use_colormaps:
                if 'GABA' in synMech['label']:
                    secSynColors[secName]['I'] = colormaps[1][secSynCount]
                else:
                    secSynColors[secName]['E'] = colormaps[0][secSynCount]
            else:
                secSynColors[secName] = {'E': synColors['E'],
                                         'I': synColors['I']}

        if bool(sec['synMechs']):
            secSynCount += 1

    return secSynColors


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
    # axs.set_xlim([0,1500])
    axs.set_title('Presynaptic Spike Trains')

    fig.savefig(os.path.join(sim_dir,f'{sim_label}-presyn_spike_trains.png'),bbox_inches='tight',dpi=300)

    return spike_trains


def plot_soma(simData, soma_name, sim_label, sim_dir):

    t = np.array(simData['t'])
    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])

    ### PLOT SOMA MEMBRANE POTENTIAL ###
    fig, axs = plt.subplots(figsize=(8,5))
    axs.plot(t, V_soma)
    axs.set_title(f'Soma Membrane Potential')
    axs.set_ylabel('Voltage (mV)')
    axs.set_xlabel('Time (ms)')
    # axs.vlines(200,-80,20,color='k')

    fig.tight_layout()
    fig.savefig(os.path.join(sim_dir,f'{sim_label}-soma_pot.png'),bbox_inches='tight',dpi=300)

    ### PLOT SOMA AND APICAL DISTAL COMPARTMENT MEMBRANE POTENTIALS FOR BACKPROPAGATION ###
    if 'V_apic_32' in simData.keys():
        fig, axs = plt.subplots(figsize=(8,5))
        axs.plot(t, V_soma, label='soma')
        axs.plot(t, np.array(simData['V_apic_32']['cell_0']), label='apic_32')
        axs.legend(loc='upper right')
        axs.set_title(f'Membrane Potentials')
        axs.set_ylabel('Voltage (mV)')
        axs.set_xlabel('Time (ms)')

        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-backpropagation.png'),bbox_inches='tight',dpi=300)


def plot_secs(simData, soma_name, spike_trains, sim_label, sim_dir, sec_syn_colors):

    sec_traces = list(simData.keys())
    [sec_traces.remove(key) for key in ['spkt', 'spkid', 't', 'V_soma', 'avgRate', '__dict__'] if key in sec_traces]
    num_secs = len(sec_traces)//3

    t = np.array(simData['t'])
    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])
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
            # for i, spk_train_name in enumerate(spk_train_names):
            #     shift = 4*i  # gid%len(spk_train_names)
            #     axs[idx_0].vlines(spike_trains[spk_train_name], (avg_V_soma/2-1)+shift, (avg_V_soma/2+2)+shift, color='royalblue', zorder=12)

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
                    # for i, spk_train_name in enumerate(spk_train_names):
                    #     shift = (4*i+1)*0.0001  # gid%len(spk_train_names)
                    #     axs[idx_1].vlines(spike_trains[spk_train_name], [((max(I_ampa)/4)-0.00005)+shift], [((max(I_ampa)/4)+0.00005)+shift], color='royalblue', zorder=13)

                    plot_spike_times = True
            if len(I_nmda) > 100:
                axs[idx_1].plot(t, I_nmda, label='NMDA', zorder=2, color='salmon', path_effects=[path_effects.SimpleLineShadow(offset=(0.5,-0.5)),path_effects.Normal()])  # 'salmon'
                
                if not plot_spike_times:
                    axs[idx_1].vlines(t_spikes, [(max(I_nmda)/2)-0.00005], [(max(I_nmda)/2)+0.00005], 'k', zorder=12)
                    # for i, spk_train_name in enumerate(spk_train_names):
                    #     shift = (4*i+1)*0.0001  # gid%len(spk_train_names)
                    #     axs[idx_1].vlines(spike_trains[spk_train_name], [((max(I_nmda)/2)-0.00005)+shift], [((max(I_nmda)/2)+0.00005)+shift], color='royalblue', zorder=13)

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

            axL.plot(t, simData[f'V_{soma_name}']['cell_0'], label='soma')
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


def plot_isolated_syn_traces(simData, soma_name, syn_secs, syns_type, num_syns, sim_label, sim_dir, output_dir, syn_colors):

    t = np.array(simData['t'])
    time_window_path = os.path.join(output_dir,'time_windows', 'eap_time_windows.pkl')
    with open(time_window_path,'rb') as fp:
        time_windows = pickle.load(fp)

    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])
    t_spikes = t[np.where(V_soma>10)]

    if len(t_spikes) > 0:
        slice_group = syns_type
    else:
        slice_group = 'soma' if 'distal' in syns_type else syns_type
    
    plot_flag = True
    try:
        slice_start = time_windows[slice_group][num_syns][0]
        slice_end = time_windows[slice_group][num_syns][1]
    except KeyError:
        plot_flag = False


    ### PLOT SYNAPSE LOCATION MEMBRANE POTENTIALS OF ISOLATED EAP ###
    if plot_flag:
        t_window = t[slice_start:slice_end]

        fig, axs = plt.subplots(figsize=(8,5))
        for syn_sec in syn_secs:
            if f'V_{syn_sec}' in simData.keys():
                V_sec = np.array(simData[f'V_{syn_sec}']['cell_0'])
                V_sec_window = V_sec[slice_start:slice_end]
            
                axs.plot(t_window, V_sec_window, color=syn_colors['E'])

        axs.set_title(f'Synapse Location Membrane Potentials')
        axs.set_ylabel('Voltage (mV)')
        axs.set_xlabel('Time (ms)')

        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-isolated_syns_pot.png'),bbox_inches='tight',dpi=300)

def plot_isolated_soma_pot(simData, soma_name, syns_type, num_syns, sim_label, sim_dir, output_dir):

    t = np.array(simData['t'])
    time_window_path = os.path.join(output_dir, 'time_windows', 'eap_time_windows.pkl')
    with open(time_window_path,'rb') as fp:
        time_windows = pickle.load(fp)

    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])
    # t_spikes = t[np.where(V_soma>10)]

    plot_flag, slice_start, slice_end, t_spike = utils.process.get_isolated_time_window(syns_type, num_syns, output_dir)
    
    # slice_groups = [syns_type, 'soma']
    
    # for slice_group in slice_groups:
    #     try:
    #         slice_start = time_windows[slice_group][num_syns][0]
    #         slice_end = time_windows[slice_group][num_syns][1]
    #         t_spike = time_windows[slice_group][num_syns][2]
    #         plot_flag = True
    #         break
    #     except KeyError:
    #         plot_flag = False

    ### PLOT SYNAPSE LOCATION MEMBRANE POTENTIALS OF ISOLATED EAP ###
    if plot_flag:
        t_window = t[slice_start:slice_end]

        fig, axs = plt.subplots(figsize=(8,5))
        V_soma_window = V_soma[slice_start:slice_end]
        axs.plot(t_window, V_soma_window, color='tab:blue')
        axs.set_title(f'Soma Membrane Potential')
        axs.set_ylabel('Voltage (mV)')
        axs.set_xlabel('Time (ms)')

        xticks = [(int(t_spike.round(0))-2)+2*i for i in range(4)]
        axs.set_xticks(xticks)

        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-isolated_soma_pot.png'),bbox_inches='tight',dpi=300)


def plot_isolated_traces(simData, soma_name, syn_secs, syns_type, num_syns, sim_label, sim_dir, output_dir, syn_colors):

    t = np.array(simData['t'])
    time_window_path = os.path.join(output_dir, 'time_windows', 'eap_time_windows.pkl')
    with open(time_window_path,'rb') as fp:
        time_windows = pickle.load(fp)

    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])
    # t_spikes = t[np.where(V_soma>10)]

    # slice_groups = [syns_type, 'soma']
    
    # for slice_group in slice_groups:
    #     try:
    #         slice_start = time_windows[slice_group][num_syns][0]
    #         slice_end = time_windows[slice_group][num_syns][1]
    #         t_spike = time_windows[slice_group][num_syns][2]
    #         plot_flag = True
    #         break
    #     except KeyError:
    #         plot_flag = False

    plot_flag, slice_start, slice_end, t_spike = utils.process.get_isolated_time_window(syns_type, num_syns, output_dir)

    ### PLOT SYNAPSE LOCATION MEMBRANE POTENTIALS OF ISOLATED EAP ###
    if plot_flag:
        t_window = t[slice_start:slice_end]

        fig, axs = plt.subplots(1, 2, figsize=(5,4))
        axs = axs.ravel()

        V_soma_window = V_soma[slice_start:slice_end]

        axs[0].plot(t_window, V_soma_window, color='tab:blue')
        for syn_sec in syn_secs:
            if f'V_{syn_sec}' in simData.keys():
                V_sec = np.array(simData[f'V_{syn_sec}']['cell_0'])
                V_sec_window = V_sec[slice_start:slice_end]
            
                axs[1].plot(t_window, V_sec_window, color=syn_colors['E'])

        axs[0].set_title(f'Soma')
        axs[0].set_ylabel('Voltage (mV)')
        axs[1].set_title('Synapse Location')

        ylims = [-85,25]
        xticks = [(int(t_spike.round(0))-2)+2*i for i in range(4)]
        for ax in axs:
            ax.set_ylim(ylims)
            ax.set_xticks(xticks)
            ax.set_xlabel('Time (ms)')

        fig.tight_layout()
        fig.savefig(os.path.join(sim_dir,f'{sim_label}-isolated_traces.png'),bbox_inches='tight',dpi=300)


def plot_isolated_LFP(simData, soma_name, syns_type, num_syns, sim_label, sim_dir, output_dir):
    t = np.array(simData['t'])
    dt = t[1] - t[0]

    V_soma = np.array(simData[f'V_{soma_name}']['cell_0'])
    t_spikes = t[np.where(V_soma>10)]
    # t_bound = 600 if len(np.where(t_spikes>600)[0])>0 else 500
    
    # t_bound = 500 if 'distal' in syns_type else 600

    plot_flag = True
    plot_spike = True

    if len(t_spikes) > 0:
        t_bound = int(100*np.floor(t_spikes[-1]/100))
        t_spike = t_spikes[np.where(t_spikes>t_bound)[0][0]]
        slice_start = int((t_spike - 2.25)/dt)
        slice_end = int((t_spike + 4.5)/dt)

        utils.process.save_eap_time(syns_type, num_syns, slice_start, slice_end, t_spike, output_dir)
    else:
        plot_flag, slice_start, slice_end, t_spike = utils.process.get_isolated_time_window(syns_type, num_syns, output_dir)


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
            if plot_spike:
                axs[ax_i].axvline(t_spike,alpha=0.4,color='k',zorder=1,linestyle='--')
            xticks = [(int(t_spike.round(0))-2)+2*i for i in range(4)]
            axs[ax_i].set_xticks(xticks)
            axs[ax_i].set_title(i)

        axs[-1].set_xlabel('Time (ms)')
        axs[-2].set_xlabel('Time (ms)')

        fig.suptitle('Normalized Amplitude')
        fig.tight_layout()

        fig.savefig(os.path.join(sim_dir,f'{sim_label}-isolated_LFP.png'),bbox_inches='tight',dpi=300)


def plot_eap_kernel(type_probes, sim_dir, sim_label):

    data_df, probes = utils.process.load_eap_probe_data(sim_dir, sim_label)
    num_probes, num_channels, _ = np.shape(probes)

    closest_probe = utils.process.find_n_closest_probes(probes, 1, max_channel_i=num_channels//2)
    closest_probe_df = data_df[data_df.probe_num.isin(closest_probe)]
    
    max_channel = utils.process.get_probe_max_channel(closest_probe_df)
    include_channels = [i for i in range(max_channel-15, max_channel+16)]

    if 'close' in type_probes:
        kernel_probes = utils.process.find_n_closest_probes(probes, num_probes//2, max_channel_i=max_channel)
    else:
        kernel_probes = utils.process.find_n_farthest_probes(probes, num_probes//2, max_channel_i=max_channel)

    kernel_probes_df = data_df[data_df.probe_num.isin(kernel_probes)]

    prep_probes_df, probe_amps = utils.process.preprocess_probes(max_channel, include_channels, kernel_probes, kernel_probes_df)
    avg_waves = utils.process.get_avgerage_waves(include_channels, prep_probes_df)

    max_channel_loc = np.argwhere(include_channels == max_channel)[0][0]
    max_avg_wave = avg_waves[max_channel_loc]

    ### Plot EAP kernel and channel amplitudes ###
    icefire_cmap = sns.color_palette("icefire_r", as_cmap=True)

    fig, axs = plt.subplots(1,2, figsize=(12,8), width_ratios=[3,1], layout='constrained')
    axs.ravel()

    im = axs[0].imshow(avg_waves, interpolation='nearest', aspect='auto', cmap=icefire_cmap, zorder=0)

    max_chan_i = np.argwhere(include_channels == max_channel)[0][0]
    max_avg_wave = avg_waves[max_chan_i]

    axs[0].plot((-max_avg_wave*29)+4, color='white', zorder=12)

    max_amp = 0
    for probe_amp in probe_amps:
        temp_max = np.max(probe_amp)
        max_amp = np.max([temp_max, max_amp])
        axs[1].plot(probe_amp, include_channels, color='grey', alpha=0.4)

    axs[0].set_yticks([(len(include_channels)-1)*i for i in [0, 0.25, 0.5, 0.75, 1]])
    axs[0].set_yticklabels([(max_channel-include_channels[0])*i for i in [10, 5, 0, -50, -10]])
    axs[0].set_ylabel(r'Distance ($\mu$m)')
    axs[0].set_xticks([])

    axs[1].set_yticks([])
    axs[1].set_ylim(include_channels[0], include_channels[-1])
    axs[1].set_xticks([10*round(max_amp/10)*i for i in [0, 0.5, 1]])
    axs[1].set_xlabel(r'Amplitude ($\mu$V)')

    cbar = fig.colorbar(im, ax=axs[0], location='bottom', label='Ve (normalized)', pad=-0.04)

    fig.suptitle(f'{type_probes.capitalize()} (n = 50)');
    fig.savefig(os.path.join(sim_dir, f'{sim_label}-EAP_kernel-{type_probes}_probes.png'), dpi=300)

