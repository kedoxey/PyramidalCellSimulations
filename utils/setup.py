import os
import pickle

import numpy as np
import pandas as pd

from neuron import h

h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')


### Functions for analyzing morphology and modifying model specifications ###

def get_compartments(hoc_fname, cell, cell_name, soma_name, group_name='all'):

    secs = np.array(list(cell['secs'].keys()))
    basal_group = list(secs[np.where(np.char.find(secs,'dend')>=0)])
    apic_group = list(secs[np.where(np.char.find(secs,'apic')>=0)])
    all_group = list(secs)
    soma = [soma_name]

    if group_name == 'basal':
        return basal_group
    elif group_name == 'apical':
        return apic_group
    elif group_name == 'apical_distal':
        return get_secs_from_dist(hoc_fname, cell_name, soma_name, 0.5, 1, secs_lim='apic')
    elif group_name == 'apical_proximal':
        return get_secs_from_dist(hoc_fname, cell_name, soma_name, 0, 0.5, secs_lim='apic')
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


def get_total_soma_distance(cell, soma_name):
    
    # soma = cell.soma_0
    soma = getattr(cell, soma_name)

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
                    if isinstance(mech['gmax'], list):
                        for i, gmax in enumerate(mech['gmax']):
                            gmax *= toggle
                            mech['gmax'][i] = gmax
                    else:
                        mech['gmax'] *= toggle
                    mech.pop('__dict__', None)
    
    return cell_params


def get_secs_from_dist(filename, cell_name, soma_name, lb, ub=1, secs_lim='all'):
    cell = get_hoc_cell(filename, cell_name)
    # soma = cell.soma_0
    soma = getattr(cell, soma_name)

    total_distance, _ = get_total_soma_distance(cell, soma_name)
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


# TODO: implement for OB mitral, tufted, and granule cell models
def update_cell_params(cell_params, cell_name, params_path):

    params = pd.read_pickle(params_path)

    model_params = params[params['cell'] == cell_name]

    for col in model_params:
        if 'cell' in col:
            continue
        
        param = col
        value = model_params[col].values[0][0]
        secs = model_params[col].values[0][1]

        temp = 5

        # for sec_name, sec_params in cell_params['secs'].items():
        #     sec_type = sec_name.split('_')[0]

        #     if sec_type in secs:
                
    return cell_params


def define_electrode_geom(num_probes, total_channels, sampling_dist, sim_label, sim_dir, model_dir):

    if sampling_dist == 'fixed_normal':
        # fixed normal distribution
        r_pos = np.random.normal(loc=20, scale=5, size=num_probes)
    
        min_radius = 10
        r_pos = [r if r > min_radius else min_radius for r in r_pos]

    elif sampling_dist == 'adjusted_uniform':
        # adjusted with detectability limit
        detect_limit = np.load(os.path.join(model_dir, 'detectability_limit.npy'))

        r_pos = np.random.uniform(low=10.,high=detect_limit,size=num_probes)

    else:
        raise Exception('Sampling distribution incorrectly specified.')

    channel_lb = total_channels // 2 if total_channels > 1 else total_channels
    channel_ub = total_channels // 2 if total_channels > 1 else total_channels
    grid_spacing = 10

    # sample angle from uniform distribution
    theta_pos = np.random.uniform(-np.pi, np.pi, size=num_probes)

    all_xs = np.multiply(r_pos, np.cos(theta_pos))
    all_zs = np.multiply(r_pos, np.sin(theta_pos))

    all_ys = np.random.uniform(low=-0.5*grid_spacing,
                               high=0.5*grid_spacing,
                               size=num_probes)
    
    rec_probes = []
    probe_list = []

    lb = int(channel_lb*grid_spacing)
    ub = int(channel_ub*grid_spacing)

    for y_center, x_center, z_center in zip(all_ys, all_xs, all_zs):

        y_lower = np.arange(y_center - lb, y_center, grid_spacing)
        y_upper = np.arange(y_center, y_center+ub, grid_spacing)

        y_channels = list(y_lower) + list(y_upper)

        probe = [[x_center, yi, z_center] for yi in y_channels]
        probe_list.append(probe)
        
        rec_probes += probe
    
    rec_electrode = rec_probes

    # save probe arrangement
    file_name = f'{sim_label}-recording_probe_locs.pkl'
    file_path = os.path.join(sim_dir, file_name)
    
    with open(file_path, 'wb') as fp:
        pickle.dump(probe_list, fp, protocol=3)

    return rec_electrode


def define_detect_lim_geom(sim_label, sim_dir):

    r_pos = np.arange(10,125,5) # 23 r-values
    thetas = np.random.uniform(-np.pi,np.pi,size=10)

    all_xs = []
    all_zs = []
    all_r_thetas = []

    for ri in r_pos:

        xi = np.multiply(ri,np.cos(thetas))
        zi = np.multiply(ri,np.sin(thetas))
        
        r_thetas = [(ri, theta) for theta in thetas]

        all_xs+=list(xi)
        all_zs+=list(zi)
        all_r_thetas.extend(r_thetas)

    all_xs = np.array(all_xs)
    all_zs = np.array(all_zs)      
        
    all_ys = np.arange(-10,15,5) # can compute SNR/amplitude decay along the Y-axis

    rec_probes = []
    probe_list = []
    r_thetas_list = []

    for xi, zi, r_thetas in zip(all_xs,all_zs, all_r_thetas):
        probe = [[xi,yi,zi] for yi in all_ys]
        probe_list.append(probe)

        r_thetas_list.extend([r_thetas for yi in all_ys])

        rec_probes += probe

                
    rec_electrode = rec_probes

    # save probe arrangement
    file_name = f'{sim_label}-recording_probe_locs.pkl'
    file_path = os.path.join(sim_dir, file_name)
    
    with open(file_path, 'wb') as fp:
        pickle.dump(probe_list, fp, protocol=3)

    return rec_electrode, r_thetas_list
