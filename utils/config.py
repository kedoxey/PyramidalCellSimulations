import os
import yaml
import shutil

cwd = os.getcwd()


### Functions for reading and writing configuration files and configuring file system ###

### Configuration files ###
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
