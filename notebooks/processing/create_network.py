import matplotlib
import os
import json
import model_helpers as mh
from pyneuroml import pynml
from pyneuroml.analysis import generate_current_vs_frequency_curve

nmldb_id = "NMLCL000073"
model_source = 'NeuroML'
model_name = "L5PC"

cwd = os.getcwd()
models_dir = os.path.join(cwd, 'models')
model_dir = os.path.join(models_dir, model_source, nmldb_id)

mh.download_from_nmldb(nmldb_id, model_source)

cell_fname = os.path.join(model_dir,model_name) + '.cell.nml'

cell_nml = pynml.read_neuroml2_file(cell_fname)

na_erev = 50

generate_current_vs_frequency_curve(cell_fname,
                                    cell_nml.id,
                                    start_amp_nA=-0.02,
                                    end_amp_nA=0.06,
                                    step_nA=0.01,
                                    pre_zero_pulse=20,
                                    post_zero_pulse=20,
                                    analysis_duration=100,
                                    temperature='34degC',
                                    plot_voltage_traces=True)

temp = 5

