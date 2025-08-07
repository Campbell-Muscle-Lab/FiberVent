# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 12:43:23 2023

@author: ken
"""

import os
import json
import copy
import shutil

from pathlib import Path

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


from ..batch import batch

def util_Frank_Starling(json_setup_file_string):
    """ Evaluates isovolumic cardiac cycle at different volumes """

    # Load the setup file
    with open(json_setup_file_string, 'r') as f:
        json_data = json.load(f)
        MyoVent_test = json_data['MyoVent_test']
    
    # Load the Frank_Starling structure
    fs = MyoVent_test['Frank_Starling']

    # Get the model and options file, setting the base directory as we go
    if not ('relative_to' in fs):
        base_dir = ''
        model_file_string = fs['model_file']
        options_file_string = fs['options_file']
    elif (fs['relative_to'] == 'this_file'):
        base_dir = Path(json_setup_file_string).parent.absolute()
        model_file_string = os.path.join(base_dir, fs['model_file'])
        options_file_string = os.path.join(base_dir, fs['options_file'])
    else:
        base_dir = fs['relative_to']
        model_file_string = os.path.join(base_dir, fs['model_file'])
        options_file_string = os.path.join(base_dir, fs['options_file'])

    # Load the base model
    with open(model_file_string, 'r') as f:
        base_model = json.load(f)
        
    # Load the base options
    with open(options_file_string, 'r') as f:
        base_options = json.load(f)
        
    # Try to clean out the sim_data folder
    try:
        print('Trying to remove %s' % os.path.join(base_dir, fs['sim_folder']))
        shutil.rmtree(os.path.join(base_dir, fs['sim_folder']),
                      ignore_errors = True)
    except OSError as e:
        print("Error: %s : %s" % (os.path.join(base_dir, fs['sim_folder']),
                                  e.strerror))
        
    # Set the sim data directories
    sim_input_dir = os.path.join(base_dir,
                                 fs['sim_folder'],
                                 'sim_input')
    if not os.path.isdir(sim_input_dir):
        os.makedirs(sim_input_dir)
    
    sim_output_dir = os.path.join(base_dir,
                                 fs['sim_folder'],
                                 'sim_output')
    if not os.path.isdir(sim_output_dir):
        os.makedirs(sim_output_dir)
        
    # Create a dictionary for the batch
    b = dict()
    
    # Copy across the exe structure
    b['MyoVent_batch'] = dict()
    b['MyoVent_batch']['MyoVentCpp_exe'] = MyoVent_test['MyoVentCpp_exe']
        
    # Cycle through the venous compliance factors, adding jobs into an array
    jobs = [];    
    for i, vcf in enumerate(fs['venous_compliance_factors']):
        
        print(vcf)
        
        new_model = copy.deepcopy(base_model)
        
        # Set all the resistances except the inflow very high
        n_compartments = len(new_model['circulation']['compartments']['resistance'])
        for j in range(n_compartments):
            if not (j==0):
                new_model['circulation']['compartments']['resistance'][j] = 1e5
        
        # Vary the venous compliance
        new_model['circulation']['compartments']['compliance'][-1] = \
            vcf * new_model['circulation']['compartments']['compliance'][-1]
        
        # Write the new model
        new_model_file_string = os.path.join(sim_input_dir,
                                             ('model_%i.json' % i))
        
        with open(new_model_file_string, 'w') as f:
            json.dump(new_model, f, indent=4)
            
        # Copy and write the options
        new_options = copy.deepcopy(base_options)
        new_options_file_string = os.path.join(sim_input_dir,
                                               ('options_%i.json' % i))
        
        with open(new_options_file_string, 'w') as f:
            json.dump(new_options, f, indent=4)
        
        # Copy the protocol
        prot = dict()
        prot['protocol'] = fs['protocol']
        
        protocol_file_string = os.path.join(sim_input_dir,
                                            ('protocol_%i.json' % i))

        with open(protocol_file_string, 'w') as f:
            json.dump(prot, f, indent=4)
            
        # Create the results file string
        results_file_string = os.path.join(sim_output_dir,
                                           ('results_%i.txt' % i))
            
        # Add in the job
        j = dict()
        j['model_file'] = new_model_file_string
        j['options_file'] = new_options_file_string
        j['protocol_file'] = protocol_file_string
        j['results_file'] = results_file_string
        
        jobs.append(j)
        
    # Add into the batch
    b['MyoVent_batch']['job'] = jobs
            
    # Write the batch
    batch_file_string = os.path.join(base_dir, fs['sim_folder'], 'batch.json')
    with open(batch_file_string, 'w') as f:
        json.dump(b, f, indent=4)
        
    # Now run the batch
    batch.run_batch(batch_file_string)
    
    # Now create the figure
    create_Frank_Starling_figure(b)
        
def create_Frank_Starling_figure(batch_data):
    """ Makes a figure from the batch data """
    
    thin_color_set = ['r', 'g']
    thick_color_set = ['r', 'g', 'b', 'm']
    
    # Create a figure
    fig = plt.figure(constrained_layout = False)
    fig.set_size_inches([14, 10])
    spec = gridspec.GridSpec(nrows = 6, ncols = 2, figure = fig,
                             wspace = 2, hspace = 1)
    
    ax=[]
    vent_sp = []
    vent_dp = []
    hsl = []
    stress_total = []
    stress_cb = []
    stress_int_pas = []
    stress_ext_pas = []
    
    
    for (i,j) in enumerate(batch_data['MyoVent_batch']['job']):
        
        results_file_string = j['results_file']
        d = pd.read_csv(results_file_string, sep='\t')
        
        # print(d.columns.to_list())
        
        if (i==0):
            ax.append(fig.add_subplot(spec[0,0]))
            ax.append(fig.add_subplot(spec[1,0]))
            ax.append(fig.add_subplot(spec[2,0]))
            ax.append(fig.add_subplot(spec[3,0]))
            ax.append(fig.add_subplot(spec[4,0]))
            ax.append(fig.add_subplot(spec[5,0]))
        
        ax[0].plot(d['time'], d['pressure_0'])
        ax[1].plot(d['time'], d['hs_length'])
        ax[2].plot(d['time'], d['hs_stress'])
        ax[3].plot(d['time'], d['memb_Ca_cytosol'])
        
        for k, str in enumerate(['myof_a_off', 'myof_a_on']):
            ax[4].plot(d['time'], d[str], color = thin_color_set[k])
            
        for k, str in enumerate(['myof_m_pop_0', 'myof_m_pop_1',
                                 'myof_m_pop_2', 'myof_m_pop_3']):
            if (str in d.columns):
                if (i==0):
                    lab = str
                else:
                    lab = None
                ax[5].plot(d['time'], d[str], color = thick_color_set[k],
                           label=lab)
        # Hold data
        hsl.append(d['hs_length'].max())
        vent_sp.append(d['pressure_0'].max())
        vent_dp.append(d.loc[171]['pressure_0'])
        stress_total.append(d['hs_stress'].max())
        stress_cb.append(d['myof_stress_cb'].max())
        stress_int_pas.append(d['myof_stress_int_pas'].max())
        stress_ext_pas.append(d['myof_stress_ext_pas'].max())


    ax[0].set_ylabel('Ventricular\npressure')

    ax[1].set_ylabel('Half-sarcomere\nlength\n(nm)')
    
    ax[2].set_ylabel('Stress')

    ax[3].set_ylabel('Ca\nconcentration\n(M)')

    ax[4].set_ylabel('Thin\npopulations')
        
    ax[5].legend(bbox_to_anchor=(1.05, 1))
    ax[5].set_ylabel('Thick\npopulations')
    ax[5].set_xlabel('Time (s)')

    # Plot values
    ax.append(fig.add_subplot(spec[0,1]))
    ax[6].plot(hsl, vent_sp, 'o-')
    ax[6].set_ylabel('Isovolumic\nVentricular\npressure\n(mmHg)')
    ax[6].set_ylim([0, np.amax(vent_sp)])
    
    ax.append(fig.add_subplot(spec[1,1]))
    ax[7].plot(hsl, stress_total, 'o-')
    ax[7].set_ylabel('Total stress')
    ax[7].set_ylim([0, np.amax(stress_total)])
    
    ax.append(fig.add_subplot(spec[2,1]))
    ax[8].plot(hsl, stress_cb, 'o-')
    ax[8].set_ylabel('cb stress')
    ax[8].set_ylim([0, np.amax(stress_total)])

    ax.append(fig.add_subplot(spec[3,1]))
    ax[9].plot(hsl, stress_int_pas, 'o-', label='Int')
    ax[9].plot(hsl, stress_ext_pas, 'o-', label='Ext')
    ax[9].set_ylabel('Pas stress')
    ax[9].legend(bbox_to_anchor=(1.05, 1))
    
    ax.append(fig.add_subplot(spec[4,1]))
    ax[10].plot(hsl, vent_dp, 'o-')
    ax[10].plot([hsl[0], hsl[-1]], [8, 8], 'k:')
    ax[10].set_ylabel('Diastolic\npressure')
    ax[10].set_xlabel('Half-sarcomere len\gth (nm)')

    fig.savefig('c:/temp/sl.png')