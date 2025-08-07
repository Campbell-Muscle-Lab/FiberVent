# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 17:22:45 2023

@author: Campbell
"""

import os
import json
import shutil
import copy
import re

from pathlib import Path

import numpy as np

from ..batch import batch

def characterize_model(json_analysis_file_string, figures_only=False):
    """ Code takes a json struct that includes a model file and runs the
        analyses described in thta file """
        
    # Check for the analysis file
    try:
        with open(json_analysis_file_string, 'r') as f:
            json_data = json.load(f)
            anal_struct = json_data['FiberVent_setup']
    except:
        print('characterize_model() problem')
        print('Invalid file: %s' % json_analysis_file_string)
        exit(1)
        
    # If there is a manipulations section, use that to create the
    # appropriate models
    if ("manipulations" in anal_struct['model']):
        json_analysis_file_string = \
            generate_model_files(json_analysis_file_string)
        
    # Loop through the characterizations
    for ch in anal_struct['characterization']:
        
        if (ch['type'] == 'isovolumic'):
            deduce_isovolumic_properties(json_analysis_file_string,
                                         ch,
                                         figures_only)
            
        if (ch['type'] == 'freeform'):
            deduce_freeform_properties(json_analysis_file_string,
                                       ch,
                                       figures_only)

def generate_model_files(json_analysis_file_string):
    """ Clones base model with modifications to facilitate comparisons """
    
    # First load the file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        model_struct = json_data['FiberVent_setup']['model']
    
    # Deduce the base model file string
    base_dir = ''
    if ('relative_to' in model_struct):
        if (model_struct['relative_to'] == 'this_file'):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = model_struct['relative_to']
    base_model_file_string = os.path.join(base_dir,
                                          model_struct['manipulations']['base_model'])
    base_model_file_string = str(Path(base_model_file_string).absolute().resolve())
    

    # Now deduce where to put the adjusted model files
    # We can use the base_dir from above
    generated_dir = os.path.join(base_dir, model_struct['manipulations']['generated_folder'])
    generated_dir = str(Path(generated_dir).absolute().resolve())
        
    # Clean the generated dir
    try:
        print('Trying to remove %s' % generated_dir)
        shutil.rmtree(generated_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s' % (generated_dir, e.strerror))
        
    if not os.path.isdir(generated_dir):
        os.makedirs(generated_dir)
        
    # Now copy the sim_options file across
    # We can still use the base dir above
    orig_options_file_string = os.path.join(base_dir, model_struct['options_file'])
    temp, file_name = os.path.split(orig_options_file_string)
    new_options_file = os.path.join(generated_dir, file_name)
    
    shutil.copy(orig_options_file_string, new_options_file)

    # Load up the base model
    with open(base_model_file_string, 'r') as f:
        base_model = json.load(f)

    # Now work out how many adjustments you need to make
    adjustments = model_struct['manipulations']['adjustments']
    if ('multipliers' in adjustments[0]):
        no_of_models = len(adjustments[0]['multipliers'])
    else:
        no_of_models = 1
    
    generated_models = []
        
    # Loop through them
    for i in range(no_of_models):
        
        # Copy the base model
        adj_model = copy.deepcopy(base_model)
                
        for (j,a) in enumerate(adjustments):
            
            if (a['class'] == 'circulation'):
                
                # Check whether a['variable'] includes numbers
                if (bool(re.search(r'\d', a['variable']))):
                    # Likely entry into an array
                    temp = re.findall(r'\d+', a['variable'])
                    cpt_number = int(temp[0])
                    var_name = a['variable'].split('_')[0]
                    
                    y = np.asarray(adj_model['FiberVent']['circulation']['compartments'][var_name],
                                   dtype = np.float32)
                    
                    base_value = y[cpt_number - 1]
                                   
                    y[cpt_number - 1] = base_value * a['multipliers'][i]
                    
                    adj_model['FiberVent']['circulation']['compartments'][var_name] = y.tolist()
                    
                else:
                    # Should be a single entry
                    base_value = adj_model['FiberVent']['circulation'][a['variable']]

                    value = base_value * a['multipliers'][i]
                    
                    adj_model['FiberVent']['circulation'][a['variable']] = value
                    
            elif (a['class'] == 'aortic_valve'):
                base_value = adj_model['FiberVent']['circulation']['ventricle']['valves']['aortic'][a['variable']]
                
                value = base_value * a['multipliers'][i]
                
                adj_model['FiberVent']['circulation']['ventricle']['valves']['aortic'][a['variable']] = value
                
            elif (a['class'] == 'mitral_valve'):
                base_value = adj_model['FiberVent']['circulation']['ventricle']['valves']['mitral'][a['variable']]
                
                value = base_value * a['multipliers'][i]
                
                adj_model['FiberVent']['circulation']['ventricle']['valves']['mitral'][a['variable']] = value
            
            elif (a['class'] == 'growth_control'):
                gc = adj_model['FiberVent']['growth']['control'][a['control_number']-1]

                base_value = gc[a['variable']]
                
                value = base_value * a['multipliers'][i]
                
                adj_model['FiberVent']['growth']['control'][a['control_number']-1][a['variable']] = \
                    value
                    
            elif (a['class'] == 'mitochondria'):
                mito_dict = adj_model['FiberVent']['circulation']['ventricle']['myocardium']['mitochondria']

                base_value = mito_dict[a['variable']]
                
                value = base_value * a['multipliers'][i]
                
                adj_model['FiberVent']['circulation']['ventricle']['myocardium']['mitochondria'][a['variable']] = \
                    value
                    
            elif ((a['class'] == 'm_kinetics') or
                    (a['class'] == 'c_kinetics')):

                kinetics_structure = adj_model['FiberVent']['circulation']['ventricle']['myocardium']['contraction'] \
                    ['model']['muscle']['half_sarcomere'][a['class']][a['isotype']-1]['state'][a['state']-1]

                # Special case for kinetics
                if ('extension' in a):
                    base_value = a['extension']
                    
                    value = base_value * a['multipliers'][i]
                    
                    kinetics_structure['extension'] = value
                else:
                    # Transition parameters
                    y = np.asarray(kinetics_structure['transition'][a['transition']-1] \
                                   ['rate_parameters'], dtype = np.float32)
                    
                    base_value = y[a['parameter_number'] - 1]
                    value = base_value * a['multipliers'][i]
                        
                    y[a['parameter_number']-1] = value
                    
                    kinetics_structure['transition'][a['transition']-1]['rate_parameters'] = \
                        y.tolist()
                        
                # Insert back into model
                adj_model['FiberVent']['circulation']['ventricle']['myocardium']['contraction'] \
                    ['model']['muscle']['half_sarcomere'][a['class']][a['isotype']-1]['state'][a['state']-1] = \
                        kinetics_structure                  
                
            elif ((a['class'] == 'titin_parameters') or
                  (a['class'] == 'extracellular_parameters')):
                
                base_value = adj_model['FiberVent']['circulation']['ventricle']['myocardium']['contraction'] \
                    ['model']['muscle']['half_sarcomere'][a['class']][a['variable']]
                
                value = base_value * a['multipliers'][i]
                
                adj_model['FiberVent']['circulation']['ventricle']['myocardium']['contraction'] \
                    ['model']['muscle']['half_sarcomere'][a['class']][a['variable']] = value
                
            else:
                base_value = adj_model[a['class']][a['variable']]
                
                value = base_value * a['multipliers'][i]
                
                if (a['output_type'] == 'int'):
                    adj_model[a['class']][a['variable']] = int(value)
                    
                if (a['output_type'] == 'float'):
                    adj_model[a['class']][a['variable']] = float(value)
                    
                # Check for NaN
                if (np.isnan(value)):
                    adj_model[a['class']][a['variable']] = 'null'
    
        # Now generate the model file string
        model_file_string = 'model_%i.json' % (i+1)
        
        # We need the full path to write it to disk
        adj_model_file_string = os.path.join(generated_dir,
                                             model_file_string)

        with open(adj_model_file_string, 'w') as f:
            json.dump(adj_model, f, indent=4)

        # Append the model files
        if not (model_struct['relative_to'] == 'this_file'):
            model_file_string = adj_model_file_string
            
        generated_models.append(model_file_string)
        
    # Update the set up file
    
    # Add in the model files
    json_data['FiberVent_setup']['model']['model_files'] = generated_models
    
    # Delete the adjustments
    del(json_data['FiberVent_setup']['model']['manipulations'])
    
    # Generate a new setup file string
    generated_setup_file_string = os.path.join(generated_dir,
                                               'generated_setup.json')
    
    # Write to file
    with open(generated_setup_file_string, 'w') as f:
        json.dump(json_data, f, indent=4)
        
    # Return the new filename
    return (generated_setup_file_string)


def deduce_freeform_properties(json_analysis_file_string,
                               characterization_struct,
                               figures_only = False):
    """ Runs a sequence of freeform simulations """
       
    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        anal_struct = json_data['FiberVent_setup']
    
    # Create a batch
    FiberV_batch = dict()
    
    # Add in the MyoVentCpp stuff
    FiberV_batch['FiberVentCpp_exe'] = anal_struct['FiberVentCpp_exe']
    
    # Deduce the base model files
    if ('relative_to' in anal_struct['model']):
        if (anal_struct['model']['relative_to'] == 'this_file'):
            base_dir = str(Path(json_analysis_file_string).parent)
        else:
            base_dir = anal_struct['model']['relative_to']
    else:
        base_dir = ''
    
    base_model_files = []
    for mf in anal_struct['model']['model_files']:
        base_model_files.append(os.path.join(base_dir, mf))
        
    # Deduce the base options file
    # Can use the base dir from above
    base_options_file = os.path.join(base_dir,
                                    anal_struct['model']['options_file'])
    
    
    # Now do stuff based on the characterization
    
    # First, work out the base dir
    
    if ('relative_to' in characterization_struct):
        if (characterization_struct['relative_to'] == 'this_file'):
            base_dir = str(Path(json_analysis_file_string).parent)
        else:
            base_dir = characterization_struct['relative_to']
    else:
        base_dir = ''
        
    # Now make the output dir
    top_data_dir = os.path.join(base_dir, characterization_struct['sim_folder'])
    
    # Now clean the output_dir and make it
    if not (figures_only):
        try:
            print('Trying to remove %s' % top_data_dir)
            shutil.rmtree(top_data_dir, ignore_errors = True)
        except OSError as e:
            print('Error: %s : %s' % (top_data_dir, e.strerror))
            
        if not os.path.isdir(top_data_dir):
            os.makedirs(top_data_dir)

    # Set some counters and arrays for the simulations
    sim_counter = 1
    job = []
                
    # Now cycle through the model files
    for (model_counter, mf) in enumerate(base_model_files):
        
        # Load the base model
        with open(mf, 'r') as f:
            base_model = json.load(f)
            
        # Load the options file
        with open(base_options_file, 'r') as f:
            base_options = json.load(f)
        
        # Now cycle through volume_factors
        for cond_ind in range(characterization_struct['no_of_conditions']):

            # Create a dictionary for the job
            j = dict()

            # Make the input folder
            sim_input_folder = os.path.join(top_data_dir,
                                            'sim_input',
                                            ('%i' % sim_counter))
            
            if not (os.path.isdir(sim_input_folder)):
                os.makedirs(sim_input_folder)

            # Copy the volume and make some adjustments
            new_model = copy.deepcopy(base_model)
                
            # If m_n is listed, set that
            if ('m_n' in characterization_struct):
                new_model['FiberVent']['circulation']['ventricle']['myocardium']['contraction'] \
                    ['model']['muscle']['half_sarcomere']['thick_structure']['m_n'] = \
                        characterization_struct['m_n']
                
            # Now write the model to file
            new_model_file_string = os.path.join(sim_input_folder,
                                                 'model.json')
            
            with open(new_model_file_string, 'w') as f:
                json.dump(new_model, f, indent=4)
                
            # Copy the options
            new_options = copy.deepcopy(base_options)
            
            new_options_file_string = os.path.join(sim_input_folder,
                                                   'options.json')
            
            with open(new_options_file_string, 'w') as f:
                json.dump(new_options, f, indent=4)
                
            # Create a protocol
            prot = dict()
            prot['protocol'] = dict()
            prot['protocol']['time_step_s'] = characterization_struct['time_step_s']
            prot['protocol']['no_of_time_steps'] = round(
                characterization_struct['sim_duration_s'] / characterization_struct['time_step_s'])
            
            # Add in an activation
            if ('activation' in characterization_struct):
                prot['activation'] = []
                for i in range(len(characterization_struct['activation'])):
                    act = copy.deepcopy(characterization_struct['activation'][i])
                    if (np.any(np.asarray(act['simulation']) == (cond_ind+1))):
                        del act['simulation']
                        prot['activation'].append(act)                
            
            # Add in the perturbation
            if ('perturbation' in characterization_struct):
                prot['perturbation'] = []
                for i in range(len(characterization_struct['perturbation'])):
                    pert = copy.deepcopy(characterization_struct['perturbation'][i])
                    if (np.any(np.asarray(pert['simulation']) == (cond_ind+1))):
                        del pert['simulation']
                        prot['perturbation'].append(pert)
            
            # Write it to file
            new_protocol_file_string = os.path.join(sim_input_folder,
                                                    'protocol.json')
            
            with open(new_protocol_file_string, 'w') as f:
                json.dump(prot, f, indent=4)
                
            # Set the output folder
            sim_output_folder = os.path.join(top_data_dir,
                                            'sim_output',
                                            ('%i' % sim_counter))
            
            if not (os.path.isdir(sim_output_folder)):
                os.makedirs(sim_output_folder)
                
            new_results_file_string = os.path.join(sim_output_folder,
                                              'sim_output.txt')
            
            # Add in an output handler
            if ('template_files' in characterization_struct):
                oh = dict()
                oh['templated_images'] = []
                for (i,template_f) in enumerate(characterization_struct['template_files']):
                    tf = dict()
                    tf['relative_to'] = 'this_file'
                    tf['template_file_string'] = os.path.join(base_dir,
                                                              template_f)
                    tf['output_file_string'] = os.path.join(sim_output_folder,
                                                            'output_%i_%i' % (sim_counter, i+1))
                    tf['output_image_formats'] = ['png', 'svg']
                    oh['templated_images'].append(tf)
                
                # Now make the file
                new_output_handler_file_string = os.path.join(sim_input_folder,
                                                              'output_handler.json')
                
                with open(new_output_handler_file_string, 'w') as f:
                    json.dump(oh, f, indent=4)
                
            # Add the job
            j['model_file'] = str(Path(new_model_file_string).resolve().absolute())
            j['options_file'] = str(Path(new_options_file_string).resolve().absolute())
            j['protocol_file'] = str(Path(new_protocol_file_string).resolve().absolute())
            j['results_file'] = str(Path(new_results_file_string).resolve().absolute())
            
            if ('template_files' in characterization_struct):
                j['output_handler_file'] = str(Path(new_output_handler_file_string).resolve().absolute())
            
            # Update the counter
            sim_counter = sim_counter + 1
            
            job.append(j)
        
    # Add the job to the batch
    FiberV_batch['job'] = job
    
    # Insert everything into a MyoVent batch
    FiberVent_batch = dict()
    FiberVent_batch['FiberVent_batch'] = FiberV_batch
    
    # Now create the fig_data that makes figures
    batch_figs = dict()
    
    batch_output_dir = str(os.path.join(top_data_dir, 'sim_output'))
    
    # Rates
    batch_figs['rates'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = batch_output_dir
    fig['output_image_file'] = os.path.join(batch_output_dir,
                                            'rates')
    fig['output_image_formats'] = ['png']

    if ('formatting' in characterization_struct):
        fig['formatting'] = characterization_struct['formatting']

    batch_figs['rates'].append(fig)
    
    # Superposed traces
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(top_data_dir,
                                         'sim_output')
    fig['output_image_file'] = os.path.join(top_data_dir,
                                            'sim_output',
                                            'superposed_traces')
    fig['output_image_formats'] = ['png']
    
    fig['formatting'] = dict()
    fig['formatting']['y_label_pad'] = 35
    fig['formatting']['tight_layout'] = True
    fig['formatting']['y_label_fontsize'] = 8
    fig['formatting']['legend_fontsize'] = 7
    fig['formatting']['legend_bbox_to_anchor'] = [1.05, 1.2]
    
    fig['layout'] = dict()
    fig['layout']['panel_height'] = 0.7
    fig['layout']['left_margin'] = 0.1
    fig['layout']['right_margin'] = 0.1
    fig['layout']['grid_hspace'] = 0.2
    
    batch_figs['superposed_traces'].append(fig)
    
    if ('espvr_start_time_s' in characterization_struct):
        # Now make a pv plot
        batch_figs['pv_loops'] = []
        fig = dict()
        fig['relative_to'] = 'False'
        fig['results_folder'] = os.path.join(top_data_dir,
                                             'sim_output')
        fig['output_image_file'] = os.path.join(top_data_dir,
                                                'sim_output',
                                                'pv_loops')
        
        fig['output_image_formats'] = ['png']
        if ('output_image_formats' in characterization_struct):
            fig['output_image_formats'] = characterization_struct['output_image_formats']
        
        if ('espvr_start_time_s' in characterization_struct):
            fig['espvr_start_time_s'] = characterization_struct['espvr_start_time_s']
    
        fig['layout'] = dict()
        fig['layout']['panel_height'] = 3.5
        fig['layout']['top_margin'] = 0.05
        fig['layout']['bottom_margin'] = 0.25
        fig['layout']['grid_wspace'] = 0.5
        fig['layout']['grid_hspace'] = 0.2
        
        fig['formatting'] = dict()
        fig['formatting']['x_label_pad'] = 10
        fig['formatting']['tight_layout'] = True
        
        batch_figs['pv_loops'].append(fig)
    
    # Add in the figures
    FiberVent_batch['FiberVent_batch']['batch_figures'] = batch_figs
        
    # Create the batch file
    batch_file_string = os.path.join(top_data_dir,
                                     'batch.json')
    
    batch_file_string = str(Path(batch_file_string).absolute().resolve())
    
    with open(batch_file_string, 'w') as f:
        json.dump(FiberVent_batch, f, indent=4)
    
    # Now run it
    batch.run_batch(batch_file_string, figures_only)
            
def deduce_isovolumic_properties(json_analysis_file_string,
                                 characterization_struct,
                                 figures_only = False):
    """ Runs a (near-)isovolumic analysis """
    
    
    
    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_data = json.load(f)
        anal_struct = json_data['FiberVent_setup']
    
    # Create a batch
    FiberV_batch = dict()
    
    # Add in the MyoVentCpp stuff
    FiberV_batch['FiberVentCpp_exe'] = anal_struct['FiberVentCpp_exe']
    
    # Deduce the base model files
    if ('relative_to' in anal_struct['model']):
        if (anal_struct['model']['relative_to'] == 'this_file'):
            base_dir = str(Path(json_analysis_file_string).parent)
        else:
            base_dir = anal_struct['model']['relative_to']
    else:
        base_dir = ''
    
    base_model_files = []
    for mf in anal_struct['model']['model_files']:
        base_model_files.append(os.path.join(base_dir, mf))
        
    # Deduce the base options file
    # Can use the base dir from above
    base_options_file = os.path.join(base_dir,
                                    anal_struct['model']['options_file'])
    
    
    # Now do stuff based on the characterization
    
    # First, work out the base dir
    
    if ('relative_to' in characterization_struct):
        if (characterization_struct['relative_to'] == 'this_file'):
            base_dir = str(Path(json_analysis_file_string).parent)
        else:
            base_dir = characterization_struct['relative_to']
    else:
        base_dir = ''
        
    # Now make the output dir
    top_data_dir = os.path.join(base_dir, characterization_struct['sim_folder'])
    
    # Now clean the output_dir and make it
    if not (figures_only):
        try:
            print('Trying to remove %s' % top_data_dir)
            shutil.rmtree(top_data_dir, ignore_errors = True)
        except OSError as e:
            print('Error: %s : %s' % (top_data_dir, e.strerror))
            
        if not os.path.isdir(top_data_dir):
            os.makedirs(top_data_dir)

    # Set some counters and arrays for the simulations
    sim_counter = 1
    job = []
                
    # Now cycle through the model files
    for (model_counter, mf) in enumerate(base_model_files):
        
        # Load the base model
        with open(mf, 'r') as f:
            base_model = json.load(f)
            
        # Load the options file
        with open(base_options_file, 'r') as f:
            base_options = json.load(f)
        
        # Now cycle through volume_factors
        for (vol_counter, slack_vol_factor) in \
            enumerate(characterization_struct['ventricular_slack_volume_factors']):
    
                # Create a dictionary for the job
                j = dict()
    
                # Make the input folder
                sim_input_folder = os.path.join(top_data_dir,
                                                'sim_input',
                                                ('%i' % sim_counter))
                
                if not (os.path.isdir(sim_input_folder)):
                    os.makedirs(sim_input_folder)
    
                # Copy the volume and make some adjustments
                new_model = copy.deepcopy(base_model)
                
                # First adjust the slack volume of the ventricle
                y = new_model['MyoVent']['circulation']['compartments']['slack_volume'][0]
                
                new_model['MyoVent']['circulation']['compartments']['slack_volume'][0] = \
                    slack_vol_factor * y
                    
                # # Now set the resistance of the valves
                # y = new_model['MyoVent']['circulation']['compartments']['resistance'][1]
                
                # new_model['MyoVent']['circulation']['compartments']['resistance'][1] = \
                #     characterization_struct['aortic_valve_resistance_factor'] * y
                    
                # y = new_model['MyoVent']['circulation']['compartments']['resistance'][0]
                
                # new_model['MyoVent']['circulation']['compartments']['resistance'][0] = \
                #     characterization_struct['mitral_valve_resistance_factor'] * y
                    
                # If m_n is listed, set that
                if ('m_n' in characterization_struct):
                    new_model['MyoVent']['circulation']['ventricle']['myocardium']['contraction'] \
                        ['model']['muscle']['half_sarcomere']['thick_structure']['m_n'] = \
                            characterization_struct['m_n']
                    
                # Now write the model to file
                new_model_file_string = os.path.join(sim_input_folder,
                                                     'model.json')
                
                with open(new_model_file_string, 'w') as f:
                    json.dump(new_model, f, indent=4)
                    
                # Copy the options
                new_options = copy.deepcopy(base_options)
                
                new_options_file_string = os.path.join(sim_input_folder,
                                                       'options.json')
                
                with open(new_options_file_string, 'w') as f:
                    json.dump(new_options, f, indent=4)
                    
                # Create a protocol
                prot = dict()
                prot['protocol'] = dict()
                prot['protocol']['time_step_s'] = characterization_struct['time_step_s']
                prot['protocol']['no_of_time_steps'] = round(
                    characterization_struct['sim_duration_s'] / characterization_struct['time_step_s'])
                
                # Add in the perturbation
                prot['perturbation'] = characterization_struct['perturbation']
                
                # Write it to file
                new_protocol_file_string = os.path.join(sim_input_folder,
                                                        'protocol.json')
                
                with open(new_protocol_file_string, 'w') as f:
                    json.dump(prot, f, indent=4)
                    
                # Set the output folder
                sim_output_folder = os.path.join(top_data_dir,
                                                'sim_output',
                                                ('%i' % sim_counter))
                
                if not (os.path.isdir(sim_output_folder)):
                    os.makedirs(sim_output_folder)
                    
                new_results_file_string = os.path.join(sim_output_folder,
                                                  'sim_output.txt')
                    
                # Add the job
                j['model_file'] = str(Path(new_model_file_string).resolve().absolute())
                j['options_file'] = str(Path(new_options_file_string).resolve().absolute())
                j['protocol_file'] = str(Path(new_protocol_file_string).resolve().absolute())
                j['results_file'] = str(Path(new_results_file_string).resolve().absolute())
                
                # Update the counter
                sim_counter = sim_counter + 1
                
                job.append(j)
        
    # Add the job to the batch
    FiberV_batch['job'] = job
    
    # Insert everything into a MyoVent batch
    FiberVent_batch = dict()
    FiberVent_batch['FiberVent_batch'] = FiberV_batch
    
    # Now create the fig_data that makes figures
    batch_figs = dict()
    
    batch_figs['superposed_traces'] = []
    fig = dict()
    fig['relative_to'] = "False"
    fig['results_folder'] = os.path.join(top_data_dir,
                                         'sim_output')
    fig['output_image_file'] = os.path.join(top_data_dir,
                                            'sim_output',
                                            'superposed_traces')
    fig['output_image_formats'] = ['png']
    
    fig['formatting'] = dict()
    fig['formatting']['y_label_pad'] = 35
    fig['formatting']['tight_layout'] = True
    fig['formatting']['y_label_fontsize'] = 8
    fig['formatting']['legend_fontsize'] = 7
    fig['formatting']['legend_bbox_to_anchor'] = [1.05, 1.2]
    
    fig['layout'] = dict()
    fig['layout']['panel_height'] = 0.7
    fig['layout']['left_margin'] = 0.1
    fig['layout']['right_margin'] = 0.1
    fig['layout']['grid_hspace'] = 0.2
    
    batch_figs['superposed_traces'].append(fig)
    
    if ('espvr_start_time_s' in characterization_struct):
        # Now make a pv plot
        batch_figs['pv_loops'] = []
        fig = dict()
        fig['relative_to'] = 'False'
        fig['results_folder'] = os.path.join(top_data_dir,
                                             'sim_output')
        fig['output_image_file'] = os.path.join(top_data_dir,
                                                'sim_output',
                                                'pv_loops')
        
        fig['output_image_formats'] = ['png']
        if ('output_image_formats' in characterization_struct):
            fig['output_image_formats'] = characterization_struct['output_image_formats']
        
        if ('espvr_start_time_s' in characterization_struct):
            fig['espvr_start_time_s'] = characterization_struct['espvr_start_time_s']
    
        fig['layout'] = dict()
        fig['layout']['panel_height'] = 3.5
        fig['layout']['top_margin'] = 0.05
        fig['layout']['bottom_margin'] = 0.25
        fig['layout']['grid_wspace'] = 0.5
        fig['layout']['grid_hspace'] = 0.2
        
        fig['formatting'] = dict()
        fig['formatting']['x_label_pad'] = 10
        fig['formatting']['tight_layout'] = True
        
        batch_figs['pv_loops'].append(fig)
    
    # Add in the figures
    FiberVent_batch['FiberVent_batch']['batch_figures'] = batch_figs
        
    # Create the batch file
    batch_file_string = os.path.join(top_data_dir,
                                     'batch.json')
    
    batch_file_string = str(Path(batch_file_string).absolute().resolve())
    
    with open(batch_file_string, 'w') as f:
        json.dump(FiberVent_batch, f, indent=4)
    
    # Now run it
    batch.run_batch(batch_file_string, figures_only)
        
