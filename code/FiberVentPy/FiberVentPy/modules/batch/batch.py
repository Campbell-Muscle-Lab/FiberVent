# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 15:28:53 2022

@author: kscamp3
"""

import os
import json

from pathlib import Path

import multiprocessing
import threading
import subprocess

from ..output_handler import output_handler as oh
from ..display import analyses

def run_batch(json_batch_file_string, figures_only=False):
    """ Runs >=1 job using multithreading """
    
    # Load the batch file
    with open(json_batch_file_string, 'r') as f:
        json_data = json.load(f)
        FiberVent_batch = json_data['FiberVent_batch']
    
    # Pull off the exe structure
    exe_structure = FiberVent_batch['FiberVentCpp_exe']
    exe_string = exe_structure['exe_file']
    if not ('relative_to' in exe_structure):
        exe_string = os.path.abspath(exe_string)
    elif (exe_structure['relative_to'] == 'this_file'):
        base_dir = Path(json_batch_file_string).parent.absolute()
        exe_string = os.path.join(base_dir, exe_string)
    else:
        base_dir = exe_structure['relative_to']
        exe_string = os.path.join(base_dir, exe_string)
    
    # Parse the job data into a list of command strings
    job_data = FiberVent_batch['job']
    
    command_strings = []
    results_file_strings = []
    
    for i,j in enumerate(job_data):
        # Build up the command string
        com_string = exe_string
        for f in ['model_file', 'options_file',
                  'protocol_file', 'results_file']:
            fs = j[f]
            if not ('relative_to' in j):
                fs = os.path.abspath(fs)
            elif (j['relative_to'] == 'this_file'):
                base_directory = Path(json_batch_file_string).parent.absolute()
                fs = os.path.join(base_directory, fs)
            else:
                base_directory = j['relative_to']
                fs = os.path.join(base_directory, fs)
            # Store the results file in case you need it to make a figure
            # using the output_handler system
            if (f == 'results_file'):
                results_file_strings.append(fs)

            com_string = '%s "%s"' % (com_string, fs)

        # Add in the system_id
        com_string = ('%s %i') % (com_string, (i+1))

        command_strings.append(com_string)

    # Check the batch to see if max threads have been specified
    if ('max_threads' in FiberVent_batch):
        requested_max_threads = FiberVent_batch['max_threads']
    else:
        requested_max_threads = float("inf")
        
    # Get max threads available
    available_threads = multiprocessing.cpu_count()-1

    # Set processes to mininmum of requested and available
    num_processes = int(min([requested_max_threads, available_threads]))

    if (figures_only == False):
        print('Running batch using %i threads' % num_processes)
                
        # Now run the batch
        my_list = command_strings
    
        threads = []
        while threads or my_list:
            if (len(threads) < num_processes) and my_list:
                t = threading.Thread(target=worker, args=[my_list.pop()])
                t.setDaemon(True)
                t.start()
                threads.append(t)
            else:
                for thread in threads:
                    if not thread.is_alive():
                        threads.remove(thread)  

    # At this point we have run all the simulations
    # Run the output handlers in parallel
   
    for i, j in enumerate(job_data):
        if ('output_handler_file' in j):
            fs = j['output_handler_file']
            if not ('relative_to' in j):
                fs = os.path.abspath(fs)
            elif (j['relative_to'] == 'this_file'):
                base_directory = Path(json_batch_file_string).parent.absolute()
                fs = os.path.join(base_directory, fs)
            else:
                base_directory = j['relative_to']
                fs = os.path.join(base_directory, fs)

            oh.output_handler(fs,
                   sim_results_file_string=results_file_strings[i])
            # try:
            #     oh.output_handler(fs,
            #                   sim_results_file_string=results_file_strings[i])
            # except:
            #     print('Could not implement output_handler for: %s' %
            #           results_file_strings[i])
            
    # Run the batch figures
    if ('batch_figures' in FiberVent_batch):
        
        if ('rates' in FiberVent_batch['batch_figures']):
            print('Now generating rates figure')
            for fig_data in FiberVent_batch['batch_figures']['rates']:
                analyses.create_rates_figure(fig_data,
                                             json_batch_file_string)
                
        if ('superposed_traces' in FiberVent_batch['batch_figures']):
            for fig_data in FiberVent_batch['batch_figures']['superposed_traces']:
                analyses.create_superposed_traces_figure(
                    fig_data,
                    json_batch_file_string)
                
        if ('pv_loops' in FiberVent_batch['batch_figures']):
            for fig_data in FiberVent_batch['batch_figures']['pv_loops']:
                analyses.create_pv_loops_figure(
                    fig_data,
                    json_batch_file_string)
    
def worker(cmd):
    subprocess.call(cmd)

   
