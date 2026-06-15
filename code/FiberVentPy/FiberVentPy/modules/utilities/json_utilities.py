import os
import json
import copy

from pathlib import Path
import shutil

def return_FiberVentCpp_exe_dict(json_analysis_file_string):
    """ Returns a dictionary for the FiberVentCpp executable information """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)
    
    # Extract the FiberCpp_exe struct
    FiberVentCpp_exe_dict = json_dict['FiberVent_setup']['FiberVentCpp_exe']

    # If we are in a relative mode, adapt the FiberVentCpp_exe for absolute paths
    # because the new simulations will be run from a different folder

    if (FiberVentCpp_exe_dict['relative_to'] == "this_file"):
        base_dir = Path(json_analysis_file_string).parent.absolute()
        FiberVentCpp_exe_dict['relative_to'] = 'False'
        FiberVentCpp_exe_dict['exe_file'] = \
            str(Path(os.path.join(base_dir, FiberVentCpp_exe_dict['exe_file'])).
                absolute().resolve())

    return FiberVentCpp_exe_dict

def return_base_dir(json_analysis_file_string,
                    dict_key,
                    append_key = [],
                    append_two_keys = [],
                    append_string = [],
                    dict_index = 0):
    """ Return the base directory for the specified dict
        If the dict is an array, use the dict_id'th member """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    # Pull off the target_dict
    target_dict = json_dict['FiberVent_setup'][dict_key]

    if (isinstance(target_dict, list)):
        target_dict = target_dict[dict_index]

    # See if we need to adjust the filesystem paths
    if ('relative_to' in target_dict):
        if (target_dict['relative_to'] == "this_file"):
            base_dir = Path(json_analysis_file_string).parent.absolute()
        else:
            base_dir = target_dict['relative_to']
        if not (append_key == []):
            base_dir = os.path.join(base_dir, target_dict[append_key])
        if not (append_two_keys == []):
            base_dir = os.path.join(base_dir, target_dict[append_two_keys[0]][append_two_keys[1]])
        if not (append_string == []):
            base_dir = os.path.join(base_dir, append_string)
    else:
        base_dir = target_dict[append_key]
        if not (append_string == []):
            base_dir = os.path.join(base_dir, append_string)

    # Return
    return base_dir

def return_model_file_strings(json_analysis_file_string):
    """ Return a list of absolute paths to model files """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    model_dict = json_dict['FiberVent_setup']['model']

    # Get the base directory
    base_dir = return_base_dir(json_analysis_file_string, 'model')

    # Get the model file strings
    model_file_strings = []
    for mfs in model_dict['model_files']:
        model_file_strings.append(
            str(Path(os.path.join(base_dir, mfs)).resolve().absolute()))

    # Return
    return model_file_strings

def return_options_file_string(json_analysis_file_string):
    """ Return an absolute path to the options file """

    # Load the analysis file
    with open(json_analysis_file_string, 'r') as f:
        json_dict = json.load(f)

    model_dict = json_dict['FiberVent_setup']['model']

    # Get the base directory
    base_dir = return_base_dir(json_analysis_file_string, 'model')

    # Get the file_string
    options_file_string = os.path.join(base_dir, model_dict['options_file'])

    # Return
    return options_file_string

def prepare_clean_dir(clean_dir):
    # Clean the dir
    try:
        print('Trying to remove: %s' % clean_dir)
        shutil.rmtree(clean_dir, ignore_errors = True)
    except OSError as e:
        print('Error: %s : %s')

    # Make the clean_dir if it does not exist
    if not os.path.exists(clean_dir):
        os.makedirs(clean_dir)