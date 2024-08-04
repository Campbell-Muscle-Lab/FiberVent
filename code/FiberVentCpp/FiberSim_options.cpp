/**
 * @file    FiberSim_options.cpp
 * @brief   Source file for the FiberSim_options class
 * @author  Ken Campbell
 */

#include <cstdio>
#include <filesystem>
#include <string>

#include "FiberSim_options.h"
#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

namespace fs = std::filesystem;

// Constructor
FiberSim_options::FiberSim_options(const rapidjson::Value& doc)
{
    // Initialise

    // Set default values
    x_pos_rel_tol = 0.001;                  /**< default value of absolute tolerance for
                                                 calculating x positions */

    x_vector_max_iterations = 100;          /**< default value of max iterations for
                                                 x vector calculation */

    hs_force_control_max_delta_hs_length = 100;
                                            /**< default value of max_delta_hs_length for
                                                 force-control for a lattice */

    adjacent_bs = 0;                        /**< default value of adjacent binding sites
                                                 available for myosin or mybpc attachment
                                                 0 restricts to nearest site */

    max_rate = 1000.0;                      /**< default value for maximum rate allowed
                                                 in calculations */

    thin_filament_sub_steps = 10;            /**< default value for number of sub_steps
                                                 used for thin filament kinetics */

    lambda_jitter = 0.0;                    /**< default value for lambda jitter */

    dump_precision = 6;                     /**< default value for dump precision */

    calculate_x_mode = 1;                   /**< default value, implies use sparse_F_and_G method
                                                 to calculate x positions */

    sprintf_s(rate_file_string, _MAX_PATH, "");
                                            /**< default value for rate file string */

    myofibril_force_tolerance = 0.001;      /**< default value for the force tolerance
                                                 for myofibril multiroot calculations */

    myofibril_max_iterations = 100;         /**< default value for the maximum number of
                                                 iterations in myofibril multiroot
                                                 calculations */

    min_hs_length = 500.0;                  /**< default value for the half-sarcomere length
                                                 at which the simulation collapses and gives
                                                 up trying to keep the calculations going */

    no_of_worker_threads = 0;               /**< default value for number of worker threads */

    // Update values from log file
    set_FiberSim_options_from_JSON_file_string(doc);
}

// Destructor
FiberSim_options::~FiberSim_options()
{
    // Tidy up

    if (log_mode > 0)
    {
        // Close log file
        fclose(log_file);
    }
}

void FiberSim_options::set_FiberSim_options_from_JSON_file_string(const rapidjson::Value& options)
{
    // Check we have entries and set them

    if (JSON_functions::is_JSON_member(options, "max_rate"))
    {
        JSON_functions::check_JSON_member_number(options, "max_rate");
        max_rate = options["max_rate"].GetDouble();
    }

    if (JSON_functions::is_JSON_member(options, "thin_filament_sub_steps"))
    {
        JSON_functions::check_JSON_member_int(options, "thin_filament_sub_steps");
        thin_filament_sub_steps = options["thin_filament_sub_steps"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "adjacent_bs"))
    {
        JSON_functions::check_JSON_member_int(options, "adjacent_bs");
        adjacent_bs = options["adjacent_bs"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "x_pos_rel_tol"))
    {
        JSON_functions::check_JSON_member_number(options, "x_pos_rel_tol");
        x_pos_rel_tol = options["x_pos_rel_tol"].GetDouble();
    }

    // In FiberSim versions > 1.2.0, x_pos_tol notation is preferred over x_pos_rel_tol

    if (JSON_functions::is_JSON_member(options, "x_pos_tol"))
    {
        JSON_functions::check_JSON_member_number(options, "x_pos_tol");
        x_pos_rel_tol = options["x_pos_tol"].GetDouble();
    }

    if (JSON_functions::is_JSON_member(options, "x_vector_max_iterations"))
    {
        JSON_functions::check_JSON_member_int(options, "x_vector_max_iterations");
        x_vector_max_iterations = options["x_vector_max_iterations"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "calculate_x_mode"))
    {
        JSON_functions::check_JSON_member_number(options, "calculate_x_mode");
        calculate_x_mode = options["calculate_x_mode"].GetInt();
    }

    if (JSON_functions::is_JSON_member(options, "hs_force_control_max_delta_hs_length"))
    {
        JSON_functions::check_JSON_member_number(options, "hs_force_control_max_delta_hs_length");
        hs_force_control_max_delta_hs_length = options["hs_force_control_max_delta_hs_length"].GetDouble();
    }

    if (JSON_functions::is_JSON_member(options, "min_hs_length"))
    {
        JSON_functions::check_JSON_member_number(options, "min_hs_length");
        min_hs_length = options["min_hs_length"].GetDouble();
    }

    // Check if lambda_jitter was specified.
    if (JSON_functions::is_JSON_member(options, "lambda_jitter"))
    {
        JSON_functions::check_JSON_member_number(options, "lambda_jitter");
        lambda_jitter = options["lambda_jitter"].GetDouble();
    }

    // Check for rand_seed - set to empty string if missing
    if (JSON_functions::is_JSON_member(options, "rand_seed"))
    {
        JSON_functions::check_JSON_member_string(options, "rand_seed");
        sprintf_s(rand_seed, _MAX_PATH, "%s", options["rand_seed"].GetString());
    }
    else
    {
        sprintf_s(rand_seed, _MAX_PATH, "");
    }

    // Check for myofibrils
    if (JSON_functions::is_JSON_member(options, "myofibrils"))
    {
        const rapidjson::Value& myofibrils = options["myofibrils"];

        printf("\n\nKen was here\n");

        JSON_functions::check_JSON_member_number(myofibrils, "force_tolerance");
        myofibril_force_tolerance = myofibrils["force_tolerance"].GetDouble();

        JSON_functions::check_JSON_member_int(myofibrils, "max_iterations");
        myofibril_max_iterations = myofibrils["max_iterations"].GetInt();

        JSON_functions::check_JSON_member_number(myofibrils, "max_delta_hs_length");
        myofibril_max_delta_hs_length = myofibrils["max_delta_hs_length"].GetDouble();
    }

    // Now check for logging
    if (JSON_functions::is_JSON_member(options, "logging"))
    {
        const rapidjson::Value& logging = options["logging"];

        JSON_functions::check_JSON_member_string(logging, "relative_to");
        sprintf_s(log_relative_to, _MAX_PATH, "%s", logging["relative_to"].GetString());

        JSON_functions::check_JSON_member_string(logging, "log_folder");
        sprintf_s(log_folder, _MAX_PATH, "%s", logging["log_folder"].GetString());
    }

    // Check for rate logging
    if (JSON_functions::is_JSON_member(options, "rate_files"))
    {
        const rapidjson::Value& rate_files = options["rate_files"];

        JSON_functions::check_JSON_member_string(rate_files, "file");
        sprintf_s(rate_file_string, _MAX_PATH, "%s", rate_files["file"].GetString());
    }

    // Now check for status files
    if (JSON_functions::is_JSON_member(options, "status_files"))
    {
        const rapidjson::Value& status_files = options["status_files"];

        JSON_functions::check_JSON_member_string(status_files, "status_folder");
        sprintf_s(status_folder, _MAX_PATH, "%s", status_files["status_folder"].GetString());

        JSON_functions::check_JSON_member_string(status_files, "time_steps");
        sprintf_s(time_steps_string, _MAX_PATH, "%s", status_files["time_steps"].GetString());

        // Check if the dump precision was specified.
        if (JSON_functions::is_JSON_member(status_files, "dump_precision"))
        {
            JSON_functions::check_JSON_member_int(status_files, "dump_precision");
            dump_precision = status_files["dump_precision"].GetInt();
        }
    }
}

void FiberSim_options::adapt_file_paths(std::string results_file_string)
{
    //! Code updates paths for logging simulations
    
    // Variables

    // Code

    fs::path results_file = results_file_string;
    fs::path results_path = results_file.parent_path();

    // Create the rate_file_string
    if (strlen(rate_file_string) > 0)
    {
        fs::path rate_path = results_path / rate_file_string;
        sprintf_s(rate_file_string, _MAX_PATH, "%s", rate_path.string().c_str());
    }

    // Now the status_files
    if (strlen(status_folder) > 0)
    {
        fs::path status_path = fs::absolute(results_path / status_folder);

        // Make sure the status folder exists
        if (fs::exists(status_path))
        {
            // Clean the directory
            printf("Cleaning status_folder: %s\n", status_path.string().c_str());
            for (auto const& dir_entry : fs::recursive_directory_iterator(status_path))
            {
                fs::remove(dir_entry);
            }
        }
        else
        {
            // Create the directory
            if (fs::create_directories(status_path))
            {
                printf("Status folder created at: %s\n", status_path.string().c_str());
            }
            else
            {
                printf("Status folder could not be created: %s\n", status_path.string().c_str());
                exit(1);
            }

        }
        // Set the status folder
        sprintf_s(status_folder, _MAX_PATH, "%s", status_path.string().c_str());

        // Parse the time_steps string
        std::string ts_string = time_steps_string;

        size_t first_sep = ts_string.find_first_of(":");
        size_t last_sep = ts_string.find_last_of(":");

        start_status_time_step = (int)std::stoi(ts_string.substr(0, first_sep));
        skip_status_time_step = (int)std::stoi(ts_string.substr((first_sep + 1), last_sep));
        stop_status_time_step = (int)std::stoi(ts_string.substr(last_sep + 1));
    }
}

void FiberSim_options::write_FiberSim_options_to_file(void)
{
    // Code writes FiberSim_options to file

    // Variables
    char output_file_string[_MAX_PATH];
    FILE* output_file;

    // Code
    sprintf_s(output_file_string, _MAX_PATH, "%s\\%s",
        log_folder, "FiberSim_options.log");

    errno_t err = fopen_s(&output_file, output_file_string, "w");
    if (err != 0)
    {
        printf("Options log file file: %s\ncould not be opened\n",
            output_file_string);
        exit(1);
    }

    fprintf_s(output_file, "log_folder: %s\n", log_folder);
    fprintf_s(output_file, "log_file_string: %s\n", log_file_string);
    fprintf_s(output_file, "max_rate: %g\n", max_rate);
    fprintf_s(output_file, "dump_precision: %i\n", dump_precision);

    // Tidy up
    fclose(output_file);
}
