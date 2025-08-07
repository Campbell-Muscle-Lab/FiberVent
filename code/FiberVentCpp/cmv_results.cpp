/**
/* @file		cmv_results.cpp
/* @brief		Source file for a cmv_results object
/* @author		Ken Campbell
*/

#include <stdio.h>
#include <iostream>
#include <filesystem>
#include <string>

#include "cmv_results.h"
#include "cmv_system.h"
#include "cmv_options.h"
#include "cmv_model.h"

#include "circulation.h"
#include "baroreflex.h"
#include "hemi_vent.h"
#include "muscle.h"
#include "FiberSim_muscle.h"
#include "FiberSim_half_sarcomere.h"

#include "gsl_vector.h"
#include "gsl_math.h"
#include "gsl_const_mksa.h"

using namespace std;
using namespace std::filesystem;

struct stats_structure {
	double mean_value;
	double min_value;
	double max_value;
	double sum;
};


// Constructor
cmv_results::cmv_results(cmv_system* set_p_parent_cmv_system, int set_no_of_time_points)
{
	// Initialise

	// Code
	printf("cmv_results constructor()\n");

	p_parent_cmv_system = set_p_parent_cmv_system;
	p_cmv_options = p_parent_cmv_system->p_cmv_options;

	no_of_defined_results_fields = 0;

	no_of_time_points = set_no_of_time_points;

	no_of_beats = 0;

	time_field_index = -1;
	new_beat_field_index = -1;
	pressure_vent_field_index = -1;
	pressure_veins_field_index = -1;
	volume_vent_field_index = -1;
	flow_mitral_valve_field_index = -1;
	flow_aortic_valve_field_index = -1;
	muscle_length_field_index = -1;
	muscle_stress_int_pas_field_index = -1;
	muscle_ATP_flux_field_index = -1;
//	myof_stress_int_pas_field_index = -1;
//	myof_ATP_flux_field_index = -1;
	vent_stroke_work_field_index = -1;
	vent_stroke_energy_used_field_index = -1;
	vent_efficiency_field_index = -1;
	vent_ejection_fraction_field_index = -1;
	vent_ATP_used_per_s_field_index = -1;
	vent_stroke_volume_field_index = -1;
	vent_cardiac_output_field_index = -1;

	// Special case
	pressure_arteries_field_index = -1;
}

// Destructor
cmv_results::~cmv_results(void)
{
	// Code

	// Delete all of the defined vectors
	for (int i = 0; i < no_of_defined_results_fields; i++)
	{
		gsl_vector_free(gsl_results_vectors[i]);
	}
}

// Other functions

void cmv_results::add_results_field(std::string field_name, double* p_double)
{
	//! Functions adds a field to the results object

	// Variables
	int new_index;			// index of new field

	bool index_set;

	// Code

	// Get the number of fields that have been defined already
	new_index = no_of_defined_results_fields;

	// Update the results_fields and the data source
	results_fields[new_index] = field_name;
	p_data_sources[new_index] = p_double;

	// Create a gsl_vector to hold the data and initialise to NaN
	gsl_results_vectors[new_index] = gsl_vector_alloc(no_of_time_points);
	gsl_vector_set_all(gsl_results_vectors[new_index], GSL_NAN);

	// Check for specific indices

	index_set = false;

	if (field_name == "time")
	{
		time_field_index = new_index;
		index_set = true;
	}

	if (field_name == "hr_new_beat")
	{
		new_beat_field_index = new_index;
		index_set = true;
	}

	if (field_name == "hr_heart_rate_bpm")
	{
		hr_heart_rate_bpm_field_index = new_index;
		index_set = true;
	}

	if (field_name == "pressure_0")
	{
		pressure_vent_field_index = new_index;
		index_set = true;
	}

	if (field_name == "volume_0")
	{
		volume_vent_field_index = new_index;
		index_set = true;
	}

	string venous_pressure = "pressure_" +
		to_string(p_parent_cmv_system->p_cmv_model->circ_no_of_compartments - 1);
	if (field_name == venous_pressure)
	{
		pressure_veins_field_index = new_index;
		index_set = true;
	}

	if (field_name == "flow_0")
	{
		flow_mitral_valve_field_index = new_index;
		index_set = true;
	}

	if (field_name == "flow_1")
	{
		flow_aortic_valve_field_index = new_index;
		index_set = true;
	}

	if (field_name == "muscle_length")
	{
		muscle_length_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_stroke_work_J")
	{
		vent_stroke_work_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_stroke_energy_used_J")
	{
		vent_stroke_energy_used_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_efficiency")
	{
		vent_efficiency_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_ejection_fraction")
	{
		vent_ejection_fraction_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_ATP_used_per_s")
	{
		vent_ATP_used_per_s_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_stroke_volume")
	{
		vent_stroke_volume_field_index = new_index;
		index_set = true;
	}

	if (field_name == "vent_cardiac_output")
	{
		vent_cardiac_output_field_index = new_index;
		index_set = true;
	}

	if (p_parent_cmv_system->p_circulation->p_baroreflex != NULL)
	{
		string b_string = "pressure_" +
			to_string(p_parent_cmv_system->p_circulation->p_baroreflex->
				baro_P_compartment);
		if (field_name == b_string)
		{
			pressure_arteries_field_index = new_index;
			index_set = true;
		}
	}

	// Update the number of defined fields
	no_of_defined_results_fields = no_of_defined_results_fields + 1;
}

void cmv_results::save_status(int t_index)
{
	//! Updates data records
	
	// Variables
	FiberSim_muscle* p_FiberSim_muscle = p_parent_cmv_system->p_circulation->p_hemi_vent->p_muscle->p_FiberSim_muscle;
	
	// Code

	// Cycle through the defined data fields updating the results vectors
	for (int i = 0; i < no_of_defined_results_fields; i++)
	{
		gsl_vector_set(gsl_results_vectors[i], t_index, *p_data_sources[i]);
	}

	// Check if we need to dump the status file
	if (p_FiberSim_muscle != NULL)
	{
		p_FiberSim_muscle->dump_hs_status(p_parent_cmv_system->sim_t_index);
	}
}

void cmv_results::calculate_sub_vector_statistics(gsl_vector* gsl_v, int start_index, int stop_index,
	stats_structure* p_stats)
{
	//! Function calculates stats
	
	// Variables
	int index;
	double value;
	double holder = 0;
	
	// Code

	p_stats->min_value = GSL_POSINF;
	p_stats->max_value = -GSL_POSINF;
	p_stats->mean_value = GSL_NAN;
	p_stats->sum = 0.0;

	if (start_index < 0)
		return;

	for (index = start_index; index <= stop_index; index++)
	{
		value = gsl_vector_get(gsl_v, index);
		holder = holder + value;
		p_stats->min_value = GSL_MIN(value, p_stats->min_value);
		p_stats->max_value = GSL_MAX(value, p_stats->max_value);
	}
	p_stats->sum = holder;
	p_stats->mean_value = holder / (double)(stop_index - start_index + 1);
}

int cmv_results::write_data_to_file(std::string output_file_string)
{
	//! Function writes data to file

	// Variables
	FILE* output_file;

	// Code
	cout << "Writing simulation results to: " << output_file_string << " points \n";

	// Make sure results directory exists
	path output_file_path(output_file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			cout << "\nError: Results folder could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	errno_t err = fopen_s(&output_file, output_file_string.c_str(), "w");
	if (err != 0)
	{
		cout << "Results file: " << output_file_string << " could not be opened\n";
		exit(1);
	}

	// Write header
	for (int i = 0; i < no_of_defined_results_fields; i++)
	{
		fprintf_s(output_file, "%s", results_fields[i].c_str());
		if (i == (no_of_defined_results_fields - 1))
			fprintf_s(output_file, "\n");
		else
			fprintf_s(output_file, "\t");
	}

	// Now data
	for (int i = 0; i < no_of_time_points; i = i + 1)
	{
		for (int j = 0; j < no_of_defined_results_fields; j++)
		{
			fprintf_s(output_file, "%g", gsl_vector_get(gsl_results_vectors[j], i));
			if (j == (no_of_defined_results_fields - 1))
				fprintf_s(output_file, "\n");
			else
				fprintf_s(output_file, "\t");
		}
	}

	cout << "Closing output_file: " << output_file_string << "\n";

	// Tidy up
	fclose(output_file);

	return(1);
}

double cmv_results::return_stroke_work(int start_t_index, int stop_t_index)
{
	//! Calculate stroke work via Shoelace formula
	
	// Variables
	double holder;

	// Code

	holder = 0.0;

	// Implement Shoelace
	for (int i = start_t_index; i <= stop_t_index; i++)
	{
		int j = i + 1;
		if (j > stop_t_index)
			j = 0;

		holder = holder +
			(gsl_vector_get(gsl_results_vectors[pressure_vent_field_index], i) +
				gsl_vector_get(gsl_results_vectors[pressure_vent_field_index], j)) *
			(gsl_vector_get(gsl_results_vectors[volume_vent_field_index], i) -
				gsl_vector_get(gsl_results_vectors[volume_vent_field_index], j));
	}

	holder = 0.5 * holder;

	// Now convert volumes to m^3
	holder = holder * 0.001;

	// Now convert pressure to Pa
	holder = holder * (0.001 * GSL_CONST_MKSA_METER_OF_MERCURY);

	return holder;
}

double cmv_results::return_energy_used(int start_t_index, int stop_t_index)
{
	//! Returns energy used as the integral of the ATP used over
	//! the cycle in moles multipled by the energy from ATP
	//!	Delta_G_ATP is set in model file as Joules / mole, 45000
	//!	https://equilibrator.weizmann.ac.il/static/classic_rxns/classic_reactions/atp.html

	// Variables

	double delta_t;
	double holder;
	double energy_used;

	// Code

	holder = 0.0;

	// Integrate energy used over time

	// Calculate the time-step
	delta_t = gsl_vector_get(gsl_results_vectors[time_field_index], stop_t_index) -
		gsl_vector_get(gsl_results_vectors[time_field_index], stop_t_index - 1);

	for (int i = start_t_index; i <= stop_t_index; i++)
	{
		holder = holder +
			gsl_vector_get(gsl_results_vectors[vent_ATP_used_per_s_field_index], i);
	}

	// Calculate energy per second
	energy_used = holder * delta_t *
		p_parent_cmv_system->p_circulation->p_hemi_vent->p_muscle->muscle_delta_G_ATP;

	return energy_used;
}

void cmv_results::backfill_beat_data(gsl_vector* gsl_v, double value, int start_t_index, int stop_t_index)
{
	//! Backfills data for a cardiac cycle
	
	// Variables

	// Code

	for (int i = start_t_index; i <= stop_t_index; i++)
	{
		gsl_vector_set(gsl_v, i, value);
	}
}
