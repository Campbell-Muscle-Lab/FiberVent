/**
/* @file		growth_control.cpp
/* @brief		Source file for a growth control object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include "math.h"

#include <iostream>
#include <regex>

#include "growth_control.h"

#include "cmv_model.h"
#include "cmv_system.h"
#include "cmv_results.h"
#include "cmv_options.h"

#include "growth.h"
#include "circulation.h"
#include "hemi_vent.h"
#include "muscle.h"
#include "heart_rate.h"
#include "membranes.h"


#include "kinetic_scheme.h"
#include "transition.h"

#include "FiberSim_muscle.h"
#include "FiberSim_half_sarcomere.h"

#include "gsl_math.h"
#include "gsl_statistics.h"
#include "gsl/gsl_fit.h"

struct cmv_model_gc_structure {
	string type;
	string level;
	string signal;
	double set_point;
	double prop_gain;
	double deriv_gain;
	double control_period_s;
	double max_rate;
};

// Constructor
growth_control::growth_control(growth* set_p_parent_growth, int set_gc_number,
	cmv_model_gc_structure* p_struct)
{
	// Initialise

	// Code
	
	// Set pointers
	p_parent_growth = set_p_parent_growth;

	p_parent_circulation = p_parent_growth->p_parent_circulation;
	p_parent_cmv_system = p_parent_circulation->p_parent_cmv_system;
	p_cmv_model = p_parent_circulation->p_cmv_model;

	// Initialise with safe options
	p_cmv_results_beat = NULL;
	p_cmv_options = NULL;

	gc_p_signal = NULL;

	// Other variables
	gc_number = set_gc_number;

	gc_output = 0.0;

	gc_prop_output = 0.0;
	gc_deriv_output = 0.0;

	gc_prop_signal = 0.0;
	gc_deriv_signal = 0.0;

	gc_signal_assigned = false;

	// Others come from the structure
	gc_type = p_struct->type;
	gc_level = p_struct->level;
	gc_signal = p_struct->signal;
	gc_prop_gain = p_struct->prop_gain;
	gc_deriv_gain = p_struct->deriv_gain;
	gc_set_point = p_struct->set_point;
	gc_control_period_s = p_struct->control_period_s;
	gc_max_rate = p_struct->max_rate;

	gc_control_points = 0;
	gc_control_t = NULL;
	gc_control_y = NULL;
}

// Destructor
growth_control::~growth_control(void)
{
	// Tidy up
	if (gc_control_t != NULL)
		free(gc_control_t);

	if (gc_control_y != NULL)
		free(gc_control_y);
}

// Other functions
void growth_control::initialise_simulation(void)
{
	//! Code initialises simulation

	// Variables
	string temp_string;
	
	// Initialise options
	p_cmv_options = p_parent_circulation->p_cmv_options;

	// Now add in the results
	p_cmv_results_beat = p_parent_circulation->p_cmv_results_beat;

	// Find the controlled variable
	set_gc_p_signal();

	// Set the number of control points
	if (gsl_isnan(gc_control_period_s))
	{
		// Means control_period was not set in file
		gc_control_points = 1;
	}
	else
	{
		gc_control_points = (int)(gc_control_period_s /
			p_parent_cmv_system->time_step_s);
	}

	// Allocate space
	gc_control_t = (double*)malloc(gc_control_points * sizeof(double));
	gc_control_y = (double*)malloc(gc_control_points * sizeof(double));

	// Initialise
	for (int i = 0 ; i < gc_control_points ; i++)
	{
		gc_control_t[i] = GSL_NAN;
		gc_control_y[i] = GSL_NAN;
	}

	// Add fields
	temp_string = "gc_" + to_string(gc_number) + "_prop_signal";
	p_cmv_results_beat->add_results_field(temp_string, &gc_prop_signal);

	temp_string = "gc_" + to_string(gc_number) + "_deriv_signal";
	p_cmv_results_beat->add_results_field(temp_string, &gc_deriv_signal);

	temp_string = "gc_" + to_string(gc_number) + "_prop_output";
	p_cmv_results_beat->add_results_field(temp_string, &gc_prop_output);

	temp_string = "gc_" + to_string(gc_number) + "_deriv_output";
	p_cmv_results_beat->add_results_field(temp_string, &gc_deriv_output);

	temp_string = "gc_" + to_string(gc_number) + "_output";
	p_cmv_results_beat->add_results_field(temp_string, &gc_output);
}

void growth_control::implement_time_step(double time_step_s, bool new_beat)
{
	//! Implements time-step
	
	// Variables
	int i;

	// Code

	// Update the arrays
	for (i = 0; i < (gc_control_points - 1); i++)
	{
		gc_control_t[i] = gc_control_t[i + 1];
		gc_control_y[i] = gc_control_y[i + 1];
	}
	gc_control_t[gc_control_points - 1] = p_parent_cmv_system->cum_time_s;

	if (gc_set_point != 0.0)
		gc_control_y[gc_control_points - 1] = (*gc_p_signal - gc_set_point) /
			gc_set_point;
	else
		gc_control_y[gc_control_points - 1] = (*gc_p_signal - gc_set_point);

	// Calculate the output
	calculate_output();
}

void growth_control::set_gc_p_signal(void)
{
	//! Code sets the pointer for gc_p_signal to the appropriate double

	// Variables

	// Code

	// Find the variable this object is controlling
	if (gc_level == "MyoSim_half_sarcomere")
	{
		printf("Growth for MyoSim not yet implemented\n");
		exit(1);
		/*
		if (gc_signal == "myof_stress_int_pas")
		{
			gc_p_signal = &p_parent_circulation->p_hemi_vent->p_hs->p_myofilaments->myof_stress_int_pas;
			gc_signal_assigned = true;
		}

		if (gc_signal == "myof_mean_stress_int_pas")
		{
			gc_p_signal = &p_parent_circulation->p_hemi_vent->p_hs->p_myofilaments->myof_mean_stress_int_pas;
			gc_signal_assigned = true;
		}
		*/
	}

	if (gc_level == "FiberSim_half_sarcomere")
	{
		if (gc_signal == "fs_stress_titin")
		{
			gc_p_signal = &p_parent_circulation->p_hemi_vent->p_muscle->p_FiberSim_muscle->p_FiberSim_hs->hs_titin_force;
			gc_signal_assigned = true;
		}
	}

	if (gc_level == "muscle")
	{
		if (gc_signal == "muscle_ATP_concentration")
		{
			gc_p_signal = &p_parent_circulation->p_hemi_vent->p_muscle->muscle_ATP_concentration;
			gc_signal_assigned = true;
		}
	}

	if (gc_level == "ventricle")
	{
		if (gc_signal == "vent_cardiac_output")
		{
			gc_p_signal = &p_parent_circulation->p_hemi_vent->vent_cardiac_output;
			gc_signal_assigned = true;
		}
	}

	if (gc_signal_assigned == false)
	{
		cout << "Growth control " << gc_level << ", " << gc_signal << " not assigned\n";
		exit(1);
	}
	else
	{
		cout << "Growth control " << gc_level << ", " << gc_signal <<
			" current_value: " << *gc_p_signal << "\n";
	}
}

void growth_control::calculate_output(void)
{
	//! Code sets the output value

	// Variables
	double c0, c1, cov00, cov10, cov11, sumsq;

	// Code

	// Calculate the prop signal from the mean of the array
	gc_prop_signal = gsl_stats_mean(gc_control_y, 1, gc_control_points);

	// Calculate the deriv signal from the slope of the array
	if (gc_control_points > 1)
	{
		gsl_fit_linear(gc_control_t, 1, gc_control_y, 1, gc_control_points,
			&c0, &c1, &cov00, &cov10, &cov11, &sumsq);

		gc_deriv_signal = c1;
	}
	else
	{
		gc_deriv_signal = 0.0;
	}

	// Check growth is active
	if (p_parent_growth->growth_active > 0.0)
	{
		gc_prop_output = gc_prop_gain * gc_prop_signal;
		gc_deriv_output = gc_deriv_gain * gc_deriv_signal;
	}
	else
	{
		gc_prop_output = 0.0;
		gc_deriv_output = 0.0;
	}

	// Add prop and deriv together
	gc_output = gc_prop_output + gc_deriv_output;

	// Check for NaNs
	if (gsl_isnan(gc_output))
	{
		gc_output = 0.0;
	}

	// Limit

	// Check the max_rate is >=0, this can cause issues if you are
	// adjusting it via the protocol
	if ((gc_max_rate < 0.0) || (gc_max_rate < 1e-10))
	{
		gc_max_rate = 0.0;
	}

	// Now limit to max rate
	if (gc_output >= 0)
	{
		gc_output = GSL_MIN(gc_output, gc_max_rate);
	}
	else
	{
		gc_output = GSL_MAX(gc_output, -gc_max_rate);
	}
}
