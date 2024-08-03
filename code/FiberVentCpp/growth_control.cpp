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
#include "gsl/gsl_fit.h"

struct cmv_model_gc_structure {
	string type;
	string level;
	string signal;
	double set_point;
	double prop_gain;
	double deriv_gain;
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
	gc_slope = 0.0;

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
	gc_max_rate = p_struct->max_rate;

	gc_deriv_points = 0;
	gc_deriv_x = NULL;
	gc_deriv_y = NULL;
}

// Destructor
growth_control::~growth_control(void)
{
	// Tidy up
	if (gc_deriv_x != NULL)
		free(gc_deriv_x);

	if (gc_deriv_y != NULL)
		free(gc_deriv_y);
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

	if (gc_type == "concentric")
	{
		// Need to calculate the deriv
		gc_deriv_points = p_cmv_options->growth_control_deriv_points;
		gc_deriv_x = (double*)malloc(gc_deriv_points * sizeof(double));
		gc_deriv_y = (double*)malloc(gc_deriv_points * sizeof(double));

		printf("\n\ngc_deriv_points: %i\n\n", gc_deriv_points);
	}

	// Set the derivitives if necessary
	if (gc_deriv_points > 0)
	{
		for (int i = 0; i < gc_deriv_points; i++)
		{
			gc_deriv_x[i] = GSL_NAN;
			gc_deriv_y[i] = GSL_NAN;
		}
	}

	// Add fields
	temp_string = "gc_" + to_string(gc_number) + "_prop_signal";
	p_cmv_results_beat->add_results_field(temp_string, &gc_prop_signal);

	temp_string = "gc_" + to_string(gc_number) + "_deriv_signal";
	p_cmv_results_beat->add_results_field(temp_string, &gc_deriv_signal);

	temp_string = "gc_" + to_string(gc_number) + "_output";
	p_cmv_results_beat->add_results_field(temp_string, &gc_output);

	temp_string = "gc_" + to_string(gc_number) + "_slope";
	p_cmv_results_beat->add_results_field(temp_string, &gc_slope);
}

void growth_control::implement_time_step(double time_step_s, bool new_beat)
{
	//! Implements time-step
	
	// Variables

	// Code

		// Update the arrays
	if (gc_deriv_points > 0)
	{
		for (int i = 0; i < (gc_deriv_points - 1); i++)
		{
			gc_deriv_x[i] = gc_deriv_x[i + 1];
			gc_deriv_y[i] = gc_deriv_y[i + 1];
		}
		gc_deriv_x[gc_deriv_points - 1] = p_parent_cmv_system->cum_time_s;
		gc_deriv_y[gc_deriv_points - 1] = *gc_p_signal;
	}

	// If we are calculating the slope, and it's a new beat, update
	if (gc_deriv_points > 0)
	{
		calculate_slope();
	}

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

	// Code

	if (p_parent_growth->growth_active > 0.0)
	{
		if (gc_set_point != 0.0)
			gc_prop_signal = p_parent_growth->gr_master_rate *
				gc_prop_gain * (*gc_p_signal - gc_set_point) / gc_set_point;
		else
			gc_prop_signal = p_parent_growth->gr_master_rate *
			gc_prop_gain * (*gc_p_signal - gc_set_point);

		gc_deriv_signal = gc_deriv_gain * gc_slope;
	}
	else
	{
		gc_prop_signal = 0.0;
		gc_deriv_signal = 0.0;
	}

	gc_output = gc_prop_signal + gc_deriv_signal;

	// Limit
	if (gc_output >= 0)
		gc_output = GSL_MIN(gc_output, gc_max_rate);
	else

		gc_output = GSL_MAX(gc_output, -gc_max_rate);
}

void growth_control::calculate_slope(void)
{
	//! Function calculates the rate of change of the control signal
	
	// Variables
	bool nan_flag = false;
	double c0, c1, cov00, cov10, cov11, sumsq;

	// Code

	// Check there are no NANs
	for (int i = 0; i < gc_deriv_points; i++)
	{
		if (gsl_isnan(gc_deriv_x[i]))
			nan_flag = true;
		if (gsl_isnan(gc_deriv_y[i]))
			nan_flag = true;
	}
	if (nan_flag)
	{
		gc_slope = 0;
		return;
	}

	// Otherwise continue
	gsl_fit_linear(gc_deriv_x, 1, gc_deriv_y, 1, gc_deriv_points,
		&c0, &c1, &cov00, &cov10, &cov11, &sumsq);

	// Set the slope
	gc_slope = c1;
}