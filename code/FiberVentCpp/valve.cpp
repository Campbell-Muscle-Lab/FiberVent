/**
/* @file		half_sarcomere.cpp
/* @brief		Source file for a half_sarcomere object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include <iostream>

#include "cmv_model.h"
#include "valve.h"
#include "hemi_vent.h"
#include "circulation.h"
#include "cmv_results.h"
#include "membranes.h"
#include "muscle.h"
#include "heart_rate.h"

#include "gsl_errno.h"
#include "gsl_odeiv2.h"
#include "gsl_math.h"
#include "gsl_roots.h"

struct cmv_model_valve_structure {
	string name;
	double mass;
	double eta;
	double k;
	double leak;
};

// Constructor
valve::valve(hemi_vent* set_p_parent_hemi_vent, cmv_model_valve_structure* set_p_cmv_model_structure)
{
	//! Constructor

	// Code

	// Set the pointers to the parent system
	p_parent_hemi_vent = set_p_parent_hemi_vent;
	p_cmv_model = p_parent_hemi_vent->p_cmv_model;

	// Set safe values
	p_cmv_options = NULL;
	p_cmv_results_beat = NULL;

	// Update from cmv_model
	p_cmv_model_valve = set_p_cmv_model_structure;

	valve_name = p_cmv_model_valve->name;
	valve_mass = p_cmv_model_valve->mass;
	valve_eta = p_cmv_model_valve->eta;
	valve_k = p_cmv_model_valve->k;
	valve_leak = p_cmv_model_valve->leak;

	// Initialise
	valve_pos = 0.0;
	valve_last_pos = 0.0;
	valve_vel = 0.0;
}

// Destructor
valve::~valve(void)
{
	//! Destructor

	// Code

	// Tidy up
}

// Other functions
void valve::initialise_simulation(void)
{
	//! Code initialises simulation
	
	// Variables
	string label_string;
	
	// Code

	// Set options from parent
	p_cmv_options = p_parent_hemi_vent->p_cmv_options;

	// Set results from parent
	p_cmv_results_beat = p_parent_hemi_vent->p_cmv_results_beat;

	// Now add the results fields
	label_string = valve_name + "_valve_pos";

	p_cmv_results_beat->add_results_field(label_string, &valve_pos);
}

// This function is not a member of the valve class but is used to interace
// with the GSL ODE system. It must appear before the class members
// that calls it, and communicates with the valve class through a pointer to
// the class object

int valve_derivs(double t, const double y[], double f[], void* params)
{
	// Function sets dV/dt for the compartments

	// Variables
	(void)(t);							// Prevents warning for unused variable

	valve* p_valve = (valve*)params;	// Pointer to valve

	circulation* p_circ = p_valve->p_parent_hemi_vent->p_parent_circulation;

	double pressure_difference = GSL_NAN;

	// Code

	// y[0] is the position of the valve
	// y[1] is the velocity

	if (p_valve->valve_name == "aortic")
	{
		pressure_difference = p_circ->circ_pressure[0] - p_circ->circ_pressure[1];
		//cout << "pd: " << pressure_difference << "\n";
	}

	if (p_valve->valve_name == "mitral")
	{
		pressure_difference = p_circ->circ_pressure[p_circ->circ_no_of_compartments - 1] - p_circ->circ_pressure[0];
		//cout << "pd: " << pressure_difference << "\n";
	}

	if (gsl_isnan(pressure_difference))
	{
		cout << "Pressure difference not defined for valve\n";
		exit(1);
	}

	f[0] = y[1];
	f[1] = (1.0 / p_valve->valve_mass) *
		((-p_valve->valve_eta * y[1]) - (p_valve->valve_k * y[0]) +
			pressure_difference);

	// Return
	return GSL_SUCCESS;
}

void valve::implement_time_step(double time_step_s)
{
	//! Implements time-step
	
	// Variables

	double eps_abs = 1e-6;
	double eps_rel = 1e-6;

	int status;

	double t_start_s = 0.0;
	double t_stop_s = time_step_s;

	double* y_calc = NULL;

	// Code

	// Allocate memory for y
	y_calc = (double*)malloc(2 * sizeof(double));

	// Fill y_calc
	y_calc[0] = valve_pos;
	y_calc[1] = valve_vel;

	gsl_odeiv2_system sys = { valve_derivs, NULL, 2, this };

	gsl_odeiv2_driver* d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
			0.5 * time_step_s, eps_abs, eps_rel);

	status = gsl_odeiv2_driver_apply(d, &t_start_s, t_stop_s, y_calc);

	gsl_odeiv2_driver_free(d);

	// Unpack
	valve_pos = y_calc[0];
	valve_vel = y_calc[1];

	// Bounds
	if (valve_pos > 1.0)
	{
		valve_pos = 1.0;
		valve_vel = 0.0;
	}

	if (valve_pos < valve_leak)
	{
		valve_pos = valve_leak;
		valve_vel = 0.0;
	}

	// Tidy up
	free(y_calc);

/*
	double pressure_difference;

	circulation* p_circ = p_parent_hemi_vent->p_parent_circulation;

	if (valve_name == "aortic")
	{
		pressure_difference = p_circ->circ_pressure[0] - p_circ->circ_pressure[1];
		//cout << "pd: " << pressure_difference << "\n";
	}

	if (valve_name == "mitral")
	{
		pressure_difference = p_circ->circ_pressure[p_circ->circ_no_of_compartments - 1] - p_circ->circ_pressure[0];
		//cout << "pd: " << pressure_difference << "\n";
	}

	double temp_x;

	temp_x = valve_last_pos;

	valve_pos = (pressure_difference + ((valve_eta * valve_last_pos) / time_step_s)) /
		(valve_k + (valve_eta / time_step_s));

	valve_last_pos = temp_x;

	// Bounds
	if (valve_pos > 1.0)
	{
		valve_pos = 1.0;
		valve_vel = 0.0;
	}

	if (valve_pos < valve_leak)
	{
		valve_pos = valve_leak;
		valve_vel = 0.0;
	}
*/
}
