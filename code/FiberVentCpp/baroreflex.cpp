/**
/* @file		baroreflex.cpp
/* @brief		Source file for a baroreflex object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include "math.h"

#include "baroreflex.h"

#include "cmv_model.h"
#include "cmv_system.h"
#include "cmv_results.h"
#include "cmv_options.h"

#include "circulation.h"
#include "reflex_control.h"

#include "gsl_math.h"

// Constructor
baroreflex::baroreflex(circulation* set_p_parent_circulation)
{
	// Initialise

	// Code
	std::cout << "baroreflex_constructor()\n";

	// Set pointers
	p_parent_circulation = set_p_parent_circulation;
	p_cmv_model = p_parent_circulation->p_cmv_model;
	p_parent_cmv_system = p_parent_circulation->p_parent_cmv_system;

	// Initialise with safe options
	p_cmv_results_beat = NULL;
	p_cmv_options = NULL;

	for (int i = 0; i < MAX_NO_OF_REFLEX_CONTROLS; i++)
	{
		p_rc[i] = NULL;
	}

	// Set safe options
	baro_active = 0.0;

	baro_A = 0.5;
	baro_B = 0.5;

	// Set from model
	baro_P_set = p_cmv_model->baro_P_set;
	baro_S = p_cmv_model->baro_S;
	baro_k_drive = p_cmv_model->baro_k_drive;
	baro_k_recov = p_cmv_model->baro_k_recov;
	baro_P_compartment = p_cmv_model->baro_P_compartment;

	// Set the controls
	no_of_reflex_controls = p_cmv_model->no_of_rc_controls;

	for (int i = 0; i < no_of_reflex_controls ; i++)
	{
		p_rc[i] = new reflex_control(this, (i + 1), p_cmv_model->p_rc[i]);
	}
}

// Destructor
baroreflex::~baroreflex(void)
{
	//! hemi_vent destructor

	// Tidy up
	for (int i = 0; i < MAX_NO_OF_REFLEX_CONTROLS; i++)
	{
		if (p_rc[i] != NULL)
			delete p_rc[i];
	}
}

// Other functions
void baroreflex::initialise_simulation(void)
{
	//! Code initialises simulation
	
	// Variables
	
	// Initialise options
	p_cmv_options = p_parent_circulation->p_cmv_options;

	// Now add in the results
	p_cmv_results_beat = p_parent_circulation->p_cmv_results_beat;

	// And now daughter objects
	for (int i = 0; i < no_of_reflex_controls; i++)
	{
		p_rc[i]->initialise_simulation();
	}

	// Add fields
	p_cmv_results_beat->add_results_field("baro_active", &baro_active);
	p_cmv_results_beat->add_results_field("baro_P_set", &baro_P_set);
	p_cmv_results_beat->add_results_field("baro_A", &baro_A);
	p_cmv_results_beat->add_results_field("baro_B", &baro_B);
}

void baroreflex::implement_time_step(double time_step_s)
{
	//! Implements time-step
	
	// Variables

	// Code

	calculate_B_a();

	calculate_B_b(time_step_s);

	// And now daughter objects
	for (int i = 0; i < no_of_reflex_controls; i++)
	{
		p_rc[i]->implement_time_step(time_step_s);
	}
}

void baroreflex::calculate_B_a(void)
{
	//! Function updates B_a
	
	// Code
	baro_A = 1.0 / (1.0 +
		exp(-baro_S * (p_parent_circulation->circ_pressure[baro_P_compartment - 1] -
			baro_P_set)));
}

void baroreflex::calculate_B_b(double time_step_s)
{
	//! Function updates B_b using an Euler step
	
	// Variables
	double delta;

	// Code

	// Calculate the delta over the time-step
	if (baro_active > 0)
	{
		if (baro_A >= 0.5)
		{
			delta = -baro_k_drive * (baro_A - 0.5) * baro_B;
		}
		else
		{
			delta = -baro_k_drive * (baro_A - 0.5) * (1.0 - baro_B);
		}
	}
	else
	{
		delta = -baro_k_recov * (baro_B - 0.5);
	}

	// Update
	baro_B = baro_B + (delta * time_step_s);

	// Limit
	baro_B = GSL_MIN(baro_B, 1.0);
	baro_B = GSL_MAX(baro_B, 0.0);
}
