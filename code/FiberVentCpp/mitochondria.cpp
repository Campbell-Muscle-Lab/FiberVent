/**
/* @file		mitochondria.cpp
/* @brief		Source file for a mitochondria object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include "mitochondria.h"

#include "muscle.h"

#include "hemi_vent.h"
#include "cmv_model.h"
#include "cmv_options.h"
#include "cmv_results.h"

// Constructor
mitochondria::mitochondria(muscle* set_p_parent_muscle)
{
	//! Constructor

	// Code
	// Set the pointer to the parent system
	p_parent_muscle = set_p_parent_muscle;
	p_cmv_model = p_parent_muscle->p_cmv_model;

	// Set other pointers safe
	p_cmv_results_beat = NULL;
	p_cmv_options = NULL;

	// Initialize
	mito_ATP_generation_rate = p_cmv_model->mito_ATP_generation_rate;

	mito_volume = 0.001 * p_parent_muscle->p_parent_hemi_vent->vent_wall_volume *
		(1.0 - p_parent_muscle->muscle_prop_fibrosis) *
		(1.0 - p_parent_muscle->muscle_prop_myofilaments);

	mito_ATP_generated_M_per_liter_per_s = 0.0;
}

// Destructor
mitochondria::~mitochondria(void)
{
	//! Destructor

	// Code
}

// Other functions
void mitochondria::initialise_simulation(void)
{
	//! Function adds data fields to main results object

	// Variables

	// Initialize

	// Set the options
	p_cmv_options = p_parent_muscle->p_cmv_options;

	// Set the pointer to the results object
	p_cmv_results_beat = p_parent_muscle->p_cmv_results_beat;

	// Now add the results fields
	
	p_cmv_results_beat->add_results_field("mito_volume", &mito_volume);
	p_cmv_results_beat->add_results_field("mito_ATP_generated_M_per_liter_per_s", &mito_ATP_generated_M_per_liter_per_s);
}

void mitochondria::implement_time_step(double time_step_s)
{
	//! Function updates membrane object by a time-step

	// Variables

	double mito_volume;			// volume of mitochondria in liters in a liter of myocardium

	// Code

	mito_volume = 0.001 *
		(1.0 - p_parent_muscle->muscle_prop_fibrosis) *
		(1.0 - p_parent_muscle->muscle_prop_myofilaments);

	mito_ATP_generated_M_per_liter_per_s = mito_volume * mito_ATP_generation_rate;
}
