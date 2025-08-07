/**
/* @file		growth.cpp
/* @brief		Source file for a growth object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include "math.h"

#include "growth.h"

#include "cmv_model.h"
#include "cmv_system.h"
#include "cmv_results.h"
#include "cmv_options.h"

#include "growth_control.h"
#include "circulation.h"
#include "hemi_vent.h"
#include "muscle.h"

#include "gsl_math.h"

// Constructor
growth::growth(circulation* set_p_parent_circulation)
{
	// Initialise

	// Code

	// Set pointers
	p_parent_circulation = set_p_parent_circulation;
	p_cmv_model = p_parent_circulation->p_cmv_model;
	p_parent_cmv_system = p_parent_circulation->p_parent_cmv_system;

	// Initialise with safe options
	p_cmv_results_beat = NULL;
	p_cmv_options = NULL;

	for (int i = 0; i < MAX_NO_OF_GROWTH_CONTROLS; i++)
	{
		p_gc[i] = NULL;
	}

	growth_active = 0.0;

	// Set from model

	gr_master_rate = p_cmv_model->gr_master_rate;

	// Set the controls
	no_of_growth_controls = p_cmv_model->no_of_gc_controls;

	for (int i = 0; i < no_of_growth_controls ; i++)
	{
		p_gc[i] = new growth_control(this, (i + 1), p_cmv_model->p_gc[i]);
	}
}

// Destructor
growth::~growth(void)
{
	//! Growth destructor

	// Tidy up
	for (int i = 0; i < MAX_NO_OF_GROWTH_CONTROLS; i++)
	{
		if (p_gc[i] != NULL)
			delete p_gc[i];
	}
}

// Other functions
void growth::initialise_simulation(void)
{
	//! Code initialises simulation
	
	// Variables
	
	// Initialise options
	p_cmv_options = p_parent_circulation->p_cmv_options;

	// Now add in the results
	p_cmv_results_beat = p_parent_circulation->p_cmv_results_beat;

	// And now daughter objects
	for (int i = 0; i < no_of_growth_controls; i++)
	{
		p_gc[i]->initialise_simulation();
	}

	// Add fields
	p_cmv_results_beat->add_results_field("growth_active", &growth_active);
}

void growth::implement_time_step(double time_step_s, bool new_beat)
{
	//! Implements time-step
	
	// Variables
	double internal_r;
	double wall_thickness;
	
	double delta_relative_wall_thickness;
	double delta_relative_n_hs;

	double delta_n_hs;
	double delta_muscle_length;

	// Code

	// Zero relative changes
	delta_relative_wall_thickness = 0.0;
	delta_relative_n_hs = 0.0;
	
	// Handle daughter objects
	for (int i = 0; i < no_of_growth_controls; i++)
	{
		p_gc[i]->implement_time_step(time_step_s, new_beat);

		if (gsl_isnan(p_gc[i]->gc_output))
			continue;

		if (p_gc[i]->gc_type == "concentric")
		{
			delta_relative_wall_thickness = delta_relative_wall_thickness +
				(time_step_s * p_gc[i]->gc_output);
		}

		if (p_gc[i]->gc_type == "eccentric")
		{
			delta_relative_n_hs = delta_relative_n_hs +
				(time_step_s * p_gc[i]->gc_output);
		}
	}

	// Apply the concentric growth
	if (fabs(delta_relative_wall_thickness) > 0.0)
	{
		internal_r = p_parent_circulation->p_hemi_vent->
			return_internal_radius_for_chamber_volume(p_parent_circulation->circ_volume[0]);
		wall_thickness = p_parent_circulation->p_hemi_vent->
			return_wall_thickness_for_chamber_volume(p_parent_circulation->circ_volume[0]);
		
		// Equation below is for hemisphere
		/*p_parent_circulation->p_hemi_vent->vent_wall_volume =
			1000 * ((2.0 / 3.0) * M_PI *
				pow((internal_r + (wall_thickness * (1.0 + delta_relative_wall_thickness))), 3.0)) -
			p_parent_circulation->circ_volume[0];
		*/

		// This is for ellipsoid
		p_parent_circulation->p_hemi_vent->vent_wall_volume =
			1000 * ((2.0 / 3.0) * M_PI * 
				pow((internal_r + (wall_thickness * (1.0 + delta_relative_wall_thickness))), 2.0) *
				(p_parent_circulation->p_hemi_vent->vent_chamber_height  + 
					(wall_thickness * (1.0 + delta_relative_wall_thickness)))) -
			p_parent_circulation->circ_volume[0];
	}

	// Now the eccentric growth
	if (fabs(delta_relative_n_hs) > 0.0)
	{
		// Work out how far half-sarcomeres move using chain rule
		delta_n_hs = delta_relative_n_hs * p_parent_circulation->p_hemi_vent->vent_n_hs;

		delta_muscle_length = -(delta_n_hs * p_parent_circulation->p_hemi_vent->p_muscle->muscle_length) /
			p_parent_circulation->p_hemi_vent->vent_n_hs;

		p_parent_circulation->p_hemi_vent->p_muscle->change_muscle_length(delta_muscle_length, time_step_s);

		// Update the wall volume
		p_parent_circulation->p_hemi_vent->vent_wall_volume =
			p_parent_circulation->p_hemi_vent->vent_wall_volume *
			(1.0 + (delta_n_hs / p_parent_circulation->p_hemi_vent->vent_n_hs));

		// Update n_hs
		p_parent_circulation->p_hemi_vent->vent_n_hs =
			p_parent_circulation->p_hemi_vent->vent_n_hs + delta_n_hs;
	}
}
