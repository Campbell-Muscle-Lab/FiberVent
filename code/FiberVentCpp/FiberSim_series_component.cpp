/**
/* @file		FiberSim_series_component.cpp
/* @brief		Source file for a FiberSim_series_component object
/* @author		Ken Campbell
*/

#include <cstdio>

#include "FiberSim_series_component.h"

#include "FiberSim_muscle.h"
#include "FiberSim_model.h"

#include "gsl_math.h"

#include "muscle.h"

#include "cmv_model.h"
#include "cmv_results.h"



// Constructor
FiberSim_series_component::FiberSim_series_component(FiberSim_muscle* set_p_parent_fs_muscle)
{
	//! Constructor

	p_parent_fs_muscle = set_p_parent_fs_muscle;

	// Update from model
	FiberSim_model* p_fs_model;
	p_fs_model = p_parent_fs_muscle->p_parent_muscle->p_cmv_model->p_fs_model;

	sc_k_stiff = p_fs_model->sc_k_stiff;
	sc_eta = p_fs_model->sc_eta;

	// Set some defaults
	sc_extension = 0.0;
	sc_last_extension = 0.0;
	sc_force = 0.0;
}

// Destructor
FiberSim_series_component::~FiberSim_series_component(void)
{
	// Destructor
}

// Other functions
void FiberSim_series_component::initialise_for_simulation(void)
{
	//! Initialises for simulation
	
	// Variables
	p_cmv_results_beat = p_parent_fs_muscle->p_cmv_results_beat;

	// Add results fields
	p_cmv_results_beat->add_results_field("fs_sc_extension", &sc_extension);
	p_cmv_results_beat->add_results_field("fs_sc_last_extension", &sc_last_extension);
	p_cmv_results_beat->add_results_field("fs_sc_force", &sc_force);

}

double FiberSim_series_component::return_series_extension(double muscle_force, double time_step_s)
{
	//! Returns the extension of the series component for a given force
	//! Assumption is that F = kx + eta*(x-x_last)/delta_t

	// Variables
	double ext;

	// Code
	if (gsl_isnan(sc_k_stiff))
	{
		printf("Error in fs_series component: sc_k_stiff is NAN\n");
		exit(1);
	}

	if (time_step_s > 0.0)
	{
		ext = (muscle_force + ((sc_eta * sc_last_extension) / time_step_s)) /
			(sc_k_stiff + (sc_eta / time_step_s));
	}
	else
	{
		ext = (muscle_force / sc_k_stiff);
	}

	return ext;
}

double FiberSim_series_component::return_series_force_for_length(double test_length, double time_step_s)
{
	//! Returns force for test_length

	// Variables
	double elastic_force;
	double viscous_force;
	double test_force;

	// Code
	elastic_force = sc_k_stiff * test_length;

	if (time_step_s > 0)
	{
		viscous_force = sc_eta * (test_length - sc_extension) / time_step_s;
	}
	else
	{
		viscous_force = 0.0;
	}

	test_force = elastic_force + viscous_force;

	return test_force;
}
