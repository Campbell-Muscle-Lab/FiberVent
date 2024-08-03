/**
/* @file		membranes.cpp
/* @brief		Source file for a membranes object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include "membranes.h"
#include "muscle.h"
#include "cmv_model.h"
#include "cmv_options.h"
#include "cmv_results.h"

#include "gsl_errno.h"
#include "gsl_odeiv2.h"


// Constructor
membranes::membranes(muscle* set_p_parent_muscle)
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
	memb_Ca_cytosol = 0.0;
	memb_Ca_sr = p_cmv_model->memb_Ca_content;
	memb_activation = 0.0;
	memb_t_open_left_s = 0.0;

	memb_t_open_s = p_cmv_model->memb_t_open_s;
	memb_k_serca = p_cmv_model->memb_k_serca;
	memb_k_leak = p_cmv_model->memb_k_leak;
	memb_k_active = p_cmv_model->memb_k_active;

	memb_J_release = 0.0;
	memb_J_uptake = 0.0;

}

// Destructor
membranes::~membranes(void)
{
	//! Destructor

	// Code
}

// Other functions
void membranes::initialise_simulation(void)
{
	//! Function adds data fields to main results object

	// Variables

	// Initialize

	// Set the options
	p_cmv_options = p_parent_muscle->p_cmv_options;

	// Set the pointer to the results object
	p_cmv_results_beat = p_parent_muscle->p_cmv_results_beat;

	// Now add the results fields
	p_cmv_results_beat->add_results_field("memb_Ca_cytosol", &memb_Ca_cytosol);
	p_cmv_results_beat->add_results_field("memb_Ca_sr", &memb_Ca_sr);
	p_cmv_results_beat->add_results_field("memb_activation", &memb_activation);
	p_cmv_results_beat->add_results_field("memb_t_open_left_s", &memb_t_open_left_s);
	p_cmv_results_beat->add_results_field("memb_k_serca", &memb_k_serca);
	p_cmv_results_beat->add_results_field("memb_k_leak", &memb_k_leak);
	p_cmv_results_beat->add_results_field("memb_k_active", &memb_k_active);
	p_cmv_results_beat->add_results_field("memb_J_release", &memb_J_release);
	p_cmv_results_beat->add_results_field("memb_J_uptake", &memb_J_uptake);
}

// This function is not a member of the membranes class but is used to interace
// with the GSL ODE system. It must appear before the membranes class members
// that calls it, and communicates with the membrane class through a pointer to
// the class object

int memb_calculate_derivs(double t, const double y[], double f[], void* params)
{
	//! Function sets derivs

	// Variables
	(void)(t);

	membranes* p_memb = (membranes *)params;

	// Code

	// Update fluxes
	p_memb->calculate_fluxes(y);

	// f[0] is rate of change of cytosol concentration
	// f[1] is rate of change of sr concentration

	f[0] = p_memb->memb_J_release - p_memb->memb_J_uptake;
	f[1] = -f[0];

	return GSL_SUCCESS;
}

void membranes::implement_time_step(double time_step_s, bool new_beat)
{
	//! Function updates membrane object by a time-step

	// Variables
	double eps_abs = 1e-6;
	double eps_rel = 1e-6;

	int status;

	double t_start_s = 0.0;
	double t_stop_s = time_step_s;

	double y[2] = { memb_Ca_cytosol, memb_Ca_sr };

	// Code
	if (new_beat)
	{
		memb_t_open_left_s = memb_t_open_s;
	}
	else
	{
		memb_t_open_left_s = memb_t_open_left_s - time_step_s;
	}

	if (memb_t_open_left_s > 0.0)
	{
		memb_activation = 1.0;
	}
	else
	{
		memb_activation = 0.0;
	}

	gsl_odeiv2_system sys = { memb_calculate_derivs, NULL, 2, this };

	gsl_odeiv2_driver* d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
			0.5*time_step_s, eps_abs, eps_rel);

	status = gsl_odeiv2_driver_apply(d, &t_start_s, t_stop_s, y);

	gsl_odeiv2_driver_free(d);

	if (status != GSL_SUCCESS)
	{
		std::cout << "Integration problem in membranes::implement_time_step\n";
	}
	else
	{
		//Unpack
		memb_Ca_cytosol = y[0];
		memb_Ca_sr = y[1];
	}
}

void membranes::calculate_fluxes(const double y[])
{
	//! Function calculates fluxes

	// Variables
	double Ca_cytosol = y[0];
	double Ca_sr = y[1];

	// Code
	memb_J_release = (memb_k_leak + (memb_activation * memb_k_active)) * Ca_sr;
	memb_J_uptake = memb_k_serca * Ca_cytosol;
}
