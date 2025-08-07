/**
/* @file		muscle.cpp
/* @brief		Source file for a muscle object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include "cmv_model.h"
#include "cmv_options.h"
#include "cmv_results.h"

#include "FiberSim_model.h"

#include "muscle.h"

#include "hemi_vent.h"
#include "membranes.h"
#include "mitochondria.h"

#include "heart_rate.h"

#include "MyoSim_muscle.h"

#include "FiberSim_muscle.h"
#include "FiberSim_half_sarcomere.h"
#include "FiberSim_series_component.h"

#include "gsl_errno.h"
#include "gsl_math.h"
#include "gsl_roots.h"
#include "gsl_const_num.h"

struct stats_structure {
	double mean_value;
	double min_value;
	double max_value;
	double sum;
};

// Constructor
muscle::muscle(hemi_vent* set_p_parent_hemi_vent)
{
	//! Constructor

	// Code
	printf("muscle constructor()\n");

	// Set the pointers to the parent system
	p_parent_hemi_vent = set_p_parent_hemi_vent;
	p_cmv_model = p_parent_hemi_vent->p_cmv_model;
	p_cmv_system = p_parent_hemi_vent->p_parent_cmv_system;
	p_cmv_options = p_parent_hemi_vent->p_cmv_options;

	// Create the daugher objects
	p_heart_rate = new heart_rate(this);
	p_membranes = new membranes(this);
	p_mitochondria = new mitochondria(this);

	// Create the correct half-sarcomere model
	if (p_cmv_model->p_fs_model == NULL)
	{
		p_MyoSim_muscle = new MyoSim_muscle(this);
		p_FiberSim_muscle = NULL;
	}
	else
	{
		// We are making a FiberSim muscle
		p_MyoSim_muscle = NULL;
		p_FiberSim_muscle = new FiberSim_muscle(this);

		muscle_length = p_FiberSim_muscle->fs_m_length;
		muscle_stress = p_FiberSim_muscle->fs_m_stress;
	}

	// Set safe values
	p_cmv_results_beat = NULL;

	// Initialise
	muscle_ATP_used_total_per_liter_per_s = 0.0;
	muscle_reference_length = p_cmv_model->mus_reference_hs_length;
	muscle_delta_G_ATP = p_cmv_model->mus_delta_G_ATP;
	muscle_ATP_concentration = p_cmv_model->mus_initial_ATP_concentration;
	muscle_prop_fibrosis = p_cmv_model->mus_prop_fibrosis;
	muscle_prop_myofilaments = p_cmv_model->mus_prop_myofilaments;
}

// Destructor
muscle::~muscle(void)
{
	//! Destructor

	// Code

	// Tidy up
	delete p_heart_rate;
	delete p_membranes;
	delete p_mitochondria;

	delete p_FiberSim_muscle;
	delete p_MyoSim_muscle;
}

// Other functions
void muscle::initialise_simulation(void)
{
	//! Code initialises simulation

	// Variables
	double slack_length;

	// Code

	// Set results from parent
	p_cmv_results_beat = p_parent_hemi_vent->p_cmv_results_beat;

	// Now add the results fields
	p_cmv_results_beat->add_results_field("muscle_length", &muscle_length);
	p_cmv_results_beat->add_results_field("muscle_stress", &muscle_stress);
	p_cmv_results_beat->add_results_field("muscle_ATP_used_myof_per_liter_per_s",
		&muscle_ATP_used_myof_per_liter_per_s);
	p_cmv_results_beat->add_results_field("muscle_ATP_used_memb_per_liter_per_s",
		&muscle_ATP_used_memb_per_liter_per_s);
	p_cmv_results_beat->add_results_field("muscle_ATP_used_total_per_liter_per_s",
		&muscle_ATP_used_total_per_liter_per_s);
	p_cmv_results_beat->add_results_field("muscle_ATP_generated_per_liter_per_s",
		&muscle_ATP_generated_per_liter_per_s);
	p_cmv_results_beat->add_results_field("muscle_ATP_concentration", &muscle_ATP_concentration);
	p_cmv_results_beat->add_results_field("muscle_prop_fibrosis", &muscle_prop_fibrosis);

	// Now initialise daughter objects
	p_heart_rate->initialise_simulation();
	p_membranes->initialise_simulation();
	p_mitochondria->initialise_simulation();

	if (p_FiberSim_muscle != NULL)
	{
		p_FiberSim_muscle->initialise_for_simulation();

		p_cmv_results_beat->add_results_field("muscle_stress_int_pas",
			&p_FiberSim_muscle->p_FiberSim_hs->hs_titin_force);
	}
	else
	{
		printf("Not yet implemented\n");
		exit(1);
		//p_MyoSim_muscle->initialise_for_simulation();
	}

	// Now adjust muscle to slack length
	slack_length = return_muscle_length_for_stress(0.0, 0.00);

	change_muscle_length(slack_length - muscle_length, 0.0);

	// Adjust the series component last length
	if (p_FiberSim_muscle != NULL)
	{
		p_FiberSim_muscle->p_FiberSim_sc->sc_last_extension =
			p_FiberSim_muscle->p_FiberSim_sc->sc_extension;
	}
	else
	{
		printf("Not yet implemented\n");
		exit(1);
	}

	printf("Muscle is slack\n");
	printf("Slack_length: %g\n", slack_length);
	printf("muscle_force: %g\n", muscle_stress);
}

bool muscle::implement_time_step(double time_step_s)
{
	//! Implements time-step

	// Variables
	bool new_beat = false;

	// Code

	// Update daughter objects
	new_beat = p_heart_rate->implement_time_step(time_step_s);
	
	p_membranes->implement_time_step(time_step_s, new_beat);

	p_mitochondria->implement_time_step(time_step_s);

	// Apply to the correct muscle
	if (p_FiberSim_muscle != NULL)
	{
		p_FiberSim_muscle->implement_time_step(time_step_s);
	}
	else
	{
		printf("Not yet implemented\n");
		exit(1);

		//p_myofilaments->implement_time_step(time_step_s);
	}

	// Update the ATP
	calculate_muscle_ATP_concentration(time_step_s);

	// Return new beat status
	return (new_beat);
}

void muscle::change_muscle_length(double delta_ml, double time_step_s)
{
	//! Changes muscle length by delta_ml

	// Code

	// Adjust the muscle length
	muscle_length = muscle_length + delta_ml;

	if (p_FiberSim_muscle != NULL)
	{
		p_FiberSim_muscle->change_muscle_length(delta_ml, time_step_s);
		muscle_stress = p_FiberSim_muscle->fs_m_stress;
	}
	else
	{
		/*
		// And the cross-bridges
		p_myofilaments->move_cb_populations(delta_hsl);

		// Now update wall stress
		p_myofilaments->calculate_stresses();

		hs_stress = p_myofilaments->myof_stress_total;
		*/
	}
}

double muscle::return_wall_stress_after_test_delta_ml(double delta_ml, double time_step_s)
{
	//! Function returns wall stress after given delta_hsl
	
	// Variables
	double wall_stress;

	// Code
	if (p_FiberSim_muscle != NULL)
	{
		wall_stress = p_FiberSim_muscle->return_wall_stress_after_test_delta_ml(delta_ml, time_step_s);
	}
	else
	{
		printf("Not yet implemented\n");
		wall_stress = 0.0;
		exit(1);
		//wall_stress = p_myofilaments->return_stress_after_delta_hsl(delta_hsl);
	}

	return wall_stress;
}

struct gsl_root_params
{
	muscle* p_muscle_temp;
	double stress_target;
	double time_step_s;
};

double muscle_stress_root_finder(double x, void* params)
{
	//! Function for the root finder
	
	// Variables
	double new_stress;
	double stress_difference;

	// Code
	
	// Unpack the pointer
	struct gsl_root_params* p = (struct gsl_root_params*)params;

	// Calculate the new stress after a length change
	new_stress = p->p_muscle_temp->return_wall_stress_after_test_delta_ml(x, p->time_step_s);

	// Calculate the difference between the new stress and the target
	stress_difference = new_stress - p->stress_target;

	// Return the difference
	return stress_difference;
}

double muscle::return_muscle_length_for_stress(double target_stress, double time_step_s)
{
	//! Code returns the muscle_length at which force is equal to the target

	// Variables
	double x_lo = -p_cmv_options->force_control_max_delta_hs_length;
	double x_hi = -x_lo;
	double r;

	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;

	gsl_function F;
	struct gsl_root_params params = { this, target_stress, time_step_s};

	int status;
	int iter = 0;
	int max_iter = 100;

	double epsabs = 0.01;
	double epsrel = 0.01;

	// Code

	F.function = &muscle_stress_root_finder;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);

		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi, epsabs, epsrel);

	} while ((status == GSL_CONTINUE) && (iter < max_iter));

	gsl_root_fsolver_free(s);

	// Calculate the length

	return (muscle_length + r);
}

void muscle::calculate_muscle_ATP_concentration(double time_step_s)
{
	//! Function updates hs ATP concentration

	// Variables
	double d_heads;					// number of heads per liter
	double cb_number_density;
	double myof_ATP_flux;

	// Code

	if (p_FiberSim_muscle != NULL)
	{
		myof_ATP_flux = p_FiberSim_muscle->p_FiberSim_hs->ATP_flux;

		cb_number_density = p_cmv_model->p_fs_model->m_filament_density *
			(p_FiberSim_muscle->p_FiberSim_hs->m_cbs_per_thick_filament);
	}
	else
	{
		// MyoSim
		printf("Not yet implemented\n");
		exit(1);

		//ATP_flux = p_MyoSim_muscle->p_MyoSim_half_sarcomere->myof_ATP_flux;
		//cb_number_density = p_MyoSim_musclemyofilaments->myof_cb_number_density;
	}

	// First calculate ATP used by myofilaments based on a liter of muscle

	// Calculate density of heads per liter
	d_heads = 0.001 *
		(1.0 - muscle_prop_fibrosis) * muscle_prop_myofilaments *
			cb_number_density *
			(1.0 / (1e-9 * muscle_reference_length));

	muscle_ATP_used_myof_per_liter_per_s = - d_heads * myof_ATP_flux /
		GSL_CONST_NUM_AVOGADRO;

	// Now calculate ATP used by SERCA
	muscle_ATP_used_memb_per_liter_per_s = -(1.0 - muscle_prop_fibrosis) *
		muscle_prop_myofilaments *
		p_membranes->memb_J_uptake *
		p_cmv_model->memb_SERCA_ATPase_excess_Ca_factor /
		p_cmv_model->memb_SERCA_Ca_per_ATP;

	// Now add to get total consumption
	muscle_ATP_used_total_per_liter_per_s = muscle_ATP_used_myof_per_liter_per_s +
		muscle_ATP_used_memb_per_liter_per_s;

	// Now save ATP generated
	muscle_ATP_generated_per_liter_per_s = p_mitochondria->mito_ATP_generated_M_per_liter_per_s;

	// Euler step
	muscle_ATP_concentration = muscle_ATP_concentration +
			((muscle_ATP_used_total_per_liter_per_s + muscle_ATP_generated_per_liter_per_s) *
				time_step_s);
}

void muscle::update_beat_metrics(void)
{
	//! Update beat metrics

	// Variables
	stats_structure* p_stats;

	// Code

	p_stats = new stats_structure;

	if (p_cmv_results_beat->muscle_length_field_index >= 0)
	{
		p_cmv_results_beat->calculate_sub_vector_statistics(
			p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->muscle_length_field_index],
			0, p_parent_hemi_vent->p_parent_cmv_system->beat_t_index,
			p_stats);
	}

	// Tidy up
	delete p_stats;
}

