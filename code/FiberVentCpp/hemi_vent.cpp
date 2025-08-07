/**
/* @file		hemi_vent.cpp
/* @brief		Source file for a hemi_vent object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include "math.h"

#include "cmv_system.h"
#include "cmv_model.h"
#include "cmv_results.h"
#include "cmv_options.h"

#include "hemi_vent.h"
#include "circulation.h"
#include "valve.h"

#include "muscle.h"

#include "heart_rate.h"

#include "FiberSim_muscle.h"
#include "FiberSim_half_sarcomere.h"
#include "FiberSim_series_component.h"

#include "gsl_errno.h"
#include "gsl_roots.h"
#include "gsl_math.h"
#include "gsl_const_mksa.h"
#include "gsl_const_num.h"

struct stats_structure {
	double mean_value;
	double min_value;
	double max_value;
	double sum;
};

// Constructor
hemi_vent::hemi_vent(circulation* set_p_parent_circulation)
{
	// Initialise

	// Code
	std::cout << "hemi_vent_constructor()\n";

	// Set pointers
	p_parent_circulation = set_p_parent_circulation;
	p_cmv_model = p_parent_circulation->p_cmv_model;
	p_parent_cmv_system = p_parent_circulation->p_parent_cmv_system;

	// Initialise with safe options
	p_cmv_results_beat = NULL;
	p_cmv_options = p_parent_circulation->p_cmv_options;

	vent_wall_density = p_cmv_model->vent_wall_density;
	vent_wall_volume = p_cmv_model->vent_wall_volume;

//	vent_chamber_height = p_cmv_model->vent_chamber_height;
	vent_z_scale = p_cmv_model->vent_z_scale;
	vent_z_exp = p_cmv_model->vent_z_exp;

	vent_stroke_work_J = GSL_NAN;
	vent_stroke_energy_used_J = GSL_NAN;
	vent_efficiency = GSL_NAN;
	vent_ejection_fraction = GSL_NAN;
	vent_ATP_used_per_s = 0.0;
	vent_cardiac_output = GSL_NAN;
	vent_stroke_volume = GSL_NAN;
	vent_cardiac_output = GSL_NAN;

	// Initialise child half-sarcomere
	p_muscle = new muscle(this);

	// Initialise aortic valve
	p_av = new valve(this, p_cmv_model->p_av);

	// Initialise mitral valve
	p_mv = new valve(this, p_cmv_model->p_mv);
}

// Destructor
hemi_vent::~hemi_vent(void)
{
	//! hemi_vent destructor

	// Tidy up
	delete p_muscle;
	delete p_av;
	delete p_mv;
}

// Other functions
void hemi_vent::initialise_simulation(void)
{
	//! Code initialises simulation

	// Variables

	// Initialise options
	p_cmv_options = p_parent_circulation->p_cmv_options;

	if (p_cmv_options->hv_thick_wall_approximation == "True")
		vent_thick_wall_multiplier = 1.0;
	else
		vent_thick_wall_multiplier = 0.0;

	// Now add in the results
	p_cmv_results_beat = p_parent_circulation->p_cmv_results_beat;

	// And now daughter objects

	p_av->initialise_simulation();
	p_mv->initialise_simulation();

	printf("hemi_vent: initialise_check\n");

	p_muscle->initialise_simulation();

	// Deduce the slack circumference of the ventricle and
	// set the number of half-sarcomeres
	// The p_hs->initialisation set hs_length so that stress was 0

	printf("Slack_m_length: %g\n", p_muscle->p_FiberSim_muscle->fs_m_length);
	printf("Slack_hs_length: %g\n", p_muscle->p_FiberSim_muscle->p_FiberSim_hs->hs_length);
	printf("Slack_sc_extension: %g\n", p_muscle->p_FiberSim_muscle->p_FiberSim_sc->sc_extension);
	//exit(1);



	vent_circumference = return_lv_circumference_for_chamber_volume(
		p_parent_circulation->circ_slack_volume[0]);

	vent_n_hs = 1e9 * vent_circumference / p_muscle->muscle_length;

	double vent_diam = vent_circumference / 3.14159;

	cout << "\n\nvent_circum: " << vent_circumference << " vent_z_scale: " << vent_z_scale <<
		" vent_volume: " << p_parent_circulation->circ_slack_volume[0] <<
		" vent_diam:" << vent_diam << "\n\n\n";

	// Add fields
	p_cmv_results_beat->add_results_field("vent_wall_volume", &vent_wall_volume);
	p_cmv_results_beat->add_results_field("vent_wall_thickness", &vent_wall_thickness);
	p_cmv_results_beat->add_results_field("vent_chamber_radius", &vent_chamber_radius);
	p_cmv_results_beat->add_results_field("vent_chamber_height", &vent_chamber_height);
	p_cmv_results_beat->add_results_field("vent_n_hs", &vent_n_hs);
	p_cmv_results_beat->add_results_field("vent_stroke_work_J", &vent_stroke_work_J);
	p_cmv_results_beat->add_results_field("vent_stroke_energy_used_J", &vent_stroke_energy_used_J);
	p_cmv_results_beat->add_results_field("vent_efficiency", &vent_efficiency);
	p_cmv_results_beat->add_results_field("vent_ejection_fraction", &vent_ejection_fraction);
	p_cmv_results_beat->add_results_field("vent_ATP_used_per_s", &vent_ATP_used_per_s);
	p_cmv_results_beat->add_results_field("vent_stroke_volume", &vent_stroke_volume);
	p_cmv_results_beat->add_results_field("vent_cardiac_output", &vent_cardiac_output);
}

bool hemi_vent::implement_time_step(double time_step_s)
{
	//! Implements time-step

	// Variables
	bool new_beat = false;

	// Code
	p_av->implement_time_step(time_step_s);
	p_mv->implement_time_step(time_step_s);

	new_beat = p_muscle->implement_time_step(time_step_s);

	// Calculate energy used per s
	calculate_vent_ATP_used_per_s();

	// Return
	return (new_beat);
}

struct gsl_thickness_root_params
{
	hemi_vent* p_hemi_vent_temp;
	double chamber_volume;
	double wall_volume_target;
};

double hemi_vent_thickness_root_finder(double x, void* params)
{
	//! Function for the root finder
	
	// Variables
	double external_volume;
	double inner_radius;
	double inner_height;
	double calculated_wall_volume;
	double volume_difference;

	// Code

	// Unpack the pointer
	struct gsl_thickness_root_params* p = (struct gsl_thickness_root_params*)(params);

	// Get the chamber dimensions

	inner_radius = p->p_hemi_vent_temp->return_internal_radius_for_chamber_volume(p->chamber_volume);

	inner_height = p->p_hemi_vent_temp->return_chamber_height(inner_radius);
	
	// Calculate external_volume
	external_volume = 1000.0 * (2.0 / 3.0) * M_PI * pow(inner_radius + x, 2.0) * (inner_height + x);

	// Estimated wall volume
	calculated_wall_volume = external_volume - p->chamber_volume;

	// How far off are we
	volume_difference = calculated_wall_volume - p->wall_volume_target;

	// Return the difference
	return volume_difference;
}

double hemi_vent::wall_thickness_root_finder(double cv)
{
	// Function uses an iterative technique to find the thickness of an elliptical ventricle

	// Variables
	double x_lo = 0;
	double x_hi = 0.05;
	double x;

	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;

	gsl_function F;
	struct gsl_thickness_root_params params = { this, cv, vent_wall_volume };

	int status;
	int iter = 0;
	int max_iter = 100;

	double epsabs = 1e-5;
	double epsrel = 1e-5;

	// Code

	F.function = &hemi_vent_thickness_root_finder;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		x = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi, epsabs, epsrel);
	} while ((status == GSL_CONTINUE) && (iter < max_iter));

	gsl_root_fsolver_free(s);

	// Return thickness
	return x;
}

double hemi_vent::return_wall_thickness_for_chamber_volume(double cv)
{
	//! Code sets object value of wall thickness
	//! Volumes are in liters, dimensions are in m
	
	// Variables

	double thickness;

	// Code

	thickness = wall_thickness_root_finder(cv);

	/*
	// Variables
	double internal_r;

	double h;

	double thickness;

	// Code
	internal_r = return_internal_radius_for_chamber_volume(cv);

	h = return_chamber_height(internal_r);

	thickness = pow((0.001 * (cv + vent_wall_volume)) / ((2.0 / 3.0) * h * M_PI), (1.0 / 2.0)) -
		internal_r;

	*/

	return thickness;
}

double hemi_vent::return_lv_circumference_for_chamber_volume(double cv)
{
	//! Code returns lv circumference for chamber volume
	//! Volumes are in liters, dimensions are in m
	//! based on 2 * pi * (internal_r + wall_thickness)

	// Variables
	double thickness;
	double lv_circum;

	// Code
	thickness = return_wall_thickness_for_chamber_volume(cv);

	lv_circum = 2.0 * M_PI *
		(return_internal_radius_for_chamber_volume(cv) + (0.5 * thickness));

	return lv_circum;
}

double hemi_vent::return_internal_radius_for_chamber_volume(double cv)
{
	//! Returns internal radius in meters for a given chamber volume in liters

	// Variables
	double r;

	double rel_hsl;

	// Code

	if (cv < 0.0)
		cv = 0.0;

	rel_hsl = (p_muscle->muscle_length / p_muscle->muscle_reference_length);

	r = pow(((3.0 * 0.001 * cv) /
		(2.0 * M_PI * pow(rel_hsl, vent_z_exp) * vent_z_scale)), (1.0 / 3.0));

	return r;
}

double hemi_vent::return_pressure_for_chamber_volume(double cv, double time_step_s)
{
	//! Code returns pressure for a given chamber volume

	// Variables
	double new_lv_circumference;
	double new_muscle_length;
	double delta_muscle_length;
	double new_stress;
	double internal_r;
	double P_in_Pascals;
	double P_in_mmHg;

	// Code
	new_lv_circumference = return_lv_circumference_for_chamber_volume(cv);

	new_muscle_length = 1.0e9 * new_lv_circumference / vent_n_hs;

	delta_muscle_length = new_muscle_length - p_muscle->muscle_length;

	// Deduce stress for the new hs_length
/*	double xt = 0.0;
	for (int i = 0; i < 10; i++)
	{
		double ns;
		ns = p_hs->return_wall_stress_after_delta_hsl(xt, time_step_s);
		printf("xt: %g    ns: %g\n", xt, ns);
		xt = xt + 1.0;
	}
*/

	new_stress = p_muscle->return_wall_stress_after_test_delta_ml(delta_muscle_length, time_step_s);

	new_stress = GSL_MAX(-1000.0, new_stress);

	internal_r = return_internal_radius_for_chamber_volume(cv);

	vent_wall_thickness = return_wall_thickness_for_chamber_volume(cv);

	// Pressure from Laplace's law
	// https://www.annalsthoracicsurgery.org/action/showPdf?pii=S0003-4975%2810%2901981-8

	if (internal_r < 1e-6)
	{
		P_in_Pascals = 0.0;
		cout << "Hemi_vent, internal_r ~= 0.0 problem\n";
		cout << "delta_muscle_length: " << delta_muscle_length << "\n";
		cout << "new_stress: " << new_stress << "\n";
		cout << "wall_thickness: " << vent_wall_thickness << "\n";
	}
	else
	{
		P_in_Pascals = (new_stress * vent_wall_thickness *
			(2.0 + (vent_thick_wall_multiplier * (vent_wall_thickness / internal_r)))) /
			internal_r;
	}

	P_in_mmHg = P_in_Pascals / (0.001 * GSL_CONST_MKSA_METER_OF_MERCURY);

	return P_in_mmHg;
}

void hemi_vent::update_chamber_volume(double new_volume, double time_step_s)
{
	//! Function updates the chamber volume

	// Variables
	double new_circumference;
	double delta_circumference;
	double delta_muscle_length;

	// Code

	// Update

	new_circumference = return_lv_circumference_for_chamber_volume(new_volume);

	delta_circumference = new_circumference - vent_circumference;

	delta_muscle_length = 1e9 * delta_circumference / vent_n_hs;

	p_muscle->change_muscle_length(delta_muscle_length, time_step_s);

	vent_circumference = new_circumference;
}

void hemi_vent::update_beat_metrics(void)
{
	//! Code updates beat metrics

	// Variables
	double cardiac_cycle_s;

	// Update beat values
	vent_stroke_work_J = p_cmv_results_beat->return_stroke_work(0, p_parent_cmv_system->beat_t_index);
	vent_stroke_energy_used_J = p_cmv_results_beat->return_energy_used(0, p_parent_cmv_system->beat_t_index);
	vent_efficiency = -vent_stroke_work_J / vent_stroke_energy_used_J;

	// Calculate the ejection fraction
	stats_structure* p_v_stats = new stats_structure;

	// Calculate stroke volume
	p_cmv_results_beat->calculate_sub_vector_statistics(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->volume_vent_field_index],
		0, p_parent_cmv_system->beat_t_index, p_v_stats);

	vent_stroke_volume = p_v_stats->max_value - p_v_stats->min_value;

	vent_ejection_fraction = vent_stroke_volume / p_v_stats->max_value;

	// Calculate period of cardiac cycle to get cardiac output
	cardiac_cycle_s = p_parent_cmv_system->cum_time_s -
		gsl_vector_get(p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->time_field_index], 0);

	if (cardiac_cycle_s > 0.0)
	{
		vent_cardiac_output = 60.0 * vent_stroke_volume / cardiac_cycle_s;
	}

	// Backfill results
	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->hr_heart_rate_bpm_field_index],
		(60.0 / cardiac_cycle_s), 0, p_parent_cmv_system->beat_t_index);
	
	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->vent_stroke_work_field_index],
		vent_stroke_work_J, 0, p_parent_cmv_system->beat_t_index);

	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->vent_stroke_energy_used_field_index],
		vent_stroke_energy_used_J, 0, p_parent_cmv_system->beat_t_index);

	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->vent_efficiency_field_index],
		vent_efficiency, 0, p_parent_cmv_system->beat_t_index);

	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->vent_ejection_fraction_field_index],
		vent_ejection_fraction, 0, p_parent_cmv_system->beat_t_index);

	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->vent_stroke_volume_field_index],
		vent_stroke_volume, 0, p_parent_cmv_system->beat_t_index);

	p_cmv_results_beat->backfill_beat_data(
		p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->vent_cardiac_output_field_index],
		vent_cardiac_output, 0, p_parent_cmv_system->beat_t_index);

	// Update hs metrics
	p_muscle->update_beat_metrics();
}

double hemi_vent::return_chamber_height(double r)
{
	//! Function returns chamber height
	
	// Variables
	
	double h;

	// Code

	h = r * vent_z_scale * pow((p_muscle->muscle_length / p_muscle->muscle_reference_length), vent_z_exp);

	return h;
}

void hemi_vent::calculate_vent_ATP_used_per_s()
{
	//! Function updates vent_ATP_used_per_s

	// Variables
	vent_ATP_used_per_s = vent_wall_volume *
		p_muscle->muscle_ATP_used_total_per_liter_per_s;
}