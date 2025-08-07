/**
/* @file		circulation.cpp
/* @brief		Source file for a circulation object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include "circulation.h"

#include "cmv_system.h"
#include "cmv_model.h"
#include "cmv_protocol.h"
#include "cmv_options.h"
#include "cmv_results.h"

#include "hemi_vent.h"
#include "muscle.h"
#include "valve.h"
#include "baroreflex.h"
#include "growth.h"

#include "gsl_math.h"
#include "gsl_errno.h"
#include "gsl_odeiv2.h"

struct stats_structure{
	double mean_value;
	double min_value;
	double max_value;
	double sum;
};

// Constructor
circulation::circulation(cmv_system* set_p_parent_cmv_system = NULL)
{
	//! Constructor

	// Code

	// Initialise

	// Set the pointer to the parent system
	p_parent_cmv_system = set_p_parent_cmv_system;
	p_cmv_model = p_parent_cmv_system->p_cmv_model;

	// Set pointers
	p_cmv_options = p_parent_cmv_system->p_cmv_options;
	p_cmv_results_beat = NULL;

	// Now initialise other objects
	circ_blood_volume = p_cmv_model->circ_blood_volume;

	circ_no_of_compartments = p_cmv_model->circ_no_of_compartments;

	// Set up arrays
	circ_resistance = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_compliance = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_slack_volume = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_inertance = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_static = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_pressure = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_volume = (double*)malloc(circ_no_of_compartments * sizeof(double));
	circ_flow = (double*)malloc(circ_no_of_compartments * sizeof(double));

	// Initialise, noting total slack_volume as we go
	// Start with the compartments at slack volume
	circ_total_slack_volume = 0.0;

	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		circ_resistance[i] = p_cmv_model->circ_resistance[i];
		circ_compliance[i] = p_cmv_model->circ_compliance[i];
		circ_slack_volume[i] = p_cmv_model->circ_slack_volume[i];
		circ_inertance[i] = p_cmv_model->circ_inertance[i];
		circ_static[i] = p_cmv_model->circ_static[i];
		circ_pressure[i] = 0.0;
		circ_volume[i] = circ_slack_volume[i];

		circ_flow[i] = 0.0;

		circ_total_slack_volume = circ_total_slack_volume +
			circ_volume[i];
	}

	// Excess blood goes in veins
	cout << "Blood volume: " << circ_blood_volume << "\n";
	cout << "Total slack volume: " << circ_total_slack_volume << "\n";

	if (circ_blood_volume < circ_total_slack_volume)
	{
		cout << "Error: blood volume is less than total slack volume\n";
		exit(1);
	}
	circ_volume[circ_no_of_compartments - 1] = circ_volume[circ_no_of_compartments - 1] +
		(circ_blood_volume - circ_total_slack_volume);

	// Make a hemi-vent object
	p_hemi_vent = new hemi_vent(this);

	// Set the pointer to the aortic valve
	p_av = p_hemi_vent->p_av;
	p_mv = p_hemi_vent->p_mv;

	// Make a new baroreflex object if it has been defined
	if (!gsl_isnan(p_cmv_model->baro_S))
	{
		p_baroreflex = new baroreflex(this);
	}
	else
	{
		p_baroreflex = NULL;
	}

	// Make a growth object if has been defined
	if (p_cmv_model->no_of_gc_controls > 0)
	{
		p_growth = new growth(this);
	}
	else
	{
		p_growth = NULL;
	}

	// Check for filling control
	if (!gsl_isnan(p_cmv_model->filling_venous_pressure))
	{
		filling_control = true;
		filling_control_venous_pressure = p_cmv_model->filling_venous_pressure;
	}
	else
	{
		filling_control = false;
		filling_control_venous_pressure = GSL_NAN;
	}

}

// Destructor
circulation::~circulation(void)
{
	//! Destructor

	// Code

	delete p_hemi_vent;

	if (p_baroreflex != NULL)
		delete p_baroreflex;

	if (p_growth != NULL)
		delete p_growth;

	free(circ_resistance);
	free(circ_compliance);
	free(circ_slack_volume);
	free(circ_inertance);
	free(circ_static);
	free(circ_pressure);
	free(circ_volume);
	free(circ_flow);
}

// Other functions

void circulation::initialise_simulation(void)
{
	//! Code initialises simulation
	
	// Variables

	// Code

	// Set the options
	p_cmv_options = p_parent_cmv_system->p_cmv_options;

	// Set the results
	p_cmv_results_beat = p_parent_cmv_system->p_cmv_results_beat;

	// Set the protocol
	p_cmv_protocol = p_parent_cmv_system->p_cmv_protocol;

	printf("circulation check 1\n");

	// Now handle daughter objects
	p_hemi_vent->initialise_simulation();

	printf("circulation check 2\n");

	if (p_baroreflex != NULL)
		p_baroreflex->initialise_simulation();

	printf("circulation check\n");

	if (p_growth != NULL)
		p_growth->initialise_simulation();

	// Add data fields
	p_cmv_results_beat->add_results_field("circ_blood_volume", &circ_blood_volume);

	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		string label = string("pressure_") + to_string(i);
		p_cmv_results_beat->add_results_field(label, &circ_pressure[i]);
	}

	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		string label = string("volume_") + to_string(i);
		p_cmv_results_beat->add_results_field(label, &circ_volume[i]);
	}

	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		string label = string("flow_") + to_string(i);
		p_cmv_results_beat->add_results_field(label, &circ_flow[i]);
	}

	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		string label = string("static_") + to_string(i);
		p_cmv_results_beat->add_results_field(label, &circ_static[i]);
	}
}

// This function is not a member of the circulation class but is used to interace
// with the GSL ODE system. It must appear before the circulation class members
// that calls it, and communicates with the circulation class through a pointer to
// the class object

int circ_vol_derivs(double t, const double y[], double f[], void* params)
{
	// Function sets dV/dt for the compartments

	// Variables
	(void)(t);						// Prevents warning for unused variable

	circulation* p_circ = (circulation*)params;
									// Pointer to circulation

	double* flow_calc = (double*)malloc(p_circ->circ_no_of_compartments * sizeof(double));
									// array of doubles to hold flows

	// Code

	// Zero flows
	for (int i = 0; i < p_circ->circ_no_of_compartments; i++)
	{
		flow_calc[i] = 0.0;
	}
	
	// Calculate the flows between the compartments
	p_circ->calculate_flows(y, flow_calc);

	// Now adjust volumes
	//! if compartments are
	//! [0] - [1] - [2] - [3] - .... [n-1]
	//! circ_flow[i] is flow from compartment [i-1] to compartment [i] through resistance i
	//! rate of change of [i] is flow[i-1] - flow[i]
	
	for (int i = 0; i < (p_circ->circ_no_of_compartments-1) ; i++)
	{
		f[i] = flow_calc[i] - flow_calc[i+1];
	}

	f[p_circ->circ_no_of_compartments-1] = flow_calc[p_circ->circ_no_of_compartments - 1] - flow_calc[0];

	// Tidy up
	free(flow_calc);

	// Return
	return GSL_SUCCESS;
}

bool circulation::implement_time_step(double time_step_s)
{
	//! Code advances by 1 time-step

	// Variables
	bool new_beat = false;

	int status;

	double t_start_s = 0.0;
	double t_stop_s = time_step_s;
	
	double eps_abs = 1e-6;
	double eps_rel = 1e-6;

	double* vol_calc = (double*)malloc(circ_no_of_compartments * sizeof(double));

	// Code

	// Update the hemi_vent object, which includes
	// updating the daughter objects
	new_beat = p_hemi_vent->implement_time_step(time_step_s);

	// Fill the volume array
	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		vol_calc[i] = circ_volume[i];
	}

	// Calculate pressures
	calculate_pressures(vol_calc, circ_pressure, time_step_s);

	// Now adjust the compartment volumes by integrating flows.
	gsl_odeiv2_system sys =
		{ circ_vol_derivs, NULL, (size_t)circ_no_of_compartments, this };

	gsl_odeiv2_driver* d =
		gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45,
			0.5 * time_step_s, eps_abs, eps_rel);

	status = gsl_odeiv2_driver_apply(d, &t_start_s, t_stop_s, vol_calc);

	gsl_odeiv2_driver_free(d);

	// Unpack the arrays
	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		circ_volume[i] = vol_calc[i];
	}

	// Update the hemi_vent with the new volume
	p_hemi_vent->update_chamber_volume(circ_volume[0], time_step_s);

	p_hemi_vent->vent_chamber_radius = 
		p_hemi_vent->return_internal_radius_for_chamber_volume(circ_volume[0]);

	p_hemi_vent->vent_chamber_height =
		p_hemi_vent->return_chamber_height(p_hemi_vent->vent_chamber_radius);

	// Make sure total volume remains constant
	double holder = 0.0;
	for (int i = 0; i < circ_no_of_compartments; i++)
	{
		holder = holder + circ_volume[i];
	}

	// Any adjustment goes in veins
	double adjustment = (circ_blood_volume - holder);
	circ_volume[circ_no_of_compartments - 1] = 
		circ_volume[circ_no_of_compartments - 1] + adjustment;

	// Warning
	if (fabs(adjustment) > 1e-4)
	{
		cout << "Blood volume mismatch\n";
		exit(1);
	}

	// Update data flows for data
	calculate_flows(circ_volume, circ_flow);

	// Update the baroreflex, which includes updating the daughter objects
	if (p_baroreflex != NULL)
	{
		p_baroreflex->baro_active = p_cmv_protocol->return_activation("baroreflex",
			p_parent_cmv_system->cum_time_s);
		p_baroreflex->implement_time_step(time_step_s);
	}

	// Update the growth, which includes updating the daughter objects
	if (p_growth != NULL)
	{
		p_growth->growth_active = p_cmv_protocol->return_activation("growth",
			p_parent_cmv_system->cum_time_s);
		p_growth->implement_time_step(time_step_s, new_beat);
	}

	// Check for active filling control
	if (filling_control)
	{
		int control_active = p_cmv_protocol->return_activation("filling_control",
			p_parent_cmv_system->cum_time_s);

		if (control_active)
		{
			int venous_compartment = p_cmv_model->circ_no_of_compartments - 1;

			double required_volume = circ_slack_volume[venous_compartment] +
				(filling_control_venous_pressure * circ_compliance[venous_compartment]);

			double delta_volume = required_volume - circ_volume[venous_compartment];

			// Adjust
			circ_volume[venous_compartment] = circ_volume[venous_compartment] + delta_volume;
			circ_blood_volume = circ_blood_volume + delta_volume;
		}
	}

	// Tidy up
	free(vol_calc);

	// Return
	return (new_beat);
}

void circulation::calculate_pressures(const double v[], double p[], double time_step_s)
{
	//! Function calculates pressures
	
	// Code
	
	p[0] = p_hemi_vent->return_pressure_for_chamber_volume(v[0], time_step_s);

	// Calculate the other pressures
	for (int i = 1; i < circ_no_of_compartments; i++)
	{
		p[i] = ((v[i] - circ_slack_volume[i]) / circ_compliance[i]) + circ_static[i];
	}
}

void circulation::calculate_flows(const double v[], double flow[])
{
	//! Function sets the values of circ_flows[] based on an
	//! array of compartment volumes
	//! if compartments are
	//! [0] - [1] - [2] - [3] - .... [n-1]
	//! circ_flow[i] is flow from compartment [i-1] to compartment [i] through resistance i

	// Variables
	double p_diff;

	double time_step_s = 0.001;

	// Code

	// Calculate the flows
	// These are flows from the aorta through to the veins
	for (int i = 2; i < circ_no_of_compartments; i++)
	{
		p_diff = (circ_pressure[i - 1] - circ_static[i-1]) - 
			(circ_pressure[i] - circ_static[i]);
		flow[i] = p_diff / circ_resistance[i];
	}

	// Special case for flow through mitral valve
	p_diff = (circ_pressure[circ_no_of_compartments - 1] - circ_static[circ_no_of_compartments -1 ]) - 
				(circ_pressure[0] - circ_static[0]);
	flow[0] = fabs(p_mv->valve_pos) * p_diff / circ_resistance[0];

	// Special case for flow through aortic valve
	p_diff = (circ_pressure[0] - circ_static[0]) - 
				(circ_pressure[1] - circ_static[1]);
	flow[1] = fabs(p_av->valve_pos) * p_diff / circ_resistance[1];
}

void circulation::update_beat_metrics(void)
{
	//! Update beat metrics in daughter objects

	// Variables
	stats_structure* p_stats;
	stats_structure* p_stats_2;

	// Code
	
	p_stats = new stats_structure;
	p_stats_2 = new stats_structure;


	if (p_cmv_results_beat->pressure_arteries_field_index >= 0)
	{
		p_cmv_results_beat->calculate_sub_vector_statistics(
			p_cmv_results_beat->gsl_results_vectors[p_cmv_results_beat->pressure_arteries_field_index],
			0, p_parent_cmv_system->beat_t_index,
			p_stats);

		cout << "Arterial pressure: " << p_stats->max_value << " / " << p_stats->min_value << "\n";
	}

	p_hemi_vent->update_beat_metrics();

	// Tidy up
	delete p_stats;
	delete p_stats_2;
}

