/**
/* @file		reflex_control.cpp
/* @brief		Source file for a reflex control object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include "math.h"

#include <iostream>
#include <regex>

#include "reflex_control.h"

#include "cmv_model.h"
#include "cmv_system.h"
#include "cmv_results.h"
#include "cmv_options.h"

#include "baroreflex.h"
#include "circulation.h"
#include "hemi_vent.h"

#include "muscle.h"
#include "heart_rate.h"
#include "membranes.h"
#include "muscle.h"

#include "kinetic_scheme.h"
#include "transition.h"

#include "FiberSim_muscle.h"
#include "FiberSim_half_sarcomere.h"
#include "FiberSim_kinetic_scheme.h"
#include "FiberSim_m_state.h"
#include "FiberSim_transition.h"

#include "gsl_math.h"

struct cmv_model_rc_structure {
	string level;
	string variable;
	double k_control;
	double k_recov;
	double para_factor;
	double symp_factor;
};

// Constructor
reflex_control::reflex_control(baroreflex* set_p_parent_baroreflex, int set_rc_number,
	cmv_model_rc_structure* p_struct)
{
	// Initialise

	// Code
	
	// Set pointers
	p_parent_baroreflex = set_p_parent_baroreflex;

	p_parent_circulation = p_parent_baroreflex->p_parent_circulation;
	p_parent_cmv_system = p_parent_circulation->p_parent_cmv_system;
	p_cmv_model = p_parent_circulation->p_cmv_model;

	// Initialise with safe options
	p_cmv_results_beat = NULL;
	p_cmv_options = NULL;

	// Other variables
	rc_number = set_rc_number;
	rc_baro_C = 0.5;

	// Others come from the structure
	rc_level = p_struct->level;
	rc_variable = p_struct->variable;
	rc_k_control = p_struct->k_control;
	rc_k_recov = p_struct->k_recov;
	rc_para_factor = p_struct->para_factor;
	rc_symp_factor = p_struct->symp_factor;
}

// Destructor
reflex_control::~reflex_control(void)
{
	// Tidy up
}

// Other functions
void reflex_control::initialise_simulation(void)
{
	//! Code initialises simulation
	
	// Variables
	bool reflex_assigned = false;
	string temp_string;
	
	// Initialise options
	p_cmv_options = p_parent_circulation->p_cmv_options;

	// Now add in the results
	p_cmv_results_beat = p_parent_circulation->p_cmv_results_beat;

	// Find the controlled variable
	set_controlled_variable();

	// Add fields
	temp_string = "rc_" + to_string(rc_number) + "_C";
	p_cmv_results_beat->add_results_field(temp_string, &rc_baro_C);

	temp_string = "rc_" + rc_level + "_" + rc_variable;
	p_cmv_results_beat->add_results_field(temp_string, p_controlled_variable);
}

void reflex_control::implement_time_step(double time_step_s)
{
	//! Implements time-step
	
	// Variables

	// Code

	// Calculate the baro_C signal, and then the output via the mapping function
	calculate_baro_C(time_step_s);

	calculate_output();

	// Now update controlled value
	*p_controlled_variable = rc_output;
}

void reflex_control::calculate_baro_C(double time_step_s)
{
	//! Function calculates the B_c signal
	
	// Variables
	double delta;

	double baro_active = p_parent_baroreflex->baro_active;
	double baro_B = p_parent_baroreflex->baro_B;

	// Code

	// Caculate the derivitive
	if (baro_active > 0.0)
	{
		if (baro_B > 0.5)
		{
			delta = rc_k_control * (baro_B - 0.5) * (1.0 - rc_baro_C);
		}
		else
		{
			delta = rc_k_control * (baro_B - 0.5) * rc_baro_C;
		}
	}
	else
	{
		delta = -rc_k_control * (rc_baro_C - 0.5);
	}

	// Update over the time-step
	rc_baro_C = rc_baro_C + (delta * time_step_s);

	// Limit
	rc_baro_C = GSL_MIN(rc_baro_C, 1.0);
	rc_baro_C = GSL_MAX(rc_baro_C, 0.0);
}

void reflex_control::extract_digits(string test_string, int digits[], int no_of_digits)
{
	//! Function fills digits array with [0-9] extracted from string
	//! See https://en.cppreference.com/w/cpp/regex/regex_iterator
	
	// Variables
	int counter;

	// Code
	regex ex("[0-9]");

	auto digits_begin = sregex_iterator(test_string.begin(), test_string.end(), ex);
	auto digits_end = sregex_iterator();

	counter = 0;
	for (regex_iterator i = digits_begin; i != digits_end; i++)
	{
		smatch match = *i;
		digits[counter] = atoi((match.str().c_str()));
		counter = counter + 1;
	}
}

void reflex_control::set_controlled_variable(void)
{
	//! Code sets the pointer for p_controlled_variable to the appropriate double

	// Variables
	bool reflex_assigned = false;

	// Code

		// Find the variable this object is controlling
	if (rc_level == "heart_rate")
	{
		if (rc_variable == "t_RR_interval_s")
		{
			p_controlled_variable = &p_parent_circulation->p_hemi_vent->p_muscle->p_heart_rate->hr_t_RR_interval_s;
			reflex_assigned = true;
		}
	}

	if (rc_level == "membranes")
	{
		if (rc_variable == "k_active")
		{
			p_controlled_variable = &p_parent_circulation->p_hemi_vent->p_muscle->p_membranes->memb_k_active;
			reflex_assigned = true;
		}

		if (rc_variable == "k_serca")
		{
			p_controlled_variable = &p_parent_circulation->p_hemi_vent->p_muscle->p_membranes->memb_k_serca;
			reflex_assigned = true;
		}
	}

	if (rc_level == "FiberSim_half_sarcomere")
	{
		FiberSim_half_sarcomere* p_FiberSim_hs = p_parent_circulation->p_hemi_vent->p_muscle->
			p_FiberSim_muscle->p_FiberSim_hs;

		if (rc_variable == "k_on")
		{
			p_controlled_variable = &p_FiberSim_hs->a_k_on;
			reflex_assigned = true;
		}

		if (rc_variable.rfind("m_state", 0) == 0)
		{
			int no_of_digits = 3;
			int digits[3];
			int state_index;
			int transition_index;
			int parameter_index;

			for (int i = 0; i < no_of_digits; i++)
				digits[i] = 0;

			extract_digits(rc_variable, digits, 3);

			state_index = digits[0] - 1;
			transition_index = digits[1] - 1;
			parameter_index = digits[2] - 1;

			// This is tricky because the variable is stored in a gsl_vector
			gsl_vector* p_gsl_v = p_FiberSim_hs->p_fs_model->p_m_scheme[0]->p_m_states[state_index]->p_transitions[transition_index]->rate_parameters;
			p_controlled_variable = gsl_vector_ptr(p_gsl_v, parameter_index);

			reflex_assigned = true;
		}
	}

	/*
	if (rc_level == "myofilaments")
	{
		if (rc_variable == "k_on")
		{
			p_controlled_variable = &p_parent_circulation->p_hemi_vent->p_hs->p_myofilaments->myof_a_k_on;
			reflex_assigned = true;
		}

		if (rc_variable.rfind("m_state", 0) == 0)
		{
			int no_of_digits = 3;
			int digits[3];
			int state_index;
			int transition_index;
			int parameter_index;

			for (int i = 0; i < no_of_digits; i++)
				digits[i] = 0;

			extract_digits(rc_variable, digits, 3);

			state_index = digits[0] - 1;
			transition_index = digits[1] - 1;
			parameter_index = digits[2] - 1;

			// This is tricky because the variable is stored in a gsl_vector
			printf("Reflex_control:: This needs to be fixed\n");
			exit(1);
			gsl_vector* p_gsl_v = p_parent_circulation->p_hemi_vent->p_hs->p_myofilaments->p_m_scheme->p_m_states[state_index]->p_transitions[transition_index]->rate_parameters;
			p_controlled_variable = p_gsl_v->data + (parameter_index) * sizeof(p_gsl_v->stride);

			reflex_assigned = true;
		}
	}
	*/

	if (rc_level == "circulation")
	{
		if (rc_variable.rfind("resistance", 0) == 0)
		{
			int no_of_digits = 1;
			int digits[1];
			int compartment_index;

			for (int i = 0; i < no_of_digits; i++)
				digits[i] = 0;

			extract_digits(rc_variable, digits, 1);

			compartment_index = digits[0] - 1;

			p_controlled_variable = &p_parent_circulation->circ_resistance[compartment_index];

			reflex_assigned = true;
		}

		if (rc_variable.rfind("compliance", 0) == 0)
		{
			int no_of_digits = 1;
			int digits[1];
			int compartment_index;

			for (int i = 0; i < no_of_digits; i++)
				digits[i] = 0;

			extract_digits(rc_variable, digits, 1);

			compartment_index = digits[0] - 1;

			p_controlled_variable = &p_parent_circulation->circ_compliance[compartment_index];

			reflex_assigned = true;
		}
	}

	// Assign
	rc_base_value = *p_controlled_variable;
	rc_symp_value = rc_symp_factor * rc_base_value;
	rc_para_value = rc_para_factor * rc_base_value;

	if (reflex_assigned == false)
	{
		cout << "Reflex control " << rc_level << ", " << rc_variable << " not assigned\n";
		exit(1);
	}
	else
	{
		cout << "Reflex control " << rc_level << ", " << rc_variable <<
			" base: " << rc_base_value << " para: " << rc_para_value << " symp: " << rc_symp_value << "\n";
	}
}

void reflex_control::calculate_output(void)
{
	//! Code sets the output value

	// Variables
	double m;

	// Code
	if (rc_baro_C >= 0.5)
	{
		m = (rc_symp_value - rc_base_value) / 0.5;
	}
	else
	{
		m = (rc_base_value - rc_para_value) / 0.5;
	}

	rc_output = rc_base_value + m * (rc_baro_C - 0.5);
}