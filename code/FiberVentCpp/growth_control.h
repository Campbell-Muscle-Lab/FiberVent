#pragma once

/**
/* @file		growth_control.h
/* @brief		Header file for a growth_control object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include "gsl_vector.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_model;
class cmv_system;
class cmv_results;
class cmv_options;

class growth;
class circulation;

struct cmv_model_gc_structure;

using namespace::std;

class growth_control
{
public:
	/**
	 * Constructor
	 */
	growth_control(growth* set_p_parent_growth, int set_rc_number,
		cmv_model_gc_structure* p_struct);

	/**
	 * Destructor
	 */
	 ~growth_control(void);

	// Variables

	growth* p_parent_growth;				/**< pointer to the parent growth object */

	cmv_system* p_parent_cmv_system;		/**< pointer to parent cmv_system */

	cmv_model* p_cmv_model;					/**< pointer to cmv_model object */

	cmv_results* p_cmv_results_beat;		/**< pointer to cmv_results object
													holding data at full time
													resolution for one beat */

	cmv_options* p_cmv_options;				/**< pointer to cmv_options object */

	circulation* p_parent_circulation;		/**< pointer the parent circulation */

	// Other variables
	int gc_number;							/**< int with the reflex control number
													subtract 1 to get the index */

	string gc_type;							/**< string with the growth control type
													eccentric, or concentric */

	string gc_level;						/**< string with the level of the controlled
													variable */

	string gc_signal;						/**< string with the name of the feedback
													signal */

	double gc_set_point;					/**< double with base value of control */

	double gc_prop_gain;					/**< double with the gain of the prop control */

	double gc_deriv_gain;					/**< double with the gain of the deriv control */

	double gc_max_rate;						/**< double with the max rate of growth */

	double gc_output;						/**< double with the output of the control */

	double* gc_p_signal;					/**< pointer to the variable driving the
													growth control */

	bool gc_signal_assigned;				/**< bool defining whether the signal
													has been defined correctly */

	double gc_control_period_s;				/**< time in s over which to average
													the control signal and its
													rate of change */

	int gc_control_points;					/**< no of time-steps in control_period_s */

	double * gc_control_t;					/**< pointer to an array of doubles
													holding control time */

	double * gc_control_y;					/**< pointer to an array of doubles
													holding control variable */

	double gc_prop_signal;					/**< double with proportional signal */

	double gc_deriv_signal;					/**< double with derivitive signal */

	double gc_prop_output;					/**< double with proportional output */

	double gc_deriv_output;					/**< double with derivitive output */

	// Other functions
	void initialise_simulation(void);

	void implement_time_step(double time_step_s, bool new_beat);

	void set_gc_p_signal(void);

//	void extract_digits(string test_string, int digits[], int no_of_digits);

//	void calculate_baro_C(double time_step_s);

	void calculate_output(void);

	void calculate_mean(void);

	void calculate_slope(void);
};
