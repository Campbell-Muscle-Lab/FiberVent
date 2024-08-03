#pragma once

/**
/* @file		reflex_control.h
/* @brief		Header file for a reflex_control object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_model;
class cmv_system;
class cmv_results;
class cmv_options;

class baroreflex;
class circulation;

struct cmv_model_rc_structure;

using namespace::std;

class reflex_control
{
public:
	/**
	 * Constructor
	 */
	reflex_control(baroreflex* set_p_parent_baroreflex, int set_rc_number,
		cmv_model_rc_structure* p_struct);

	/**
	 * Destructor
	 */
	 ~reflex_control(void);

	// Variables

	baroreflex* p_parent_baroreflex;		/**< pointer to the parent baroreflex */

	cmv_system* p_parent_cmv_system;		/**< pointer to parent cmv_system */

	cmv_model* p_cmv_model;					/**< pointer to cmv_model object */

	cmv_results* p_cmv_results_beat;		/**< pointer to cmv_results object
													holding data at full time
													resolution for a single beat */

	cmv_options* p_cmv_options;				/**< pointer to cmv_options object */

	circulation* p_parent_circulation;		/**< pointer the parent circulation */

	// Other variables
	int rc_number;							/**< int with the reflex control number
													subtract 1 to get the index */

	string rc_level;						/**< string with the level of the controlled
													variable */

	string rc_variable;						/**< string with the name of the controlled
													variable */

	double rc_base_value;					/**< double with base value of control */

	double rc_para_factor;					/**< double */
	double rc_symp_factor;					/**< double */

	double rc_para_value;					/**< double with control value under max
													parasympathetic stimulation
													= para_factor * base_value */

	double rc_symp_value;					/**< double with control value under max
													sympathetic stimulation
													= symp_factor * base_value */

	double rc_baro_C;						/**< double with baro_C signal */

	double rc_k_control;

	double rc_k_recov;

	double rc_output;						/**< double with the current control value */

	double* p_controlled_variable;			/**< double to the variable managed
													by the reflex control */

	// Other functions
	void initialise_simulation(void);

	void implement_time_step(double time_step_s);

	void set_controlled_variable(void);

	void extract_digits(string test_string, int digits[], int no_of_digits);

	void calculate_baro_C(double time_step_s);

	void calculate_output(void);
};
