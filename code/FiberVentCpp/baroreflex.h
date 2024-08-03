#pragma once

/**
/* @file		baroreflex.h
/* @brief		Header file for a baroreflex object
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

class circulation;
class reflex_control;

class baroreflex
{
public:
	/**
	 * Constructor
	 */
	baroreflex(circulation* set_p_parent_circulation);

	/**
	 * Destructor
	 */
	 ~baroreflex(void);

	 // Variables

	cmv_system* p_parent_cmv_system;		/**< pointer to parent cmv_system */

	cmv_model* p_cmv_model;					/**< pointer to cmv_model object */

	cmv_results* p_cmv_results_beat;		/**< pointer to cmv_results object
													holding data at full time
													resolution for one beat */

	cmv_options* p_cmv_options;				/**< pointer to cmv_options object */

	circulation* p_parent_circulation;		/**< pointer the parent circulation */

	reflex_control* p_rc[MAX_NO_OF_REFLEX_CONTROLS];
											/**< array of pointers to reflex controls */

	int no_of_reflex_controls;				/**< int with the number of reflex controls */

	// Other variables
	double baro_active;						/**< double defining whether the reflex is
													1.0, active
													0.0, inactive */

	double baro_A;							/**< double for the baroreflex afferent */

	double baro_B;							/**< double for the baroreflex balance */

	double baro_P_set;						/**< double for the current pressure setpoint */

	double baro_S;							/**< double for the S variable */

	double baro_k_drive;					/**< double for the baro_k_drive parameter */
	
	double baro_k_recov;					/**< double for the baro_k_recovery parameter */

	int baro_P_compartment;					/**< int for the compartment number used to
													calculate the B_a signal, subtract 1
													to get the compartment index */

	// Other functions
	void initialise_simulation(void);

	void implement_time_step(double time_step_s);

	void calculate_B_a(void);

	void calculate_B_b(double time_step_s);
};
