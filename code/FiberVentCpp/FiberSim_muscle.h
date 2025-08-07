#pragma once

/**
/* @file		FiberSim_muscle.h
/* @brief		Header file for a FiberSim_muscle object
/* @author		Ken Campbell
*/

#include <iostream>

#include "gsl_vector.h"

// Forward declaration

class cmv_system;
class cmv_model;
class cmv_results;
class cmv_options;

class muscle;

class FiberSim_half_sarcomere;
class FiberSim_series_component;
class FiberSim_options;

class FiberSim_muscle
{
public:
	// Variables
	muscle* p_parent_muscle;							/**< Pointer to the parent muscle */

	cmv_results* p_cmv_results_beat;					/**< Pointer to a cmv_results object */

	FiberSim_half_sarcomere* p_FiberSim_hs;				/**< Pointer to a FiberSim half-sarcomere object */

	FiberSim_series_component* p_FiberSim_sc;			/**< Pointer to a FiberSim series component */

	FiberSim_options* p_FiberSim_options;				/**< Pointer to a FiberSim options object */

	double fs_m_length;									/**< double holding the length of the muscle */

	double fs_m_stress;									/**< double holding the stress in the muscle */

	int dump_status_counter;							/**< integer used to handle when to dump status files */

	// Functions

	/**
	* Constructor
	*/
	FiberSim_muscle(muscle* set_p_parent_muscle);

	/**
	* Destructur
	*/
	~FiberSim_muscle(void);

	// Other functions

	void initialise_for_simulation(void);

	void implement_time_step(double time_step_s);		/**< Runs sarcomere kinetics */

	int change_muscle_length(double delta_ml, double time_step_s);
														/** Changes muscle length */

	size_t worker_length_control_myofibril_with_series_compliance(
		const gsl_vector* x, void* p, gsl_vector* f);

	double return_muscle_length_for_force(double target_force, double time_step_s);
														/** Returns the muscle length for
															a given force */

	double calculate_delta_ml_for_force(double target_force, double time_step_s);

	double return_wall_stress_after_test_delta_ml(double delta_ml, double time_step_s);

	void dump_hs_status(int time_step);

	void write_rates_file(void);
};