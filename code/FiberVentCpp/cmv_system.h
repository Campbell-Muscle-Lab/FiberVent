#pragma once

/**
/* @file		cmv_system.h
/* @brief		Header file for a cmv_system object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include <string>

// Forward declarations
class cmv_model;
class cmv_options;
class cmv_protocol;
class cmv_results;
class circulation;
class hemi_vent;

using namespace std;

class cmv_system
{
public:
	/**
	 * Constructor
	 */
	cmv_system(string JSON_model_file_string, string options_file_string, int system_id);
		/**

	* Destructor
	*/
	~cmv_system(void);

	// Variables
	cmv_model* p_cmv_model;					/**< Pointer to cmv_model object */

	cmv_options* p_cmv_options;				/**< Poniter to cmv_options object */

	cmv_protocol* p_cmv_protocol;			/**< Pointer to a cmv_protocol object */

	cmv_results* p_cmv_results_summary;	/**< Pointer to cmv_results holding
													down-sampled data for the
													simulation */

	cmv_results* p_cmv_results_beat;		/**< Pointer to cmv_results holding
													data for a beat */

	circulation* p_circulation;				/**< Pointer to a circulation */

	double time_step_s;						/**< double holding the simuation
													time_step in s */

	int sim_t_index;						/**< integer holding index in the simulation */

	int beat_t_index;						/**< integer holding index in the
													beat_results object */

	int summary_t_index;					/**< integer holding index in the
													summary results object */

	double cum_time_s;						/**< double, with system time in s */

	double last_write_s;					/**< system time in s at which results file
													was last updated. Used for progressive
													file writes */

	int system_id;

	// Functions

	/**
	/* function ensures p_clone has same fields as p_source where
	* p_clone and p_source are both cmv_results objects
	*/
	void clone_results_fields(cmv_results* p_source, cmv_results* p_clone);

	/**
	/* function runs a simulation
	*/
	void run_simulation(string protocol_file_string,
		string results_file_string);

	void add_fields_to_cmv_results_beat();

	bool implement_time_step(double time_step_s);

	void update_beat_metrics();

	void update_cmv_results_summary();

	bool sim_time_dumps_to_summary(double sim_time);
};