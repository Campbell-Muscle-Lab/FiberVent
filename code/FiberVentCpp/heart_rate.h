#pragma once

/**
/* @file		hemi_vent.h
/* @brief		Header file for a hemi_vent object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_system;
class cmv_model;
class cmv_results;
class muscle;


class heart_rate
{
public:
	/**
	 * Constructor
	 */
	heart_rate(muscle* set_p_parent_muscle);

	/**
	 * Destructor
	 */
	 ~heart_rate(void);

	 // Variables

	cmv_system* p_parent_cmv_system;		/**< pointer to parent cmv_system */

	cmv_model* p_cmv_model;					/**< pointer to cmv_model object */

	cmv_results* p_cmv_results_beat;		/**< pointer to cmv_results object
													holding data at full time
													resolution for one beat */

	muscle* p_parent_muscle;				/**< pointer to parent muscle */

	double hr_new_beat;						/**< double, set to 0.0 normally but to
													1.0 on a new beat */

	double hr_t_RR_interval_s;				/**< double, with RR interval in seconds */

	double hr_t_countdown_s;				/**< double, with t to next heart-beat
													in s */

	double hr_heart_rate_bpm;				/**< double with heart rate in bpm */

	// Other functions
	
	/**
	/* function prepares for simulation
	*/
	void initialise_simulation();

	/**
	* function manages the hr_activation variable which determines the
	* initiation of a heart-beat
	*/
	bool implement_time_step(double time_step);
};
