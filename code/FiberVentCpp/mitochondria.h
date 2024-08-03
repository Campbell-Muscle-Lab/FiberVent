#pragma once

/**
/* @file		mitchondria.h
/* @brief		Header file for a mitochondria object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class muscle;
class cmv_model;
class cmv_options;
class cmv_results;

class mitochondria
{
public:
	/**
	 * Constructor
	 */
	mitochondria(muscle* set_parent_muscle);

	/**
	* Destructor
	*/
	~mitochondria(void);

	// Variables

	muscle* p_parent_muscle;			/**< pointer to parent half-sarcomere */

	cmv_model* p_cmv_model;				/**< pointer to cmv_model object */

	cmv_options* p_cmv_options;			/**< pointer to cmv_options object */

	cmv_results* p_cmv_results_beat;	/**< pointer to cmv_results object
												holding data at full time
												resolution for one beat */

	double mito_ATP_generation_rate;	/**< double with rate mitochondria
												generate ATP, in M s^-1 m^-3 */

	double mito_volume;					/**< double holding the mitochondrial
												volume in m^3 */

	double mito_ATP_generated_M_per_liter_per_s;
										/**< double holding the moles of ATP
												generated per liter of tissue
												per s */

	/**
	/* function adds data fields and vectors to the results objet
	*/
	void initialise_simulation(void);

	/**
	/* function updates state variables
	*/
	void implement_time_step(double time_step_s);
};
