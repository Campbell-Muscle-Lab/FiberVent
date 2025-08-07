#pragma once

/**
/* @file		membranes.h
/* @brief		Header file for a membranes object
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

class membranes
{
public:
	/**
	 * Constructor
	 */
	membranes(muscle* set_p_parent_muscle);

	/**
	* Destructor
	*/
	~membranes(void);

	// Variables

	muscle* p_parent_muscle;			/**< pointer to parent half-sarcomere */

	cmv_model* p_cmv_model;				/**< pointer to cmv_model object */

	cmv_options* p_cmv_options;			/**< pointer to cmv_options object */

	cmv_results* p_cmv_results_beat;	/**< pointer to cmv_results object
												holding data at full time
												resolution for one beat */

	double memb_Ca_cytosol;				/**< double with cytosolic Ca_conc in M */

	double memb_Ca_sr;					/**< double with sarcoplasmc reticular
												Ca concentration in M */

	double memb_activation;				/**< double with membrane activation
												1.0 when active
												0.0 when inactive */

	double memb_t_open_s;				/**< double with active duration in s */

	double memb_t_open_left_s;			/**< double with active duration left
												in s */

	double memb_k_serca;				/**< double describing rate at which
												SERCA pumps Ca in to
												sarcoplasmic reticulum
												in M^-1 s^-1 */

	double memb_k_leak;					/**< double desribing rate at which
												Ca leaks from sarcoplasmic reticulum
												in M s^-1 */

	double memb_k_active;				/**< double describing rate at which
												Ca is released from sarcoplasmic reticulum
												when the membrane is active
												in M s^-1 */

	double memb_J_release;				/**< double describing Ca release flux
												 in M s^-1 */

	double memb_J_uptake;				/**< double describing Ca uptake flux
												 in M s^-1 */


	/**
	/* function adds data fields and vectors to the results objet
	*/
	void initialise_simulation(void);

	/**
	/* function updates state variables
	*/
	void implement_time_step(double time_step_s, bool new_beat);

	/**
	/*
	* function calculates derivs
	*/
	//int calculate_derivs(double t, const double y[], double f[], void* params);

	/**
	/*
	* function calculates fluxes based on input concentrations
	*/
	void calculate_fluxes(const double y[]);

};
