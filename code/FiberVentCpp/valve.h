#pragma once

/**
/* @file		half_sarcomere.h
/* @brief		Header file for a half_sarcomere object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_system;
class cmv_model;
class cmv_results;
class cmv_options;

class hemi_vent;

struct cmv_model_valve_structure;

class valve
{
public:
	/**
	* Constructor
	*/
	valve(hemi_vent* set_p_parent_hemi_vent, cmv_model_valve_structure* set_p_structure);

	/**
	* Destructor
	*/
	~valve(void);

	// Variables
	hemi_vent* p_parent_hemi_vent;
	
	cmv_model* p_cmv_model;							/**< Pointer to the cmv_model object */

	cmv_results* p_cmv_results_beat;				/**< Pointer to cmv_results object
															holding data at full time
															resolution for one beat */

	cmv_options* p_cmv_options;						/**< Pointer to cmv_options */

	cmv_model_valve_structure* p_cmv_model_valve;	/**< Pointer to the structure in the
															model for the valve */

	double valve_pos;								/**< Double holding valve status
															1.0 = open
															0.0 = closed */

	double valve_last_pos;							/**< Double holding valve status
															on last iteration */
	
	double valve_vel;								/**< Double holding valve velocity */

	string valve_name;								/**< string holding valve name */

	double valve_mass;								/**< Double holding the valve mass */

	double valve_eta;								/**< Double holding the valve resistance */

	double valve_k;									/**< Double holding the valve stiffness */

	double valve_leak;								/**< double holding the valve leak
															0 if doesn't leak
															<0 if it does */

	void initialise_simulation(void);
	
	void implement_time_step(double time_step_s);
};