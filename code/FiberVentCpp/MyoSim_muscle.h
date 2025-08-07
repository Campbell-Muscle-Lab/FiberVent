#pragma once

/**
/* @file		MyoSim_muscle.h
/* @brief		Header file for a MyoSim_muscle object
/* @author		Ken Campbell
*/

#include <iostream>

// Forward declaration

class cmv_system;
class cmv_model;
class cmv_results;
class cmv_options;

class hemi_vent;
class muscle;

class MyoSim_half_sarcomere;
class MyoSim_series_component;

class MyoSim_muscle
{
public:
	// Variables
	muscle* p_parent_muscle;							/**< Pointer to a parent muscle */

	MyoSim_half_sarcomere* p_MyoSim_hs;					/**< Pointer to a MyoSim half-sarcomere object */

	MyoSim_series_component* p_MyoSim_sc;				/**< Pointer to a MyoSim series component */

	// Functions

	/**
	* Constructor
	*/
	MyoSim_muscle(muscle* set_p_parent_muscle);

	/**
	* Destructur
	*/
	~MyoSim_muscle(void);

	// Others
	void implement_time_step(double time_step_s);
};