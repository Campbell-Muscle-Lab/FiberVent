#pragma once

/**
/* @file		MyoSim_series_component.h
/* @brief		Header file for a MyoSim_series_component object
/* @author		Ken Campbell
*/

#include <iostream>

// Forward declaration
class MyoSim_muscle;

class MyoSim_series_component;
{
public:

	// Variables
	
	MyoSim_muscle* p_parent_myo_muscle;						/**< Pointer to the parent MyoSim muscle */

	double sc_extension;									/**< double holding sc_length in nm */

	double sc_k_stiff;										/**< double holding sc_stiffness in N M^-1 */

	double sc_force;										/**< double holding the sc_force in N m^-2 */

	// Functions

	/**
	* Constructor
	*/
	MyoSim_series_component(MyoSim_muscle* set_p_parent_myo_muscle);

	/**
	* Destructur
	*/
	~MyoSim_series_component(void);
};