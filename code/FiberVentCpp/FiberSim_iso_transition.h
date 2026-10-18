#pragma once

/**
* @file		FiberSim_transition.h
* @brief	Header file for the transition class
* @author	Ken Campbell
*/

#include "FiberSim_half_sarcomere.h"
#include "rapidjson/document.h"
#include "JSON_functions.h"
#include "FiberSim_options.h"

#include "gsl_vector.h"

// Forward declaration
class FiberSim_half_sarcomere;
class FiberSim_iso_type;

class FiberSim_iso_transition
{
public:

	// Variables

	FiberSim_iso_type* p_parent_iso_type;		/**< pointer to the parent isotype */

	FiberSim_options* p_fs_options;	/**< pointer to the FiberSim options */

	int new_iso_type;				/**< integer defining the new isotype */

	char rate_type[_MAX_PATH];		/**< char array defining the transition type */

	gsl_vector* rate_parameters;	/**< gsl_vector holding parameter variables */

	// Functions

	// Constructor

	/**
	* Normal constructor called with an entry into a JSON document
	*/
	FiberSim_iso_transition(const rapidjson::Value& iso_tr, FiberSim_iso_type* p_iso_type);

	/**
	* Constructor that sets default values
	* Called when the transition is empty
	*/
	FiberSim_iso_transition();

	/**
	* Destructor
	*/
	~FiberSim_iso_transition(void);

	// Functions

	/**
	* double calculate_rate(double x)
	* @param x double defining the cb_x position
	* @return the rate in units of s^-1
	*/
	double calculate_rate(double x, FiberSim_half_sarcomere* p_hs);
};
