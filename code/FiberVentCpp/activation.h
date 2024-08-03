#pragma once

/**
/* @file		activation.h
/* @brief		Header file for an activation object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_protocol;

using namespace::std;

class activation
{
public:
	/**
	 * Constructor
	 */
	activation(string set_activation_type, double set_t_start_s, double t_stop_s);

	/**
	 * Destructor
	 */
	 ~activation(void);

	 // Variables
	 string activation_type;				/**< string holding the activation type, typically
													baroreflex
													growth
													filling_pressure
													*/

	 double t_start_s;						/**< double defining when the activation starts */

	 double t_stop_s;						/**< double definining when the activation stops */

	 // Functions

	 int return_status(string test_type, double t_test_s);
};
