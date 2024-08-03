#pragma once

/**
/* @file		growth_activation.h
/* @brief		Header file for a baroreflex object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_protocol;

class growth_activation
{
public:
	/**
	 * Constructor
	 */
	growth_activation(double set_t_start_s, double t_stop_s);

	/**
	 * Destructor
	 */
	 ~growth_activation(void);

	 // Variables
	 double t_start_s;						/**< double defining when the activation starts */

	 double t_stop_s;						/**< double definining when the activation stops */

	 // Functions

	 int return_status(double t_test_s);
};
