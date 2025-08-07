#pragma once

/**
/* @file		cmv_protocol.h
/* @brief		Header file for a cmv_protocol object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include <iostream>
#include <string>

// Definitions for JSON parsing
#ifndef _RAPIDJSON_DOCUMENT
#define _RAPIDJSON_DOCUMENT
#include "rapidjson/document.h"
#endif

#include "cmv_system.h"

#include "global_definitions.h"

using namespace std;

class activation;
class perturbation;

class cmv_protocol
{
public:
	/**
	 * Constructor
	 */
	cmv_protocol(cmv_system* set_p_cmv_system, string protocol_file_string);

	/**
	* Destructor
	*/
	~cmv_protocol(void);

	// Variables
	cmv_system* p_cmv_system;				/**< pointer to the parent system */

	string protocol_file_string;			/**< string for the protocol file */

	double time_step_s;						/**< double holding time_step in s */
	
	int no_of_time_steps;					/**< int holding number of time-steps */

	int no_of_activations;					/**< int holding the number of activations */

	int no_of_perturbations;				/**< int holding the number of perturbations */

	activation* p_activation[MAX_NO_OF_ACTIVATIONS];
											/**< an array of pointers to activation
													objects */

	perturbation* p_perturbation[MAX_NO_OF_PERTURBATIONS];
											/**< an array of pointers to perturbation
													objects */

	// Functions

	/**
	/* Function initialises protocol object from file
	*/
	void initialise_protocol_from_JSON_file(string JSON_protocol_file_string);

	double return_activation(string activation_type, double time_s);

	void impose_perturbations(double sim_time_s);
};