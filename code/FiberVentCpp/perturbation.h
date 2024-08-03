#pragma once

/**
/* @file		perturbation.h
/* @brief		Header file for the perturbation class
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>
#include <regex>

#include "global_definitions.h"

// Forward declararations
class cmv_protocol;

struct perturbation_struct;

using namespace::std;

class perturbation
{
public:
	/**
	 * Constructor
	 */
	perturbation(cmv_protocol* set_p_cmv_protocol, perturbation_struct* set_perturbuation_data);

	/**
	 * Destructor
	 */
	 ~perturbation(void);

	 // Variables
	 cmv_protocol* p_cmv_protocol;			/**< pointer to parent cmv_protocol */

	 string class_name;						/**< string holding the perturbation class
													circulation, myofilaments, baroreflex etc. */

	 string variable;						/**< the variable that is perturbed */

	 double t_start_s;						/**< double defining when the perturbation starts */

	 double t_stop_s;						/**< double definining when the perturbation stops */

	 double total_change;					/**< double defining the total change */

	 double prot_time_step_s;				/**< double with simulation time_step in s */

	 double increment;						/**< change per time-step */

	 // Functions

	 int return_status(string test_type, double t_test_s);

	 void impose(double time_step_s);

	 void extract_digits(string test_string, int digits[], int no_of_digits);

	 void make_adjustment(double* p_double, double increment);


};
