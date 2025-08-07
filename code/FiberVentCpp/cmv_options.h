#pragma once

/**
/* @file		cmv_options.h
/* @brief		Header file for a cmv_options object
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

// Forward declarations
class FiberSim_options;
class MyoSim_options;


using namespace std;

class cmv_options
{
public:
	/**
	 * Constructor
	 */
	cmv_options(string set_options_file_string);

	/**
	* Destructor
	*/
	~cmv_options(void);

	// Variables
	string options_file_string;				/**< string for the protocol file */

	double bin_min;							/**< double with minimum cross-bridge
													bin position */

	double bin_max;							/**< double with maximum cross-bridge
													bin position */

	double bin_width;						/**< double with width of cross-bridge
													bin */

	double max_rate;						/**< double with maximum rate for a
													cross-bridge rate in s^-1 */

	string rates_dump_relative_to;			/**< string defining path type
													for rates_dump file */

	string rates_dump_file_string;			/**< string with rates file */

	string cb_dump_relative_to;				/**< string defining path type
													for cb_dump file */

	string cb_dump_file_string;				/**< string with cb dump file */

	double force_control_max_delta_hs_length;
											/**< double with maximum length change for the
													force root finder */

	string hv_thick_wall_approximation;		/**< string defining thick wall approximation
													If True, use thick-wall approximation
													otherwise, use thin-wall */

	double beat_length_s;					/**< double defining the length in s
													of a cmv_results object
													for a beat */

	int beat_length_points;					/**< int defining number of time-points
													in the cmv_results object for
													a beat */

	double summary_time_step_s;				/**< double defining the time-step
													in seconds between points in
													the summary output */

	int summary_points;						/**< int defining the number of
													time-points in the
													cmv_results object for the
													summary output */

	double write_every_s;					/**< double defining how often to write the
													output file - useful to evaluate
													long simulations */

	FiberSim_options* p_FiberSim_options;	/**< pointer to a FiberSim options object */

	MyoSim_options* p_MyoSim_options;		/**< pointer to a MyoSim options object */

	/**
	/* Function initialises protocol object from file
	*/
	void initialise_options_from_JSON_file(string JSON_options_file_string);

};
