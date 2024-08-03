#pragma once

/**
/* @file		circulation.h
/* @brief		Header file for a circulation object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_system;
class cmv_model;
class cmv_protocol;
class cmv_options;
class cmv_results;

class hemi_vent;
class valve;
class baroreflex;
class growth;

class circulation
{
public:
	/**
	 * Constructor
	 */
	circulation(cmv_system* set_p_parent_cmv_system);

	/**
	* Destructor
	*/
	~circulation(void);

	// Variables

	cmv_system* p_parent_cmv_system;					/**< Pointer to the parent_cmv_system */

	cmv_model* p_cmv_model;								/**< Pointer to the cmv_model */

	cmv_protocol* p_cmv_protocol;						/**< Pointer to the cmv protocol object */

	cmv_options* p_cmv_options;							/**< Pointer to cmv_options */

	cmv_results*  p_cmv_results_beat;					/**< Pointer to cmv_results object
																holding data at full
																time resolution for the last
																beat */

	hemi_vent* p_hemi_vent;								/**< Pointer to a hemi_vent object */

	valve* p_av;										/**< Pointer to the aortic valve */

	valve* p_mv;										/**< Pointer to the mitral valve */

	baroreflex* p_baroreflex;							/**< Pointer to the baroreflex object */

	growth* p_growth;									/**< Pointer to the growth object */

	bool filling_control;								/**< bool defining whether filling control
																has been defined */

	double filling_control_venous_pressure;				/**< double with venous pressure when filling
																is being controlled */

	double circ_blood_volume;							/**< double holding total blood volume
																in liters */

	int circ_no_of_compartments;						/**< integer holding number of
																compartments */

	double* circ_resistance;							/**< Pointer to array of doubles
																holding resistances for
																each compartment */

	double* circ_compliance;							/**< Pointer to array of doubles
																holding compliances for
																each compartment */

	double* circ_slack_volume;							/**< Pointer to array of doubles
																holding slack_volumes for
																each compartment */

	double* circ_inertance;								/**< Pointer to array of doubles
																holding inertances for
																each compartment */

	double* circ_pressure;								/**< Pointer to array of doubles
																holding pressure in each
																compartment */

	double* circ_static;								/**< Pointer to array of doubles
																holding static pressure in
																each compartment */

	double* circ_volume;								/**< Pointer to array of doubles
																holding pressure in each
																compartment */

	double* circ_flow;									/**< Pointer to array of doubles
																holding flows between
																compartments */

	double* circ_last_flow;								/**< Pointer to array of doubles
																holding flows between
																compartments on last
																time step */


	double circ_total_slack_volume;						/**< double holding total slack
																volume in liters */

	// Functions

	void initialise_simulation(void);

	bool implement_time_step(double time_step_s);

	void calculate_pressures(const double v[], double p[], double time_step_s);

	void calculate_flows(const double v[], double flow[]);

	void update_beat_metrics(void);
};