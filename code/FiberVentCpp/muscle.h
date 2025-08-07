#pragma once

/**
/* @file		muscle.h
/* @brief		Header file for a muscle object
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

class membranes;
class heart_rate;
class mitochondria;

class MyoSim_muscle;
class FiberSim_muscle;


class muscle
{
public:
	/**
	* Constructor
	*/
	muscle(hemi_vent* set_p_parent_hemi_vent);

	/**
	* Destructor
	*/
	~muscle(void);

	// Variables
	hemi_vent* p_parent_hemi_vent;
	
	cmv_model* p_cmv_model;							/**< Pointer to the cmv_model object */

	cmv_results* p_cmv_results_beat;				/**< Pointer to cmv_results object
															holding data at full time
															resolution for one beat */

	cmv_options* p_cmv_options;						/**< Pointer to cmv_options */

	cmv_system* p_cmv_system;						/**< Pointer to the cmv system */

	heart_rate* p_heart_rate;						/**< Pointer to heart-rate object */

	membranes* p_membranes;							/**< Pointer to the daughter
															membranes object */

	mitochondria* p_mitochondria;					/**< Pointer to the daughter
															mitochondria object */

	double muscle_prop_fibrosis;					/**< Double with the proportion of the
															hs cross-sectional area
															that is fibrosis */

	double muscle_prop_myofilaments;				/**< Double with the proportion of the
															hs non-fibrotic area that is
															that is myofilaments, remainder
															is mitochondria */

	double muscle_length;
	double muscle_reference_length;
	double muscle_stress;

	double muscle_ATP_used_myof_per_liter_per_s;	/**< double with ATP used by myofilaments
															in moles per liter of muscle
															per s */

	double muscle_ATP_used_memb_per_liter_per_s;	/**< double with ATP used by membranes
															in moles per liter of muscle
															per s */

	double muscle_ATP_used_total_per_liter_per_s;	/**< double with total ATP used in moles
															per liter of muscle per s */

	double muscle_ATP_generated_per_liter_per_s;	/**< double with ATP generated in moles
															per liter of muscle per s */

	double muscle_ATP_concentration;				/**< double with ATP concentration */

	double muscle_delta_G_ATP;						/**< double with energy in Joules
															per mole of ATP */

	FiberSim_muscle* p_FiberSim_muscle;				/**< pointer to a FiberSim muscle
															object */

	MyoSim_muscle* p_MyoSim_muscle;					/**< pointer to a MyoSim muscle
															object */

	/**
	/* function adds data fields and vectors to the results object
	*/
	
	void initialise_simulation(void);
	
	bool implement_time_step(double time_step_s);

	void change_muscle_length(double delta_ml, double time_step_s);

	double return_wall_stress_after_test_delta_ml(double delta_ml, double time_step_s);

	double return_muscle_length_for_stress(double target_stress, double time_step_s);

	void calculate_muscle_ATP_concentration(double time_step);
	
	void update_beat_metrics(void);
};