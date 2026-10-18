#pragma once

/**
* @file		FiberSim_iso_switch_scheme.h
* @brief	header file for the iso_switch_scheme class
* @author	Ken Campbell
*/

#include "rapidjson/document.h"
#include "JSON_functions.h"

#include "gsl_vector.h"

#include "FiberSim_model.h"
#include "FiberSim_options.h"

#include "global_definitions.h"

// Forward declarations
class FiberSim_half_sarcomere;
class FiberSim_iso_type;

class FiberSim_iso_scheme
{
public:

	// Variables

	FiberSim_model* p_fs_model;			/**< pointer to a FiberSim_model object */

	FiberSim_options* p_fs_options;		/**< pointer to a FiberSim_options objects */

	int no_of_isotypes;					/**< int defining the number of isotypes */

	FiberSim_iso_type* p_iso_types[MAX_NO_OF_ISOTYPES];
	/**< pointer to an array of isotype objects */

	int max_no_of_transitions;			/**< int defining the maximum number of transitions
											 from a state */

	gsl_vector* transition_probs;		/**< gsl_vector holding transition probabilities */

	// Functions

	/**
	* Constructor
	* takes a FiberSim_model and parses it to give the iso_switch scheme
	*/
	FiberSim_iso_scheme(const rapidjson::Value& iso_sch, FiberSim_model* set_p_fs_model,
		FiberSim_options* set_p_fs_options);

	/**
	* Destructor
	*/
	~FiberSim_iso_scheme(void);
};
