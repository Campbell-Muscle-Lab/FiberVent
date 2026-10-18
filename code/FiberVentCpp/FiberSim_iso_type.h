#pragma once

/**
* @file		FiberSim_iso_type.h
* @brief	Header file for the iso_type class
* @author	Ken Campbell
*/

#include "gsl_vector.h"
#include "global_definitions.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"


class FiberSim_iso_scheme;
class FiberSim_iso_transition;

class FiberSim_iso_type
{
public:

	// Variables

	FiberSim_iso_scheme* p_parent_iso_scheme;
	/**< pointer to the parent iso scheme */

	int iso_number;					/**< integer defining the iso number */

	int no_of_transitions;			/**< integer defining the number of transitions
										 from the iso_typee */

	FiberSim_iso_transition* p_iso_transitions[MAX_NO_OF_TRANSITIONS];

	/**
	* Constructor
	*/
	FiberSim_iso_type(const rapidjson::Value& iso_ty, FiberSim_iso_scheme* set_p_parent_scheme);

	/**
	* Destuctor
	*/
	~FiberSim_iso_type(void);

};

