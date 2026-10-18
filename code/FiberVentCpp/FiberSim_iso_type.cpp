/**
* @file		FiberSim_iso_type.cpp
* @brief	Source file for the iso_type class
* @author	Ken Campbell
*/

#include <cstdio>

#include "FiberSim_iso_type.h"
#include "FiberSim_iso_scheme.h"
#include "FiberSim_iso_transition.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"

// Constructor
FiberSim_iso_type::FiberSim_iso_type(const rapidjson::Value& iso_ty, FiberSim_iso_scheme* set_p_parent_scheme)
{
	p_parent_iso_scheme = set_p_parent_scheme;

	JSON_functions::check_JSON_member_int(iso_ty, "number");
	iso_number = iso_ty["number"].GetInt();

	// Pull array of transitions
	JSON_functions::check_JSON_member_array(iso_ty, "transition");
	const rapidjson::Value& trans = iso_ty["transition"];

	for (int i = 0; i < (int)trans.Size(); i++)
	{
		p_iso_transitions[i] = new FiberSim_iso_transition(trans[i], this);
	}

	// Fill in gaps if not all transitions are set
	if (trans.Size() < (size_t)(p_parent_iso_scheme->max_no_of_transitions))
	{
		for (int i = trans.Size(); i < p_parent_iso_scheme->max_no_of_transitions; i++)
		{
			p_iso_transitions[i] = new FiberSim_iso_transition();
		}
	}

	no_of_transitions = trans.Size();
}

// Destructor
FiberSim_iso_type::~FiberSim_iso_type(void)
{
	// Tidy up
	for (int i = 0; i < p_parent_iso_scheme->max_no_of_transitions; i++)
		delete p_iso_transitions[i];

}
