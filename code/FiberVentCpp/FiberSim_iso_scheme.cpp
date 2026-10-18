/**
* @file		FiberSim_iso_scheme.cpp
* @brief	Source file for the iso_scheme class
* @author	Ken Campbell
*/

#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "FiberSim_model.h"
#include "FiberSim_options.h"
#include "FiberSim_iso_scheme.h"
#include "FiberSim_iso_type.h"
#include "FiberSim_iso_transition.h"

#include "gsl_vector.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"

using namespace std::filesystem;

// Constructor
FiberSim_iso_scheme::FiberSim_iso_scheme(const rapidjson::Value& iso_sch,
	FiberSim_model* set_p_fs_model, FiberSim_options* set_p_fs_options)
{
	// Initialise

	// Set the pointer to the model
	p_fs_model = set_p_fs_model;

	// Set the pointer to the options
	p_fs_options = set_p_fs_options;

	// Check the isotypes
	JSON_functions::check_JSON_member_array(iso_sch, "type");
	const rapidjson::Value& types = iso_sch["type"];

	// Deduce the max number of transitions - we need this to make the types
	max_no_of_transitions = 0;
	for (rapidjson::SizeType i = 0; i < types.Size(); i++)
	{
		const rapidjson::Value& trans = types[i]["transition"];
		max_no_of_transitions = GSL_MAX(max_no_of_transitions, (int)trans.Size());
	}

	// Allocate the transition_probs_vector
	transition_probs = gsl_vector_alloc(max_no_of_transitions);

	// Now cycle through making the isotypes
	for (rapidjson::SizeType i = 0; i < types.Size(); i++)
	{
		p_iso_types[i] = new FiberSim_iso_type(types[i], this);
	}

	no_of_isotypes = types.Size();
}

// Destructor
FiberSim_iso_scheme::~FiberSim_iso_scheme(void)
{
	//! Destructor

	// Tidy up

	// Delete the iso_types
	for (int type_counter = 0; type_counter < no_of_isotypes; type_counter++)
	{
		delete p_iso_types[type_counter];
	}

	// Delete the transition probs
	gsl_vector_free(transition_probs);
}

// Functions
