/**
/* @file		activation.cpp
/* @brief		Source file for an activation object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "activation.h"

using namespace::std;

// Constructor
activation::activation(string set_activation_type, double set_t_start_s, double set_t_stop_s)
{
	// Initialise

	// Code
	activation_type = set_activation_type;
	t_start_s = set_t_start_s;
	t_stop_s = set_t_stop_s;
}

// Destructor
activation::~activation(void)
{
	// Destructor
	// 
	// Tidy up
}

// Other functions
int activation::return_status(string test_type, double t_test_s)
{
	// Returns 1 if activation type matches test and is active, 0 otherwise

	if (test_type != activation_type)
		return 0;

	if ((t_test_s >= t_start_s) && (t_test_s <= t_stop_s))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
