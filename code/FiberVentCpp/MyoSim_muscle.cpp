/**
/* @file		MyoSim_muscle.cpp
/* @brief		Source file for a MyoSim_muscle object
/* @author		Ken Campbell
*/

#include <cstdio>

#include "muscle.h"

#include "MyoSim_muscle.h"

// Constructor
MyoSim_muscle::MyoSim_muscle(muscle* set_p_parent_muscle)
{
	//! Constructor
	
	// Set pointer to parent
	p_parent_muscle = set_p_parent_muscle;
}

// Destructor
MyoSim_muscle::~MyoSim_muscle(void)
{
	//! Destructor
}

// Other functions
void MyoSim_muscle::implement_time_step(double time_step_s)
{
	//! 
}