/**
/* @file		MyoSim_options.cpp
/* @brief		Source file for the MyoSim_options class
/* @author		Ken Campbell
*/

#include "stdio.h"

#include "MyoSim_options.h"

#include "cmv_options.h"

// Constructor
MyoSim_options::MyoSim_options(cmv_options* set_p_cmv_options)
{
	//! Constructor
	
	p_cmv_options = set_p_cmv_options;
}

// Destructor
MyoSim_options::~MyoSim_options()
{

}