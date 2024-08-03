#pragma once

/**
/* @file		MyoSim_options.h
/* @brief		Header file for a MyoSim_options object
/* @author		Ken Campbell
*/

#include <iostream>

// Forward declarations
class cmv_options;

class MyoSim_options
{
public:
	/**
	* Constructor
	*/
	MyoSim_options(cmv_options* set_p_cmv_options);

	/**
	* Destructor
	*/
	~MyoSim_options(void);

	// Variables
	cmv_options* p_cmv_options;
};
