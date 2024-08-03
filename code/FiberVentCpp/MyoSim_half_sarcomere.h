#pragma once

/**
/* @file		MyoSim_half_sarcomere.h
/* @breif		Header file for a MyoSim half-sarcomere object
/* @author		Ken Campbell
*/

class half_sarcomere;

class MyoSim_half_sarcomere
{
public:
	/**
	* Constructor
	*/
	MyoSim_half_sarcomere(half_sarcomere* set_p_hs);

	/**
	* Destructor
	*/
	~MyoSim_half_sarcomere(void);

	// Variables
	half_sarcomere* p_hs;				/**< pointer to a half_sarcomere */
};