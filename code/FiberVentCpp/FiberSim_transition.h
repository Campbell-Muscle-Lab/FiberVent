#pragma once

/**
* @file		transition.h
* @brief	Header file for the transition class
* @author	Ken Campbell
*/

//#include "FiberSim_half_sarcomere.h"
//#include "FiberSim_m_state.h"

#include "rapidjson/document.h"
#include "JSON_functions.h"

#include "gsl_vector.h"

// Forward declaration
class FiberSim_half_sarcomere;
class FiberSim_m_state;

class FiberSim_transition
{
public:

	// Variables

	FiberSim_m_state* p_parent_m_state;		/**< pointer to parent m_state */

	int new_state;					/**< integer defining the new state */

	char rate_type[_MAX_PATH];		/**< char array defining the transition type */

	char transition_type;			/**< char defining 'a' attachment, 'd' detachment, 'n' neutral */

	int uses_ATP;					/**< integer, defining how many molecules of ATP the
											transition uses */

	gsl_vector* rate_parameters;	/**< gsl_vector holding parameter variables */

	// Functions

	// Constructor

	/**
	* Normal constructor called with an entry into a JSON document
	*/
	FiberSim_transition(const rapidjson::Value& m_ks, FiberSim_m_state* set_p_m_state);

	/**
	* Constructor that sets default values
	* Called when the transition is empty
	*/
	FiberSim_transition();

	/**
	* Destructor
	*/
	~FiberSim_transition(void);

	// Functions

	/**
	* double calculate_rate(double x)
	* @param x double defining the cb_x position - the bs_x position
	* @return the rate in units of s^-1
	*/
	double calculate_rate(double x, double x_ext, double node_force,
							int mybpc_state, int mybpc_iso,
							short int active_neigh,
							FiberSim_half_sarcomere* p_fs_hs);
};
