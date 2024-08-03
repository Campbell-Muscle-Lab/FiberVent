#pragma once

/**
* @file		kinetic_scheme.h
* @brief	header file for the kinetic_scheme class
* @author	Ken Campbell
*/

#include <iostream>

#include "rapidjson/document.h"
#include "JSON_functions.h"

#include "gsl_vector.h"

#include "global_definitions.h"

// Forward declarations
class cmv_model;
class cmv_options;
class m_state;
class myofilaments;

using namespace std;

class kinetic_scheme
{
public:

	// Variables

	cmv_model* p_cmv_model;					/**< pointer to parent model */

	cmv_options* p_cmv_options;				/**< pointer to a cmv_options object */

	myofilaments* p_parent_myofilaments;	/**< pointer to parent myofilaments object */

	int no_of_states;						/**< int defining the number of different states */

	int no_of_detached_states;				/**< int defining the number of detached states */

	int no_of_attached_states;				/**< int defining the number of attached states */

	int max_no_of_transitions;				/**< int defining the maximum number of transitions from a state */

	int first_DRX_state;					/**< int defining first DRX state */

	int last_ATP_requiring_transition;		/**< int with the index of the last transition
													that needs ATP */

	m_state* p_m_states[MAX_NO_OF_KINETIC_STATES];
											/**< pointer to an array of m_state objects */

	// Functions

	/**
	* Constructor
	* takes a cmv_model and parses it to give the kinetic scheme
	*/
	kinetic_scheme(const rapidjson::Value& m_ks, cmv_model* set_p_cmv_model);

	/**
	* Destructor
	*/
	~kinetic_scheme(void);

	/**
	* void set_transition_types(void)
	* loops through transitions from each state and sets the type, 'a' for attach,
	* 'd' for detach, and 'n' for neither
	* needs to be run after all the states are known so can't be included in the
	* transition constructor as these are built in sequence with each state
	* @return void
	*/
	void set_transition_types(void);

	/**
	* updates 
	*/
	void initialise_simulation(myofilaments* set_p_parent_myofilaments);

	/**
	* void write_kinetic_scheme_to_file(char output_file_string)
	* writes kinetic_scheme to specified file in JSON format
	* @return void
	*/
	void write_kinetic_scheme_to_file(char output_file_string[]);

	/**
	* void write_rate_functions_to_file()
	* as a JSON file structure
	* @return void
	*/
	void write_rate_functions_to_file();

};