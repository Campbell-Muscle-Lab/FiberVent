/**
* @file		kinetic_scheme.cpp
* @brief	Source file for the kinetic_scheme class
* @author	Ken Campbell
*/

#include <cstdio>
#include <stdio.h>
#include <iostream>
#include <filesystem>

#include "kinetic_scheme.h"

#include "cmv_model.h"
#include "cmv_options.h"
#include "m_state.h"
#include "transition.h"

#include "JSON_functions.h"
#include "rapidjson\document.h"

using namespace std::filesystem;
using namespace std;

// Constructor
kinetic_scheme::kinetic_scheme(const rapidjson::Value& m_ks, cmv_model* set_p_cmv_model)
{
	//! Constructor for kinetic scheme

	// Code

	// Initialise
	p_cmv_model = set_p_cmv_model;

	// Set other options safely
	p_cmv_options = NULL;
	p_parent_myofilaments = NULL;

	// Pull no_of_states
	JSON_functions::check_JSON_member_int(m_ks, "no_of_states");
	no_of_states = m_ks["no_of_states"].GetInt();

	JSON_functions::check_JSON_member_int(m_ks, "max_no_of_transitions");
	max_no_of_transitions = m_ks["max_no_of_transitions"].GetInt();

	// Pull array
	JSON_functions::check_JSON_member_array(m_ks, "scheme");
	const rapidjson::Value& scheme = m_ks["scheme"];

	// Set the numbers of detached and attached states as we
	// are initialising

	no_of_detached_states = 0;
	no_of_attached_states = 0;

	// Note the first DRX state as we go
	first_DRX_state = -1;

	for (rapidjson::SizeType i = 0; i < scheme.Size(); i++)
	{
		p_m_states[i] = new m_state(scheme[i], this);

		if (first_DRX_state < 0)
		{
			if (p_m_states[i]->state_type == 'D')
				first_DRX_state = (int)(i + 1);
		}

		if (p_m_states[i]->state_type == 'A')
		{
			no_of_attached_states = no_of_attached_states + 1;
		}
		else
		{
			no_of_detached_states = no_of_detached_states + 1;
		}
	}

	// Now that we know the state properties, set the transition types
	set_transition_types();
}

// Destructor
kinetic_scheme::~kinetic_scheme(void)
{
	//! Destructor

	// Tidy up

	// Delete the kinetic states
	for (int state_counter = 0; state_counter < no_of_states;
		state_counter++)
	{
		delete p_m_states[state_counter];
	}
}

// Functions
void kinetic_scheme::set_transition_types(void)
{
	//! Cycle through the states, identifying the transition type for each one
	
	// Variables
	int t_index = 0;

	// Code
	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		char current_state_type = p_m_states[state_counter]->state_type;

		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			int new_state = p_m_states[state_counter]->p_transitions[t_counter]->new_state;

			if (new_state == 0)
			{
				// Transition is not allowed - skip out
				continue;
			}

			char new_state_type = p_m_states[new_state - 1]->state_type;

			if (current_state_type == 'S')
			{
				if (new_state_type == 'S')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'D')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'A')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'a';
				}
			}

			if (current_state_type == 'D')
			{
				if (new_state_type == 'S')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'D')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
				if (new_state_type == 'A')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'a';
				}
			}

			if (current_state_type == 'A')
			{
				if (new_state_type == 'S')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'd';
				}
				if (new_state_type == 'D')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'd';
				}
				if (new_state_type == 'A')
				{
					p_m_states[state_counter]->p_transitions[t_counter]->transition_type = 'n';
				}
			}
		
			// Prepare for next transition
			t_index = t_index + 1;
		}
	}
}

void kinetic_scheme::initialise_simulation(myofilaments* set_p_parent_myofilaments)
{
	//! Updates the scheme ready for a simulation

	// Code

/*
	// Updates the parent
	p_parent_myofilaments = set_p_parent_myofilaments;

	// Update the options
	p_cmv_options = p_parent_myofilaments->p_cmv_options;
	
	// Code
	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			int new_state = p_m_states[state_counter]->p_transitions[t_counter]->new_state;

			if (new_state == 0)
			{
				// Transition is not allowed - skip out
				continue;
			}
			p_m_states[state_counter]->p_transitions[t_counter]->p_cmv_options =
				p_cmv_options;
		}
	}

	if (p_cmv_options->rates_dump_file_string != "")
	{
		write_rate_functions_to_file();
	}
*/
}

void kinetic_scheme::write_rate_functions_to_file(void)
{
	/*
	//! Writes rate functions to output file based on data in p_cmv_options

	// Variables
	int counter = 0;			// integer counter

	m_state* p_m_state;			// pointer to an m_state
	transition* p_trans;		// pointer to a transition

	int new_state;				// integer for new state

	FILE* output_file;

	path options_file_path;		// path for options file
	path output_file_path;		// path for output file

	string output_file_string;	// string version of output file

	// Code

	// First, make sure we need to write the file, return if file_string is empty
	if (p_cmv_options->rates_dump_file_string.empty())
		return;

	// Now adjust the rates file
	if (p_cmv_options->rates_dump_relative_to != "")
	{
		path base_dir;

		if (p_cmv_options->rates_dump_relative_to == "this_file")
		{
			options_file_path = path(p_cmv_options->options_file_string);
			base_dir = options_file_path.parent_path();
		}
		else
		{
			base_dir = path(p_cmv_options->rates_dump_relative_to);
		}

		output_file_path = base_dir / p_cmv_options->rates_dump_file_string;
		output_file_string = output_file_path.string();
	}
	else
	{
		output_file_string = p_cmv_options->rates_dump_file_string;
	}

	// Make sure directory exists
	output_file_path = absolute(path(output_file_string));

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: Folder for rates file could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	errno_t err = fopen_s(&output_file, output_file_string.c_str(), "w");

	if (err != 0)
	{
		printf("write_rate_functions_to_file(): %s\ncould not be opened\n",
			output_file_string.c_str());
		exit(1);
	}
	else
	{
		cout << "Writing rate functions to: " << output_file_string << "\n";
	}

	// Write the JSON bracket
	fprintf_s(output_file, "{\n\t\"scheme\":\n\"\n");

	// Cycle through transitions and rates writing the column headers
	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		p_m_state = p_m_states[state_counter];

		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			p_trans = p_m_state->p_transitions[t_counter];
			new_state = p_trans->new_state;

			if (new_state > 0)
			{
				// It's a transition
				counter = counter + 1;
				if (counter == 1)
					fprintf_s(output_file, "x\tr_%i", counter);
				else
					fprintf_s(output_file, "\tr_%i", counter);
			}
		}
	}
	fprintf_s(output_file, "\n");

	// Cycle through x values and bins
	for (double x = p_cmv_options->bin_min ; x <= p_cmv_options->bin_max ;
						x = x + p_cmv_options->bin_width)
	{
		fprintf_s(output_file, "%8g", x);

		for (int state_counter = 0; state_counter < no_of_states; state_counter++)
		{
			p_m_state = p_m_states[state_counter];

			for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
			{
				p_trans = p_m_state->p_transitions[t_counter];
				new_state = p_trans->new_state;

				if (new_state > 0)
				{
					// It's a transition
					double x_ext = p_m_state->extension;
					double rate = p_trans->calculate_rate(x, x_ext, 0, p_parent_myofilaments->p_parent_hs->hs_length);

					fprintf_s(output_file, "\t%8g", rate);
				}
			}
		}
		fprintf_s(output_file, "\n");
	}

	// Close the JSON bracket
	fprintf_s(output_file, "\"\n}\n");

	fclose(output_file);
	*/
}

void kinetic_scheme::write_kinetic_scheme_to_file(char output_file_string[])
{
	//! Writes kinetics scheme to output file

	// Variables
	FILE* output_file;

	// Code

	// Make sure directory exists
	path output_file_path(output_file_string);

	if (!(is_directory(output_file_path.parent_path())))
	{
		if (create_directories(output_file_path.parent_path()))
		{
			std::cout << "\nCreating folder: " << output_file_path.string() << "\n";
		}
		else
		{
			std::cout << "\nError: Results folder could not be created: " <<
				output_file_path.parent_path().string() << "\n";
			exit(1);
		}
	}

	// Check file can be opened, abort if not
	errno_t err = fopen_s(&output_file, output_file_string, "w");
	if (err != 0)
	{
		printf("write_kinetic_scheme_to_file(): %s\ncould not be opened\n",
			output_file_string);
		exit(1);
	}

	// Kinetic scheme information
	fprintf_s(output_file, "{\n");
	fprintf_s(output_file, "\t\"m_kinetics\": {\n");
	fprintf_s(output_file, "\t\t\"no_of_states\": %i,\n", no_of_states);
	fprintf_s(output_file, "\t\t\"max_no_of_transitions\": %i,\n", max_no_of_transitions);
	fprintf_s(output_file, "\t\t\"scheme\": [\n");

	for (int state_counter = 0; state_counter < no_of_states; state_counter++)
	{
		fprintf_s(output_file, "\t\t{\n");
		fprintf_s(output_file, "\t\t\t\"number\": %i,\n", p_m_states[state_counter]->state_number);
		fprintf_s(output_file, "\t\t\t\"type\": \"%c\",\n", p_m_states[state_counter]->state_type);
		fprintf_s(output_file, "\t\t\t\"extension\": \"%g\",\n", p_m_states[state_counter]->extension);
		fprintf_s(output_file, "\t\t\t\"transition\":\n");
		fprintf_s(output_file, "\t\t\t[\n");

		for (int t_counter = 0; t_counter < max_no_of_transitions; t_counter++)
		{
			fprintf_s(output_file, "\t\t\t\t{\n");
			fprintf_s(output_file, "\t\t\t\t\t\"new_state\": %i,\n",
				p_m_states[state_counter]->p_transitions[t_counter]->new_state);
			fprintf_s(output_file, "\t\t\t\t\t\"transition_type\": \"%c\",\n",
				p_m_states[state_counter]->p_transitions[t_counter]->transition_type);
			fprintf_s(output_file, "\t\t\t\t\t\"rate_type\": \"%s\",\n",
				p_m_states[state_counter]->p_transitions[t_counter]->rate_type);
			fprintf_s(output_file, "\t\t\t\t\t\"rate_parameters\": [");

			for (int p_counter = 0; p_counter < MAX_NO_OF_RATE_PARAMETERS; p_counter++)
			{
				fprintf_s(output_file, "%g",
					gsl_vector_get(p_m_states[state_counter]->p_transitions[t_counter]->rate_parameters,
						p_counter));
				if (p_counter == (MAX_NO_OF_RATE_PARAMETERS - 1))
					fprintf_s(output_file, "]\n");
				else
					fprintf_s(output_file, ", ");
			}
			fprintf_s(output_file, "\t\t\t\t}");

			if (t_counter == (max_no_of_transitions - 1))
				fprintf_s(output_file, "\n");
			else
				fprintf_s(output_file, ",\n");
		}
		fprintf_s(output_file, "\t\t\t]\n");

		fprintf_s(output_file, "\t\t}");
		if (state_counter == (no_of_states - 1))
			fprintf_s(output_file, "\n");
		else
			fprintf_s(output_file, ",\n");
	}
	fprintf_s(output_file, "\t}\n");
	fprintf_s(output_file, "}\n");

	// Tidy up
	fclose(output_file);
}
