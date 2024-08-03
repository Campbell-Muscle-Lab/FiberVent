/**
/* @file		cmv_protocol.cpp
/* @brief		Source file for a cmv_protocol object
/* @author		Ken Campbell
*/

#include <iostream>
#include <filesystem>

#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

#include "cmv_system.h"
#include "cmv_protocol.h"
#include "activation.h"
#include "perturbation.h"

#include "gsl_math.h"

using namespace std::filesystem;
using namespace std;

struct perturbation_struct {
	string class_name;
	string variable;
	double t_start_s;
	double t_stop_s;
	double total_change;
};

// Constructor
cmv_protocol::cmv_protocol(cmv_system* set_p_cmv_system, string set_protocol_file_string)
{
	// Initialise

	// Code

	// Initialise variables
	p_cmv_system = set_p_cmv_system;

	protocol_file_string = set_protocol_file_string;
	time_step_s = 0.0;
	no_of_time_steps = 0;

	no_of_activations = 0;
	no_of_perturbations = 0;

	// Null the activation pointers
	for (int i = 0; i < MAX_NO_OF_ACTIVATIONS; i++)
	{
		p_activation[i] = NULL;
	}

	// Null the perturbation pointers
	for (int i = 0; i < MAX_NO_OF_PERTURBATIONS; i++)
	{
		p_perturbation[i] = NULL;
	}

	// Now update from file
	initialise_protocol_from_JSON_file(protocol_file_string);
}

// Destructor
cmv_protocol::~cmv_protocol(void)
{
	// Initialise

	// Code

	// Tidy up
	for (int i = 0; i < MAX_NO_OF_ACTIVATIONS; i++)
	{
		if (p_activation[i] != NULL)
			delete p_activation[i];
	}

	for (int i = 0; i < MAX_NO_OF_PERTURBATIONS; i++)
	{
		if (p_perturbation[i] != NULL)
			delete p_perturbation[i];
	}
}

void cmv_protocol::initialise_protocol_from_JSON_file(string protocol_file_string)
{
	//! Code initialises a protocol from file

	// Variables
	errno_t file_error;
	FILE* fp;
	char readBuffer[65536];

	string activation_type;
	double t_start_s;
	double t_stop_s;

	// Code
	file_error = fopen_s(&fp, protocol_file_string.c_str(), "rb");
	if (file_error != 0)
	{
		cout << "Error opening protocol file: " << protocol_file_string;
		exit(1);
	}

	rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

	rapidjson::Document doc;
	doc.ParseStream(is);

	fclose(fp);

	// Now trying to read file
	cout << "Parsing protocol file: " << protocol_file_string << "\n";

	// Load the protocol object
	JSON_functions::check_JSON_member_object(doc, "protocol");
	const rapidjson::Value& prot = doc["protocol"];

	// Load protocol variables
	JSON_functions::check_JSON_member_number(prot, "time_step_s");
	time_step_s = prot["time_step_s"].GetDouble();

	JSON_functions::check_JSON_member_int(prot, "no_of_time_steps");
	no_of_time_steps = prot["no_of_time_steps"].GetInt();

	// Check for activations
	if (JSON_functions::check_JSON_member_exists(doc, "activation"))
	{
		JSON_functions::check_JSON_member_array(doc, "activation");
		const rapidjson::Value& act = doc["activation"];

		// Pull aray
		for (rapidjson::SizeType i = 0; i < act.Size(); i++)
		{
			const rapidjson::Value& temp = act[i];

			JSON_functions::check_JSON_member_string(temp, "type");
			activation_type = temp["type"].GetString();
			
			JSON_functions::check_JSON_member_number(temp, "t_start_s");
			t_start_s = temp["t_start_s"].GetDouble();

			JSON_functions::check_JSON_member_number(temp, "t_stop_s");
			t_stop_s = temp["t_stop_s"].GetDouble();

			p_activation[i] = new activation(activation_type, t_start_s, t_stop_s);

			no_of_activations = no_of_activations + 1;
		}
	}

	// Check for perturbations
	if (JSON_functions::check_JSON_member_exists(doc, "perturbation"))
	{
		JSON_functions::check_JSON_member_array(doc, "perturbation");
		const rapidjson::Value& pert = doc["perturbation"];

		perturbation_struct* p_struct = new perturbation_struct;

		// Pull array
		for (rapidjson::SizeType i = 0; i < pert.Size(); i++)
		{
			const rapidjson::Value& temp = pert[i];

			JSON_functions::check_JSON_member_string(temp, "class");
			p_struct->class_name = temp["class"].GetString();

			JSON_functions::check_JSON_member_string(temp, "variable");
			p_struct->variable = temp["variable"].GetString();

			JSON_functions::check_JSON_member_number(temp, "t_start_s");
			p_struct->t_start_s = temp["t_start_s"].GetDouble();

			JSON_functions::check_JSON_member_number(temp, "t_stop_s");
			p_struct->t_stop_s = temp["t_stop_s"].GetDouble();

			JSON_functions::check_JSON_member_number(temp, "total_change");
			p_struct->total_change = temp["total_change"].GetDouble();

			p_perturbation[i] = new perturbation(this, p_struct);

			no_of_perturbations = no_of_perturbations + 1;
		}

		delete p_struct;
	}
}

double cmv_protocol::return_activation(string activation_type, double time_s)
{
	//! Function returns the activation status for a given type
	
	// Variables
	double activation = 0.0;

	// Code
	for (int i = 0; i < no_of_activations; i++)
	{
		activation = activation + p_activation[i]->return_status(activation_type, time_s);
	}

	activation = GSL_MIN(1.0, activation);

	return (activation);
}

void cmv_protocol::impose_perturbations(double sim_time_s)
{
	//! Function cycles through perturbations, applying changes as required
	
	// Variables

	// Code
	for (int i = 0; i < no_of_perturbations; i++)
	{
		p_perturbation[i]->impose(sim_time_s);
	}
}

