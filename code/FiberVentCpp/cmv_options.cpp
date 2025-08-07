/**
/* @file		cmv_options.cpp
/* @brief		Source file for a cmv_options object
/* @author		Ken Campbell
*/

#include <iostream>
#include <filesystem>

#include "JSON_functions.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"

#include "cmv_options.h"
#include "FiberSim_options.h"
#include "MyoSim_options.h"

#include "gsl_math.h"

using namespace std::filesystem;
using namespace std;

// Constructor
cmv_options::cmv_options(string set_options_file_string)
{
	// Initialise

	// Code

	// Initialise variables
	options_file_string = set_options_file_string;

	// Set pointers
	p_FiberSim_options = NULL;
	p_MyoSim_options = NULL;

	// Set defaults
	write_every_s = GSL_NAN;

	// Now update from file
	initialise_options_from_JSON_file(options_file_string);
}

// Destructor
cmv_options::~cmv_options(void)
{
	// Initialise

	// Code
	delete p_FiberSim_options;
	delete p_MyoSim_options;
}

void cmv_options::initialise_options_from_JSON_file(string options_file_string)
{
	//! Code initialises options from file

	// Variables
	errno_t file_error;
	FILE* fp;
	char readBuffer[65536];

	// Code
	file_error = fopen_s(&fp, options_file_string.c_str(), "rb");
	if (file_error != 0)
	{
		cout << "Error opening options file: " << options_file_string;
		exit(1);
	}

	rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

	rapidjson::Document doc;
	doc.ParseStream(is);

	fclose(fp);

	// Now trying to read file
	cout << "Parsing options file: " << options_file_string << "\n";

	JSON_functions::check_JSON_member_object(doc, "FiberVent_options");
	const rapidjson::Value& opt = doc["FiberVent_options"];

	// Check for muscle
	if (JSON_functions::check_JSON_member_exists(opt, "muscle"))
	{
		const rapidjson::Value& mus = opt["muscle"];

		// Load protocol variables
		JSON_functions::check_JSON_member_number(mus, "force_control_max_delta_hs_length");
		force_control_max_delta_hs_length = mus["force_control_max_delta_hs_length"].GetDouble();
	}

	// Check for hemi_vent
	if (JSON_functions::check_JSON_member_exists(opt, "hemi_vent"))
	{
		const rapidjson::Value& hv = opt["hemi_vent"];

		// Load protocol variables
		JSON_functions::check_JSON_member_string(hv, "thick_wall_approximation");
		hv_thick_wall_approximation = hv["thick_wall_approximation"].GetString();
	}

	// Check for results
	if (JSON_functions::check_JSON_member_exists(opt, "results"))
	{
		const rapidjson::Value& res = opt["results"];

		JSON_functions::check_JSON_member_number(res, "beat_length_s");
		beat_length_s = res["beat_length_s"].GetDouble();

		JSON_functions::check_JSON_member_number(res, "summary_time_step_s");
		summary_time_step_s = res["summary_time_step_s"].GetDouble();

		if (JSON_functions::check_JSON_member_exists(res, "write_every_s"))
		{
			write_every_s = res["write_every_s"].GetDouble();
		}
	}

	// FiberSim
	if (JSON_functions::check_JSON_member_exists(opt, "FiberSim"))
	{
		const rapidjson::Value& fs = opt["FiberSim"];

		p_FiberSim_options = new FiberSim_options(fs);
	}

	// MyoSim
	if (JSON_functions::check_JSON_member_exists(opt, "MyoSim"))
	{
		const rapidjson::Value& myo = opt["MyoSim"];

		// Load protocol variables
		JSON_functions::check_JSON_member_number(myo, "bin_min");
		bin_min = myo["bin_min"].GetDouble();

		JSON_functions::check_JSON_member_number(myo, "bin_max");
		bin_max = myo["bin_max"].GetDouble();

		JSON_functions::check_JSON_member_number(myo, "bin_width");
		bin_width = myo["bin_width"].GetDouble();

		JSON_functions::check_JSON_member_number(myo, "max_rate");
		max_rate = myo["max_rate"].GetDouble();

		// Check for rates dump
		if (JSON_functions::check_JSON_member_exists(myo, "rates_dump"))
		{
			const rapidjson::Value& rd = myo["rates_dump"];

			if (JSON_functions::check_JSON_member_exists(rd, "relative_to"))
			{
				rates_dump_relative_to = rd["relative_to"].GetString();
			}
			else
			{
				rates_dump_relative_to = "";
			}

			JSON_functions::check_JSON_member_string(rd, "file_string");
			rates_dump_file_string = rd["file_string"].GetString();
		}
		else
		{
			rates_dump_relative_to = "";
			rates_dump_file_string = "";
		}

		// Check for cb dump
		if (JSON_functions::check_JSON_member_exists(myo, "cb_dump"))
		{
			const rapidjson::Value& cd = myo["cb_dump"];

			if (JSON_functions::check_JSON_member_exists(cd, "relative_to"))
			{
				cb_dump_relative_to = cd["relative_to"].GetString();
			}
			else
			{
				cb_dump_relative_to = "";
			}

			JSON_functions::check_JSON_member_string(cd, "file_string");
			cb_dump_file_string = cd["file_string"].GetString();
		}
		else
		{
			cb_dump_relative_to = "";
			cb_dump_file_string = "";
		}
	}
}
