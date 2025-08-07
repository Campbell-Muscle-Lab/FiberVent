#pragma once

/**
/* @file		cmv_results.h
/* @brief		Header file for a cmv_results object
/* @author		Ken Campbell
*/

#include "stdio.h"
#include <iostream>
#include <string>

#include "global_definitions.h"
#include "cmv_system.h"

#include "gsl_vector.h"

using namespace std;

struct stats_structure;

class cmv_system;
class cmv_options;

class cmv_results
{
public:
	/**
	 * Constructor
	 */
	cmv_results(cmv_system* set_p_parent_cmv_system, int set_no_of_time_points);

	/**
	* Destructor
	*/
	~cmv_results(void);

	// Variables
	cmv_system* p_parent_cmv_system;		/**< pointer to parent cmv_system */

	cmv_options* p_cmv_options;				/**< pointer to cmv options */

	std::string results_fields[MAX_NO_OF_RESULT_FIELDS];
											/**< array of strings defining the
													data fields in a results object */

	gsl_vector* gsl_results_vectors[MAX_NO_OF_RESULT_FIELDS];
											/**< array of gsl_vectors holding data */

	int no_of_defined_results_fields;		/**< integer defining the number of
													results fields that have been defined */

	double* p_data_sources[MAX_NO_OF_RESULT_FIELDS];
											/**< array of pointers to doubles  holding
													the sources of the data */

	int no_of_time_points;					/**< integer defining the number of
													time-points in a result file */

	int no_of_beats;						/**< integer counting number of beats
													written to record */

	int time_field_index;					/**< integer holding the index for the
													time field */

	int new_beat_field_index;				/**< integer holding the index for the
													new_beat field */

	int hr_heart_rate_bpm_field_index;		/**< integer holding the index for the
													heart_rate_bpm */

	int pressure_vent_field_index;			/**< integer holding the index for the
													ventricular pressure */

	int volume_vent_field_index;			/**< integer holding the index for the
													ventricular volume */

	int pressure_arteries_field_index;		/**< integer holding the index for the
													arterial pressure */

	int pressure_veins_field_index;			/**< integer holding the index for the
													venous pressure */

	int flow_mitral_valve_field_index;		/**< integer holding the index for flow
													through the mitral valve */

	int flow_aortic_valve_field_index;		/**< integer holding the index for flow
													through the aortic valve */

	int muscle_length_field_index;			/**< integer holding the index for the
													muscle length */


	int muscle_stress_int_pas_field_index;	/**< integer holding the index for the
													internal passive stress */

	int muscle_ATP_flux_field_index;		/**< integer holding the index for the
													muscle ATPase field */

	int vent_stroke_work_field_index;		/**< integer holding the index for the
													vent stroke work field */

	int vent_stroke_energy_used_field_index;		/**< integer holding the index for the
													vent energy used field */

	int vent_efficiency_field_index;		/**< integer holding the index for the
													vent efficiency field */

	int vent_ejection_fraction_field_index;	/**< integer holding the index for the
													vent ejection fraction field */

	int vent_ATP_used_per_s_field_index;	/**< integer holding the index for the
													vent ATP used per s field */

	int vent_stroke_volume_field_index;		/**< integer holding the index for the
													vent stroke volume field */

	int vent_cardiac_output_field_index;	/**< integer holding the index for the
													vent cardiac output volume field */

	// Functions

	void add_results_field(std::string field_name, double* p_double);
											/**< function adds a double to the results
													object */

	void save_status(int t_index);

	int write_data_to_file(string output_file_string);
											/**< write data to file */

	//void calculate_beat_metrics(int t_beat_index);

	void calculate_sub_vector_statistics(gsl_vector* gsl_v,
		int start_index, int stop_index,
		stats_structure* p_stats_structure);

	double return_stroke_work(int start_t_index, int stop_t_index);

	double return_energy_used(int start_t_index, int stop_t_index);

	void backfill_beat_data(gsl_vector* gsl_v, double value,
		int start_t_index, int stop_t_index);

};