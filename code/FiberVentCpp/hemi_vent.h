#pragma once

/**
/* @file		hemi_vent.h
/* @brief		Header file for a hemi_vent object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>

#include "global_definitions.h"

// Forward declararations
class cmv_model;
class cmv_system;
class cmv_results;
class cmv_options;

class circulation;

class valve;
class muscle;



class hemi_vent
{
public:
	/**
	 * Constructor
	 */
	hemi_vent(circulation* set_p_parent_circulation);

	/**
	 * Destructor
	 */
	 ~hemi_vent(void);

	 // Variables

	cmv_system* p_parent_cmv_system;		/**< pointer to parent cmv_system */

	cmv_model* p_cmv_model;					/**< pointer to cmv_model object */

	cmv_results* p_cmv_results_beat;		/**< pointer to cmv_results object
													holding data at full time
													resolution for the last
													beat */

	cmv_options* p_cmv_options;				/**< pointer to cmv_options object */

	circulation* p_parent_circulation;		/**< pointer the parent circulation */

	valve* p_av;							/**< pointer to the aortic valve */

	valve* p_mv;							/**< pointer to the mitral valve */

	muscle* p_muscle;						/**< pointer to child muscle */

	double vent_wall_density;				/**< double with wall density in kg m^-3 */

	double vent_wall_volume;				/**< double with wall volume in liters */

	double vent_chamber_height;				/**< double with internal base to apex distance
													in meters calculated so that
													h = radius * vent_z_scale * 
														(hsl/hsl_ref)^vent_z_exp
											*/

	double vent_z_scale;					/**< double for chamber height as above */

	double vent_z_exp;						/**< double for chamber height as above */

	double vent_chamber_volume;				/**< double with chamner volume in liters */

	double vent_chamber_pressure;			/**< double with pressure in the ventricle in mm Hg */

	double vent_chamber_radius;				/**< double with chamber inner radius in m */

	double vent_wall_thickness;				/**< double with wall thickness in m */

	double vent_n_hs;						/**< double with the number of half-sarcomeres
													around the ventricular circumference */

	double vent_circumference;				/**< double with circumference of ventricle
													in m */

	double vent_thick_wall_multiplier;		/**< double with
													1.0 if using thick wall approximation
													0.0 if not */

	double vent_stroke_work_J;				/**< double with stroke work in J for a cardiac cycle */

	double vent_stroke_energy_used_J;		/**< double with energy_used in J for a cardiac cycle */

	double vent_efficiency;					/**< double efficiency for a cardiac cycle */

	double vent_ejection_fraction;			/**< double with ejection fraction for a cardiac cycle */

	double vent_ATP_used_per_s;				/**< double with ATP used per s by the ventricle */

	double vent_stroke_volume;				/**< double with vent stroke volume
													in liters */

	double vent_cardiac_output;				/**< double with vent cardiac output
													in liter per minute */

	// Other functions
	void initialise_simulation(void);

	bool implement_time_step(double time_step_s);

	double return_wall_thickness_for_chamber_volume(double cv);

	double wall_thickness_root_finder(double cv);

	double return_internal_radius_for_chamber_volume(double cv);

	double return_chamber_height(double r);

	double return_lv_circumference_for_chamber_volume(double cv);

	double return_pressure_for_chamber_volume(double cv, double time_step_s);

	void update_chamber_volume(double new_volume, double time_step_s);

	void update_beat_metrics(void);

	void calculate_vent_ATP_used_per_s(void);
};
