#pragma once

/**
/* @file		cmv_model.h
/* @brief		Header file for a cmv_model object
/* @author		Ken Campbell
*/

#include "stdio.h"

#include <iostream>
#include <string>

// Definitions for JSON parsing
#ifndef _RAPIDJSON_DOCUMENT
#define _RAPIDJSON_DOCUMENT
#include "rapidjson/document.h"
#endif

#include "global_definitions.h"

// Forward declararations
class half_sarcomere;
class membranes;
class heart_rate;
class kinetic_scheme;

class FiberSim_model;

struct cmv_model_valve_structure;
struct cmv_model_rc_structure;
struct cmv_model_gc_structure;

using namespace std;

class cmv_model
{
public:
	/**
	* Constructor
	*/
	cmv_model(string JSON_model_file_string);

	/**
	* Destructor
	*/
	~cmv_model(void);

	// Variables

	string version_string;				/**< string holding the model version */

	double temperature_K;				/**< double with temperature in K */

	// Circulation
	int circ_no_of_compartments;		/**< integer holding number of compartments */

	double circ_blood_volume;			/**< double holding total blood volume */

	double* circ_resistance;			/**< pointer to array of doubles with
												resistances to individual compartments */

	double* circ_compliance;			/**< pointer to array of doubles with
												compliance of individual compartments */

	double* circ_static;				/**< pointer to array of doubles with
												static pressure in each compartment */

	double* circ_slack_volume;			/**< pointer to array of doubles with
												slack volume of individual compartments */

	double* circ_inertance;				/**< pointer to array of doubles holding
												inertances of individual compartments */

	// Baroreflex
	double baro_P_set;					/**< double with baroreflex set point in mmHg */

	double baro_S;						/**< double with baroreflex S parameter (mmHg)^-1 */

	double baro_k_drive;				/**< double with baroreflex k_drive (s)^-1 */

	double baro_k_recov;				/**< double with baroreflex k_recovery (s)^-1 */

	int baro_P_compartment;				/**< int for the compartment number used to
													calculate the B_a signal, subtract 1
													to get the compartment index */

	cmv_model_rc_structure* p_rc[MAX_NO_OF_REFLEX_CONTROLS];
										/**< array of pointers to reflex controls */

	int no_of_rc_controls;				/**< integer with the number of rc_controls */

	// Growth
	double gr_master_rate;				/**< double with master rate for growth */

	double gr_shrink_eccentric_rate;
										/**< double with eccentric shrinkage rate */

	double gr_shrink_concentric_rate;
										/**< double with concentric shrinkage rate */

	cmv_model_gc_structure* p_gc[MAX_NO_OF_GROWTH_CONTROLS];
										/**< array of pointers to growth controls */

	int no_of_gc_controls;				/**< integer with the number of growth controls */

	// Filling_control

	double filling_venous_pressure;		/**< double holding the venous pressure when filling
												control is activated */

	// Ventricle
	double vent_wall_density;			/**< double with wall density in kg m^-3 */

	double vent_wall_volume;			/**< double with wall volume in liters */

	//double vent_chamber_height;			/**< double with base to apex height in m */
	
	double vent_z_scale;				/**< doubles such that chamber height is
												calculated as
												h = r * vent_z_scale * (hsl / hsl_ref)^vent_z_exp
											*/

	double vent_z_exp;				/**< as above */

	cmv_model_valve_structure* p_av;	/**< pointer to aortic valve structure */

	cmv_model_valve_structure* p_mv;	/**< pointer to mitral valve structure */

	// Heart rate
	double hr_t_RR_interval_s;			/**< double with RR interval in s */

	// Muscle
	double mus_reference_hs_length;		/**< double with hs length used for normalization */

	double mus_initial_ATP_concentration;
										/**< double with initial ATP concentration */

	double mus_delta_G_ATP;				/**< double with energy in Joules
												per mole of ATP */

	double mus_prop_fibrosis;			/**< double with the proportion of the hs
												cross-sectional area that is fibrosis */

	double mus_prop_myofilaments;		/**< Double with the proportion of the
												hs non-fibrotic area that is
												that is myofilaments, remainder
												is mitochondria */
	// Membranes
	double memb_Ca_content;				/**< double with total Ca_concentration
												in M */
	double memb_t_open_s;				/**< double with open duration in s */

	double memb_k_serca;				/**< double describing rate at which
												SERCA pumps Ca in to
												sarcoplasmic reticulum
												in M^-1 s^-1 */

	double memb_k_leak;					/**< double desribing rate at which
												Ca leaks from sarcoplasmic reticulum
												in M s^-1 */

	double memb_k_active;				/**< double describing rate at which
												Ca is released from sarcoplasmic reticulum
												when the release channels are open
												in M s^-1 */

	double memb_SERCA_Ca_per_ATP;		/**< double describing the number of Ca ions
												moved by SERCA for each ATP molecule
												Should be fixed at 2.0 */

	double memb_SERCA_ATPase_excess_Ca_factor;
										/**< double describing the additional Ca that
												must be pumped by SERCA for a given
												change in cytosolic Ca */

	// Myofilaments

	FiberSim_model* p_fs_model;			/**< pointer to a FiberSim model */

	double myof_cb_number_density;		/**< double describing the number of
												cross-bridges in a half-sarcomere
												with a cross-sectional area of 1 m^2 */

	double myof_prop_fibrosis;			/**< double describing the proportional
												cross-sectional area occupied by
												fibrosis */

	double myof_prop_myofilaments;		/**< double describing the proportional
												area of myofibrils occupied by
												myofilaments - the remainder is
												mitochondria */


	double myof_int_pas_sigma;			/**< double with sigma for intracellular
												passive stress exponential
												in N m^-2 */

	double myof_int_pas_L;				/**< double with curvature for intracellular
												passive stress exponential
												exponential in nm */

	double myof_int_pas_slack_hsl;		/**< double with slack length for
												intraceullar passive stress
												exponential in nm */

	double myof_ext_pas_sigma;			/**< double with sigma for extracellular
												passive stress exponential
												in N m^-2 */

	double myof_ext_pas_L;				/**< double with curvature for extracellular
												passive stress exponential
												exponential in nm */

	double myof_ext_pas_slack_hsl;		/**< double with slack length for
												extraceullar passive stress
												exponential in nm */

	double myof_fil_compliance_factor;	/**< double holding the myofilament compliance
													factor. When the half-sarcomere changes
													length by x, cross-bridges are displaced
													by x * this variable */

	double myof_k_cb;					/**< double describing cross-bridge stiffness
												in N m^-1 */

	kinetic_scheme* p_m_scheme;			/**< pointer to the kinetic scheme
												for myosin */

	double myof_a_k_on;					/**< double describing actin k_on */

	double myof_a_k_off;				/**< double describing actin k_off */

	double myof_a_k_coop;				/**< double describing actin k_coop */

	double myof_thick_fil_length;		/**< double describing thick filament length */

	double myof_bare_zone_length;		/**< double describing bare zone length */

	double myof_thin_fil_length;		/**< double describing thin filament length */

	// Mitochondria
	double mito_ATP_generation_rate;	/**< double with rate mitochondria
												generate ATP, in M s^-1 m^-3 */

	// Other functions

	/**
	/* Function initialises a model object from file
	*/
	void initialise_model_from_JSON_file(string JSON_model_file_string);
};