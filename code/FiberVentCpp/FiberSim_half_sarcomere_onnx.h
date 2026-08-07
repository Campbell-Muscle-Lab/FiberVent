#pragma once
#include <iostream>
#include <string>

#include "global_definitions.h"

#include "gsl_vector.h"
#include "gsl_matrix.h"
#include "gsl_math.h"

class FiberSim_half_sarcomere;

#include <onnxruntime_cxx_api.h>

class FiberSim_half_sarcomere_onnx
{
public:
	/**
	* Constructor
	*/
	FiberSim_half_sarcomere_onnx(
		FiberSim_half_sarcomere* set_p_FiberSim_half_sarcomere, Ort::Session* set_p_session);

	/**
	* Destructor
	*/
	~FiberSim_half_sarcomere_onnx(void);

	// Variables
	FiberSim_half_sarcomere* p_parent_hs;		/** pointer to parent half-sarcomere */

	Ort::Session* p_session;					/** pointer to an ONNX session */

	Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(
		OrtArenaAllocator,
		OrtMemTypeDefault);;					/** memory info for ONNX prediction */

	Ort::Value input_tensor_from_gsl_matrix_float;
												/** input tensor */

	size_t no_of_predictor_variables;			/** number of predictor variables */
	size_t no_of_predictor_time_points;			/** number of time-points for prediction */

	size_t no_of_predicted_variables;			/** number of predicted variables */
	size_t no_of_predicted_time_points;			/** number of time-points that were predicted */

	std::vector<const char*> input_names;		/** input names for model */
	std::vector<const char*> output_names;		/** output names for the model */

	std::vector<std::string> inputNamesString;	/** helper variables to get names in the right */
	std::vector<std::string> outputNamesString;	/** format for the ONNX runtime */

	gsl_matrix_float* gsl_predictor_matrix;		/** pointer to a gsl_matrix that holds the
													data for predictions */

	gsl_matrix_float* gsl_predicted_matrix;		/** pointer to a gsl_matrix that holds the
													predictions */

	size_t idx_for_hs_1_force_prediction;		/** index for hs_1_force prediction */

	// Functions
	void update_predictor_matrix(std::string update_mode);
												/**< Updates the predictor matrix */

	void update_predicted_matrix(void);

	void print_gsl_float_matrix(gsl_matrix_float* m);
};