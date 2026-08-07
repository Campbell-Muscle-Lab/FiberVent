#include "stdio.h"
#include <iostream>
#include <string>

#include "global_definitions.h"

#include "gsl_vector.h"
#include "gsl_matrix.h"
#include "gsl_math.h"

#include "FiberSim_model.h"
#include "cmv_model.h"
#include "muscle.h"
#include "FiberSim_muscle.h"
#include "FiberSim_half_sarcomere.h"

#include "FiberSim_half_sarcomere_onnx.h"

#include <onnxruntime_cxx_api.h>

// Constructor
FiberSim_half_sarcomere_onnx::FiberSim_half_sarcomere_onnx(
    FiberSim_half_sarcomere* set_p_FiberSim_half_sarcomere, Ort::Session* set_p_session)
{
    // Constructor

    // Set variables
    p_parent_hs = set_p_FiberSim_half_sarcomere;
    p_session = set_p_session;

    // Get the shape of the inputs
    auto input_type_info = p_session->GetInputTypeInfo(0);
    auto input_tensor_info = input_type_info.GetTensorTypeAndShapeInfo();
    auto input_shape = input_tensor_info.GetShape();

    no_of_predictor_variables = (size_t)(input_shape[1]);
    no_of_predictor_time_points = (size_t)(input_shape[2]);

    printf("no_of_predictor_variables: %zu\n", no_of_predictor_variables);
    printf("no_of_predictor_time_points: %zu\n", no_of_predictor_time_points);

    // Reserve and allocate gsl_matrix
    gsl_predictor_matrix = gsl_matrix_float_alloc(
        no_of_predictor_variables, no_of_predictor_time_points);

    gsl_matrix_float_set_all(gsl_predictor_matrix, GSL_NAN);

    // Create the input tensor
    std::vector<int64_t> tensor_shape = { 1,
        static_cast<int64_t>(gsl_predictor_matrix->size1),
        static_cast<int64_t>(gsl_predictor_matrix->size2) };

    input_tensor_from_gsl_matrix_float = Ort::Value::CreateTensor<float>(
        memory_info,
        gsl_predictor_matrix->data, // pointer to GSL data
        gsl_predictor_matrix->size1 * gsl_predictor_matrix->size2, // number of elements
        tensor_shape.data(),
        tensor_shape.size());

    printf("test mult: %zu %zu\n", gsl_predictor_matrix->size1, gsl_predictor_matrix->size2);

    std::cout << gsl_predictor_matrix->size1 << " "
        << gsl_predictor_matrix->size2 << " "
        << gsl_predictor_matrix->tda << std::endl;

    // Get the shape of the outputs
    auto output_type_info = p_session->GetOutputTypeInfo(0);
    auto output_tensor_info = output_type_info.GetTensorTypeAndShapeInfo();
    std::vector<int64_t> output_shape = output_tensor_info.GetShape();

    no_of_predicted_time_points = (size_t)(output_shape[0]);
    no_of_predicted_variables = (size_t)(output_shape[1]);


    printf("no_of_predicted_variables: %zu\n", no_of_predicted_variables);
    printf("no_of_predicted_time_points: %zu\n", no_of_predicted_time_points);

    // Reserve and allocate gsl_matrix
    gsl_predicted_matrix = gsl_matrix_float_alloc(
        no_of_predicted_variables, no_of_predicted_time_points);

    gsl_matrix_float_set_all(gsl_predicted_matrix, GSL_NAN);

    //// Set the input and output names

    Ort::AllocatorWithDefaultOptions allocator;

    for (size_t i = 0; i < p_session->GetInputCount(); i++)
    {
        auto nameAllocated = p_session->GetInputNameAllocated(i, allocator);
        inputNamesString.push_back(nameAllocated.get());
    }
    for (const auto& name : inputNamesString)
    {
        input_names.push_back(name.c_str());
    }

    for (size_t i = 0; i < p_session->GetOutputCount(); i++)
    {
        auto nameAllocated = p_session->GetOutputNameAllocated(i, allocator);
        outputNamesString.push_back(nameAllocated.get());
    }
    for (const auto& name : outputNamesString)
    {
        output_names.push_back(name.c_str());
    }

    // Deduce an index that will be useful later
    std::vector<string>* p_target_fields = &(p_parent_hs->p_parent_fs_muscle->p_parent_muscle->p_cmv_model->p_fs_model->onnx_target_fields);
    
    for (size_t i = 0; i < no_of_predicted_variables; i++)
    {
        if ((*p_target_fields)[i] == "hs_1_force")
        {
            idx_for_hs_1_force_prediction = i;
        }
    }
    printf("Target_index for hs_1_force: %zu\n", idx_for_hs_1_force_prediction);
}

FiberSim_half_sarcomere_onnx::~FiberSim_half_sarcomere_onnx(void)
{
    // Destructor

    // Free matrix
    gsl_matrix_float_free(gsl_predictor_matrix);
    gsl_matrix_float_free(gsl_predicted_matrix);

    printf("FS_ML_predictors::destructor\n");
}

// Other functions
void FiberSim_half_sarcomere_onnx::update_predictor_matrix(std::string update_mode)
{
    //! Updates gsl_predictor_matrix
    //! If initializing == true, fill the whole matrix
    //! else, just shuffle the values up a row, and replace the bottom row
    
    // Variables
    std::vector<string>* p_input_fields = &(p_parent_hs->p_parent_fs_muscle->p_parent_muscle->p_cmv_model->p_fs_model->onnx_input_fields);
                                            // pointer to the input fields

    double* p_double;                       // pointer to a double

    // Code

    // Cycle through the columns, trying to find the right pointer
    for (size_t i = 0; i < no_of_predictor_variables; i++)
    {
        p_double = NULL;

        if ((*p_input_fields)[i] == "hs_1_length")
        {
            p_double = &p_parent_hs->hs_length;
        }
        else if ((*p_input_fields)[i] == "hs_1_pCa")
        {
            p_double = &p_parent_hs->pCa;
        }
        else if ((*p_input_fields)[i] == "hs_1_force")
        {
            p_double = &p_parent_hs->hs_force;
        }

        if (p_double == NULL)
        {
            printf("Error: pointer in FiberSim_half_sracomere_onnx::update_predictor_matrix() was not intialized");
            printf("Error with: %s\n", (*p_input_fields)[i].c_str());
            printf("Now exiting\n\n");
        }

        // Otherwise fill the matrix
        if (update_mode == "initialize")
        {
            // Fill the whole time_sequence
            for (size_t j = 0; j < no_of_predictor_time_points; j++)
            {
                gsl_matrix_float_set(gsl_predictor_matrix, i, j, (float)(*p_double));
            }
        }
        else if (update_mode == "update")
        {
            // Add the new point
            gsl_matrix_float_set(gsl_predictor_matrix, i, (no_of_predictor_time_points - 1),
                (float)*p_double);
        }
        else if (update_mode == "shift")
        {
            // Shift the row one place to the left
            for (size_t j = 1; j < no_of_predictor_time_points; j++)
            {
                gsl_matrix_float_set(gsl_predictor_matrix, i, j - 1,
                    gsl_matrix_float_get(gsl_predictor_matrix, i, j));
            }
        }
    }
}

void FiberSim_half_sarcomere_onnx::update_predicted_matrix(void)
{
    // Returns predictions as a gsl_vector

    auto output_tensors = p_session->Run(
        Ort::RunOptions{ nullptr },
        input_names.data(),
        &input_tensor_from_gsl_matrix_float,
        1,
        output_names.data(),
        output_names.size());

    auto& output_tensor = output_tensors[0];
    auto output_shape = output_tensor.GetTensorTypeAndShapeInfo().GetShape();

    float* output_data = output_tensor.GetTensorMutableData<float>();

    for (size_t i = 0; i < (size_t)output_shape[0]; i++)
    {
        for (size_t j = 0; j < (size_t)output_shape[1]; j++)
        {
            size_t idx = (i * output_shape[1]) + j;

            gsl_matrix_float_set(
                gsl_predicted_matrix,
                i, j,
                output_data[idx]);
        }
    }
}

void FiberSim_half_sarcomere_onnx::print_gsl_float_matrix(gsl_matrix_float* m)
{
    printf("***\n");
    for (size_t r = 0; r < m->size1; r++)
    {
        for (size_t c = 0; c < m->size2; c++)
        {
            printf("%g\t\t", gsl_matrix_float_get(m, r, c));
        }
        printf("\n");
    }
    printf("\n\n");
}
