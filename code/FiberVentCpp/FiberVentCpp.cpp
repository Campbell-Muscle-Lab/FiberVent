/**
 * @file    FiberVentCpp.cpp
 * @brief   Core file for the FiberVentCpp model
 * @author  Ken Campbell
 */

#include <stdio.h>
#include <string>
#include <iostream>

// Includes
#include "cmv_system.h"

// Pointers
cmv_system* p_cmv_system;       /**< pointer to a cmv_system */

using namespace std;

int main(int argc, char* argv[])
{
    /**
    Main function
    + the entry point for FiberVentCpp
    */
    
    // Variables
    string model_file_string;
    string options_file_string;
    string protocol_file_string;
    string results_file_string;
    string system_id;

    // Set inputs
    model_file_string = argv[1];
    options_file_string = argv[2];
    protocol_file_string = argv[3];
    results_file_string = argv[4];
    system_id = argv[5];

    // Initialize
    p_cmv_system = new cmv_system(model_file_string, options_file_string, stoi(system_id));

    p_cmv_system->run_simulation(protocol_file_string, results_file_string);

    // Tidy up
    delete p_cmv_system;

    // Close
    printf("Closing FiberVentCpp\n");

    return(1);

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
