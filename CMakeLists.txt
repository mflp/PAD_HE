# Copyright (c) Microsoft Corporation. All rights reserved.
# Licensed under the MIT license.

cmake_minimum_required(VERSION 3.12)
project(SEALExamples VERSION 3.6.2 LANGUAGES CXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )

# If not called from root CMakeLists.txt
if(NOT DEFINED SEAL_BUILD_EXAMPLES)
    set(SEAL_BUILD_EXAMPLES ON)

    # Import Microsoft SEAL
    find_package(SEAL 3.6.2 EXACT REQUIRED)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()

if(SEAL_BUILD_EXAMPLES)
    add_executable(sealexamples)
    target_sources(sealexamples
        PRIVATE
            ${CMAKE_CURRENT_LIST_DIR}/examples.cpp        # The class for the console menu and guiding the user to the preferred application  
            ${CMAKE_CURRENT_LIST_DIR}/1_bfv_basics.cpp    # DO-NOT-CONSIDER-WITHIN-THE-SCOPE-OF-THE-APPLICATION
			${CMAKE_CURRENT_LIST_DIR}/encodersplain.cpp   # The main application class of the crypto application
            ${CMAKE_CURRENT_LIST_DIR}/3_levels.cpp        # DO-NOT-CONSIDER-WITHIN-THE-SCOPE-OF-THE-APPLICATION 
            ${CMAKE_CURRENT_LIST_DIR}/4_ckks_basics.cpp   # DO-NOT-CONSIDER-WITHIN-THE-SCOPE-OF-THE-APPLICATION
            ${CMAKE_CURRENT_LIST_DIR}/5_rotation.cpp      # DO-NOT-CONSIDER-WITHIN-THE-SCOPE-OF-THE-APPLICATION
            ${CMAKE_CURRENT_LIST_DIR}/6_serialization.cpp # DO-NOT-CONSIDER-WITHIN-THE-SCOPE-OF-THE-APPLICATION
            ${CMAKE_CURRENT_LIST_DIR}/7_performance.cpp   # DO-NOT-CONSIDER-WITHIN-THE-SCOPE-OF-THE-APPLICATION
			${CMAKE_CURRENT_LIST_DIR}/rawplain.cpp        	   # The class which reads and stores the plain matrices (1)				
			${CMAKE_CURRENT_LIST_DIR}/generateplaintextMRP.cpp # The functions used to prepare plain matrices in MRP format (2.0)
			${CMAKE_CURRENT_LIST_DIR}/prepareVecMatMRP.cpp	   # The class which prepares the read matrices in MRP format   (2.1) 	
			${CMAKE_CURRENT_LIST_DIR}/initializationPLCP.cpp   # The class containing the functions which convert the prepared matrices to SEAL Plaintext and Ciphertext ob.(3)			
			${CMAKE_CURRENT_LIST_DIR}/plcpOperations.cpp       # The class containing the functions which does Ciphertext-Ciphertext and Plaintext-Ciphertext arithmetic and algebraic operations (4.0)  
			${CMAKE_CURRENT_LIST_DIR}/encryptedAppx.cpp        # The class containing the functions which does the Chebyshev approximation with different assumptions (4.1)
			${CMAKE_CURRENT_LIST_DIR}/secretShare.cpp          # The class containing the functions which does the secret sharing (4.2)
			${CMAKE_CURRENT_LIST_DIR}/applyPLCPSimulation.cpp  # The class containing the crypto application functions for each targeted equation (5)
			${CMAKE_CURRENT_LIST_DIR}/printCont.cpp            # The class containing the printing functions for control purposes (6)
    )

    if(TARGET SEAL::seal)
        target_link_libraries(sealexamples PRIVATE SEAL::seal) # original line
    elseif(TARGET SEAL::seal_shared)
        target_link_libraries (sealexamples PRIVATE SEAL::seal_shared) # original line
    else()
        message(FATAL_ERROR "Cannot find target SEAL::seal or SEAL::seal_shared")
    endif()
endif()
