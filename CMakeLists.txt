cmake_minimum_required(VERSION 3.31)
project(bertini_real VERSION 1.9.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


configure_file(brconfig.h.in config.h)







# Append the cmake/ directory to the CMAKE_MODULE_PATH
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# sadly, these don't have cmake-provided find calls.
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)



# i wrote this one, and my code for finding it probably sucks
find_package(bertini1 REQUIRED)




# Find required packages
find_package(MPI REQUIRED)
find_package(Boost 1.53 REQUIRED
    COMPONENTS
        filesystem
        timer
    )




set(CMAKE_CXX_FLAGS ${Bertini1_CPP_FLAGS})

include(CMakePrintHelpers)
cmake_print_variables(GMP_INCLUDES)
cmake_print_variables(MPFR_INCLUDES)
cmake_print_variables(Bertini1_INCLUDES)
cmake_print_variables(Bertini1_LIBRARIES)
message("include dir is " ${CMAKE_CURRENT_SOURCE_DIR}/include)



# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include

    "${PROJECT_BINARY_DIR}" 

    ${Boost_INCLUDE_DIRS}


    ${GMP_INCLUDES}
    ${MPFR_INCLUDE_DIR}

    ${Bertini1_INCLUDES}
)



# Source files
file(GLOB_RECURSE HDR_FILES include/*.hpp)
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# Executable
add_executable(bertini_real ${SRC_FILES} ${HDR_FILES})
target_link_libraries(bertini_real
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
    ${Bertini1_LIBRARIES}
)




# Installation
install(TARGETS bertini_real DESTINATION bin)
