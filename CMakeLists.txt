cmake_minimum_required(VERSION 3.31)
project(bertini_real VERSION 1.9.0 LANGUAGES C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


configure_file(brconfig.h.in config.h)







# Append the cmake/ directory to the CMAKE_MODULE_PATH
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")




FIND_PACKAGE(FLEX REQUIRED)
FIND_PACKAGE(BISON REQUIRED)

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

    ${MPI_INCLUDE_PATH}

    ${GMP_INCLUDES}
    ${MPFR_INCLUDE_DIR}

    ${Bertini1_INCLUDES}
)



# Source files

INCLUDE(files.cmake)




# https://cmake.org/cmake/help/latest/module/FindFLEX.html


flex_target(MyScanner  src/io/partitionParse.l ${CMAKE_CURRENT_BINARY_DIR}/partitionParse.yy.c
            COMPILE_FLAGS "-s --prefix=partitionParse")

#bison_target(MyParser ${CMAKE_CURRENT_BINARY_DIR}/partitionParse.yy.c ${CMAKE_CURRENT_BINARY_DIR}/parser.cpp)
#add_flex_bison_dependency(MyScanner MyParser)





# Executable
add_executable(bertini_real 
                 ${common_headers} 
                 ${common_src}
                 ${br_headers}
                 ${br_src}
                 ${FLEX_MyScanner_OUTPUTS}
                 )

add_executable(sampler 
                 ${common_headers} 
                 ${common_src}
                 ${sampler_headers}
                 ${sampler_src}
                 ${FLEX_MyScanner_OUTPUTS}
                 )


target_link_libraries(bertini_real
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
    ${Bertini1_LIBRARIES}
)



target_link_libraries(sampler
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
    ${Bertini1_LIBRARIES}
)



# Installation
install(TARGETS bertini_real sampler DESTINATION bin)
