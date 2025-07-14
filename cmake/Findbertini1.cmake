if(WIN32)
    find_path(Bertini1_INCLUDE_DIR
            NAMES bertini.h
            PATHS "$ENV{BERTINI1_DIR}/Include")
    find_library(Bertini1_LIBRARY
            NAMES bertini-parallel
            PATHS "$ENV{BERTINI1_DIR}/Lib")
elseif(UNIX)
    find_path(Bertini1_INCLUDE_DIR
            NAMES bertini.h
            PATHS "$ENV{BERTINI1_DIR}/include" ${INCLUDE_INSTALL_DIR}
            REQUIRED)
    find_library(Bertini1_LIBRARY
            NAMES bertini-parallel
            PATHS "$ENV{BERTINI1_DIR}/lib" ${LIB_INSTALL_DIR}
            REQUIRED)
endif()
set(bertini1_INCLUDE_DIR ${Bertini1_INCLUDE_DIR})
set(Bertini1_LIBRARIES ${Bertini1_LIBRARY})
set(Bertini1_CPP_FLAGS "-D_HAVE_MPI")

