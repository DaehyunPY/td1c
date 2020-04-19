# Find package Intel MKL
# https://software.intel.com/mkl
#
# This module will set the following variables in your project:
#   MKL_FOUND
#   MKL_ROOT

if(NOT MKL_FOUND)
    set(MKL_FOUND FALSE)

    foreach(DIR
            "$ENV{MKLROOT}"
            "$ENV{HOME}/intel/mkl"
            "/opt/intel/mkl")
        if((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
            message("-- Found MKL: ${DIR}")
            set(MKL_FOUND TRUE)
            set(MKL_ROOT "${DIR}")
            break()
        endif((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
    endforeach(DIR)

    if(NOT MKL_FOUND)
        message(
            FATAL_ERROR
            "Package MKL not found. Please set specify library location \
            by setting system environment variable MKLROOT."
        )
    endif(NOT MKL_FOUND)
endif(NOT MKL_FOUND)
