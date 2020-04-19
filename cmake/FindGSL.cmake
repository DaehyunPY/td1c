# Find GSL package
# https://www.gnu.org/software/gsl
#
# This module will set the following variables in your project:
#   GSL_FOUND
#   GSL_ROOT

if(NOT GSL_FOUND)
    set(GSL_FOUND FALSE)

    foreach(DIR
            "$ENV{HOME}/.local/opt/gsl"
            "$ENV{HOME}/.linuxbrew/opt/gsl"
            "/home/linuxbrew/.linuxbrew/opt/gsl")
        if((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
            message("-- Found GSL: ${DIR}")
            set(GSL_FOUND TRUE)
            set(GSL_ROOT "${DIR}")
            break()
        endif((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
    endforeach(DIR)

    if(NOT GSL_FOUND)
        message(
            FATAL_ERROR
            "Package GSL not found."
        )
    endif(NOT GSL_FOUND)
endif(NOT GSL_FOUND)
