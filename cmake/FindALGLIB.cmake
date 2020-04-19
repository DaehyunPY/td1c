# Find ALGLIB package
# https://www.alglib.net
#
# This module will set the following variables in your project:
#   ALGLIB_FOUND
#   ALGLIB_ROOT

if(NOT ALGLIB_FOUND)
    set(ALGLIB_FOUND FALSE)

    foreach(DIR
            "$ENV{HOME}/.local/opt/alglib")
        if((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
            message("-- Found ALGLIB: ${DIR}")
            set(ALGLIB_FOUND TRUE)
            set(ALGLIB_ROOT "${DIR}")
            break()
        endif((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
    endforeach(DIR)

    if(NOT ALGLIB_FOUND)
        message(
            FATAL_ERROR
            "Package ALGLIB not found."
        )
    endif(NOT ALGLIB_FOUND)
endif(NOT ALGLIB_FOUND)
