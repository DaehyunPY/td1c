# Find package SHTOOLS
# https://shtools.oca.eu/shtools/public
#
# This module will set the following variables in your project:
#   SHTOOLS_FOUND
#   SHTOOLS_ROOT

if(NOT SHTOOLS_FOUND)
    set(SHTOOLS_FOUND FALSE)

    foreach(DIR
            "$ENV{HOME}/.local/opt/shtools"
            "$ENV{HOME}/.linuxbrew/opt/shtools"
            "/home/linuxbrew/.linuxbrew/opt/shtools")
        if((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
            message("-- Found SHTOOLS: ${DIR}")
            set(SHTOOLS_FOUND TRUE)
            set(SHTOOLS_ROOT "${DIR}")
            break()
        endif((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
    endforeach(DIR)

    if(NOT SHTOOLS_FOUND)
        message(
            FATAL_ERROR
            "Package SHTOOLS not found."
        )
    endif(NOT SHTOOLS_FOUND)
endif(NOT SHTOOLS_FOUND)
