# Find FGSL package
# https://github.com/reinh-bader/fgsl
#
# This module will set the following variables in your project:
#   FGSL_FOUND
#   FGSL_ROOT

if(NOT FGSL_FOUND)
    set(FGSL_FOUND FALSE)

    foreach(DIR
            "$ENV{HOME}/.local/opt/fgsl"
            "$ENV{HOME}/.linuxbrew/opt/fgsl"
            "/home/linuxbrew/.linuxbrew/opt/fgsl")
        if((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
            message("-- Found FGSL: ${DIR}")
            set(FGSL_FOUND TRUE)
            set(FGSL_ROOT "${DIR}")
            break()
        endif((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
    endforeach(DIR)

    if(NOT FGSL_FOUND)
        message(
            FATAL_ERROR
            "Package FGSL not found."
        )
    endif(NOT FGSL_FOUND)
endif(NOT FGSL_FOUND)
