# Find Libxc package
# https://www.tddft.org/programs/libxc
#
# This module will set the following variables in your project:
#   Libxc_FOUND
#   Libxc_ROOT

if(NOT Libxc_FOUND)
    set(Libxc_FOUND FALSE)

    foreach(DIR
            "$ENV{HOME}/.local/opt/libxc@4.3"
            "$ENV{HOME}/.local/opt/libxc@4"
            "$ENV{HOME}/.local/opt/libxc"
            "$ENV{HOME}/.linuxbrew/.local/opt/libxc@4.3"
            "$ENV{HOME}/.linuxbrew/.local/opt/libxc@4"
            "$ENV{HOME}/.linuxbrew/.local/opt/libxc"
            "/home/linuxbrew/.linuxbrew/.local/opt/libxc@4.3"
            "/home/linuxbrew/.linuxbrew/.local/opt/libxc@4"
            "/home/linuxbrew/.linuxbrew/.local/opt/libxc")
        if((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
            message("-- Found Libxc: ${DIR}")
            set(Libxc_FOUND TRUE)
            set(Libxc_ROOT "${DIR}")
            break()
        endif((NOT "${DIR}" STREQUAL "") AND (EXISTS "${DIR}"))
    endforeach(DIR)

    if(NOT Libxc_FOUND)
        message(
            FATAL_ERROR
            "Package Libxc not found."
        )
    endif(NOT Libxc_FOUND)
endif(NOT Libxc_FOUND)
