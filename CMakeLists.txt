project (AutoRender)

cmake_minimum_required (VERSION 2.8)

include (CMakeDependentOption)

# Use folders
set_property(GLOBAL PROPERTY USE_FOLDERS On)


# Macro to set output directories on all configurations
macro (set_output_directories OUTPUT_PATH)
    foreach (TYPE ${ARGN})
        set (CMAKE_${TYPE}_OUTPUT_DIRECTORY ${OUTPUT_PATH})
        foreach (CONFIG ${CMAKE_CONFIGURATION_TYPES})
            string (TOUPPER ${CONFIG} CONFIG)
            set (CMAKE_${TYPE}_OUTPUT_DIRECTORY_${CONFIG} ${OUTPUT_PATH})
        endforeach ()
    endforeach ()
endmacro ()

add_subdirectory (Source)

