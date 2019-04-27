cmake_minimum_required(VERSION 3.1)

set (TARGET_NAME AutoRender)

file(GLOB_RECURSE _HEADERS  *.h *.hpp)

file(GLOB_RECURSE _SCRS  *.c *.cpp *.cxx *.cc)

add_executable (${TARGET_NAME} ${_SCRS} ${_HEADERS})

set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE cxx)