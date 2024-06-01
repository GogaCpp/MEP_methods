cmake_minimum_required(VERSION 3.22)
project(MEP_methods)

add_subdirectory(mep)

add_executable(example exempleCircle.cpp)

target_link_libraries(example
    PUBLIC
    mep)


