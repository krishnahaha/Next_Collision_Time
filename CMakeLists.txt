cmake_minimum_required(VERSION 3.14)
project(EDMD)

SET(CMAKE_C_COMPILER g++)

set(CMAKE_C_FLAGS "-o3")

add_executable(edmdMain
        ./functions/hardSphere.cpp
        edmdMain.cpp)
