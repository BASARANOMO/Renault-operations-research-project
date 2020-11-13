cmake_minimum_required(VERSION 3.3)

project(dispatch-solution-checker C)

option(CXX "enable C++ compilation" ON)
if(CXX)
    enable_language(CXX)
endif()

# set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE "Release")
# endif()
# message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

add_executable(dispatch Dispatch.cpp)
target_sources(dispatch PUBLIC
    Fournisseur.cpp
    Instance.cpp
    Route.cpp
    Usine.cpp
    Solution.cpp
)

