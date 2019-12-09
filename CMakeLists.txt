cmake_minimum_required(VERSION 3.15)
project(MaxService)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/gurobi752/linux64/include)

link_directories(/opt/gurobi752/linux64/lib)

set(SOURCE_FILES main.cpp src/Graph.cpp headers/Graph.h headers/Include.h src/Model.cpp headers/Model.h src/Arc.cpp headers/Arc.h)

add_executable(MaxService ${SOURCE_FILES})

target_link_libraries(MaxService gurobi_c++ gurobi90 pthread m)
