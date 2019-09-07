cmake_minimum_required(VERSION 3.14)
project(MRP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD")


include_directories("/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include")
include_directories("/opt/ibm/ILOG/CPLEX_Studio1261/concert/include")

link_directories("/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic")
link_directories("/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic")

add_executable(MRP main.cpp src/Data.cpp headers/Data.h headers/Include.h src/Model.cpp headers/Model.h src/Arc.cpp headers/Arc.h)

target_link_libraries(MRP concert ilocplex cplex m pthread)