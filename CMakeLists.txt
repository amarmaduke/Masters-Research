cmake_minimum_required(VERSION 2.6)

project(sunsol)

include_directories("${PROJECT_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

add_executable(main main.cpp json/json.cpp json/JSON_checker.cpp serial/force.cpp)

set(EXTRA_LIBS sundials_cvode sundials_nvecserial m)

target_link_libraries(main ${EXTRA_LIBS})
