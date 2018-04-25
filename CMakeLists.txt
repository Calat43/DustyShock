cmake_minimum_required(VERSION 3.7)
project(DustyShock)

set(CMAKE_CXX_STANDARD 11)
# show all warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Wextra")

set(SOURCE_FILES
        main.cpp
        common_use.cpp
        common_use.hpp
        gas_shock.cpp
        gas_shock.hpp
        explicit_shock.cpp
        explicit_shock.hpp
        cells_shock.cpp
        cells_shock.hpp)
add_executable(DustyShock ${SOURCE_FILES})
target_link_libraries(DustyShock m)