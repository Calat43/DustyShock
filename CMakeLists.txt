cmake_minimum_required(VERSION 3.7)
project(DustyShock)

set(CMAKE_C_STANDARD 11)
# show all warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -Wextra")

set(SOURCE_FILES
        main.c
        common_use.c
        common_use.h
        gas_shock.c
        gas_shock.h
        explicit_shock.c
        explicit_shock.h cells_shock.c cells_shock.h)
add_executable(DustyShock ${SOURCE_FILES})
target_link_libraries(DustyShock m)