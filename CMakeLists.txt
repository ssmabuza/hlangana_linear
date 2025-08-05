project(hlangana)

cmake_minimum_required(VERSION 4)

set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(LINKER_LANGUAGE CXX)

set(Trilinos_DIR "$ENV{SOFTWARE}/trilinos/lib/cmake/Trilinos")
message(STATUS "Looking for Trilinos")
find_package(Trilinos CONFIG REQUIRED)
if( Trilinos_FOUND )
  message(STATUS "Found Trilinos!")
endif()

# Some serial impl in the hydrofem library
set(HYDROFEM_DIR "$ENV{PROJECTS}/research/cfd/")

# set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lmpi -L${PETSCFEM_LIBS} -lparmetis -lmetis -lboost_system -lboost_filesystem -lboost_program_options -lpetsc -Wall -Wextra")
set(CMAKE_CXX_FLAGS "-O3 -lmpi -L${PETSCFEM_LIBS} -lparmetis -lmetis -lboost_program_options")

set(CORE "./core")
set(DISC "./disc")
set(MESH "./mesh")
set(LINEAR "./linear")
set(PROJECT_HEADERS "./")

include_directories(${CORE} ${DISC} ${LINEAR} ${PROJECT_HEADERS} ${MESH} ${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

file(GLOB SOURCES "${CORE}/*.cpp"
                  "${DISC}/*.cpp"
                  "${MESH}/*.cpp"
                  "${LINEAR}/*.cpp"
                  "${PROJECT_HEADERS}*.cpp")

add_executable(hlangana ${SOURCES})

target_link_libraries(hlangana LINK_PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} parmetis metis boost_program_options)

install(TARGETS hlangana RUNTIME DESTINATION bin)
