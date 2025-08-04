project(scalar)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(LINKER_LANGUAGE CXX)

#set(Trilinos_DIR "/home/sibusiso/Workspace/project_trilinos/lib/cmake/Trilinos")
set(Trilinos_DIR "/home/smabuza/Workspace/project_trilinos/lib/cmake/Trilinos")
message(STATUS "Looking for Trilinos")
find_package(Trilinos CONFIG REQUIRED)
if( Trilinos_FOUND )
  message(STATUS "Found Trilinos!")
endif()

set(PETSCFEM_LIBS "/home/smabuza/Workspace/usr/lib")
# set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lmpi -L${PETSCFEM_LIBS} -lparmetis -lmetis -lboost_system -lboost_filesystem -lboost_program_options -lpetsc -Wall -Wextra")
set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -lmpi -L${PETSCFEM_LIBS} -lparmetis -lmetis -lboost_system -lboost_filesystem -lboost_program_options")


set(LINEAR_SOLVERS "../IML/")
set(DEPS_HEADERS "../")
set(PROJECT_HEADERS "./")
set(PETSCFEM_HEADERS "/home/smabuza/Workspace/usr/include")
set(TVMET_HEADERS "${PETSCFEM_HEADERS}/tvmet")
set(MPI_STUFF "/usr/include/openmpi")

include_directories(${MPI_STUFF} ${DEPS_HEADERS} ${LINEAR_SOLVERS} ${PROJECT_HEADERS} ${PETSCFEM_HEADERS} ${TVMET_HEADERS} ${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

file(GLOB SOURCES "${DEPS_HEADERS}*.cpp" "${DEPS_HEADERS}*.hpp" "${PROJECT_HEADERS}*.cpp" "${PROJECT_HEADERS}*.hpp")

add_executable(scalar ${SOURCES})

target_link_libraries(scalar LINK_PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} parmetis metis boost_system boost_filesystem boost_program_options)

install(TARGETS scalar RUNTIME DESTINATION bin)
