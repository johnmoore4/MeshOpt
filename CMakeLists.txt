cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MeshOpt CXX)


set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;*.la;.so")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(Boost_USE_STATIC_LIBS   ON)
FIND_PACKAGE(Boost COMPONENTS program_options serialization mpi chrono iostreams REQUIRED)


SET(CMAKE_CXX_FLAGS "--std=c++11 -O3 ${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "")

macro(append_src FILES FILES2)
  foreach(FILE ${${FILES}})
    list(APPEND LIST ${CMAKE_CURRENT_SOURCE_DIR}/${FILE})
  endforeach(FILE)
  set(SOURCES ${SOURCES};${LIST} PARENT_SCOPE)
  foreach(FILE ${${FILES2}})
    list(APPEND LIST2 ${CMAKE_CURRENT_SOURCE_DIR}/${FILE})
  endforeach(FILE)
  set(HEADERS ${HEADERS};${LIST2} PARENT_SCOPE) 
  set(INCDIRS ${INCDIRS};${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
endmacro(append_src)



macro(move_src_up_level)
  set(SOURCES ${SOURCES} PARENT_SCOPE)
  set(HEADERS ${HEADERS} PARENT_SCOPE)
  set(INCDIRS ${INCDIRS} PARENT_SCOPE)
endmacro(move_src_up_level)

INCLUDE_DIRECTORIES(${OCC_ROOT}/inc)
INCLUDE_DIRECTORIES(${OCC_ROOT})
INCLUDE_DIRECTORIES(${LBFGS_ROOT}/include)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${ARMA_ROOT}/include)
INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES(${ARMA_ROOT}/lib)
LINK_DIRECTORIES(${OCC_ROOT}/lib)
LINK_DIRECTORIES(${LBFGS_ROOT}/lib)
LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

#SET(BLAS_LIBS mkl_intel_lp64 mkl_core mkl_sequential)
#SET(BLAS_LIB_DIRS /home/john/local/intel/mkl/lib/intel64 /home/john/local/intel/lib/intel64)

message("BLAS_LIBS: !" ${BLAS_LIBS})


LINK_DIRECTORIES(${BLAS_LIB_DIRS})

FILE(GLOB OCC_LIBS "${OCC_ROOT}/lib/*.a")

SET(MeshOpt_LIBS armadillo ${Boost_LIBRARIES} Common Mesh FiniteElement Optimization ${LBFGS_ROOT}/lib/liblbfgs.a ${BLAS_LIBS})

#ADD_SUBDIRECTORY(Common)
#ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../Solvers 
#  ${CMAKE_CURRENT_BINARY_DIR}/Solvers)

ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../Optimization 
  ${CMAKE_CURRENT_BINARY_DIR}/Optimization)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../FiniteElement 
  ${CMAKE_CURRENT_BINARY_DIR}/FiniteElement)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../Mesh
  ${CMAKE_CURRENT_BINARY_DIR}/Mesh)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../ParameterList
  ${CMAKE_CURRENT_BINARY_DIR}/ParameterList)
#ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/../Trilinos 
#  ${CMAKE_CURRENT_BINARY_DIR}/Trilinos)

ADD_SUBDIRECTORY(Common)
ADD_SUBDIRECTORY(BL)
ADD_SUBDIRECTORY(Optimization)
ADD_SUBDIRECTORY(HighOrder)


#ADD_SUBDIRECTORY(Debugging)

#move_src_up_level()

message("MeshOpt_LIBS: ${MeshOpt_LIBS}")
message("Mesh_INC: " ${Mesh_INC})
message("sources: " ${SOURCES})
message("Headers: " ${HEADERS})

add_executable(MeshOpt ${SOURCES} ${HEADERS})
#target_include_directories(MeshOpt PUBLIC ${INCDIRS} ${Mesh_INC} ${Solvers_INC} ${Optimization_INC} ${FiniteElement_INC})
target_include_directories(MeshOpt PUBLIC ${INCDIRS} ${Mesh_INC} ${Optimization_INC} ${FiniteElement_INC})

#target_link_libraries(MeshOpt ${MeshOpt_LIBS})

target_link_libraries(MeshOpt -Wl,--start-group  ${MeshOpt_LIBS} ${OCC_LIBS} ${OPENMPI_LIBS} -Wl,--end-group -Wl,--whole-archive /usr/lib/libibverbs.a -Wl,--no-whole-archive pthread m dl z)

install(TARGETS MeshOpt DESTINATION bin)


#install(FILES ${HEADERS} DESTINATION include)


