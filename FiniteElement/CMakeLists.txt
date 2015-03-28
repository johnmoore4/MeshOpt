ADD_SUBDIRECTORY(ShapeFunctions)
ADD_SUBDIRECTORY(InterpolationPoints)
ADD_SUBDIRECTORY(Polynomials)
ADD_SUBDIRECTORY(Quadrature)
ADD_SUBDIRECTORY(Math)

add_library(FiniteElement ${SOURCES} ${HEADERS})
target_include_directories(FiniteElement PUBLIC ${INCDIRS})

set(FiniteElement_INC ${INCDIRS} PARENT_SCOPE)


install(TARGETS FiniteElement DESTINATION lib)
#install(FILES ${HEADERS} DESTINATION include)

#move_src_up_level()
