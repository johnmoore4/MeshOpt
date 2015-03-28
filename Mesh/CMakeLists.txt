ADD_SUBDIRECTORY(Elements)
ADD_SUBDIRECTORY(Nodes)
ADD_SUBDIRECTORY(Src)
ADD_SUBDIRECTORY(Geometry)
ADD_SUBDIRECTORY(IO)

#move_src_up_level()


add_library(Mesh ${SOURCES} ${HEADERS})
target_include_directories(Mesh PUBLIC ${INCDIRS} ${FiniteElement_INC} ${Common_INC} ${TrilinosCompiles_INC})
set(Mesh_INC ${INCDIRS} PARENT_SCOPE)

install(TARGETS Mesh DESTINATION lib)
#install(FILES ${HEADERS} DESTINATION include)
