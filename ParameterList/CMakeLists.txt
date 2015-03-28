ADD_SUBDIRECTORY(Utility)

add_library(Common ${SOURCES} ${HEADERS})
target_include_directories(Common PUBLIC ${INCDIRS})

set(Common_INC ${INCDIRS} PARENT_SCOPE)

install(TARGETS Common DESTINATION lib)
#install(FILES ${HEADERS} DESTINATION include)


#move_src_up_level()