

SET(S
  #hmesh.cpp
  main.cpp
  HModel.cpp
  MeshOptConfig.cpp
  CaseParameterList.cpp
  BLParameterList.cpp
  HighOrderParameterList.cpp
)

SET(H
  #GlobalDefines.h
  #TrilinosTypedefs.h
  MeshOptConfig.h
  CaseParameterList.h
  BLParameterList.h
  HighOrderParameterList.h
  HModel.h
)

#append_src(${CMAKE_CURRENT_SOURCE_DIR} S H)
append_src(S H)