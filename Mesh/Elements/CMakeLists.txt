set(S
  Mel.cpp
  El1D.cpp
  El2D.cpp
  El3D.cpp
  ActiveMEl.cpp
  CompEl.cpp
  OptEl.cpp
#  OptElParam.cpp
  OptElDev.cpp
  OptElManager.cpp
  elementFactory.cpp
  ElementAnalyzer.cpp
  ElementAnalyzerManager.cpp
  MEl_to_ActiveEl.cpp
)

set(H
  Mel.h
  El1D.h
  El2D.h
  El3D.h
  ActiveMEl.h
  CompEl.h
  OptEl.h
#  OptElParam.h
  OptElManager.h
  elementFactory.h
  ElementAnalyzer.h
  ElementAnalyzerManager.h
  CEl.h
  MEl_to_ActiveEl.h
  MeshTypedefs.h
)

append_src(S H)
