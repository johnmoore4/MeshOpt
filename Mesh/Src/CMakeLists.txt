SET(S
  #generateHOMesh.cpp
  getMeshConnectivity.cpp
  MeshContainer.cpp
  ChildGenerator.cpp
#  MeshDistributor.cpp
#  MeshDistributorHDG.cpp
#  MeshDistributorDG.cpp
  HDGMesh.cpp
  MeshAnalyzer.cpp
  ElementShapeFunctionManager.cpp
#  HighOrderMeshGenerator.cpp
  ParametricCoordinateEvaluator.cpp
  NodeReParameterizer.cpp
)

SET(H
  MeshConnectivity.h
  MeshContainer.h
  ChildGenerator.h
#  MeshDistributor.h
  LocalMesh.h
  DGMesh.h
#  MeshDistributorHDG.h
#  MeshDistributorDG.h
  MeshAnalyzer.h
  ElementShapeFunctionManager.h
#  HighOrderMeshGenerator.h
  ParametricCoordinateEvaluator.h
  NodeReParameterizer.h
)

append_src(S H)
