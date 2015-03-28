SET(S
  BoundaryLayer.cpp
  OptimizeBL.cpp
  OptimizeBLSubmesh.cpp
  NodeParamsOnGeometry.cpp
  BoundaryLayerSubdivider.cpp
  BLMesh.cpp
  DefaultNodeSpacingEvaluator.cpp
)

SET(H
  BoundaryLayer.h
  BoundaryLayerSubdivider.h
  BLMesh.h
  NodeSpacingEvaluator.h
  DefaultNodeSpacingEvaluator.h
)


append_src(S H)
