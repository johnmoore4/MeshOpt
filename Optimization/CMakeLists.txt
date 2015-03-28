SET(S
  MeshOptimizer.cpp
#  GlobalMeshOptimizer.cpp
  SubmeshOptimizer.cpp
  LocalMeshMeritEvaluator.cpp
  SubmeshMeritEvaluator.cpp
)

SET(H
  MeshOptimizer.h
#  GlobalMeshOptimizer.h
  SubmeshOptimizer.h
  LocalMeshMeritEvaluator.h
  SubmeshMeritEvaluator.h
)

append_src(S H)
