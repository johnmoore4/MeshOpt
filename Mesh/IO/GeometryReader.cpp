#include "GeometryReader.h"
#include "GeometryContainer.h"


void GeometryReader::ReadGeometry(GeometryContainer& geom,std::string filename){

  STEPControl_Reader step_reader;
  step_reader.ReadFile(filename.c_str());
  step_reader.NbRootsForTransfer();
  step_reader.TransferRoots();
  geom.shape = step_reader.OneShape();
  BRepTools::Clean(geom.shape);
 
  TopTools_IndexedDataMapOfShapeListOfShape 
    EdgeFaceMap, VertexEdgeMap, VertexFaceMap;

  TopExp::MapShapesAndAncestors(geom.shape, TopAbs_VERTEX, TopAbs_EDGE, 
				VertexEdgeMap);
  TopExp::MapShapesAndAncestors(geom.shape, TopAbs_VERTEX, TopAbs_FACE, 
  				  VertexFaceMap);  
  TopExp::MapShapesAndAncestors(geom.shape, TopAbs_EDGE, TopAbs_FACE, 
				EdgeFaceMap);

  TopExp::MapShapes(geom.shape,TopAbs_FACE,geom.fmap);
  TopExp::MapShapes(geom.shape,TopAbs_EDGE,geom.emap);
  TopExp::MapShapes(geom.shape,TopAbs_VERTEX,geom.vmap);

  std::vector<myVertex>& myVertices = geom.getVerticesNC();
  std::vector<myEdge>& myEdges = geom.getEdgesNC();
  std::vector<myFace>& myFaces = geom.getFacesNC();
  for(int i=1; i<=geom.fmap.Extent(); i++){
    TopoDS_Face fc = TopoDS::Face(geom.fmap(i));
    myFaces.push_back(myFace(fc));

  }

  std::map<int,std::set<int> >& Edge2FacesMap = geom.getEdge2FacesMapNC();
  for(int i=1; i<=geom.emap.Extent(); i++){
    TopoDS_Edge eg = TopoDS::Edge(geom.emap.FindKey(i));
    myEdges.push_back(myEdge(eg));
    const TopTools_ListOfShape& listOfShapes = EdgeFaceMap.FindFromKey(eg);
    for(TopTools_ListIteratorOfListOfShape iterator(listOfShapes); 
	iterator.More(); iterator.Next()){
      TopoDS_Shape shape = iterator.Value();
      int index = geom.fmap.FindIndex(shape);
      Edge2FacesMap[i-1].insert(index-1);
    }
 
  }
  cout << "Edge2FaceMap size: " << Edge2FacesMap.size() << endl;
  std::map<int,std::set<int> >& Vertex2FacesMap = geom.getVertex2FacesMapNC();
  std::map<int,std::set<int> >& Vertex2EdgesMap = geom.getVertex2EdgesMap();

  for(int i=1; i<=geom.vmap.Extent(); i++){
    TopoDS_Vertex vx = TopoDS::Vertex(geom.vmap.FindKey(i));
    myVertices.push_back(myVertex(vx));
    // Build vertex->face map
    {
      const TopTools_ListOfShape& listOfShapes = VertexFaceMap.FindFromKey(vx);
      int cnt=0;
      for(TopTools_ListIteratorOfListOfShape iterator(listOfShapes); 
	  iterator.More(); iterator.Next()){
	TopoDS_Shape shape = iterator.Value();
	int index = geom.fmap.FindIndex(shape);
	Vertex2FacesMap[i-1].insert(index-1);
      }
    }
    // Build vertex->edge map
    {
      const TopTools_ListOfShape& listOfShapes = VertexEdgeMap.FindFromKey(vx);
      int cnt=0;
      for(TopTools_ListIteratorOfListOfShape iterator(listOfShapes); 
	  iterator.More(); iterator.Next()){
	TopoDS_Shape shape = iterator.Value();
	int index = geom.emap.FindIndex(shape);
	Vertex2EdgesMap[i-1].insert(index-1);
      }
    }
  }
   cout << "Vertex2FaceMap size: " << Vertex2FacesMap.size() << endl;


}
