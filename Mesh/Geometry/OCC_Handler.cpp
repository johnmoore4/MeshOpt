#include "OCC_Handler.h"
#include "myVertex.h"
#include "myEdge.h"
#include "myFace.h"


using namespace std;
void OCC_Handler::readSTEP(const std::string name){

  STEPControl_Reader step_reader;
  step_reader.ReadFile(name.c_str());
  cout << "after reading file" << endl;
  step_reader.NbRootsForTransfer();
  cout << "after Nbroots transfer" << endl;
  step_reader.TransferRoots();
  cout << "after transfer roots" << endl;
  _shape = step_reader.OneShape();
  cout << "after one shape" << endl;
  BRepTools::Clean(_shape);
  cout << "after clean" << endl;

  makeLists();

}

void OCC_Handler::makeLists(){
 
  

  //TopExp::MapShapesAndAncestors(_shape, TopAbs_VERTEX, TopAbs_EDGE, 
  //				  VertexEdgeMap);
  //TopExp::MapShapesAndAncestors(_shape, TopAbs_VERTEX, TopAbs_FACE, 
  //				  VertexFaceMap);  
  //TopExp::MapShapesAndAncestors(_shape, TopAbs_EDGE, TopAbs_FACE, EdgeFaceMap);
  /*
  cout << "number of faces: " << fmap.Extent() << endl;
 
  for(int i=1; i<=fmap.Extent(); i++){
    TopoDS_Face fc = TopoDS::Face(fmap(i));
    myFaces.push_back(face_ptr(new myFace(i-1,fc)));

  }
  cout << "after pushing back my faces" << endl;
  for(int i=1; i<=emap.Extent(); i++){
    TopoDS_Edge eg = TopoDS::Edge(emap.FindKey(i));
    vector<face_ptr> facelist;
    //myEdges.push_back( edge_ptr(new myEdge(i-1,eg,facelist)) );
    myEdges.push_back( edge_ptr(new myEdge(i-1,eg)) );
  }

  cout << "after pushing back my edges" << endl;

  for(int i=1; i<=vmap.Extent(); i++){
    TopoDS_Vertex vx = TopoDS::Vertex(vmap.FindKey(i));
    vector<edge_ptr> edgelist;
    vector<face_ptr> facelist;
    //myVertices.push_back(vertex_ptr(new myVertex(i-1,vx,edgelist,facelist)));
    myVertices.push_back(vertex_ptr(new myVertex(i-1,vx)));
  }
  cout << "after pushing back my vertices" << endl;
  */

  /*
  for(int i=1; i<=EdgeFaceMap.Extent(); i++){
    TopoDS_Edge eg = TopoDS::Edge(EdgeFaceMap.FindKey(i));
    const TopTools_ListOfShape& ListOfFaces = EdgeFaceMap.FindFromIndex(i);
    vector<face_ptr> facelist;
    TopTools_ListIteratorOfListOfShape it;
    for(it.Initialize(ListOfFaces); it.More(); it.Next()){
      int find = fmap.FindIndex(it.Value())-1;
      facelist.push_back(myFaces[find]);
    }
    myEdges.push_back( edge_ptr(new myEdge(i-1,eg,facelist)) );
  }
  */

  /*
  for(int i=1; i<=VertexEdgeMap.Extent(); i++){
    const TopTools_ListOfShape& ListOfEdges = VertexEdgeMap.FindFromIndex(i);
    const TopTools_ListOfShape& ListOfFaces = VertexFaceMap.FindFromIndex(i);
    TopoDS_Vertex vx = TopoDS::Vertex(VertexEdgeMap.FindKey(i));

    vector<face_ptr> facelist;
    TopTools_ListIteratorOfListOfShape it;
    cout << "B" << endl;
    for(it.Initialize(ListOfFaces); it.More(); it.Next()){
      int find = fmap.FindIndex(it.Value())-1;
      facelist.push_back(myFaces[find]);
    }

    vector<edge_ptr> edgelist;

 
    myVertices.push_back(vertex_ptr(new myVertex(i-1,vx,edgelist,facelist)));
    
  }
  

  */

  

  /*
  
  TopExp_Explorer exp0, exp1, exp2, exp3, exp4, exp5;
 

  // Free Shells
  for(exp1.Init(_shape, TopAbs_SHELL); exp1.More(); exp1.Next()){
    //cout << "free shell: " << endl;
    TopoDS_Shape shell = exp1.Current();

    TopExp::MapShapesAndAncestors(shell, TopAbs_VERTEX, TopAbs_EDGE, 
				  VertexEdgeMap);
    TopExp::MapShapesAndAncestors(shell, TopAbs_EDGE, TopAbs_FACE, EdgeFaceMap);

    if(shmap.FindIndex(shell) < 1){
      shmap.Add(shell);

      int fccnt=0;
      for(exp2.Init(shell, TopAbs_FACE); exp2.More(); exp2.Next()){
	fccnt++;
        TopoDS_Face face = TopoDS::Face(exp2.Current());
        if(fmap.FindIndex(face) < 1){
          fmap.Add(face);
	  //faces.push_back(face);
	  //surfaces.push_back(BRepAdaptor_Surface(faces.back()));
	  //myFaces.push_back(myFace(faces.back(),surfaces.back()));
	  //myFaces.push_back(myFace(face,BRepAdaptor_Surface(face)));

	  int wirecnt=0;
          for(exp3.Init(exp2.Current(), TopAbs_WIRE); exp3.More(); exp3.Next()){
	    // cout << "wire: " << wirecnt << endl;
	    wirecnt++;
	    TopoDS_Wire wire = TopoDS::Wire(exp3.Current());
            if(wmap.FindIndex(wire) < 1){
              wmap.Add(wire);

              for(exp4.Init(exp3.Current(), TopAbs_EDGE); exp4.More(); exp4.Next()){
                TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
                if(emap.FindIndex(edge) < 1){
                  emap.Add(edge);
		  edges.push_back(edge);
		  curves.push_back(BRepAdaptor_Curve(edges.back()));
		  myEdges.push_back(myEdge(edges.back(),curves.back()));

                  for(exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next()){
                    TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                    if(vmap.FindIndex(vertex) < 1){
                      vmap.Add(vertex);
		      vertices.push_back(vertex);
		    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  */

  //cout << "N edges: " << myEdges.size() << endl;
  //cout << "N faces: " << myFaces.size() << endl;
  

  /*

for(int i = 1; i <= EdgeFaceMap.Extent(); i++) {
  const TopoDS_Edge& E = TopoDS::Edge(EdgeFaceMap.FindKey(i));

  const TopTools_ListOfShape& ListOfFaces = EdgeFaceMap.FindFromIndex(i);

 

  TopTools_ListIteratorOfListOfShape it;

  for(it.Initialize(ListOfFaces); it.More(); it.Next()){
    const TopoDS_Face *test = &TopoDS::Face(it.Value());
    bool found = false;
    for(int k=0; k<myFaces.size(); k++){
      if(myFaces[k].getTopoDS_Face().IsEqual(it.Value())){
	myEdges[i-1].addFace(&myFaces[k]);
	found = true;
      }
    }
    if(found == false) cout << "did not find a face for edge " << i << endl;
  }


}
  */

  /*
  // Free Faces
  for(exp2.Init(_shape, TopAbs_FACE, TopAbs_SHELL); exp2.More(); exp2.Next()){
    cout << "free face: "  << endl;

    TopoDS_Face face = TopoDS::Face(exp2.Current());
    if(fmap.FindIndex(face) < 1){
      fmap.Add(face);

      for(exp3.Init(exp2.Current(), TopAbs_WIRE); exp3.More(); exp3.Next()){
        TopoDS_Wire wire = TopoDS::Wire(exp3.Current());
        if(wmap.FindIndex(wire) < 1){
          wmap.Add(wire);

          for(exp4.Init(exp3.Current(), TopAbs_EDGE); exp4.More(); exp4.Next()){
            TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
            if(emap.FindIndex(edge) < 1){
              emap.Add(edge);

              for(exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next()){
                TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
                if(vmap.FindIndex(vertex) < 1)
                  vmap.Add(vertex);
              }
            }
          }
        }
      }
    }
  }

  // Free Wires
  for(exp3.Init(_shape, TopAbs_WIRE, TopAbs_FACE); exp3.More(); exp3.Next()){
    cout << "in free wires" << endl;
    TopoDS_Wire wire = TopoDS::Wire(exp3.Current());
    if(wmap.FindIndex(wire) < 1){
      wmap.Add(wire);

      for(exp4.Init(exp3.Current(), TopAbs_EDGE); exp4.More(); exp4.Next()){
        TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
        if(emap.FindIndex(edge) < 1){
          emap.Add(edge);

          for(exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next()){
            TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
            if(vmap.FindIndex(vertex) < 1)
              vmap.Add(vertex);
          }
        }
      }
    }
  }

  // Free Edges
  for(exp4.Init(_shape, TopAbs_EDGE, TopAbs_WIRE); exp4.More(); exp4.Next()){
    
    cout << "in free edges" << endl;
    TopoDS_Edge edge = TopoDS::Edge(exp4.Current());
    if(emap.FindIndex(edge) < 1){
      emap.Add(edge);

      for(exp5.Init(exp4.Current(), TopAbs_VERTEX); exp5.More(); exp5.Next()){
        TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
        if(vmap.FindIndex(vertex) < 1)
          vmap.Add(vertex);
      }
    }
  }

  // Free Vertices
  for(exp5.Init(_shape, TopAbs_VERTEX, TopAbs_EDGE); exp5.More(); exp5.Next()){
    cout << "in free vertices" << endl;

    TopoDS_Vertex vertex = TopoDS::Vertex(exp5.Current());
    if(vmap.FindIndex(vertex) < 1)
      vmap.Add(vertex);
  }

  */

}




/*
arma::vec myEdge::edgeParamsFromPoint(const arma::vec &pt){

  gp_Pnt pnt(pt(0),pt(1),pt(2));
  GeomAPI_ProjectPointOnCurve proj(pnt, _curveHandle, umin, umax);
  if(!proj.NbPoints()) cout << "Curve projection failed!" << endl;
  arma::vec params(1);
  params(0) = proj.LowerDistanceParameter();
    
  return params;

}
*/

/*
arma::vec myFace::faceParamsFromPoint(const arma::vec &pt){

  gp_Pnt pnt(pt(0),pt(1),pt(2));

  GeomAPI_ProjectPointOnSurf proj(pnt, _surfaceHandle, umin, umax,vmin,vmax);
  if(!proj.NbPoints()) cout << "Curve projection failed!" << endl;
  arma::vec params(2);
  proj.LowerDistanceParameters(params(0),params(1));


  return params;
}
*/
