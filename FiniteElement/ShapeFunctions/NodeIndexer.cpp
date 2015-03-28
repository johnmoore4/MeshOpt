#include "NodeIndexer.h"
#include <iostream>
#include <cmath>
#include <map>

void PopulateGMSHTypes(int mytype,
		       std::vector<std::pair<int,int> >& GMSH2Type, 
		       NodeIndexer& node_indexer){

  const std::vector<int>& types = node_indexer.getGMSHTypes();
  for(int i=0; i<types.size(); i++){
    GMSH2Type[types[i]] = std::make_pair(mytype,i);
  }

}

int NodeIndexer::Initialize(int p){
  SetOrder(p);
  ComputeLocalFaceIndices();
}

int NodeIndexer::ComputeLocalFaceIndices(){
  std::map<indtype,indtype> local;
  local_face_indices.resize(NumChildren());
  rev_local_face_indices.clear();
  int cnt = 0;
  for(int ch = 0; ch < NumChildren(); ch++){
    const std::vector<indtype>& children = getChildNodes(ch);
    local_face_indices[ch].resize(children.size());
    for(int n = 0; n < children.size(); n++){

      auto it = local.find(children[n]);
      if(it == local.end()){
	it = local.insert(it,std::make_pair(children[n],cnt++));
	//it = local[children[n]] = cnt++;
	rev_local_face_indices.push_back(children[n]);
      }
      local_face_indices[ch][n] = it->second;
      
    }
  }
}
const std::vector<indtype> NodeIndexer::getExteriorNodes() const{
   const std::vector<indtype>& interior = getInteriorNodes();
   
    std::vector<bool> is_interior(Ndof(),false);
    for(int i=0; i<interior.size(); i++){
      is_interior[interior[i]] = true;
    }

    std::vector<indtype> exterior(Ndof()-interior.size());
    for(int i=0, cnt = 0; i < is_interior.size(); i++){
      if(!is_interior[i]) exterior[cnt++] = i;
    }

    //for(int i=0; i<exterior.size(); i++){
    //  std::cout << exterior[i] << " ";
    //}
    //std::cout << std::endl;
    return exterior;

}

NodeIndexFactory::NodeIndexFactory(){
  GMSH2Type.resize(120);


  {
    PointNodeIndexer Point(1);
    PopulateGMSHTypes(0,GMSH2Type,Point);
  }
  {
    LineNodeIndexer Line(1);
    PopulateGMSHTypes(1,GMSH2Type,Line);
  }
  {
    TriNodeIndexer Tri(1);
    PopulateGMSHTypes(2,GMSH2Type,Tri);
  }
  {
    QuadNodeIndexer Quad(1);
    PopulateGMSHTypes(3,GMSH2Type,Quad);
  }
  {
    TetNodeIndexer Tet(1);
    PopulateGMSHTypes(4,GMSH2Type,Tet);
  }
  {
    HexNodeIndexer Hex(1);
    PopulateGMSHTypes(5,GMSH2Type,Hex);
  }
  {
    PrismNodeIndexer Prism(1);
    PopulateGMSHTypes(6,GMSH2Type,Prism);
  }
  {
    PyramidNodeIndexer Pyramid(1);
    PopulateGMSHTypes(7,GMSH2Type,Pyramid);
  }
 
}

const int NodeIndexFactory::TypeFromGMSHType(const int gmtype){
  return GMSH2Type[gmtype].first;
}
const int NodeIndexFactory::OrderFromGMSHType(const int gmtype){
  return GMSH2Type[gmtype].second;
}

const std::shared_ptr<NodeIndexer> 
NodeIndexFactory::NodeIndexerFromType(int type, int order){
  switch(type){
  case 0:
    return std::make_shared<PointNodeIndexer>(order);
  case 1:
    return std::make_shared<LineNodeIndexer>(order);
  case 2:
    return std::make_shared<TriNodeIndexer>(order);
  case 3:
    return std::make_shared<QuadNodeIndexer>(order);
  case 4:
    return std::make_shared<TetNodeIndexer>(order);
  case 5:
    return std::make_shared<HexNodeIndexer>(order);
  case 6:
    return std::make_shared<PrismNodeIndexer>(order);
  case 8:
    return std::make_shared<PyramidNodeIndexer>(order);
  default:
    throw std::runtime_error("Invalid element type > 7!");
  }
}

const NodeIndexer* NodeIndexFactory::
getNodeIndexer(const int type, const int p){
  
  auto& ni = indexers[type][p];
  if(!ni){
    ni = NodeIndexerFromType(type,p);
    ni->Initialize(p);
  }
  return ni.get();
  
}

/*
const NodeIndexer* NodeIndexFactory::
getNodeIndexer(const int type, const int p){
  switch(type){
  case 0:
    {
      if(!Point){
	Point = std::unique_ptr<PointNodeIndexer>
	  (new PointNodeIndexer(p));
	//Point->Initialize(p);
      }
      return Point.get();
    }
  case 1:
    {
      if(!Line){
	Line = std::unique_ptr<LineNodeIndexer>(new LineNodeIndexer(p));
	//Line->Initialize(p);
      }
      if(Line->getP() != p) Line->SetOrder(p);
      return Line.get();
    }
  case 2:
    {
      if(!Tri){
	Tri = std::unique_ptr<TriNodeIndexer>(new TriNodeIndexer(p));
	//Tri->Initialize(p);
      }
      if(Tri->getP() != p) Tri->SetOrder(p);
      return Tri.get();
    }
  case 3:
    {
      if(!Quad){
	Quad = std::unique_ptr<QuadNodeIndexer>(new QuadNodeIndexer(p));
	//Quad->Initialize(p);
      }
      return Quad.get();
    }
  case 4:
    {
      if(!Tet){
	Tet = std::unique_ptr<TetNodeIndexer>(new TetNodeIndexer(p));
	//Tet->Initialize(p);
      }
      return Tet.get();
    }
  case 5:
    {
      if(!Hex){
	Hex = std::unique_ptr<HexNodeIndexer>(new HexNodeIndexer(p));
	//Hex->Initialize(p);
      }
      return Hex.get();
    }
  case 6:
    {
      if(!Prism){
	Prism =std::unique_ptr<PrismNodeIndexer>
	  (new PrismNodeIndexer(p));
	//Prism->Initialize(p);
      }
      return Prism.get();
    }
  case 7:
    {
      if(!Pyramid){
	Pyramid = std::unique_ptr<PyramidNodeIndexer>
	  (new PyramidNodeIndexer(p));
	//Pyramid->Initialize(p);
      }
      return Pyramid.get();
    }
  }
}
*/

template<class T>
void recursiveGMSHIndex(NodeIndexer* current, std::vector<indtype>& gmsh_index){

  int cnt = 0;

  using std::cout;
  using std::endl;

  int p = current->getP();
  gmsh_index.resize(current->Ndof(p));

  const std::vector<indtype>& cornern = current->getCornerNodes();
  const std::vector<indtype>& gmsh_cn = current->gmshCornerNodes();
  
  //for(int i=0; i<cornern.size(); i++)gmsh_index[cornern[i]] = cnt++;
  for(int i=0; i<cornern.size(); i++)gmsh_index[cnt++] = cornern[gmsh_cn[i]]; 

  for(int i=0; i<current->NumChildren(); i++){
    const std::vector<indtype>& edgen = current->getChildNodes(i);
    for(int j=1; j<edgen.size()-1; j++){
      //gmsh_index[edgen[j]] = cnt++;
      gmsh_index[cnt++] = edgen[j];
    }
  }
 

  const std::vector<indtype>& intnodes = current->getInteriorNodes();
 
 

  if(intnodes.size() == 1){
    //gmsh_index[intnodes[0]] = cnt++;
    gmsh_index[cnt++] = intnodes[0];
  }
  else if(intnodes.size() > 1){
    int ord;
 
    if(current->getMyType() == 2) ord = p-3;
    else if(current->getMyType() == 3) ord = p-2;
    else throw std::runtime_error("Only recursive Tri and Quad supported!");
    
    T newel(ord);
    
    const std::vector<indtype>& temp = newel.getGMSHIndex();
    for(int i=0; i<intnodes.size(); i++){
      gmsh_index[cnt++] = intnodes[temp[i]];
    }
  }

  
}


void PointNodeIndexer::SetOrder(const int p){

  order = p;

  gmsh_index = {0};
    
  gmsh_types = {15,15,15,15,15,15,15,15,15,15,15};
}

void LineNodeIndexer::SetOrder(const int p){
  

  // Set Object order
  order = p;

  // GMSH type
  gmsh_types = {1,1,8,26,27,28,62,63,64,65,66};

  // Corner nodes
  corner_nodes = {0,p};

  gmsh_corner_nodes = {0,1};

  // Interior Nodes
  const int ndof = Ndof(p);
  interior_nodes.resize(std::max(p-1,0));
  for(int i=0; i<p-1; i++){
    interior_nodes[i] = i+1;
  }

  // Child nodes
  child_nodes.resize(2);
  for(int i=0; i<2; i++) child_nodes[i].resize(1);
  child_nodes[0][0] = 0;
  child_nodes[1][0] = p;
 

  // Reversed nodes
  reversed_nodes.resize(ndof);
  for(int i=0; i<=p; i++) reversed_nodes[i] = p-i;


  // Child nodes
  

  // Oriented nodes
  oriented_nodes.resize(1);
  oriented_nodes[0].resize(p+1);
  for(int i=0; i<=p; i++) oriented_nodes[0][i] = i;

  // GMSH index
  gmsh_index.resize(p+1);
  gmsh_index[0] = 0;
  gmsh_index[1] = p;
  for(int i=1; i<p; i++){
    gmsh_index[1+i] = i;
  }


  // GMSh element type and index
  switch(p){
  case 1:
    gmsh_type = 1;
    break;
  case 2:
    gmsh_type = 8;
    break;
  case 3:
    gmsh_type = 26;
    break;
  case 4:
    gmsh_type =  27;
    break;
  case 5:
    gmsh_type = 28;
    break;
  default:
    gmsh_type =  56 + p;
  }

  
}

void TriNodeIndexer::SetOrder(const int p){
 

  // Set Object order
  order = p;

  my_type = 2;

  // GMSH type
  gmsh_types = {2,2,9,21,23,25,42,43,44,45};

  // Corner Nodes
  corner_nodes = {0,p,Ndof(p)-1};

  gmsh_corner_nodes = {0,1,2,3};

  // Interior Nodes
  const int ndof = Ndof(p);
  const int nint = Ndof(p-3);
  interior_nodes.resize(std::max(nint,0));

  
  for(int i=1, cnt=0; i<p; i++){
    for(int j=1; j<(p-i); j++,cnt++){
      interior_nodes[cnt] = ndof - Ndof(p-i)+j;
    }
  }
  

  // Reversed nodes
  reversed_nodes.resize(ndof);
  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<=(p-j); i++,cnt++){
      reversed_nodes[cnt] = ndof-Ndof(p-j) + p-j-i;
    }
  }
  

  
  // Child nodes
  child_nodes.resize(3);
  for(int i=0; i<3; i++) child_nodes[i].resize(p+1);

  for(int i=0; i<=p; i++){
    child_nodes[0][i] = i;
    child_nodes[1][i] = ndof - Ndof(p-i-1)-1;
    child_nodes[2][i] = ndof - Ndof(i);
  }

  // Oriented nodes
  oriented_nodes.resize(3);
  for(int i=0; i<3; i++) oriented_nodes[i].resize(ndof);
  for(auto i=0; i<ndof; i++) oriented_nodes[0][i] = i;
  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<=(p-j); i++,cnt++){
      oriented_nodes[1][cnt] = ndof - Ndof(p-i-1)-1 -j;
      oriented_nodes[2][cnt] = ndof - Ndof(i+j) + j;
    }
  }



  // GMSH element type and index
  switch(p){
  case 1:
    gmsh_type = 2;
    break;
  case 2:
    gmsh_type = 9;
    break;
  case 3:
    gmsh_type = 21;
    break;
  case 4:
    gmsh_type =  23;
    break;
  case 5:
    gmsh_type = 25;
    break;
  default:
    gmsh_type =  36 + p;
  }
  
  // GMSH index
 


  recursiveGMSHIndex<TriNodeIndexer>(this,gmsh_index);


  //gmsh_index = getCornerNodes();
  //gmsh_index = {0,1,2};

}

void QuadNodeIndexer::SetOrder(const int p){
  // Set Object order
  order = p;

  my_type = 3;

  // GMSH type
  gmsh_types = {3,3,10,36,37,38,47,48,49,50,51};

  // Corner nodes
  corner_nodes = {0, p,Ndof(p)-p-1, Ndof(p)-1};

  gmsh_corner_nodes = {0,1,3,2};
  
  // Interior Nodes
  interior_nodes.resize(std::max(Ndof(p-2),0));
  for(int i=1, cnt=0; i<p; i++){
    for(int j=1; j<p; j++,cnt++){
      interior_nodes[cnt] = i*(p+1)+j;
    }
  }

  // Reversed nodes
  reversed_nodes.resize(Ndof(p));
  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<=p; i++,cnt++){
      reversed_nodes[cnt] = j*(p+1) + (p-i);
    }
  }

  // Child nodes
  child_nodes.resize(4);
  for(int i=0; i<4; i++) child_nodes[i].resize(p+1);
  for(int i=0; i<=p; i++){
    child_nodes[0][i] = i;
    child_nodes[1][i] = (i+1)*(p+1)-1;
    child_nodes[2][i] = Ndof(p)-i-1;
    child_nodes[3][i] = (p-i)*(p+1);
  }

  // Oriented nodes
  oriented_nodes.resize(4);
  for(int i=0; i<4; i++) oriented_nodes[i].resize(Ndof(p));
  for(int j=0,cnt=0; j<=p; j++){
    for(int i=0; i<=p; i++,cnt++){
      oriented_nodes[0][cnt] = j*(p+1)+i;
      oriented_nodes[1][cnt] = i*(p+1)+(p-j);
      oriented_nodes[2][cnt] = (p-i) + (p-j)*(p+1);
      oriented_nodes[3][cnt] = (p-i)*(p+1)+j;
    }
  }

  // GMSH element type
  switch(p){
  case 1:
    gmsh_type = 3;
    break;
  case 2:
    gmsh_type = 10;
    break;
  case 3:
    gmsh_type = 36;
    break;
  case 4:
    gmsh_type = 37;
    break;
  case 5:
    gmsh_type = 38;
    break;
  default:
    gmsh_type = 41 + p;
  } 

  // GMSH index

  recursiveGMSHIndex<QuadNodeIndexer>(this,gmsh_index);

  //gmsh_index = getCornerNodes();
  //gmsh_index = {0,1,3,2};
}

void TetNodeIndexer::SetOrder(const int p){
  // Set Object order
  order = p;

  // GMSH type
  gmsh_types = {4,4,11,29,30,31,71,72,73,74,75};

  TriNodeIndexer tri(p);

  // Corner nodes
  corner_nodes = {0, p, tri.Ndof(p)-1, Ndof(p)-1};

  // Interior Nodes
  int nint = Ndof(p-4);
  interior_nodes.resize(std::max(nint,0));

  for(int k=1, cnt=0; k<p; k++){
    for(int j=1; j<(p-k); j++){
      for(int i=1; i<(p-j-k); i++,cnt++){
	interior_nodes[cnt] = 
	  i + tri.Ndof(p-k)-tri.Ndof(p-j-k) + Ndof(p)-Ndof(p-k);
      }
    }
  }

  // Reversed nodes

  // Child nodes

  child_nodes.resize(4);
  for(int i=0; i<4; i++) child_nodes[i].resize(tri.Ndof(p));
  const std::vector<indtype>& tri_reversed = tri.getReversedNodes();

  for(int j=0, cnt=0; j<=p; j++){
    for(int i=0; i<=(p-j); i++,cnt++){
      child_nodes[0][cnt] = tri.Ndof(p)-tri.Ndof(p-j) + i;
      child_nodes[1][cnt] = Ndof(p)-Ndof(p-j) + i;
      child_nodes[2][cnt] = Ndof(p)-Ndof(p-j) + tri.Ndof(p-j)-tri.Ndof(p-j-i);
      child_nodes[3][cnt] = Ndof(p)-Ndof(p-j) + tri.Ndof(p-j) -
	tri.Ndof(p-j-i-1) -1;
    }
  }
  std::vector<indtype> temp = child_nodes[0];
  for(int i=0; i<tri.Ndof(p); i++){
    child_nodes[0][i] = temp[tri_reversed[i]];
  }
  temp = child_nodes[2];
  for(int i=0; i<tri.Ndof(p); i++){
    child_nodes[2][i] = temp[tri_reversed[i]];
  }


  // Oriented nodes

  // GMSH element type
  switch(p){
  case 1:
    gmsh_type = 4;
    break;
  case 2:
    gmsh_type = 11;
    break;
  case 3:
    gmsh_type = 29;
    break;
  case 4:
    gmsh_type = 30;
    break;
  case 5:
    gmsh_type = 31;
    break;
  default:
    gmsh_type = 65 + p;
  }

  // GMSH index
  switch(p){

  case 1:
    gmsh_index = {0, 1, 2, 3};
    break;
  case 2:
    gmsh_index = {0, 2, 5, 9, 1, 4, 3, 6, 8, 7};
    break;
  case 3:
    gmsh_index = {0, 3, 9, 19, 1, 2, 6, 8, 7, 4, 16, 10, 18, 15, 17, 12, 5, 11, 13, 14};
    break;
  case 4:
    gmsh_index = {0, 4, 14, 34, 1, 2, 3, 8, 11, 13, 12, 9, 5, 31, 25, 15, 33, 30, 24, 32, 27, 18, 6, 10, 7, 16, 17, 26, 19, 28, 22, 29, 21, 23, 20};
    break;
  case 5:
    gmsh_index = {0, 5, 20, 55, 1, 2, 3, 4, 10, 14, 17, 19, 18, 15, 11, 6, 52, 46, 36, 21, 54, 51, 45, 35, 53, 48, 39, 25, 7, 16, 9, 12, 13, 8, 22, 24, 47, 23, 38, 37, 26, 49, 33, 40, 43, 30, 50, 29, 34, 42, 32, 44, 27, 28, 31, 41};
    break;
  case 6:
    gmsh_index = {0, 6, 27, 83, 1, 2, 3, 4, 5, 12, 17, 45, 24, 26, 25, 22, 18, 13, 7, 80, 74, 64, 49, 28, 82, 79, 73, 63, 48, 81, 76, 67, 53, 33, 8, 23, 11, 14, 19, 20, 16, 10, 9, 15, 29, 32, 75, 30, 31, 52, 66, 65, 50, 51, 34, 77, 46, 54, 68, 71, 61, 43, 39, 58, 78, 38, 47, 70, 57, 42, 62, 21, 72, 60, 35, 37, 44, 69, 36, 41, 40, 55, 59, 56};
    break;
  case 7:
    gmsh_index = {0, 7, 35, 119, 1, 2, 3, 4, 5, 6, 14, 20, 25, 29, 32, 34, 33, 30, 26, 21, 15, 8, 116, 110, 100, 85, 64, 36, 118, 115, 109, 99, 84, 63, 117, 112, 103, 89, 69, 42, 9, 31, 13, 16, 22, 27, 28, 24, 19, 12, 11, 10, 17, 23, 18, 37, 41, 111, 38, 39, 40, 68, 88, 102, 101, 86, 65, 66, 67, 87, 43, 113, 61, 70, 90, 104, 107, 97, 82, 58, 54, 49, 75, 94, 79, 114, 48, 62, 106, 93, 74, 53, 57, 60, 83, 98, 108, 96, 78, 81, 44, 47, 59, 105, 45, 46, 52, 56, 55, 50, 91, 71, 95, 80, 92, 73, 51, 72, 76, 77};
    break;
  case 8:
    gmsh_index = {0, 8, 44, 164, 1, 2, 3, 4, 5, 6, 7, 16, 23, 29, 34, 38, 41, 43, 42, 39, 35, 30, 24, 17, 9, 161, 155, 145, 130, 109, 81, 45, 163, 160, 154, 144, 129, 108, 80, 162, 157, 148, 134, 114, 87, 52, 10, 40, 15, 18, 25, 31, 36, 37, 33, 28, 22, 14, 13, 12, 11, 19, 32, 21, 26, 27, 20, 46, 51, 156, 47, 48, 49, 50, 86, 113, 133, 147, 146, 131, 110, 82, 83, 85, 132, 84, 112, 111, 53, 158, 78, 88, 115, 135, 149, 152, 142, 127, 106, 75, 71, 66, 60, 94, 139, 103, 120, 124, 99, 159, 59, 79, 151, 138, 119, 93, 65, 70, 74, 77, 107, 128, 143, 153, 141, 98, 105, 123, 102, 126, 54, 58, 76, 150, 55, 56, 57, 64, 69, 73, 72, 67, 61, 136, 116, 89, 140, 125, 104, 137, 118, 92, 62, 68, 63, 90, 91, 117, 95, 121, 100, 122, 97, 101, 96};
    break;
  case 9:
    gmsh_index = {0, 9, 54, 219, 1, 2, 3, 4, 5, 6, 7, 8, 18, 26, 33, 39, 44, 48, 51, 53, 52, 49, 45, 40, 34, 27, 19, 10, 216, 210, 200, 185, 164, 136, 100, 55, 218, 215, 209, 199, 184, 163, 135, 99, 217, 212, 203, 189, 169, 142, 107, 63, 11, 50, 17, 20, 28, 35, 41, 46, 47, 43, 38, 32, 25, 16, 15, 14, 13, 12, 21, 42, 24, 29, 36, 37, 31, 23, 22, 30, 56, 62, 211, 57, 58, 59, 60, 61, 106, 141, 168, 188, 202, 201, 186, 165, 137, 101, 102, 105, 187, 103, 104, 140, 167, 166, 138, 139, 64, 213, 97, 108, 143, 170, 190, 204, 207, 197, 182, 161, 133, 94, 90, 85, 79, 72, 115, 194, 130, 149, 175, 179, 158, 126, 121, 154, 214, 71, 98, 206, 193, 174, 148, 114, 78, 84, 89, 93, 96, 134, 162, 183, 198, 208, 196, 120, 132, 178, 153, 125, 129, 160, 181, 157, 65, 70, 95, 205, 66, 67, 68, 69, 77, 83, 88, 92, 91, 86, 80, 73, 191, 171, 144, 109, 195, 180, 159, 131, 192, 173, 147, 113, 74, 87, 76, 81, 82, 75, 110, 112, 172, 111, 146, 145, 116, 176, 127, 150, 155, 122, 177, 119, 128, 152, 124, 156, 117, 118, 123, 151};
    break;
  case 10:
    gmsh_index = {0, 10, 65, 285, 122, 2, 3, 4, 5, 6, 73, 74, 9, 8, 137, 144, 185, 50, 85, 224, 160, 64, 63, 60, 163, 199, 227, 248, 263, 260, 256, 282, 276, 266, 251, 230, 202, 166, 121, 66, 284, 281, 275, 265, 250, 229, 201, 165, 120, 167, 68, 69, 70, 71, 127, 128, 129, 75, 280, 61, 19, 103, 40, 30, 162, 52, 62, 58, 54, 49, 43, 36, 28, 29, 17, 16, 240, 142, 191, 262, 59, 37, 86, 159, 47, 55, 155, 44, 26, 35, 184, 33, 41, 149, 24, 20, 88, 123, 124, 125, 126, 72, 7, 173, 208, 235, 255, 269, 278, 283, 277, 267, 252, 231, 203, 168, 172, 268, 169, 170, 171, 207, 234, 254, 253, 232, 204, 205, 206, 233, 270, 279, 118, 154, 195, 32, 95, 147, 189, 21, 249, 192, 196, 111, 56, 119, 164, 45, 228, 39, 273, 146, 264, 51, 188, 221, 246, 225, 151, 156, 106, 38, 101, 152, 193, 100, 271, 84, 115, 67, 259, 80, 214, 180, 83, 92, 99, 141, 110, 114, 117, 116, 200, 197, 245, 94, 236, 243, 91, 113, 132, 89, 213, 14, 104, 109, 198, 241, 261, 223, 97, 22, 12, 18, 57, 257, 176, 79, 90, 81, 82, 27, 219, 42, 48, 53, 112, 46, 93, 102, 274, 237, 210, 175, 131, 77, 209, 31, 220, 157, 161, 1, 239, 15, 179, 136, 23, 108, 134, 139, 138, 222, 105, 98, 218, 96, 211, 135, 272, 133, 143, 186, 212, 258, 238, 177, 11, 174, 107, 76, 130, 226, 187, 145, 153, 194, 217, 178, 158, 78, 25, 34, 247, 215, 242, 190, 87, 150, 181, 183, 244, 148, 182, 140, 216, 13};
    break;

  }
}

void HexNodeIndexer::SetOrder(const int p){
  order = p;
  gmsh_types = {5,5,12,92,93};

  QuadNodeIndexer quad(p);

  const std::vector<indtype>& qc = quad.getCornerNodes();
  const int quad_dof = quad.Ndof(p);

  corner_nodes.resize(8);
  for(int i = 0; i < 4; i++){
    corner_nodes[i] = qc[i];
    corner_nodes[i+4] = p*quad.Ndof(p) + qc[i];
  }

  const int nint = std::max(std::pow(p-2,3),0.0);

  const std::vector<indtype>& qi = quad.getInteriorNodes();
  child_nodes.resize(8);
  for(int i = 0; i < 6; i++){
    child_nodes[i].resize(quad.Ndof(p));
  }
  for(int i = 0, cnt = 0; i < p+1; i++){
    for(int j = 0; j < p+1; j++, cnt++){
      child_nodes[0][cnt] = p-j + (p+1)*i; 
      child_nodes[1][cnt] = j + (p+1)*i + p*quad_dof;
      child_nodes[2][cnt] = j + i*quad_dof;
      child_nodes[3][cnt] = p + (p+1)*j + i*quad_dof;
      child_nodes[4][cnt] = p-j + p*(p+1) + i*quad_dof;
      child_nodes[5][cnt] = (p-j)*(p+1) + i*quad_dof;
    }
  }
  /*
  std::cout << "hex child nodes: " << std::endl;
  for(int i = 0; i < 6; i++){
    for(int j = 0; j < pow(p+1,2); j++){
      std::cout << child_nodes[i][j] << " ";
    }
    std::cout << std::endl;
  }
  */
  gmsh_type = gmsh_types[order];

  switch(p){
  case 1:
    gmsh_index = {0, 1, 3, 2, 4, 5, 7, 6};
    break;
  case 2:
    gmsh_index = {0,8,1,9,20,11,3,13,2,10,21,12,22,26,23,15,24,14,4,16,5,17,25,18,7,19,6};
    break;
  }
}

void PrismNodeIndexer::SetOrder(const int p){

  // Set Object order
  order = p;

  // GMSH type
  gmsh_types = {6,6,13,90,91,106,107,108,109,110,111};

  TriNodeIndexer tri(p);
  QuadNodeIndexer quad(p);

  // Corner nodes
  corner_nodes = {0, p, tri.Ndof(p)-1, p*tri.Ndof(p), p*tri.Ndof(p)+p,
		  tri.Ndof(p)*(p+1)-1};

 

  // Interior Nodes
  const int nint = (p-1)*tri.Ndof(p-3);
  interior_nodes.resize(std::max(nint,0));
  const std::vector<indtype> triint = tri.getInteriorNodes();
  for(int i=1, cnt=0; i<p; i++){
    for(int j=0; j<tri.Ndof(p-3); j++, cnt++){
      interior_nodes[cnt] = i*tri.Ndof(p) + triint[j];
    }
  }

 // Reversed nodes

  // Child nodes
  child_nodes.resize(6);
  for(int i=0; i<2; i++) child_nodes[i].resize(tri.Ndof(p));
  for(int i=2; i<6; i++) child_nodes[i].resize(quad.Ndof(p));

  const std::vector<indtype> trirev = tri.getReversedNodes();
  
  for(int i=0; i<tri.Ndof(p); i++){
    child_nodes[0][i] = trirev[i];
    child_nodes[1][i] = i + p*tri.Ndof(p);
  }
  for(int k=0; k<3; k++){
    const std::vector<indtype>& trichild = tri.getChildNodes(k);
    for(int i=0,cnt=0; i<=p; i++){
      for(int j=0; j<=p; j++,cnt++){
	child_nodes[k+2][cnt] = trichild[j] + i*tri.Ndof(p);
      }
    }
  }

  

  // GMSH element type
  switch(p){
  case 1:
    gmsh_type = 6;
    break;
  case 2:
    gmsh_type = 13;
    break;
  case 3:
    gmsh_type = 90;
    break;
  case 4:
    gmsh_type = 91;
    break;
  default:
    gmsh_type = 101 + p;
  }

  // GMSH index
  switch(p){
  case 1:
   gmsh_index = {0, 1, 2, 3, 4, 5};
   break;
  case 2:
    gmsh_index = {0, 2, 5, 12, 14, 17, 1, 3, 6, 4, 8, 11, 13, 15, 16, 7, 9, 10};
    break;
  case 3:
    gmsh_index = {0, 3, 9, 30, 33, 26, 1, 2, 4, 7, 10, 20, 6, 8, 13, 23, 19, 29, 31, 32, 34, 22, 36, 38, 5, 35, 11, 12, 37, 21, 14, 24, 27, 17, 16, 18, 28, 39, 15, 25};
    break;
  case 4:
    gmsh_index = {0, 4, 14, 60, 64, 74, 1, 2, 3, 5, 9, 12, 15, 30, 45, 8, 11, 13, 19, 34, 49, 29, 44, 59, 61, 62, 63, 65, 69, 72, 68, 71, 73, 6, 10, 7, 66, 67, 70, 16, 18, 48, 46, 17, 33, 47, 31, 32, 20, 50, 57, 27, 35, 54, 42, 24, 39, 23, 28, 58, 53, 26, 43, 56, 38, 41, 21, 51, 36, 22, 52, 37, 25, 55, 40};
    break;
  case 5:
    gmsh_index = {0, 5, 20, 105, 110, 125, 1, 2, 37, 4, 6, 11, 15, 18, 21, 42, 63, 84, 10, 14, 17, 19, 26, 47, 68, 89, 41, 62, 83, 104, 106, 107, 108, 109, 111, 116, 120, 123, 115, 119, 122, 98, 7, 16, 9, 12, 13, 8, 112, 114, 93, 113, 118, 117, 22, 25, 88, 85, 23, 24, 46, 67, 87, 86, 64, 43, 44, 45, 66, 65, 27, 90, 102, 39, 48, 69, 95, 99, 81, 60, 36, 32, 53, 74, 78, 57, 31, 40, 103, 94, 35, 38, 61, 82, 101, 124, 73, 52, 56, 59, 80, 77, 28, 91, 49, 70, 30, 121, 51, 72, 3, 100, 58, 79, 29, 92, 50, 71, 34, 97, 55, 76, 33, 96, 54, 75};
    break;
  case 6:
    gmsh_index = {0, 6, 27, 168, 174, 129, 1, 2, 47, 4, 5, 7, 13, 18, 22, 25, 28, 56, 84, 112, 140, 12, 17, 21, 24, 26, 34, 62, 90, 118, 146, 55, 83, 111, 139, 136, 169, 170, 171, 172, 173, 175, 181, 186, 116, 123, 180, 185, 189, 192, 194, 8, 23, 11, 14, 19, 20, 54, 51, 9, 15, 176, 179, 191, 177, 178, 184, 151, 187, 182, 144, 29, 33, 145, 141, 30, 31, 32, 138, 89, 117, 183, 143, 142, 113, 85, 57, 58, 60, 190, 153, 59, 88, 115, 86, 87, 35, 147, 132, 53, 63, 91, 119, 114, 121, 127, 137, 109, 81, 50, 46, 41, 69, 125, 134, 78, 97, 130, 106, 74, 102, 40, 16, 166, 152, 45, 49, 52, 82, 110, 61, 164, 161, 157, 124, 96, 68, 73, 80, 167, 195, 77, 108, 133, 101, 105, 36, 148, 64, 92, 120, 39, 188, 67, 95, 193, 10, 163, 79, 107, 135, 37, 149, 65, 93, 158, 38, 150, 66, 94, 122, 44, 156, 72, 100, 128, 48, 160, 76, 104, 165, 3, 159, 75, 103, 131, 42, 154, 70, 98, 126, 43, 155, 71, 99, 162};
    break;
  case 7:
    gmsh_index = {0, 110, 35, 115, 259, 205, 1, 2, 3, 4, 5, 6, 8, 15, 21, 26, 30, 33, 36, 72, 108, 144, 92, 216, 14, 20, 25, 29, 32, 34, 43, 79, 252, 151, 187, 223, 71, 107, 143, 91, 215, 114, 253, 254, 255, 256, 257, 258, 260, 267, 273, 278, 282, 285, 266, 272, 277, 281, 241, 286, 9, 31, 13, 16, 22, 27, 28, 70, 19, 12, 11, 10, 17, 23, 18, 261, 265, 283, 262, 263, 264, 271, 276, 280, 279, 274, 268, 269, 270, 275, 37, 42, 222, 217, 38, 39, 40, 41, 78, 251, 150, 186, 221, 220, 219, 218, 181, 50, 109, 73, 74, 175, 185, 182, 75, 174, 113, 149, 184, 232, 146, 7, 111, 112, 148, 147, 44, 224, 249, 69, 80, 116, 152, 188, 231, 237, 242, 246, 213, 84, 141, 105, 66, 62, 57, 51, 87, 195, 210, 102, 123, 159, 201, 206, 76, 138, 98, 93, 129, 165, 170, 134, 145, 24, 250, 230, 56, 61, 65, 68, 106, 142, 85, 214, 248, 245, 284, 236, 194, 158, 122, 86, 180, 104, 212, 200, 97, 101, 140, 176, 209, 287, 164, 128, 133, 137, 173, 169, 45, 225, 81, 117, 153, 189, 49, 229, 178, 121, 157, 193, 67, 247, 103, 139, 77, 211, 46, 226, 82, 118, 154, 190, 47, 227, 83, 119, 155, 238, 48, 228, 177, 120, 156, 192, 55, 235, 179, 127, 163, 199, 60, 240, 96, 132, 168, 204, 64, 244, 100, 136, 172, 208, 63, 198, 99, 135, 171, 207, 58, 191, 94, 130, 166, 202, 52, 183, 88, 124, 160, 196, 53, 233, 89, 125, 161, 197, 54, 234, 90, 126, 162, 243, 59, 239, 95, 131, 167, 203};
    break;
  case 8:
    gmsh_index = {0, 8, 44, 254, 368, 404, 1, 2, 3, 4, 5, 6, 7, 9, 17, 24, 30, 35, 39, 42, 45, 90, 135, 180, 225, 270, 142, 16, 23, 29, 34, 38, 41, 43, 53, 280, 143, 361, 233, 278, 323, 89, 134, 179, 51, 269, 266, 359, 188, 362, 363, 364, 365, 366, 367, 369, 377, 384, 390, 395, 399, 402, 376, 383, 389, 394, 344, 401, 299, 10, 40, 15, 18, 25, 31, 36, 37, 33, 28, 22, 14, 13, 12, 11, 19, 32, 21, 26, 27, 20, 370, 375, 292, 371, 372, 373, 374, 382, 388, 337, 397, 284, 275, 385, 378, 379, 381, 392, 320, 329, 386, 46, 52, 322, 316, 47, 48, 49, 50, 118, 97, 315, 358, 232, 277, 321, 380, 319, 318, 317, 158, 226, 60, 136, 91, 92, 96, 276, 272, 93, 215, 95, 312, 186, 231, 391, 274, 273, 106, 182, 137, 138, 140, 230, 228, 139, 185, 345, 183, 184, 54, 151, 357, 87, 99, 144, 189, 234, 279, 332, 339, 229, 350, 354, 262, 267, 111, 177, 132, 84, 80, 75, 69, 62, 107, 227, 309, 129, 152, 197, 242, 236, 244, 305, 264, 103, 174, 125, 120, 114, 159, 249, 260, 170, 204, 255, 94, 165, 210, 61, 88, 187, 331, 190, 74, 79, 83, 86, 133, 178, 223, 157, 313, 356, 353, 349, 398, 338, 286, 241, 196, 324, 287, 113, 131, 311, 293, 119, 181, 128, 176, 221, 150, 308, 304, 403, 248, 203, 271, 164, 173, 263, 360, 169, 218, 259, 209, 214, 55, 325, 100, 145, 68, 235, 98, 124, 330, 105, 314, 195, 240, 285, 85, 303, 130, 175, 220, 265, 310, 56, 326, 101, 146, 191, 294, 281, 57, 327, 102, 147, 192, 237, 282, 58, 328, 219, 148, 193, 238, 283, 59, 387, 104, 149, 194, 239, 396, 67, 393, 112, 268, 202, 247, 400, 73, 343, 224, 163, 208, 253, 298, 78, 348, 123, 168, 213, 258, 355, 82, 352, 127, 172, 217, 141, 307, 81, 297, 126, 171, 216, 261, 306, 76, 290, 121, 166, 211, 256, 301, 70, 340, 115, 160, 205, 250, 295, 63, 333, 108, 153, 198, 243, 288, 64, 334, 109, 154, 199, 300, 289, 66, 336, 222, 156, 201, 246, 291, 77, 347, 122, 167, 212, 257, 302, 65, 335, 110, 155, 200, 245, 346, 72, 342, 117, 162, 207, 252, 351, 71, 341, 116, 161, 206, 251, 296};
    break;
case 9:
  gmsh_index = {0, 249, 54, 495, 504, 496, 1, 2, 3, 84, 5, 6, 7, 8, 10, 19, 27, 34, 40, 45, 49, 52, 55, 110, 165, 220, 275, 189, 385, 440, 18, 26, 113, 39, 44, 48, 51, 53, 64, 119, 397, 442, 497, 339, 394, 449, 109, 164, 219, 221, 329, 384, 381, 374, 363, 284, 498, 499, 500, 501, 502, 503, 505, 229, 522, 529, 535, 336, 346, 355, 513, 521, 528, 534, 539, 543, 546, 548, 11, 50, 17, 20, 28, 35, 41, 46, 47, 43, 38, 32, 25, 16, 15, 91, 13, 12, 21, 42, 24, 29, 36, 37, 31, 97, 22, 30, 506, 512, 545, 507, 508, 509, 510, 511, 520, 527, 533, 538, 410, 541, 536, 530, 523, 515, 516, 519, 401, 517, 446, 526, 532, 391, 524, 525, 56, 63, 448, 228, 57, 58, 59, 60, 61, 62, 341, 253, 436, 491, 338, 393, 447, 518, 445, 444, 443, 514, 386, 331, 276, 274, 166, 111, 112, 340, 392, 459, 33, 260, 115, 116, 172, 435, 282, 337, 531, 390, 461, 388, 404, 204, 222, 167, 168, 317, 540, 405, 169, 316, 226, 281, 335, 474, 128, 223, 224, 225, 280, 279, 65, 237, 368, 107, 120, 175, 230, 285, 197, 395, 246, 174, 334, 344, 353, 361, 377, 382, 327, 141, 217, 162, 104, 100, 95, 89, 82, 74, 129, 332, 372, 159, 184, 158, 294, 277, 342, 351, 359, 366, 379, 324, 133, 214, 155, 150, 144, 137, 192, 287, 375, 210, 247, 302, 364, 370, 320, 124, 205, 199, 254, 309, 315, 114, 303, 108, 493, 458, 81, 88, 94, 99, 103, 106, 163, 218, 273, 117, 383, 378, 429, 488, 484, 479, 473, 466, 403, 348, 293, 451, 183, 278, 357, 161, 494, 411, 143, 149, 154, 239, 216, 271, 190, 245, 433, 283, 424, 418, 356, 301, 387, 191, 118, 213, 450, 549, 127, 209, 268, 323, 441, 369, 308, 173, 259, 264, 319, 314, 66, 238, 121, 176, 231, 286, 198, 396, 72, 457, 349, 182, 438, 292, 347, 402, 105, 490, 160, 215, 270, 325, 380, 373, 67, 452, 122, 177, 232, 136, 412, 467, 68, 453, 123, 178, 233, 288, 343, 398, 69, 454, 265, 179, 234, 289, 480, 469, 70, 455, 125, 321, 235, 290, 345, 400, 71, 456, 126, 322, 437, 291, 544, 537, 80, 465, 135, 326, 439, 300, 547, 542, 87, 472, 142, 328, 252, 307, 362, 417, 93, 478, 148, 203, 258, 313, 492, 487, 98, 483, 153, 208, 263, 318, 227, 428, 102, 423, 157, 212, 267, 181, 236, 432, 101, 486, 156, 211, 266, 180, 376, 431, 96, 481, 151, 206, 261, 170, 371, 426, 90, 475, 145, 200, 255, 310, 365, 420, 83, 468, 138, 193, 248, 73, 358, 413, 75, 460, 130, 185, 240, 295, 350, 333, 76, 389, 131, 186, 241, 296, 419, 406, 79, 464, 134, 330, 244, 299, 354, 409, 23, 416, 152, 207, 262, 171, 434, 427, 77, 462, 132, 187, 242, 297, 352, 407, 78, 463, 269, 188, 243, 298, 485, 476, 86, 471, 272, 196, 251, 306, 489, 482, 92, 477, 147, 202, 257, 312, 367, 422, 14, 408, 146, 201, 256, 311, 430, 421, 4, 399, 139, 194, 9, 304, 425, 414, 85, 470, 140, 195, 250, 305, 360, 415};
  break;
  }

}

void PyramidNodeIndexer::SetOrder(const int p){
  order = p;
  
  gmsh_types = {7,7,14,118,119};
  gmsh_type = gmsh_types[p];
  
  QuadNodeIndexer quad(p);
  const std::vector<indtype>& qc = quad.getCornerNodes();
  corner_nodes.resize(5);
  for(int i = 0 ; i < 4; i++) corner_nodes[i] = qc[i];
  corner_nodes[4] = Ndof(p)-1;

  // Interior nodes
  {
    int dofcnt = pow(p+1,2);
    for(int k = 1; k < p; k++){
      const int pk = p-k;
      for(int i = 1; i < pk; i++){
	for(int j = 0; j < pk; j++){
	  interior_nodes.push_back(i*pk+1+j+dofcnt);
	}
      } 
      dofcnt+= pow(pk+1,2);
    }
  }

  // Child nodes
  TriNodeIndexer tri(p);
  child_nodes.resize(5);
  child_nodes[0].resize(quad.Ndof(p));
  //for(int i = 1; i < 5; i++) child_nodes[i].resize(tri.Ndof(p));
  
  for(int i = 0, cnt = 0; i < p+1; i++){
    for(int j = 0; j < p+1; j++, cnt++){
      child_nodes[0][cnt] = p-j + (p+1)*i; 
    }
  }
  
  {
    int dofcnt = 0;
    for(int k = 0; k < p; k++){
      const int pk = p-k;
      QuadNodeIndexer quadk(pk);
      const int nc = quadk.NumChildren();
      
      for(int c = 0; c < nc; c++){
	const std::vector<indtype>& children = quadk.getChildNodes(c);
	for(int i = 0; i < children.size(); i++){
	  indtype temp = children[i] + dofcnt;
	  child_nodes[1+c].push_back(temp);
	}
      }
      dofcnt+= quadk.Ndof(pk);
    }
    for(int i = 0; i < 4; i++){
      child_nodes[1+i].push_back(dofcnt); // Include the tip node
    }
  }
  //child_nodes[0] = {0,1,2,3};

  

  switch(p){
  case 1:
    gmsh_index = {0,1,3,2,4};
    break;
  case 2:
    std::vector<int> temp = {0,5,1,6,13,8,3,10,2,7,9,12,11,4};
    gmsh_index.resize(temp.size());
    for(int i = 0; i < temp.size(); i++){
      gmsh_index[temp[i]] = i;
    }
      //gmsh_index = {0,5,1,6,13,8,3,10,2,7,9,12,11,4};
    break;
  }


}
