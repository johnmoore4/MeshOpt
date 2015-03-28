//#define ARMA_NO_DEBUG
#include "Mel.h"
#include "Mnode.h"
#include "El1D.h"
#include "El2D.h"
#include "El3D.h"
#include "NodeIndexer.h"

//#include "shapeFunctionHandler.h"
#include <limits>
#include <functional>
#include <algorithm>  
#include <vector>
#include "assert.h"
#include "myEdge.h"
#include "myFace.h"


MEl::MEl(std::vector<gind>& nodes, unsigned char order, unsigned char bc_tag_t,
    short int geo_entity_t,bool geo_type_t):
  porder(order), bc_tag(bc_tag_t), geo_entity(geo_entity_t), geo_type(geo_type_t){
  _nodes = std::unique_ptr<gind[]>(new gind[nodes.size()]);
  for(int i=0; i<nodes.size(); i++){
    _nodes[i] = nodes[i];
  }
}
  


signed char ElementOrientation(const unique_element_ptr& el1,
			       const unique_element_ptr& el2,
			       NodeIndexFactory& indexfactory){
  const NodeIndexer* indexer = 
    indexfactory.getNodeIndexer(el1->getElementType(),el1->getOrder());

  const std::vector<indtype>& cnodes = indexer->getCornerNodes();
  const std::vector<indtype>& gmsh_ind = indexer->gmshCornerNodes();
  const std::vector<indtype>& reversed = indexer->getReversedNodes();
  
  const int ncn = cnodes.size();
  const int Nn = el1->NumNodes();
  
  const gind* nodes1 = el1->getNodes();
  const gind* nodes2 = el2->getNodes();
  const gind fn1 = nodes1[cnodes[gmsh_ind[0]]];
  const gind fn2 = nodes1[cnodes[gmsh_ind[1]]];
  /*
  std::cout << "nodes: " << std::endl;
  for(int i=0; i<Nn; i++){
    std::cout << nodes1[i] << " " << nodes2[i] << std::endl;
  }
  */
  for(int i=0; i<indexer->NumChildren(); i++){
    const std::vector<indtype>& or_pos = indexer->getOrientedNodes(i+1);
    const std::vector<indtype>& or_neg = indexer->getOrientedNodes(-(i+1));
    int cnt_pos=0, cnt_neg = 0;
    //std::cout << "orientations: " << std::endl;
    for(int n=0; n<Nn; n++){
      //std::cout << or_pos[n] << " " << or_neg[n] << std::endl;
      //std::cout << nodes1[n] << " " << nodes2[or_neg[n]] << std::endl;
      if(nodes1[n] == nodes2[or_pos[n]]) cnt_pos++;
      if(nodes1[n] == nodes2[or_neg[n]]) cnt_neg++;
    }
    //std::cout << "cnt_pos: " << cnt_pos << " cnt_neg: " << cnt_neg << std::endl;
    if(cnt_pos == Nn) return i+1;
    else if(cnt_neg == Nn) return -(i+1);

  }
  throw std::runtime_error("Orientation could not be found!");
  /*
  unsigned char orient;
  for(int i=0; i<indexer->NumChildren(); i++){
    int idx = i+1;
    //if(i == indexer->NumChildren()-1) idx = 0;
    const gind na = nodes2[cnodes[gmsh_ind[i]]];
    const gind nb = nodes2[cnodes[gmsh_ind[idx]]];
    //std::cout << na << " " << nb << std::endl;
    if(na == fn1 && nb == fn2){
      orient = (i+1);
      //std::cout << "in pos orient" << std::endl;
    }
    else if(na == fn2 && nb == fn1){
      orient =  -(i+1);
      // std::cout << "in neg orient" << std::endl;
    }
  }
  //std::cout << "orient: " << int(orient) << std::endl;
  return orient;
  */
  // cout << "Error! Elements are not orientable!" << endl;
  

  /*
  const gind* nodes1 = el1->getNodes();
  const gind* nodes2 = el2->getNodes();
  const gind tofind[] = {nodes1[cnodes0[0]],nodes1[cnodes1[0]]};

  
  for(int i=0; i<indexer->NumChildren(); i++){
    const std::vector<indtype>& cn0 = indexer->getChildNodes(i);
    int idx = i+1;
    if(i == indexer->NumChildren()-1) idx = 0;
    const std::vector<indtype>& cn1 = indexer->getChildNodes(idx);
    const gind na = nodes2[cn0[0]];
    const gind nb = nodes2[cn1[0]];
    if(na == tofind[0] && nb == tofind[1]) return (i+1);
    else if(na == tofind[1] && nb == tofind[0]) return -(i+1);
  }
  cout << "Error! Elements are not orientable!" << endl;
  */
  /*
  const int ncn1 = el1->numCornerNodes();
  const int ncn2 = el2->numCornerNodes();
  int cn1[ncn1];
  int cn2[ncn2];
  
  
  
  //const int* cn1 = el1->getCornerNodes(indexer);
  //const int* cn2 = el2->getCornerNodes();

  const int tofind[] = {cn1[0],cn1[1]};

  
  for(int i=0; i<indexer->NumBoundary(); i++){
    const int* bns = indexer->ChildCornerNodes(i);
    if(cn2[bns[0]] == tofind[0] && cn2[bns[1]] == tofind[1]) return i+1;
    else if(cn2[bns[0]] == tofind[1] && cn2[bns[1]] == tofind[0]) return -(i+1);
  }
  cout << "Error! Elements are not orientable!" << endl;
  */

}



/*
template<class T1, class T2>
class lw_vector{
private:

public:
  T2 N;
  T1* content;
  lw_vector():N(0){}
  lw_vector(T2 Nt): N(Nt){
    content = new T1[N];
  }
  ~lw_vector(){
    delete[] content;
  }
  lw_vector& operator=(lw_vector temp){
    std::swap(N,temp.N);
    std::swap(content,temp.content);
    return *this;
  }
  lw_vector(const lw_vector& other):N(other.N), content(new T1[other.N]){
    std::copy(other.content,other.content+N,content);
  }
  inline T1& operator[](T2 ind) const { return content[ind]; }
};
*/





/*
std::unique_ptr<gind[]> MEl::computeSortedNodes()const{

  const int N = numCornerNodes();
  std::vector<gind> ndind(N);
  const gind* cornern = getCornerNodes();
  for(auto i=0; i<N; i++){
    ndind[i] = cornern[i];
  }
  std::sort(ndind.begin(),ndind.end());

  std::unique_ptr<gind[]> sorted = std::unique_ptr<gind[]>(new gind[N]);
  for(auto i=0; i<N; i++){
    sorted[i] = ndind[i];
    //cout << sorted[i] << " ";
  }
  //cout << endl;
  return sorted;
  //for(auto i=0; i<N; i++) sorted_nodes[i] = ndind[i].second;
}
*/


bool ElementLess::operator()(const unique_element_ptr& el1, 
			      const unique_element_ptr& el2) const {
  const int nn1 = el1->numCornerNodes();
  const int nn2 = el2->numCornerNodes();
  if(nn2 > nn1) return true;
  else if(nn2 < nn1) return false;
  
  gind sn1[8], sn2[8];
  el1->ComputeSortedCornerNodes(sn1);
  el2->ComputeSortedCornerNodes(sn2);
  
  for(int i = 0; i < nn1; i++){
    if(sn1[i] < sn2[i]) return true;
    else if(sn1[i] > sn2[i]) return false;
  }
  return false;
}



bool ElementEqual::operator()(const unique_element_ptr& el1, 
			      const unique_element_ptr& el2) const {


  if(el1->numCornerNodes() != el2->numCornerNodes()){
    return false;
  }
  if(el1->getElementType() != el2->getElementType()){
    return false;
  }
  const int nn = el1->numCornerNodes();
  
  gind sn1[8], sn2[8];
  el1->ComputeSortedCornerNodes(sn1);
  el2->ComputeSortedCornerNodes(sn2);
  
  for(int i=0; i<nn; i++) if(sn1[i] != sn2[i]) return false;

  /*
  const unsigned char* sn1 = el1->getSortedCornerNodes();
  const gind* cn1 = el1->getNodes();

  const unsigned char* sn2 = el2->getSortedCornerNodes();
  const gind* cn2 = el2->getNodes();

  for(int i=0; i<nn; i++)if(cn2[sn2[i]] != cn1[sn1[i]])return false; 
  */

  return true;
  
}




std::size_t ElementHasher::operator()(const unique_element_ptr& el) const{
  //using boost::hash_value;
  using boost::hash_combine;
  std::hash<gind> hash_fn;
  std::size_t seed = 0.0;

  const int ncn = el->numCornerNodes();
  gind snodes[8];
  el->ComputeSortedCornerNodes(snodes);

  for(int i=0; i<ncn; i++){;
    hash_combine(seed,hash_fn(snodes[i]));
  }

  /*
  const gind* nodes = el->getNodes();

  const unsigned char* sn = el->getSortedCornerNodes();

  for(int i=0; i<ncn; i++){;
    hash_combine(seed,hash_fn(nodes[sn[i]]));
  }
  */

  return seed;
  
}

/*
inline short int const KroneckerDelta(const short int i, const short int j){
  return (i==j) ? 1 : 0; 
}
*/

//inline double const Dirac(const int i, const int j){ return (i==j) ? 1. : 0.; }
template <typename T>
double Sign(T a) {return (a >= 0.0) ? 1.0 : -1.0;}



using namespace std;
using namespace arma;



typedef std::pair<int,unsigned char> pair_type;
bool pair_greater(const pair_type a, const pair_type b){
  if(a.first > b.first) return true;
  else return false;
}

void MEl::ComputeSortedCornerNodes(gind* sorted){
  const int N = numCornerNodes();

  int cni[8];
  getCornerNodes(cni);
  for(int i=0; i<N; i++){
    sorted[i] = _nodes[cni[i]];
  }
 
  std::sort(sorted,sorted+N);
  for(int i=0; i<N; i++){
    //std::cout << sorted[i] << " ";
  }
  //std::cout << std::endl;
}

void MEl::sortNodes(unsigned char* sorted_nodes, const NodeIndexer* indexer){
  const int N = numCornerNodes();
  std::vector<std::pair<int,unsigned char> > ndind(N);
  //const gind* cornern = getCornerNodes();
  const gind* nodes = getNodes();
  const std::vector<indtype>& cnind = indexer->getCornerNodes();
  
  for(auto i=0; i<N; i++){
    //std::cout << cnind[i] << " ";
    ndind[i].first = nodes[cnind[i]];
    ndind[i].second = cnind[i];
    //ndind[i].first = cornern[i];
    //ndind[i].second = i;
  }
  //std::cout << std::endl;

  std::sort(ndind.begin(),ndind.end(),pair_greater);
  for(auto i=0; i<N; i++) sorted_nodes[i] = ndind[i].second;

  /*
  const arma::uvec& cornern = _sf->getCornerNodes();

  //const unsigned short int nn = _nodes.size();
  const unsigned short int nn = cornern.size();
  
  std::vector<pair_type> ndind(nn);
  
  for(auto i=0; i<nn; i++){
    ndind[i].first = _nodes[cornern[i]]->nd();
    ndind[i].second = i;
  }
  
  std::sort(ndind.begin(),ndind.end(),pair_greater);
  for(int i=0; i<nn; i++) _sorted_nodes[i] = ndind[i].second;
  */
}

/*
MEl::MEl(int el, unsigned char bc_tag, std::vector<node_ptr>& nodes, ShapeFunction* sf): 
  _el(el), _bc_tag(bc_tag), _nodes(nodes), _sf(sf), _curved(false){
  sortNodes();
  //_sorted_corners = new unsigned char[10];

}
*/


//arma::mat MEl::nodes() const{
  /*
  const int nn = _nodes.size();
  arma::mat temp(3,nn);
  //for(int i=0; i<nn; i++) {temp.unsafe_col(i) = _nodes[i]->getXYZ();}
  for(int i=0; i<nn; i++) {temp.unsafe_col(i) = _nodes[i]->xyzvec3();}
  return temp; 
  */
//}

//const arma::mat MEl::linearNodes() const{
  /*
  const arma::uvec& linear_nodes = _sf->getLinearNodes();
  const int nn = linear_nodes.n_elem;
  arma::mat temp(3,nn);
  for(int i=0; i<nn; i++){
    temp.unsafe_col(i) = _nodes[linear_nodes(i)]->xyzvec3();
  }
  return temp;
  */
//}






/*
std::vector<int> MEl::getSortedNodes() const {

  const int nn = _nodes.size();
  std::vector<int> sorted_nodes(nn);
  for(int i=0; i<nn; i++) sorted_nodes[i] = _nodes[i]->nd();

  std::sort(sorted_nodes.begin(),sorted_nodes.end());

  return sorted_nodes;
}
*/







void MEl::setOrder(int p){
  /*
  _porder = p;

  const int nn2 = _sf->Ndof();
 
  // Get interior nodes
  std::vector<node_ptr> intnodes = projectPhysicalNodes();

  // re-size nodes
  _nodes.resize(nn2);

  
  // Move corner nodes
  const arma::uvec& cn_new = _sf->getCornerNodes();
  const arma::uvec& cn_orig = _sf->getCornerNodesMap(); 
  const int ncorn = cn_new.n_elem;
  std::vector<node_ptr> temp(ncorn);
  for(auto i=0; i<ncorn; i++){
    temp[i] = _nodes[cn_orig(i)];
  }
  for(auto i=0; i<ncorn; i++){
    _nodes[cn_new(i)] = temp[i];
  }
  
  //getChildNodes();

  // Get interior nodes
  //std::vector<node_ptr> intnodes = projectPhysicalNodes();

  // Copy over interior nodes
  const arma::uvec& in_new = _sf->getInteriorNodes();
  for(auto i=0; i<intnodes.size(); i++){
    _nodes[in_new(i)] = intnodes[i];
  }
  */

  
}


//const std::vector<element_ptr> MEl::generateChildren(){
  /*
  int N = _sf->Nchildren();
  _children.resize(N);
 

  for(auto i=0; i<N; i++){
    const arma::uvec& child_nodes = _sf->getBoundaryNodes(i);
    const int nnc = child_nodes.size();
    ShapeFunction* child_sf = _sf->getChildShapeFunction(i);
    std::vector<node_ptr> nodes(nnc);
    for(auto j=0; j<nnc; j++) nodes[j] = _nodes[child_nodes(j) ];
    if(dim() == 2){
      myFace* face = static_cast<myFace*>(getPhysical());
     _children[i] = make_shared<El1D>(El1D(0,0,nodes,child_sf,0,face));
    }
    else if(dim() == 3){
      _children[i] = make_shared<El2D>(El2D(0,0,nodes,child_sf));
    }
    _children[i]->setParent(i,this);
    _child_orientation[i] = 1; 

  }
  return _children;
  */
//}

/*
const std::vector<element_ptr> El3D::generateChildren() {
  //int El3D::generateChildren(){
  int N = _sf->Nchildren();
  _children.resize(N);
  for(auto i=0; i<N; i++){
    const arma::uvec& child_nodes = _sf->getBoundaryNodes(i);
    const int nnc = child_nodes.size();
    ShapeFunction* child_sf = _sf->getChildShapeFunction(i);
    std::vector<node_ptr> nodes(nnc);
    for(auto j=0; j<nnc; j++) nodes[j] = _nodes[child_nodes(j) ];
    //_children[i] = make_shared<El2D>(El2D(1,_bc_tag,nodes,child_sf));
    //_children[i]->setParent(i,this);
    _child_orientation[i] = 1; 

  }
  return _children;

  //return 1;
}
*/

// Returns true if element1 has the same node ordering of element 2
//signed char elOrientation(const element_ptr &el1, const element_ptr &el2){
  /*
  const std::vector<node_ptr>& nodes1 = el1->nodes_vector();
  const std::vector<node_ptr>& nodes2 = el2->nodes_vector();
  const ShapeFunction* sftemp = el1->getShapeFunction();

  const arma::umat& perm = sftemp->getPerm();
  for(int i=0; i<perm.n_cols; i++){
    if(nodes1[perm(0,i)]->nd() == nodes2[0]->nd() && nodes1[perm(1,i)]->nd() == nodes2[1]->nd()){
      return i+1;
    }
    else if(nodes1[perm(1,i)]->nd() == nodes2[0]->nd() && nodes1[perm(0,i)]->nd() == nodes2[1]->nd()){
      return -(i+1);
    }
  }
  cout << "Invalind node comparison in elOrientation!" << endl;
  for(int i=0; i<nodes1.size(); i++){
    cout << nodes1[i]->nd() << " " << nodes2[i]->nd() << endl;
  }
  */
//}

/*
bool elGreater::operator()(const element_ptr &el1, const element_ptr &el2) const{
 
  const int nn1 = el1->numCornerNodes();
  const int nn2 = el2->numCornerNodes();
  if(nn2 > nn1) return true;
  else if(nn2 < nn1) return false;
  else{

    const std::vector<int> nodes1 = el1->getSortedCornerNodes();
    const std::vector<int> nodes2 = el2->getSortedCornerNodes();

    //const unsigned char* el1_nodes = el1->getSortedNodes2();
    //const unsigned char* el2_nodes = el2->getSortedNodes2(); 

    //const std::vector<node_ptr> & nodes1 = el1->nodes_vector();
    //const std::vector<node_ptr> & nodes2 = el2->nodes_vector();
    for(int i=0; i<nn1; i++){
      //const int en1 = el1_nodes[i];
      //const int en2 = el2_nodes[i];
      if(nodes2[i] > nodes1[i]) return true;

      //if(el2_nodes[i] > el1_nodes[i]) return true;	
      //else if(el2_nodes[i] < el1_nodes[i]) return false;
    }
  }
  //if(el2->hasPhysical() < el1->hasPhysical()) return true;
 
  return false;
  
}
*/

/*
bool MEl::isGreater(const element_ptr &el1, const element_ptr &el2){
 
  const int nn1 = el1->numNodes();
  const int nn2 = el2->numNodes();
  if(nn2 > nn1) return true;
  else if(nn2 < nn1) return false;
  else{

    const unsigned char* el1_nodes = el1->getSortedNodes2();
    const unsigned char* el2_nodes = el2->getSortedNodes2(); 

    for(int i=0; i<nn1; i++){
      if(el2_nodes[i] > el1_nodes[i]) return true;	
      else if(el2_nodes[i] < el1_nodes[i]) return false;
    }
  }
  if(el2->hasPhysical() < el1->hasPhysical()) return true;
  return false;
  
}
  

bool MEl::isEqual(const element_ptr &el1, const element_ptr &el2){
  const int nn1 = el1->numNodes();
  const int nn2 = el2->numNodes();
  if(nn1 != nn2) return false;  


  const unsigned char* el1_nodes = el1->getSortedNodes2();
  const unsigned char* el2_nodes = el2->getSortedNodes2(); 
  

  for(int i=0; i<nn1; i++){
    if(el2_nodes[i] != el1_nodes[i]){
      return false;
    }
  }
  return true;
    
}
*/

arma::mat nodeParamsOnFace(const std::vector<node_ptr> &nodes,
			   const myFace* face){
  /*
  const int nn = nodes.size();
  arma::mat uv(2,nn);
  for(int i=0; i<nn; i++){
    
    if(nodes[i]->dim() == 0){ 
      MVertNode* temp = static_cast<MVertNode*>(nodes[i].get());
      uv.col(i) = face->paramsOnFace(temp->getVertex());
 
    }
    else if(nodes[i]->dim() == 1){
      MEdgeNode* temp = static_cast<MEdgeNode*>(nodes[i].get());
      const myEdge* test = temp->getEdge();
      double param = temp->getParam();

      uv.col(i) = face->paramsOnFace(temp->getEdge(),temp->getParam());
    }
    else if(nodes[i]->dim() == 2){

      MFaceNode* temp = static_cast<MFaceNode*>(nodes[i].get());
      if(!temp) cout << "note temp is null!" << endl;
      uv.col(i) = temp->getParams();
    }
    
    
  }
  if(!is_finite(uv)){
    cout << "uv is not finite in node params on face!" << endl;
  }
  return uv;
  */
}



std::vector<node_ptr> projectPointsOnFace(const arma::mat xyz0, 
					  const arma::mat& uv0,
					  const myFace* face){

  /*
  const int nn = xyz0.n_cols;
  std::vector<node_ptr> nodes(nn);

  for(int i=0; i<nn; i++){
    arma::vec2 uv = face->faceParamsFromPoint(xyz0.colptr(i),uv0.colptr(i));
    arma::vec3 xyz = face->params2xyz(uv(0),uv(1));
    if(!is_finite(xyz)){
      
      cout << "xyz is not finite in projectPointsOnFace!" << endl;
      //cout << xyz0 << endl;
      cout << uv0 << endl;
      //cout << uv << endl;
      //cout << xyz << endl;
      //assert(1);
    }
    nodes[i] = make_shared<MFaceNode>
      (MFaceNode(0,xyz.memptr(),2,0,uv.memptr(),face));
  }

  return nodes;
  */
}

std::vector<node_ptr> createFreeNodes(arma::mat& xyz){
  /*
  const int nn = xyz.n_cols;
  std::vector<node_ptr> nodes(nn);
  for(auto i=0; i<nn; i++){
    nodes[i] = make_shared<MVolumeNode>(MVolumeNode(0,xyz.colptr(i),0));
  }
  return nodes;
  */
}






double myFastInverse(const arma::mat& A, arma::mat& B){
  const unsigned char dim = A.n_rows;
  double det;
  if(dim == 1){
    det = A[0];
    B[0] = 1.0/det;
    return det;
  }
  if(dim == 2){
    det = A[0]*A[3] - A[1]*A[2];
    B[0] = A[3];
    B[3] = A[0];
    B[1] = -A[1];
    B[2] = -A[2];
    for(unsigned char i=0; i<4; i++) B[i]/=det;
    return det;
  }
  else if(dim == 3){
    det = A[0]*A[4]*A[8] - A[0]*A[5]*A[7] - A[1]*A[3]*A[8] + 
      A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[2]*A[4]*A[6];
    
    B[0] = A[4]*A[8] - A[5]*A[7];
    B[3] = A[5]*A[6] - A[3]*A[8];
    B[6] = A[3]*A[7] - A[4]*A[6];
    B[1] = A[2]*A[7] - A[1]*A[8];
    B[4] = A[0]*A[8] - A[2]*A[6];
    B[7] = A[1]*A[6] - A[0]*A[7];
    B[2] = A[1]*A[5] - A[2]*A[4];
    B[5] = A[2]*A[3] - A[0]*A[5];
    B[8] = A[0]*A[4] - A[1]*A[3];
    for(unsigned char i=0; i<9; i++) B[i]/=det;
    return det;
  }
  else{
    assert("dimension sould always be <= 3");
  }
}



/*
const arma::mat El2D::computeJacobian(vec3* grad) const{
  const arma::uvec& cornern = _sf->getLinearNodes();
  vec3 r0 = _nodes[cornern[0]]->xyzvec3();
  vec3 r1 = _nodes[cornern[1]]->xyzvec3()-r0;
  vec3 r2 = _nodes[cornern[2]]->xyzvec3()-r0;
  vec3 n = cross(r1,r2);

  

  mat::fixed<3,2> temp;
  temp.unsafe_col(0) = grad[0];
  temp.unsafe_col(1) = grad[1];
 
  
  //temp.unsafe_col(0)/= norm(temp.unsafe_col(0),2.0);
  //temp.unsafe_col(1)/= norm(temp.unsafe_col(1),2.0);
  vec3 nn = cross(temp.unsafe_col(0),temp.unsafe_col(1));
  if(dot(nn,n) < 0.0){
    //nn = -nn;
    cout << "normal is reversed" << endl;
    temp.unsafe_col(1) = cross(temp.unsafe_col(0),nn);
  }
  else{
  //temp.unsafe_col(1) = cross(temp.unsafe_col(0),nn);
  temp.unsafe_col(1) = cross(nn,temp.unsafe_col(0));
  }
  temp.unsafe_col(0)/=norm(temp.unsafe_col(0),2.0);
  temp.unsafe_col(1)/=norm(temp.unsafe_col(1),2.0);
  mat::fixed<3,2> rawgrads;
  rawgrads.col(0) = grad[0];
  rawgrads.col(1) = grad[1];

  mat22 S = temp.t()*rawgrads;
  //if(dot(nn,n) < 0.0){
  //  S(0)*= -1;
  //}
 
  return S;

 
}

const arma::mat El3D::computeJacobian(vec3* grad) const{
  mat33 J;
  for(int i=0; i<3; i++){
    J.unsafe_col(i) = grad[i];
    J.unsafe_col(i) = grad[i];
  }
  return J;
}



const arma::vec El2D::computeGradDetJacobian(arma::vec3* grad) const{
  arma::vec6 Grad;

  const double A0 = grad[0][0];
  const double A1 = grad[0][1];
  const double A2 = grad[0][2];
  const double A3 = grad[1][0];
  const double A4 = grad[1][1];
  const double A5 = grad[1][2];

  const double A02 = A0*A0;
  const double A12 = A1*A1;
  const double A22 = A2*A2;
  const double A32 = A3*A3;
  const double A42 = A4*A4;
  const double A52 = A5*A5;
 

  double denom = A02*A42 + A02*A52 - 2.*A0*A1*A3*A4 - 2.*A0*A2*A3*A5 + A12*A32 + A12*A52 - 2.*A1*A2*A4*A5 + A22*A32 + A22*A42;
  denom = sqrt(denom);

  Grad[0] = (A0*A42 - A1*A3*A4 + A0*A52 - A2*A3*A5)/denom;
  Grad[1] = (A1*A32 - A0*A4*A3 + A1*A52 - A2*A4*A5)/denom;
  Grad[2] = (A2*A32 - A0*A5*A3 + A2*A42 - A1*A5*A4)/denom;

  Grad[3] = (A3*A12 - A0*A4*A1 + A3*A22 - A0*A5*A2)/denom;
  Grad[4] = (A4*A02 - A1*A3*A0 + A4*A22 - A1*A5*A2)/denom;
  Grad[5] = (A5*A02 - A2*A3*A0 + A5*A12 - A2*A4*A1)/denom;

  return Grad;
}

const arma::vec El3D::computeGradDetJacobian(arma::vec3* grad) const{
  arma::mat33 A, invA;
  for(int i=0; i<3; i++) A.unsafe_col(i) = grad[i];
  myFastInverse(A,invA);
  invA*= det(A);
  return vectorise(invA);
}



const arma::mat El2D::computeHessDetJacobian(arma::vec3* grad) const{
  arma::mat66 Hess;

  const double A0 = grad[0][0];
  const double A1 = grad[0][1];
  const double A2 = grad[0][2];
  const double A3 = grad[1][0];
  const double A4 = grad[1][1];
  const double A5 = grad[1][2];

  const double A02 = A0*A0;
  const double A12 = A1*A1;
  const double A22 = A2*A2;
  const double A32 = A3*A3;
  const double A42 = A4*A4;
  const double A52 = A5*A5;
 
  
  double denom = A02*A42 + A02*A52 - 2.*A0*A1*A3*A4 - 2.*A0*A2*A3*A5 + A12*A32 + A12*A52 - 2.*A1*A2*A4*A5 + A22*A32 + A22*A42;
  double denom12 = sqrt(denom);
  denom = pow(denom,1.5);
 


  Hess(0,0) = pow(A1*A5 - A2*A4,2)*(A32 + A42 + A52)/denom;
  Hess(0,1) = -(A0*A5 - A2*A3)*(A1*A5 - A2*A4)*(A32 + A42 + A52)/denom;
  Hess(0,2) = (A0*A4 - A1*A3)*(A1*A5 - A2*A4)*(A32 + A42 + A52)/denom;
  Hess(0,3) = -pow(A1*A5 - A2*A4,2)*(A0*A3 + A1*A4 + A2*A5)/denom;
  Hess(0,4) = (2.*A0*A4 - A1*A3)/denom12 - ((A0*A42 - A1*A3*A4 + A0*A52 - A2*A3*A5)*(2*A4*A02 - 2*A1*A3*A0 + 2*A4*A22 - 2*A1*A5*A2))/(2.*denom);
  Hess(0,5) = (2.*A0*A5 - A2*A3)/denom12 - ((A0*A42 - A1*A3*A4 + A0*A52 - A2*A3*A5)*(2*A5*A02 - 2*A2*A3*A0 + 2*A5*A12 - 2*A2*A4*A1))/(2.*denom);

  Hess(1,0) = -(A0*A5 - A2*A3)*(A1*A5 - A2*A4)*(A32 + A42 + A52)/denom;
  Hess(1,1) = pow(A0*A5 - A2*A3,2)*(A32 + A42 + A52)/denom;
  Hess(1,2) = -(A0*A4 - A1*A3)*(A0*A5 - A2*A3)*(A32 + A42 + A52)/denom;
  Hess(1,3) = -(A0*A4 - 2*A1*A3)/denom12 -((2*A3*A12 - 2*A0*A4*A1 + 2*A3*A22 - 2*A0*A5*A2)*(A1*A32 - A0*A4*A3 + A1*A52 - A2*A4*A5))/(2.*denom);
  Hess(1,4) = -pow(A0*A5 - A2*A3,2)*(A0*A3 + A1*A4 + A2*A5)/denom;
  Hess(1,5) = (2*A1*A5 - A2*A4)/denom12 - ((A1*A32 - A0*A4*A3 + A1*A52 - A2*A4*A5)*(2*A5*A02 - 2*A2*A3*A0 + 2*A5*A12 - 2*A2*A4*A1))/(2.*denom);

  Hess(2,0) = (A0*A4 - A1*A3)*(A1*A5 - A2*A4)*(A32 + A42 + A52)/denom;
  Hess(2,1) = -(A0*A4 - A1*A3)*(A0*A5 - A2*A3)*(A32 + A42 + A52)/denom;
  Hess(2,2) = pow(A0*A4 - A1*A3,2)*(A32 + A42 + A52)/denom;
  Hess(2,3) = -(A0*A5 - 2*A2*A3)/denom12 - ((2*A3*A12 - 2*A0*A4*A1 + 2*A3*A22 - 2*A0*A5*A2)*(A2*A32 - A0*A5*A3 + A2*A42 - A1*A5*A4))/(2.*denom);
  Hess(2,4) = -(A1*A5 - 2*A2*A4)/denom12 - ((2*A4*A02 - 2*A1*A3*A0 + 2*A4*A22 - 2*A1*A5*A2)*(A2*A32 - A0*A5*A3 + A2*A42 - A1*A5*A4))/(2.*denom);
  Hess(2,5) = -pow(A0*A4 - A1*A3,2)*(A0*A3 + A1*A4 + A2*A5)/denom;

  Hess(3,0) = -pow(A1*A5 - A2*A4,2)*(A0*A3 + A1*A4 + A2*A5)/denom;
  Hess(3,1) = -(A0*A4 - 2*A1*A3)/denom12 - ((A3*A12 - A0*A4*A1 + A3*A22 - A0*A5*A2)*(2*A1*A32 - 2*A0*A4*A3 + 2*A1*A52 - 2*A2*A4*A5))/(2.*denom);
  Hess(3,2) = -(A0*A5 - 2*A2*A3)/denom12 - ((A3*A12 - A0*A4*A1 + A3*A22 - A0*A5*A2)*(2*A2*A32 - 2*A0*A5*A3 + 2*A2*A42 - 2*A1*A5*A4))/(2.*denom);
  Hess(3,3) = pow(A1*A5 - A2*A4,2)*(A02 + A12 + A22)/denom;
  Hess(3,4) = -(A0*A5 - A2*A3)*(A1*A5 - A2*A4)*(A02 + A12 + A22)/denom;
  Hess(3,5) = (A0*A4 - A1*A3)*(A1*A5 - A2*A4)*(A02 + A12 + A22)/denom;

  Hess(4,0) = (2*A0*A4 - A1*A3)/denom12 - ((A4*A02 - A1*A3*A0 + A4*A22 - A1*A5*A2)*(2*A0*A42 - 2*A1*A3*A4 + 2*A0*A52 - 2*A2*A3*A5))/(2.*denom);
  Hess(4,1) = -pow(A0*A5 - A2*A3,2)*(A0*A3 + A1*A4 + A2*A5)/denom;
  Hess(4,2) = -(A1*A5 - 2*A2*A4)/denom12 - ((A4*A02 - A1*A3*A0 + A4*A22 - A1*A5*A2)*(2*A2*A32 - 2*A0*A5*A3 + 2*A2*A42 - 2*A1*A5*A4))/(2.*denom);
  Hess(4,3) = -(A0*A5 - A2*A3)*(A1*A5 - A2*A4)*(A02 + A12 + A22)/denom;
  Hess(4,4) = pow(A0*A5 - A2*A3,2)*(A02 + A12 + A22)/denom;
  Hess(4,5) = -(A0*A4 - A1*A3)*(A0*A5 - A2*A3)*(A02 + A12 + A22)/denom;

  Hess(5,0) = (2*A0*A5 - A2*A3)/denom12 - ((A5*A02 - A2*A3*A0 + A5*A12 - A2*A4*A1)*(2*A0*A42 - 2*A1*A3*A4 + 2*A0*A52 - 2*A2*A3*A5))/(2.*denom);
  Hess(5,1) = (2*A1*A5 - A2*A4)/denom12 - ((A5*A02 - A2*A3*A0 + A5*A12 - A2*A4*A1)*(2*A1*A32 - 2*A0*A4*A3 + 2*A1*A52 - 2*A2*A4*A5))/(2.*denom);
  Hess(5,2) = -pow(A0*A4 - A1*A3,2)*(A0*A3 + A1*A4 + A2*A5)/denom;
  Hess(5,3) = (A0*A4 - A1*A3)*(A1*A5 - A2*A4)*(A02 + A12 + A22)/denom;
  Hess(5,4) = -(A0*A4 - A1*A3)*(A0*A5 - A2*A3)*(A02 + A12 + A22)/denom;
  Hess(5,5) = pow(A0*A4 - A1*A3,2)*(A02 + A12 + A22)/denom;

}

const arma::mat El3D::computeHessDetJacobian(arma::vec3* grad) const{
  vec9 Grad;
  for(int i=0,cnt=0; i<3; i++){
    for(int j=0; j<3; j++,cnt++) Grad[j] = grad[i][j];
  }
  mat99 Hess;
  Hess(0,4) = Grad[8];  Hess(0,8) = Grad[4]; 
  Hess(0,5) = -Grad[7]; Hess(0,7) = -Grad[5];

  Hess(1,5) = Grad[6]; Hess[1,6] = Grad[5];
  Hess(1,3) = -Grad[8]; Hess(1,8) = -Grad[3];

  Hess(2,3) = Grad[7]; Hess(2,7) = Grad[3];
  Hess(2,4) = -Grad[6]; Hess(2,6) = -Grad[4];

  Hess(3,2) = Grad[7]; Hess(3,7) = Grad[2];
  Hess(3,1) = -Grad[8]; Hess(3,8) = -Grad[1];

  Hess(4,0) = Grad[8]; Hess(4,8) = Grad[0];
  Hess(4,2) = -Grad[6]; Hess(4,6) = -Grad[2];

  Hess(5,1) = Grad[6]; Hess(5,6) = Grad[1];
  Hess(5,0) = -Grad[7]; Hess(5,7) = -Grad[0];

  Hess(6,1) = Grad[5]; Hess(6,5) = Grad[1];
  Hess(6,2) = -Grad[4]; Hess(6,4) = -Grad[2];

  Hess(7,2) = Grad[3]; Hess(6,3) = Grad[2];
  Hess(7,0) = -Grad[5]; Hess(7,5) = -Grad[0];

  Hess(8,0) = Grad[4]; Hess(8,4) = -Grad[0];
  Hess(8,1) = -Grad[3]; Hess(8,3) = -Grad[1];


  return Hess;
}



const arma::vec El2D::computeGradFroNorm(arma::vec3* grad) const{

}

const arma::vec El3D::computeGradFroNorm(arma::vec3* grad) const{
  mat33 J;
  for(int i=0; i<3; i++) J.unsafe_col(i) = grad[i];
  double fronorm = norm(J,"fro");
  return vectorise(J)/fronorm;
}



const arma::mat El2D::computeHessFroNorm(arma::vec3* grad) const{

}

const arma::mat El3D::computeHessFroNorm(arma::vec3* grad) const{
  mat33 J;
  for(int i=0; i<3; i++) J.unsafe_col(i) = grad[i];
  double fronorm = norm(J,"fro");
  double fronorm3 = 1.0/pow(fronorm,3);

  mat99 Hess;
  for(int i=0; i<9; i++){
    for(int j=0; j<9; j++){
      Hess(i,j) = -J[i]*J[j]*fronorm3;
      if(i==j) Hess(i,j)+= 1.0/fronorm;
    }
  }

  return Hess;
}

*/

/*
const double MEl::computeMeritFunction() {
  const int ndof = _nodes.size();

  //const arma::uvec cornern = _sf->getLinearNodes();
  const arma::mat linnodes = linearNodes();
  arma::mat nodesideal;
  bool corner_match = true;
  for(int i=0; i<linnodes.n_elem; i++){
    if(std::abs(linnodes[i]-nodesideal[i]) > 1.0e-12){
      corner_match = false;
      break;
    }
  }
  
  if(corner_match){
    if(!is_curved()){
      //_distortion = 1.0;
      return 1.0;
    }
  }



  const double eps = 2.0*std::numeric_limits<double>::epsilon();
  const int ng =_sf->sf.n_cols;
  const int dm = dim();
  
  
  
  const mat nodesel = nodes();
  arma::mat phys_deriv[3];
  
  for(int i=0; i<dm; i++){
    phys_deriv[i] = nodesel*_sf->sfderiv.slice(i);
  }

  
  //arma::mat nodesideal;
  arma::mat ideal_deriv[3];
  for(int i=0; i<dm; i++){
    ideal_deriv[i] = nodesideal*_sf->sflinderiv.slice(i);
  }


  double element_area=0.0;
  double eta_shap=0.0;
  mat S(dm,dm), SI(dm,dm);
  mat invSI(dm,dm);
  
  
  
  for(int i=0; i<ng; i++){

    vec3 phys_gpderivs[3];
    vec3 ideal_gpderivs[3];
    for(int j=0; j<dm; j++){
      phys_gpderivs[j] = phys_deriv[j].unsafe_col(i);
      ideal_gpderivs[j] = ideal_deriv[j].unsafe_col(i);
    }

    S = computeJacobian(phys_gpderivs);
 
    SI = computeJacobian(ideal_gpderivs);

    //double dj = std::abs(det(S));
    double dj = std::abs(det(S));
    element_area+= dj*_sf->gw(i);
    if(dj < 0.0) cout << "dj: " << endl;

    invSI = inv(SI);
    //myFastInverse(SI,invSI);

    S = invSI*S;
 
    
    //double detS = det(S);
    double detS = std::abs(det(S));
    if(detS < 0.0) cout << detS << endl;

    
    double delta;
    if(detS < eps) delta = sqrt(eps*(eps-detS));
    else delta = 0.0;
    
    double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));
    
    if(sigma < 0.01) cout << sigma << endl;
    double fronorm = norm(S,"fro");

    
    double denom;
    if(dm == 1){
      denom = sigma*sigma;
    }
    else if(dm == 2){
      denom = 2.*sigma;
    }
    else{
      double cbs = cbrt(sigma);
      denom = 3.*cbs*cbs;
    }

    double eta = fronorm*fronorm/denom;

    eta_shap+= eta*eta*dj*_sf->gw(i);
    
  }
  
  double distortion = sqrt(eta_shap/abs(element_area));
  //_distortion = distortion;
 

  return distortion;
  
}

arma::vec MEl::computeGradMeritFunction() const {

  const double eps = 2.0*std::numeric_limits<double>::epsilon();
  const int ng =_sf->sf.n_cols;
  const int dm = dim();
  
  
  
  const mat nodesel = nodes();
  arma::mat phys_deriv[3];
  
  for(int i=0; i<dm; i++){
    phys_deriv[i] = nodesel*_sf->sfderiv.slice(i);
  }

  arma::mat nodesideal;

  arma::mat ideal_deriv[3];
  for(int i=0; i<dm; i++){
    ideal_deriv[i] = nodesideal*_sf->sflinderiv.slice(i);
  }


  double element_area=0.0;
  double eta_shap=0.0;
  mat S(dm,dm), SI(dm,dm);
  mat invSI(dm,dm);
  
  mat ddj_dS;
  
  for(int i=0; i<ng; i++){
    vec3 phys_gpderivs[3];
    vec3 ideal_gpderivs[3];
    for(int j=0; j<dm; j++){
      phys_gpderivs[j] = phys_deriv[j].unsafe_col(i);
      ideal_gpderivs[j] = ideal_deriv[j].unsafe_col(i);
    }

    S = computeJacobian(phys_gpderivs);
    SI = computeJacobian(ideal_gpderivs);

    double detS = det(S);
    double dj = std::abs(detS);
    
    ddj_dS = Sign(detS)*inv(S);

    element_area+= dj*_sf->gw(i);

    myFastInverse(SI,invSI);

    S = invSI*S;

  }


}
*/



/*
bool El2D::isGreater(const element2D_ptr &el1, const element2D_ptr &el2){
 
  const int nn1 = el1->numNodes();
  const int nn2 = el2->numNodes();
  if(nn2 > nn1) return true;
  else if(nn2 < nn1) return false;
  else{
    const vector<node_ptr>& el1_nodes = el1->getSortedNodes2();
    const vector<node_ptr>& el2_nodes = el2->getSortedNodes2(); 
    for(int i=0; i<el1_nodes.size(); i++){
      if(el2_nodes[i] > el1_nodes[i]) return true;	
      else if(el2_nodes[i] < el1_nodes[i]) return false;
    }
  }
  if(el2->hasFace() < el1->hasFace()) return true;
  return false;
  
}
  

bool El2D::isEqual(const element2D_ptr &el1, const element2D_ptr &el2){
  if(el2->numNodes() != el1->numNodes()) return false;

  const vector<node_ptr>& el1_nodes = el1->getSortedNodes2();
  const vector<node_ptr>& el2_nodes = el2->getSortedNodes2(); 
    
  for(int i=0; i<el1_nodes.size(); i++){
    if(el2_nodes[i] != el1_nodes[i]) return false;
  }
  return true;
    
}
*/

/*
using namespace arma;
int MTri::writeTecNodes(ofstream &out){

  const int ndof = _nodes.size();
  const int npn = _sf->sf_plot.n_cols;
  mat elnodes(3,ndof);
  for(int i=0; i<ndof; i++){
    elnodes.col(i) = _nodes[i]->getXYZ();
  }

  mat plotnodes = elnodes*_sf->sf_plot;

  for(int i=0; i<npn; i++){
    for(int j=0; j<3; j++) out << plotnodes(j,i) << " ";
    out << 1.0 << endl;
  }

  if(_trino == 300){
    //plotnodes.save("plotnodes",raw_ascii);
    //_sf->plot_connectivity.save("connectivity",raw_ascii);
  }

}

int MTri::writeTecEls(ofstream &out){
  
  
  int npref = _sf->sf_plot.n_cols;
  const int ntl = _sf->plot_connectivity.n_rows;

  for(int el=0; el<ntl; el++){
    for(int j=0; j<3; j++) out << _trino*npref + _sf->plot_connectivity(el,j) + 
			     1 << " ";
    out << endl;
  }
  
}

int MTet::writeTecNodes(ofstream &out){
 
  const int ndof = _nodes.size();
  const int npn = _sf->sf_plot.n_cols;
  mat elnodes(3,ndof);
  for(int i=0; i<ndof; i++){
    elnodes.col(i) = _nodes[i]->getXYZ();
  }
  mat plotnodes = elnodes*_sf->sf_plot;

  for(int i=0; i<npn; i++){
    for(int j=0; j<3; j++) out << plotnodes(j,i) << " ";
    out << 1.0 << endl;
  }

  if(_tetno == 0) plotnodes.save("plotnodes",raw_ascii);
}

int MTet::writeTecEls(ofstream &out){
  
  int npref = _sf->sf_plot.n_cols;
  const int ntl = _sf->plot_connectivity.n_rows;

  for(int el=0; el<ntl; el++){
    for(int j=0; j<_sf->plot_connectivity.n_cols; j++){
      out << _tetno*npref + _sf->plot_connectivity(el,j) + 1 << " ";
    }
    out << endl;
  }
  //if(_tetno == 0) _sf->plot_connectivity.save("plot_conn",raw_ascii);
}


int MTri::generateChildren(const sf_type &sfh){

  std::vector<node_ptr> nodes(2);
  
  for(int i=0; i<3; i++){
    for(int j=0; j<2; j++) nodes[j] = _nodes[_sf->child_nodes(j,i)];
    element1D_ptr temp = make_shared<El1D>(El1D(1,1,nodes,_bc_tag,&sfh->Line,0,_face));
    temp->setParent(i,this);
    //_children[i] = temp.get();
    _children.push_back(temp);
    _child_orientation[i] = true;
  }

  return 1;
}


int MTet::generateChildren(const sf_type &sfh){
  //std::vector<element2D_ptr> children;
 
  std::vector<node_ptr> nodes(3);
  //_child_orientation.resize(4);
  for(int i=0; i<4; i++){
    for(int j=0; j<3; j++) nodes[j] = _nodes[_sf->child_nodes(j,i)];
    element2D_ptr temp = make_shared<MTri>(MTri(1,1,nodes,_bc_tag,&sfh->Tri));
    //temp->setParent(0,i,this);
    temp->setParent(i,this);
    _children.push_back(temp);
    _child_orientation[i] = true;
  }
  //return children;
  return 1;
}

*/

/*
MTri::MTri(int el, int trino, std::vector<node_ptr> nodes, unsigned char
	   bc_tag, ShapeFunction* sf, const myFace* face):
  El2D(el,bc_tag,nodes,sf,face), _trino(trino){
 
}
*/

/*

arma::vec MTri::getGradMerit(int node) const{
  int nnode = _nodes.size();
  if(_gradMeritAnalytical.size() != _nodes.size()){
    cout << "size of grad merit does nt equal size of nodes!" << endl;
  }
  for(int i=0; i<nnode; i++){
    if(_nodes[i]->nd() == node){
      return _gradMeritAnalytical[i];
      //return _gradMerit[i];
    }
  }
 
  cout << "node " << node << " not found in element " << _el << "!\n";

  
}

arma::mat MTri::getHessMerit(int node) const{
  int nnode = _nodes.size();
  for(int i=0; i<nnode; i++){
    if(_nodes[i]->nd() == node){
      return _HessianMerit[i];
    }
  }

  cout << "node " << node << " not found in element " << _el << "!\n";

  
}

arma::vec MTri::gradMeritFD(node_ptr node, double step){
  
  vec pcoords0 = node->getParams();

  int npar = pcoords0.n_elem;
  vec gradMerit(npar);
  //double merit0 = computeMeritFunction();
  for(int i=0; i<npar; i++){
    vec pcoordsFD = pcoords0;
    pcoordsFD(i)+= step;
    node->updateParams(pcoordsFD);
    double meritPlus = computeMeritFunction();
    
    pcoordsFD = pcoords0;
    pcoordsFD(i)-= step;
    node->updateParams(pcoordsFD);
    double meritMinus = computeMeritFunction();

    gradMerit(i) = (meritPlus-meritMinus)/step;
  }

  node->updateParams(pcoords0);

  return gradMerit;
  
}

arma::mat MTri::HessianMeritFD(node_ptr node, double step){
  
  vec pcoords0 = node->getParams();
  int npar = pcoords0.n_elem;
 
  
  mat HessMerit(npar,npar);
  for(int i=0; i<npar; i++){
    vec pcoordsFD = pcoords0;
    pcoordsFD(i)+= step;
    node->updateParams(pcoordsFD);
    vec gradGradPlus = gradMeritFD(node,1.0e-8);
    pcoordsFD = pcoords0;
    pcoordsFD(i)-= step;
    node->updateParams(pcoordsFD);
    vec gradGradMinus = gradMeritFD(node,1.0e-8);
    HessMerit.col(i) = (gradGradPlus-gradGradMinus)/step;
  }

  node->updateParams(pcoords0);

  return HessMerit;
  
}


int MTri::computeGradMerit(void){

  
   computeGradDistortion();
  _gradMerit.resize(_nodes.size());
  _HessianMerit.resize(_nodes.size());
  for(int nd=0; nd<_nodes.size(); nd++){
    if(_nodes[nd]->dim() > 0){
      vec pcoords0 = _nodes[nd]->getParams();
      int npar = pcoords0.n_elem;
      //_gradMerit[nd].resize(npar);
      //_HessianMerit[nd].resize(npar,npar);
      //_gradMerit[nd] = gradMeritFD(_nodes[nd],1.0e-8);
      _gradMerit[nd] = _gradMeritAnalytical[nd];
      //_HessianMerit[nd] = HessianMeritFD(_nodes[nd],1.0e-4);
      // vec deb = inv(_HessianMerit[nd])*_gradMerit[nd];
      //if(!is_finite(deb)) cout << "product is not finite!" << endl;
      //if(norm(deb,2.0) > 1.0e2) cout << norm(deb,2.0) << endl;
 
    }
  }
  return 1;
  

}

void MTri::computeGradMeritNode(int ndg, arma::vec &gradM, arma::mat &HessM){
  int nd=-1;
  for(int i=0; i<_nodes.size(); i++){
    if(_nodes[i]->nd() == ndg){
      nd = i;
      break;
    }
  }

  if(nd == -1) cout << "node was not found!" << endl;
 

  wall_clock timer;
  if(_nodes[nd]->dim() > 0){   
    computeGradDistortion();
    gradM = _gradMeritAnalytical[nd];
  }
 
 

}

double MTri::computeMeritFunction(void){
  // double merit_function;
  //double eps = std::numeric_limits<double>::epsilon();
  //double delta = 10.0*eps;
  
  //double delta = 0.0;
  //double distortion_measure = computeDistortion(delta);
 

  //merit_function = 0.5*pow(distortion_measure-1.0,2.0);
  
  //_merit = merit_function;
  
  return computeDistortion();
}

double MTri::computeDistortion() {

  
  const double eps = 2.0*std::numeric_limits<double>::epsilon();
  const int ndof = _nodes.size();
  const int ng =_sf->sf2d_xi.n_cols;

  mat nodesel(3,ndof);
  for(int i=0; i<ndof; i++) nodesel.col(i) = _nodes[i]->getXYZ();
  
  mat xyz_xi = nodesel*_sf->sf2d_xi;
  mat xyz_et = nodesel*_sf->sf2d_et;


  double element_area=0.0;
  double eta_shap=0.0;
  mat22 S;
  for(int i=0; i<ng; i++){
    vec3 normal = cross(xyz_xi.unsafe_col(i),xyz_et.unsafe_col(i));
    double dj = norm(normal,2.0);
    normal/=dj;
    element_area+= dj*_sf->gw2d(i);
    S.col(0) = pseudoInv*xyz_xi.unsafe_col(i);
    S.col(1) = pseudoInv*xyz_et.unsafe_col(i);
  
    double detS = S(0,0)*S(1,1) - S(0,1)*S(1,0);

    double delta;
    if(detS < eps) delta = sqrt(eps*(eps-detS));
    else delta = 0.0;

    double sigma = 0.5*(detS + sqrt(detS*detS + 4.0*delta*delta));
    double eta = pow(norm(S,"fro"),2.0)/(2.0*sigma);
    eta_shap+= pow(eta,2.0)*dj*_sf->gw2d(i);
  }

  _distortion = sqrt(eta_shap/abs(element_area));
  _merit = 0.5*pow(_distortion-1.0,2);


  return _merit;
  
}
*/

//int MTri::initialize(){};

/*
int MTri::initialize(){
    const int ndof = _nodes.size();
    vec3 r0 = _nodes[0]->getXYZ();
    vec3 r1 = _nodes[_porder]->getXYZ()-r0;
    vec3 r2 = _nodes[ndof-1]->getXYZ()-r0;
    vec3 l = r1;
    vec3 n = cross(r1,r2);
    vec3 m = cross(n,l);
    area_ideal = norm(l,2.0);

    l/=norm(l,2.0);
    m/=norm(m,2.0);
    n/=norm(n,2.0);
    pseudoInv.row(0) = l.t();
    pseudoInv.row(1) = m.t();

    mat33 R;
    R.col(0) = l;
    R.col(1) = m;
    R.col(2) = n;
    mat33 dx;
    dx.zeros();
    dx.col(1) = r1;
    dx.col(2) = r2;

    mat33 planepts = R.t()*dx;
  
    mat22 W;
    W(0,0) = planepts(0,1); W(0,1) = planepts(0,2);
    W(1,0) = planepts(1,1); W(1,1) = planepts(1,2);

  
    mat22 invW = inv(W);
    pseudoInv = invW*pseudoInv;


}
*/

/*

int MTri::computeGradDistortion(){
 
  
  const double eps = 2.0*std::numeric_limits<double>::epsilon();
    const int ndof = _nodes.size();
 

    mat nodesel(3,ndof);
    for(int i=0; i<ndof; i++) nodesel.col(i) = _nodes[i]->getXYZ();
  
    mat xyz_xi = nodesel*_sf->sf2d_xi;
    mat xyz_et = nodesel*_sf->sf2d_et;

  

    wall_clock timer;
    timer.tic();

    mat::fixed<2,3> gradDJ;
    mat::fixed<2,3> gradDetS;
    mat::fixed<2,3> gradSigma;
    mat::fixed<2,3> gradNorm;
    mat::fixed<2,3> gradEta;
    mat::fixed<2,3> gradEtaShape;
    //gradEtaShape.zeros();

    mat grad_eta_accum(ndof,3,fill::zeros);
    mat grad_area_accum(ndof,3,fill::zeros);
    vec3 nrm;
    vec3 ons;
    ons.ones();

    const int ng = xyz_xi.n_cols;
    mat22 S;
    mat22 invS;
    double eta_shape_accum = 0.0;
    double area_accum = 0.0;
    bool isneg = false;
    mat temp(ndof,2);
    for(int i=0; i<ng; i++){
      nrm = cross(xyz_xi.col(i),xyz_et.col(i));
      double dj = norm(nrm,2.0);
      S.col(0) = pseudoInv*xyz_xi.col(i);
      S.col(1) = pseudoInv*xyz_et.col(i);
      

      double detS = S(0,0)*S(1,1) - S(0,1)*S(1,0);
      invS(0,0) = S(1,1);
      invS(1,1) = S(0,0);
      invS(0,1) = -S(0,1);
      invS(1,0)*= -S(1,0);
      invS/=detS;

      double delta;
      if(detS < eps) delta = sqrt(eps*(eps-detS));
      else delta = 0.0;

      double sigma = 0.5*(detS + sqrt(pow(detS,2.0) + 4.0*pow(delta,2.0)));

      double fronorm = norm(S,"fro");
      double eta = pow(norm(S,"fro"),2.0)/(2.0*abs(sigma));

      if(detS < 0.0) isneg = true;

      eta_shape_accum+= pow(eta,2.0)*dj*_sf->gw2d(i); 
      area_accum+= dj*_sf->gw2d(i);

      for(int j=0; j<3; j++){
	ons.zeros();
	ons(j) = 1.0;
	gradDJ(0,j) = sum(cross(ons,xyz_et.col(i))%nrm)/dj;
	gradDJ(1,j) = sum(cross(xyz_xi.col(i),ons)%nrm)/dj;
      }
      

      temp.col(0) = _sf->sf2d_xi.col(i);
      temp.col(1) = _sf->sf2d_et.col(i);
      
      gradDetS = detS*invS*pseudoInv;
      gradNorm = S.t()*pseudoInv/fronorm;

      gradSigma = 
	0.5*(gradDetS + detS*gradDetS/(sqrt(pow(detS,2.0) + 
					    4.0*pow(delta,2.0))));

      gradEta = fronorm/abs(detS)*gradNorm - 
	pow(fronorm,2.0) /(2.0*pow(abs(sigma),2.))*Sign(sigma)*gradSigma;

      gradEtaShape= (2.0*eta*dj*gradEta + pow(eta,2.0)*gradDJ);
      

      grad_eta_accum+= temp*gradEtaShape*_sf->gw2d(i);
      //grad_area_accum+= temp*gradDJ*_sf->gw2d(i);


    }
    
    
    //double distortion = sqrt(eta_shape_accum/area_ideal);
    //cout << "distortion: " << distortion << endl;

    //mat gradEtaEl = grad_eta_accum/(2*area_ideal*distortion);
    
 

    
     
    double distortion = sqrt(eta_shape_accum/abs(area_accum));
    mat gradEtaEl = (grad_eta_accum/abs(area_accum) - 
		     eta_shape_accum/pow(abs(area_accum),2.0)*
		     Sign(area_accum)*grad_area_accum)/distortion/2.0;

    //cout << "analytical time: " << timer.toc() << endl;
   _merit = 0.5*pow(distortion-1.0,2.0);
    
    vector<node_ptr>::iterator nd;
    _gradMeritAnalytical.resize(ndof);
    int it = 0;
    for(nd=_nodes.begin(); nd!= _nodes.end(); ++nd){
      mat xyz_params = (*nd)->D1();
      vec grad_merit = (distortion-1)*(gradEtaEl.row(it)*xyz_params).t();
      _gradMeritAnalytical[it] = grad_merit;
      it++;
    }

  
  return 1;
}

double MTri::HessMerit(){

  
  const int debel = 1492;
  const double eps = 2.0*std::numeric_limits<double>::epsilon();
  const int ndof = _nodes.size();
  const int ng = _sf->sf2d.n_cols;


  mat nodesel(3,ndof);
  for(int i=0; i<ndof; i++) nodesel.col(i) = _nodes[i]->getXYZ();
  
  mat xyz_xi = nodesel*_sf->sf2d_xi;
  mat xyz_et = nodesel*_sf->sf2d_et;

  mat22 S, invS;
  for(int gp=0; gp<ng; gp++){
    vec3 nrm = cross(xyz_xi.col(gp),xyz_et.col(gp));
    double detJ = norm(nrm,2.0);
    S.col(0) = pseudoInv*xyz_xi.col(gp);
    S.col(1) = pseudoInv*xyz_et.col(gp);
    double detS = S(0,0)*S(1,1) - S(0,1)*S(1,0);
    invS(0,0) = S(1,1);
    invS(1,1) = S(0,0);
    invS(0,1) = -S(0,1);
    invS(1,0)*= -S(1,0);
    invS/=detS;

    double delta;
    if(detS < eps) delta = sqrt(eps*(eps-detS));
    else delta = 0.0;
    double denom = sqrt(pow(detS,2.0) + 4.0*pow(delta,2.0));
    double sigma = 0.5*(detS + denom);
    //double FrobNorm = norm(S,"fro");
    double FrobNorm2 = pow(norm(S,"fro"),2.0);
    double eta = pow(norm(S,"fro"),2.0)/(2.0*abs(sigma));

    // Gradient of detS
    mat::fixed<2,3> grad_detS = detS*invS*pseudoInv;
    mat::fixed<2,3> grad_FrobNorm2 = 2.0*S.t()*pseudoInv;

    mat::fixed<2,3> grad_sigma = 
      0.5*(grad_detS + detS*grad_detS/(sqrt(pow(detS,2.0) + 
					    4.0*pow(delta,2.0))));

    mat::fixed<2,3> grad_eta = 0.5*(grad_FrobNorm2/sigma - 
				    FrobNorm2/pow(sigma,2.0)*grad_sigma);
	


    // Gradient of the surface scaling
    mat::fixed<2,3> grad_detJ;
    vec3 direc;
    direc.zeros();
    for(int j=0; j<3; j++){
      direc(j) = 1.0;
      grad_detJ(0,j) = sum(nrm%cross(direc,xyz_et.col(gp)))/detJ;
      grad_detJ(1,j) = sum(nrm%cross(xyz_xi.col(gp),direc))/detJ;
    }

    mat::fixed<2,3> grad_eta2detJ= (2.0*eta*detJ*grad_eta + 
				    pow(eta,2.0)*grad_detJ);
      

    rowvec3 p_xi_x_p_et = cross(pseudoInv.row(0),pseudoInv.row(1));

    double hess_detS[3][3][2];
    for(int i=0; i<3; i++){
      vec3 dir;
      dir.zeros();
      dir(i) = 1.0;
      vec3 temp = cross(dir,p_xi_x_p_et);
      for(int j=0; j<3; j++){
	hess_detS[i][j][0] = temp(j);
	hess_detS[i][j][1] = -temp(j);
      }
    }

    // Now, for the Hessians
    for(int i=0; i<2; i++){
      for(int j=0; j<3; j++){
	for(int k=0; k<=i; k++){
	  for(int l=0; l<=j; l++){
	    double hess_sigma = 0.5*(KroneckerDelta(i,k)*hess_detS[j][l][i]*
				     (1.0+detS/denom) -
				     grad_detS(i,j)*grad_detS(k,l)*
				     pow(detS,2)/pow(denom,3));
	    double hess_FrobNorm2 = 2*KroneckerDelta(i,k)*pseudoInv(i,j)*
	      pseudoInv(k,l);
	    double hess_eta =
	      0.5*( hess_FrobNorm2/sigma - 1.0/pow(sigma,2.0)*
		    grad_FrobNorm2(i,j)*grad_sigma(k,l) - 
		    (1.0/pow(sigma,2.0)*grad_sigma(i,j)*grad_FrobNorm2(k,l)-
		     2.0*FrobNorm2/pow(sigma,3.0)*grad_sigma(i,j)*
		     grad_sigma(k,l) + FrobNorm2/pow(sigma,2)* hess_sigma) );
									

	  }
	}

      }
    }

  }
  

  return 1.0;
}

int MTri::writeQuality(ofstream &out) const{

  //out << _distortion << " "  << 1.0/_distortion << endl;

}

*/
