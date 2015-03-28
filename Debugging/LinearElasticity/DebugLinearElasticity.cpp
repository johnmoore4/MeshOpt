#include "HModel.h"
#include "ElasticitySystem.h"
#include "MeshDistributor.h"


int HModel::DebugLinearElasticity(){
  using std::cout;
  using std::endl;

 
  /* 
  std::map<gind,arma::vec> Dirichlet;

 
  const double h = 0.3;

  const node_map& nodes = mesh_container.getNodes();
  const element_set& subelements = mesh_container.getSubElements();
 
  arma::vec3 np = {0.0, 0.0, 1.0};
  for(auto el = subelements.begin(); el != subelements.end(); ++el){
    const int ncn = (*el)->numCornerNodes();
    const gind* cn = (*el)->getCornerNodes();
    arma::vec3 r = nodes.at(cn[1])->getXYZ() - nodes.at(cn[0])->getXYZ();
    arma::vec3 n = cross(r,np);
    arma::vec2 temp;
    temp(0) = n(0);
    temp(1) = n(1);

    for(int i=0; i<2; i++){
      Dirichlet[cn[i]]+= -temp;
    }
  }


  
  for(auto bc = Dirichlet.begin(); bc != Dirichlet.end(); ++bc){
    const double nbc = norm(bc->second,2.0);
    if(nbc > 1.0e-10){
      bc->second/=(nbc/h);
    }
  }
  */

  MeshDistributor mesh_distributor;
  LocalMesh local_mesh = mesh_distributor.CreateDistributedMesh(mesh_container);

  //std::map<gind,arma::vec> Forces;
  //ElasticitySystem elasticity(mesh_container,sf_factory,Dirichlet,Forces);
  //elasticity.Solve();


  /*
  arma::wall_clock timer;
  const int N = 1;
  for(int i=0; i<N; i++){
    evaluator.SetUpProblem(0.03/N);

    timer.tic();
    //evaluator.Assemble(5.87e10,1.73e10);
    //evaluator.Assemble(1.0,1.0);
    cout << "Assemble time: " << timer.toc() << endl;

    timer.tic();
    evaluator.Solve();
    cout << "Solve time: " << timer.toc() << endl;

    //evaluator.UpdateMesh();
  }
  */

}
