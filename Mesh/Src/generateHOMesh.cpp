#include "HModel.h"
//#include "shapeFunctionHandler.h"
#include "MeshConnectivity.h"
#include "BoundaryLayer.h"

using namespace std;

using arma::wall_clock;

int HModel::MeshHighOrder(int porder){
  /*
  // Generate Boundary Layer
  BoundaryLayerGenerator bl_generator(mesh_container,geometry_container,
				      sf_factory);
  bl_generator.GenerateBL(0.05, 1);



  // Generate Connectivity


  bl_generator.OptimizeBL();
  cout << "after optimize BL" << endl;

  */

  ComputeMeshConnectivity(mesh_container);
  cout << "after generate mesh connectivity" << endl;

  // Generate HO mesh

  // Optimize HO mesh

}


int HModel::generateHOMesh(int porder){

  //typedef std::vector<element3D_ptr>::iterator el3dit;
  //typedef std::vector<element2D_ptr>::iterator el2dit;

  //cout << "porder: " << porder << endl;
  //cout << "_porder: " << _porder << endl;


  // Now, all elements have primary sf of order porder
  //sf_handler->initialize(porder,_porder,_optimal,0,mesh_dim);
  //wall_clock timer;
  //timer.tic();

  /*

  cout << "before 1D set order" << endl;
  // Mesh 1D interior nodes
  for(auto el1d=elements1D.begin(); el1d != elements1D.end(); ++el1d){
    (*el1d)->setOrder(porder);
  }

  // Associate 1D child nodes with 2D element
  //for(auto el2d=elements2D.begin(); el2d != elements2D.end(); ++el2d){
  //  (*el2d)->getChildNodes();
  //}  

  cout << "before 2D set order" << endl;
  // Mesh 2D interior nodes
  for(auto el2d=elements2D.begin(); el2d != elements2D.end(); ++el2d){
    (*el2d)->setOrder(porder);
  }

  // Associate 2D child nodes with 3D element
  //for(auto el3d=elements3D.begin(); el3d != elements3D.end(); ++el3d){
  //  (*el3d)->getChildNodes();
  //} 

  // Mesh 3D interior nodes
  for(auto el3d=elements3D.begin(); el3d != elements3D.end(); ++el3d){
    (*el3d)->setOrder(porder);
  }
  */
  /*
  cout << "before main interior nodes meshing..." << endl;
  for(auto el=mesh_elements.begin(); el != mesh_elements.end(); ++el){
    (*el)->setOrder(porder);
  }

  for(int dm = mesh_dim-1; dm > 0; dm--){
    cout << "before dim: " << dm << " interior nodes meshing..." << endl;
    for(auto el=subelements.begin(); el != subelements.end(); ++el){
      if((*el)->dim() == dm) (*el)->setOrder(porder);
    }
  }
  cout << "after all meshing" << endl;
  
  if(mesh_dim == 3){
    for(auto el=subelements.begin(); el != subelements.end(); ++el){
      if((*el)->dim() == 2) (*el)->getChildNodes();
    }
  }
  
  for(auto el=mesh_elements.begin(); el != mesh_elements.end(); ++el){
    (*el)->getChildNodes();
  }
  */

  /*
  // Mesh 3D interior nodes
  for(auto el3d=elements3D.begin(); el3d != elements3D.end(); ++el3d){
    (*el3d)->setOrder(porder);
  }
  
  cout << "before meshing 2d interior nodes" << endl;
  // Mesh 2D interior nodes
  for(auto el2d=elements2D.begin(); el2d != elements2D.end(); ++el2d){
    (*el2d)->setOrder(porder);
  }

  cout << "after meshing 2d interior nodes" << endl;

  // Mesh 1D interior nodes
  for(auto el1d=elements1D.begin(); el1d != elements1D.end(); ++el1d){
    (*el1d)->setOrder(porder);
  }
  cout << "after meshing 1d interior nodes" << endl;

  // Associate 1D child nodes with 2D element
  for(auto el2d=elements2D.begin(); el2d != elements2D.end(); ++el2d){
    (*el2d)->getChildNodes();
  }  

  // Associate 2D child nodes with 3D element
  for(auto el3d=elements3D.begin(); el3d != elements3D.end(); ++el3d){
    (*el3d)->getChildNodes();
   } 
  */

  //_porder = porder;

  //cout << "generateHOMesh time: " << timer.toc() << endl;

  return 1;
}
