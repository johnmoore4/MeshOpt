#include "HModel.h"
#include <iostream>
#include <memory>
#include <fstream>
//#include <omp.h>
#include "OCC_Handler.h"

//#include <Tpetra_DefaultPlatform.hpp>
//#include <Teuchos_GlobalMPISession.hpp>

using namespace std;

int main(int argc, char ** argv){

  //Teuchos::GlobalMPISession mpiSession (&argc, &argv);



  //unique_ptr<HModel> hmodel(new HModel);  
  HModel hmodel;

  //omp_set_num_threads(1);

  hmodel._optimal = false;
  //hmodel._optimal = true;

 


  try{
    //hmodel->parseCommandLineOptions(argc,argv);
    hmodel.parseCommandLineOptions(argc,argv);    
  }
  catch(std::exception& e){
    cout << "std::runtime error thrown: " << e.what() << endl;
    return 1;
  }

 
  //hmodel->readOCC();
  hmodel.ReadGeometry();
  cout << "after reading occ" << endl;
  try{
    //hmodel->readGMSH();
    hmodel.readMesh();
  }
  catch(std::exception& e){
    cout << "std::runtime error thrown: " << e.what() << endl;
    return 1;
  }
  cout << "after reading GMSH" << endl;
  
  hmodel.generateBL();

  hmodel.PrepareMesh();

  hmodel.MeshHighOrder();

    
  hmodel.WriteMesh("blmesh.msh","GMSH");
  
 


  return 0;
}


