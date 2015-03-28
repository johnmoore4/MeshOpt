#include <iostream>
#include <string>
#include <fstream>
#include <boost/program_options.hpp>
#include "HModel.h"
#include "MeshReader.h"
#include "GeometryReader.h"
#include "MeshWriter.h"
#include "BoundaryLayer.h"
#include "MeshConnectivity.h"
//#include "HighOrderMeshGenerator.h"
#include "HighOrderMeshGeneratorAndOptimizer.h"
#include "MeshOptConfig.h"
#include "CaseParameterList.h"
#include "HighOrderParameterList.h"

//#include "OCC_Handler.h"
#include "Mnode.h"
#include "Mel.h"
#include "assert.h"

using namespace std;
namespace po = boost::program_options;

HModel::HModel(void){
  
  //mesh_container = std::unique_ptr<MeshContainer>(new MeshContainer);
  mesh_reader = std::unique_ptr<GMSHReader>(new GMSHReader);
  //occ_handler = std::unique_ptr<OCC_Handler>(new OCC_Handler);
  geometry_reader = std::unique_ptr<GeometryReader>(new GeometryReader);
  mesh_opt_config = std::make_shared<MeshOptConfig>();
  mesh_opt_config->Read();
}

int HModel::readMesh(void){
  //std::cout << "gmsh file name: " << mesh_opt_config->getCaseParameters()->gmshFileName << std::endl;
  mesh_reader->ReadMesh(mesh_container,
			mesh_opt_config->getCaseParameters()->gmshFileName);

  mesh_reader->SetGeoTags(geometry_container);
  fixFaceOrientations(mesh_container);

  return 1;
}

HModel::~HModel(){}

int HModel::WriteMesh(std::string filename, std::string format){
  if(format=="GMSH"){
    GMSHWriter writer(mesh_container);
    writer.Save(filename,"ascii");
  }
  else{
    assert("Only GMSH ouptut supported!");
  }
  return 1;
}

int HModel::ReadGeometry(void){
  geometry_reader->ReadGeometry(geometry_container,
				mesh_opt_config->getCaseParameters()->
				stepFileName);

}
int HModel::parseCommandLineOptions(int argc, char ** argv){
  /*
    try {

        po::options_description desc("Allowed options");
        desc.add_options()
	  ("help", "produce help message")
	  ("mesh-file", po::value< string >(), "mesh name")
	  ("geo-file", po::value< string >(), "gometry name")
        ;

	po::positional_options_description p;
	p.add("mesh-file", 1);
	p.add("geo-file", 1);

        po::variables_map vm;        
	po::store(po::command_line_parser(argc, argv).
          options(desc).positional(p).run(), vm);
        //po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

	if(vm.count("mesh-file")){
	  string filename = vm["mesh-file"].as< string >();
	  int size = filename.size();
	  string ext = filename.substr(size-4,4);
	  if(strcmp(ext.c_str(),".msh")){
	    cout << "\"" << ext << "\" is not a valid mesh name extension!\n";
	    throw 10;
	  }
	  else{
	    cout << "mesh-file is: " << filename << endl;
	    gmsh_file_name = filename;
	  }
	}
	else{
	  cout << "mesh-file has not been set!" << endl;
	  throw 10;
	}


	if(vm.count("geo-file")){
	  string filename = vm["geo-file"].as< string >();
	  int size = filename.size();
	  string ext = filename.substr(size-4,4);
	  if(strcmp(ext.c_str(),".stp")){
	    cout << "\"" <<  ext << "\" is not a valid mesh name extension!\n";
	    throw 10;
	  }
	  else{
	    cout << "geo-file is: " << filename << endl;
	    step_file_name = filename;
	  }
	}
	else{
	  cout << "geo-file has not been set!" << endl;
	  throw 10;
	}


    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception thrown in reading command-line input!\n";
    }

  */
    return 1;
}

int HModel::PrepareMesh(){
    ComputeMeshConnectivity(mesh_container);
}

int HModel::MeshHighOrder(){
  HighOrderMeshGeneratorAndOptimizer 
    ho_mesh_generator(mesh_container,geometry_container);
  
  int porder = mesh_opt_config->getHighOrderParameters()->order;

  high_order_mesh = ho_mesh_generator.GenerateAndOptimize(porder);
  bl_generator->SubDivideBL();
  
  //HighOrderMeshGenerator ho_mesh_generator;
  //high_order_mesh = ho_mesh_generator.generateHighOrderMeshRecursive(porder);
}

int HModel::generateBL(){

  bl_generator =
    std::make_shared<BoundaryLayerGenerator>(mesh_container,
					     geometry_container,
					     sf_factory,
					     mesh_opt_config->
					     getBLParameters());
	
  /*				     
  BoundaryLayerGenerator bl_generator(mesh_container,geometry_container,
				      sf_factory,
				      mesh_opt_config->getBLParameters());
  */

  bl_generator->GenerateBL();
  std::cout << "After generate BL" << std::endl;

  bl_generator->OptimizeBL();


}


