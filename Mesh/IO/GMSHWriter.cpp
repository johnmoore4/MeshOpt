#include "MeshWriter.h"
#include "NodeIndexer.h"
#include "ActiveMEl.h"
#include "OptEl.h"
#include "ShapeFunctionMatrices.h"

#include <iomanip>

void GMSHWriter::Save(std::string filename, std::string format){
  using namespace std;

  std::string name_plus_ext = filename+".msh";

  cout << "Writing Mesh" << endl;
  NodeIndexFactory index_factory;

  ofstream out(name_plus_ext.c_str());
  out.precision(16);

  out << "$MeshFormat" << endl;
  out << 2.2 << " " <<  0 << " " << 8 << endl;
  out << "$EndMeshFormat" << endl;

  node_map& nodes = mesh.getNodesNC();
  const element_set& elements = mesh.getElements();
  const element_set& subelements = mesh.getSubElements();
  const element_set& subsubelements = mesh.getSubSubElements();

 

  // Write the nodes
  out << "$ParametricNodes" << endl;
  out << nodes.size() << endl;
  
  for(auto nd=nodes.begin(); nd != nodes.end(); ++nd){
    arma::vec3 xyz = (*nd).second->getXYZ();
    out << (*nd).first + 1 << " ";
    for(int j=0; j<3; j++){
      out << std::setprecision(16) << xyz[j] << " ";
    }
    int type = nd->second->getType();
    //out << type << " " << type;
    out << type << " " << nd->second->getGeoEntity()+1;
    if(type != 3){
      for(int i = 0; i < type; i++){
	out << " " << nd->second->getParametricCoords()[i];
      }
    }
    out << endl;
  }
  out << "$EndParametricNodes" << endl;

  //const ideal_map& idealElements = mesh.getIdealElements();

  ShapeFunctionMatricesFactory sf_factory;


  //const double plot_threshold = 0.01;

  /*
  int qplotcnt = 0;
  int cnt = 0;
  std::vector<double> quality(elements.size());
  
  for(auto el=elements.begin(); el != elements.end(); ++el,cnt++){

    const ShapeFunctionMatrices* sf = 
	sf_factory.getShapeFunction((*el)->getElementType(),
				    (*el)->getOrder(),0);

    ActiveMEl activeEl((*el).get(),index_factory,nodes,sf);

    arma::mat ideal;
    
    auto it = idealElements.find((*el).get());
    if(it != idealElements.end()) ideal = it->second;
    
    double distortion;
    if((*el)->getDim() == 2){
      OptEl2D optel(activeEl,ideal);
      distortion = optel.computeDistortion();

    }
    else if((*el)->getDim() == 3){
      OptEl3D optel(activeEl,ideal);
      distortion = optel.computeDistortion();
    }
    quality[cnt] = 1.0/distortion;
    if(quality[cnt] < plot_threshold) qplotcnt++;
 
  }
  */




  int towrite=elements.size();
  //int towrite = qplotcnt;

  
  for(auto el=subelements.begin(); el!= subelements.end(); ++el){
    if((*el)->getBCTag() != 0){
      //if((*el)->hasGeoEntity()){
      towrite++;
    }
  }
  
  /*
  for(auto el=subsubelements.begin(); el!= subsubelements.end(); ++el){
    if((*el)->hasGeoEntity()) towrite++;
  }
  */
 
  //ShapeFunctionMatricesFactory sf_factory;
  cout << "before writing elements" << endl;

  out << "$Elements" << endl;
  out << towrite << endl;
  {
    int cnt = 0;
    for(auto el=elements.begin(); el != elements.end(); ++el,cnt++){
      ActiveMEl activeEl((*el).get(),index_factory,nodes,NULL);
      activeEl.writeGMSH(out,cnt);
      
    }
    
    for(auto el=subelements.begin(); el!= subelements.end(); ++el,cnt++){
      if((*el)->getBCTag() != 0){
	//if((*el)->hasGeoEntity()){
	ActiveMEl activeEl((*el).get(),index_factory,nodes,NULL);
	activeEl.writeGMSH(out,cnt);
      }
    }
    /*
    for(auto el=subsubelements.begin(); el!= subsubelements.end(); ++el,cnt++){
      if((*el)->hasGeoEntity()){
	ActiveMEl activeEl((*el).get(),index_factory,nodes,NULL);
	activeEl.writeGMSH(out,cnt);
      }
    }
    */
  }
  out << "$EndElements" << endl;
  /*
  cout << "Before writing merits" << endl;

  const std::map<MEl*,double>& all_merits = mesh.getAllMerits();

  const double plot_threshold = 0.9;
  int qplotcnt=0;
 
  for(auto el = all_merits.begin(); el != all_merits.end(); ++el){
    //std::cout << "merit: " << el->second << endl;
    if(el->second > plot_threshold) qplotcnt++;
  }
  //const int qplotcnt = all_merits.size();
 
  out << "$ElementData" << endl;
  out << 1 << endl;
  out << "\"quality measure\"" << endl;
  out << 1 << endl;
  out << 0.0 << endl;
  out << 4 << endl;
  out << 0 << endl;
  out << 1 << endl;
  out << qplotcnt << endl;
  out << 0 << endl;
 
  {
    int cnt=0;
    for(auto el = elements.begin(); el != elements.end(); ++el, ++cnt){
      auto it = all_merits.find(el->get());
      if(it != all_merits.end() && it->second > plot_threshold){
	out << cnt + 1 << " " << 1.0/it->second << endl;
      }
    }
    for(auto el = subelements.begin(); el != subelements.end(); ++el, ++cnt){
      auto it = all_merits.find(el->get());
      if(it != all_merits.end() && it->second > plot_threshold){
	out << cnt + 1 << " " << 1.0/it->second << endl;
      }
    }

  }

  out << "$EndElementData" << endl;
  */

  out.close();
}
