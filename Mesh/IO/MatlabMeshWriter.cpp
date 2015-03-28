#include "MatlabMeshWriter.h"
#include "MeshContainer.h"
#include <fstream>
#include <iomanip>

void MatlabMeshWriter::Write(std::string filename){
  const element_set& subelements = mesh.getSubElements();
  std::map<const MEl*,gind> SEptr2Indx;
  {
    int cnt = 1;
    for(auto it = subelements.begin(); it != subelements.end(); ++it,++cnt){
      SEptr2Indx[it->get()] = cnt;
    }
  }

  std::ofstream out, outf, outdg, outp, outt;
  out.open("t2f");
  outf.open("f");
  outdg.open("dgnodes");
  outp.open("p");
  outt.open("t");

  const element_set& elements = mesh.getElements();
  const node_map& nodes = mesh.getNodes();

  {
    int cnt = 0;
    for(auto el = elements.begin(); el != elements.end(); ++el){
      const MEl* curr = el->get();
      const int nc = curr->NumChildren();
      for(int c = 0; c < nc; c++){
	const gind child = SEptr2Indx[curr->getChild(c)];
	int orient = curr->getChildOrientation(c);
	out << child*orient << " ";
      }
      out << std::endl;

      const gind* ni = curr->getNodes();
      for(int n = 0; n < curr->NumNodes(); n++){
	const MNode* currnd = nodes.at(ni[n]).get();
	const arma::vec3& xyz = currnd->xyzvec3();
	//outdg << xyz(0) << " " << xyz(1) << std::endl;
	outt << ni[n]+1 << " ";
      }
      outt << std::endl;

      for(int d = 0; d < 2; d++){
	for(int n = 0; n < curr->NumNodes(); n++){
	  const MNode* currnd = nodes.at(ni[n]).get();
	  const arma::vec3& xyz = currnd->xyzvec3();
	  outdg << std::setprecision(16) << xyz(d) <<  std::endl;
	}
      }

    }
  }

  {
    for(auto el = subelements.begin(); el != subelements.end(); ++el){
      const MEl* curr = el->get();
      int cn[8];
      curr->getCornerNodes(cn);
      const gind* nds =  curr->getNodes();
      for(int i = 0; i < curr->numCornerNodes(); i++){
	outf << nds[cn[i]]+1 << " ";
      }
      outf << "0 ";
      int tag = curr->getBCTag();
      if(tag == 0) outf << 1 << std::endl;
      else outf << -tag << std::endl;
      //outf << -curr->getBCTag() << " ";
      //outf << std::endl;
    }
  }
  for(auto nd = nodes.begin(); nd != nodes.end(); ++nd){
    const arma::vec3& xyz = nd->second->xyzvec3();
    outp << std::setprecision(16) << xyz(0) << " " << xyz(1) << std::endl;
    //std::cout << nd->first << std::endl;
  }
  out.close();
  outf.close();
  outdg.close();
  outp.close();
  outt.close();
}
