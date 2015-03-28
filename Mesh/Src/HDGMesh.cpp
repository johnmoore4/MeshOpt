#include "HDGMesh.h"
#include <fstream>
//#include "StaticTeuchosComm.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

int HDGMesh::Write(std::string format){
  //int rank = TrilinosExtern::Comm->getRank();
  int rank = 1;
  std::ofstream ofs("mesh_part"+std::to_string(rank));
  boost::archive::binary_oarchive oa(ofs);
  oa << *this;
}

int HDGMesh::Read(std::string format){
  int rank =1 ;
  //int rank = TrilinosExtern::Comm->getRank();
  std::ifstream ifs("mesh_part"+std::to_string(rank));
  boost::archive::binary_iarchive ia(ifs);
  ia >> *this;

}
