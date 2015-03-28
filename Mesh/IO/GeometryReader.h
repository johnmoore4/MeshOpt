#pragma once
#include <string>
#include <map>

class GeometryContainer;

class GeometryReader{
 private:

 protected:

 public:
  void ReadGeometry(GeometryContainer& geom, std::string filename);

};

/*
class OCCReader: pubic GeometryReader{
 private:

 protected:

 public:
  void ReadGeometry(GeometryContainer& geom, std::string filename){}
}
*/
