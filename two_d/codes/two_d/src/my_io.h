#ifndef IO_HPP
#define IO_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Eigen>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace my_io{
  
  class pyWriter{
    public:
      pyWriter() {}
      pyWriter(std::shared_ptr<const lf::mesh::Mesh>); 
      void writePointData(std::string, Eigen::VectorXd&); 
    private:
      std::shared_ptr<const lf::mesh::Mesh> mesh_p_; 
  };

  pyWriter::pyWriter(std::shared_ptr<const lf::mesh::Mesh> mesh_p):
    mesh_p_(mesh_p) {}; 

  void pyWriter::writePointData(std::string filename, Eigen::VectorXd& datav){
    std::ofstream file(filename); 
    for (const lf::mesh::Entity &ent: mesh_p_->Entities(2)){ 
      const lf::geometry::Geometry* geo_p = ent.Geometry(); 
      Eigen::MatrixXd node = lf::geometry::Corners(*geo_p);
      file << node.transpose() <<  "\t" << 
			datav[mesh_p_->Index(ent)] << std::endl;  
    }
    file.close();
  }
} // namespace closure
#endif
