#include <iostream>

#include "Util/IO/io.hpp"


#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>

#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/node.h"
#include "libmesh/point.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/unv_io.h"

int main (int argc, char ** argv)
{
    BeatIt::printBanner(std::cout);

    using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);
      libMesh::Mesh mesh(init.comm());
      std::string mesh_file = "heart.e";
      mesh.read(&mesh_file[0]);


      std::map<unsigned int, std::string> old_organ_map;
      old_organ_map[1] = "RA";
      old_organ_map[2] = "LA";
      old_organ_map[3] = "H";

      std::map<unsigned int, std::string> organ_map;
      organ_map[1] = "right_atria";
      organ_map[1] = "eustachian_valve";
      organ_map[1] = "right_atrial_appendage";
      organ_map[1] = "right_atrial_roof_muscle";
      organ_map[1] = "inferior_vena_cava";
      organ_map[1] = "superior_vena_cava";
      organ_map[1] = "inter_atrial_band";
      organ_map[1] = "aorta";
      organ_map[1] = "tricuspid_valve_ring";
      organ_map[1] = "mitral_valve_ring";
      organ_map[1] = "left_atrial_appendage";
      organ_map[1] = "pulmonary_veins";
      organ_map[1] = "pulmonary_artery";
      organ_map[2] = "left_atria";
      organ_map[3] = "left_ventricle";
      organ_map[3] = "right_ventricle";


      std::vector<std::ifstream *> csv_files;
      csv_files.push_back( new std::ifstream ("inter_atrial_band.csv") ); //0
      csv_files.push_back( new std::ifstream ("EV.csv") ); //1
      csv_files.push_back( new std::ifstream ("inferior_vena_cava.csv") ); //2
      csv_files.push_back( new std::ifstream ("superior_vena_cava.csv") ); //3
      csv_files.push_back( new std::ifstream ("right_atria_appendage.csv") );//4
      csv_files.push_back( new std::ifstream ("right_atria_roof_muscle.csv") );//5
      csv_files.push_back( new std::ifstream ("tricuspid_valve_ring.csv") );//6
      csv_files.push_back( new std::ifstream ("aorta.csv") );//7
      csv_files.push_back( new std::ifstream ("mitral_valve_ring.csv") );//8
      csv_files.push_back( new std::ifstream ("pulmonary_veins.csv") );//9
      csv_files.push_back( new std::ifstream ("right_ventricle.csv") );//10
      csv_files.push_back( new std::ifstream ("pulmonary_artery.csv") );//11
      csv_files.push_back( new std::ifstream ("left_atrial_appendage.csv") );//12
      csv_files.push_back( new std::ifstream ("cresta_terminalis.csv") );//12


      std::string line;

      int objID, elID, elID2;
      char v, v2;
      int c = 0;
      int n = 0;

      bool trick = false;
      for(int fn = 0; fn < csv_files.size(); fn++)
      {
          int ID = 4+fn;
          std::cout << "file: " << fn << std::endl;

          if(csv_files[fn]->is_open())
          {
              getline (*csv_files[fn],line);

                while ( getline (*csv_files[fn],line) )
                {
                    std::istringstream ss(line);
                    //line will have
                    ss >> objID >> v >> elID;
//                    if(fn == 10) elID = elID2;
                    //std::cout << elID << std::endl;
  //                  if(elID == 16896) trick = true;
    //                if(fn == 10 && trick == true) elID = objID;
                    auto * elem = mesh.query_elem(elID-1);
                    elem->subdomain_id() = ID;
                }
          }
      }

      std::cout << "Done" << std::endl;

      for(auto && file : csv_files) file->close();
      for(auto && file : csv_files) delete file;

      mesh.subdomain_name(1) = "right_atria";
      mesh.subdomain_name(2) = "left_atria";
      mesh.subdomain_name(3) = "left_ventricle";

      mesh.subdomain_name(4) = "inter_atrial_band";
      mesh.subdomain_name(5) = "eustachian_valve";

      mesh.subdomain_name(6) = "inferior_vena_cava";
      mesh.subdomain_name(7) = "superior_vena_cava";

      mesh.subdomain_name(8) = "right_atria_appendage";
      mesh.subdomain_name(9) = "right_atria_roof_muscle";

      mesh.subdomain_name(10) = "tricuspid_valve_ring";

      mesh.subdomain_name(11) = "aorta";

      mesh.subdomain_name(12) = "mitral_valve_ring";
      mesh.subdomain_name(13) = "pulmonary_veins";
      mesh.subdomain_name(14) = "right_ventricle";
      mesh.subdomain_name(15) = "pulmonary_artery";
      mesh.subdomain_name(16) = "left_atrial_appendage";
      mesh.subdomain_name(17) = "cresta_terminalis";
      std::cout << "Exporting" << std::endl;
      libMesh::ExodusII_IO(mesh).write("heart_with_regions.e");


    return 0;
}
