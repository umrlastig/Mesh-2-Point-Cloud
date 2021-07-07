#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/Color.h>

#include <iostream>
#include <fstream>
#include <array>

namespace OpenMesh {
	typedef TriMesh_ArrayKernelT<>	MyMesh;
}

namespace CGAL {
	typedef Simple_cartesian<double>			Kernel;
	typedef Point_set_3<Kernel::Point_3>		Point_set;
	typedef Point_set::Property_map<Color>	Color_map;
}

int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " Mesh_File PC_File" << std::endl;
		return 1;
	}

	std::cout << "Convert " << argv[1] << " into " << argv[2] << std::endl;

	OpenMesh::MyMesh mesh;
	
	mesh.request_halfedge_texcoords2D();
	mesh.request_face_texture_index();

	OpenMesh::IO::Options ropt;
	ropt += OpenMesh::IO::Options::FaceTexCoord;

	if ( ! OpenMesh::IO::read_mesh(mesh, argv[1], ropt)) {
		std::cerr << "Error loading mesh from file " << argv[1] << std::endl;
		return 1;
	}

	OpenMesh::MPropHandleT< std::map< int, std::string > > property;
	mesh.get_property_handle(property, "TextureMapping");
	auto texture_files = mesh.property(property);

	CGAL::Point_set point_set;
	CGAL::Color_map color;

	boost::tie (color, boost::tuples::ignore) = point_set.add_property_map<CGAL::Color>("color", CGAL::Color(0, 0, 0));
	auto point = point_set.insert (CGAL::Kernel::Point_3(0., 0., 0.));
	color[*point] = CGAL::Color(255, 255, 255);

	std::ofstream out_file (argv[2]);
	CGAL::write_ply_points_with_properties(out_file, point_set,
										CGAL::make_ply_point_writer (point_set.point_map()),
										std::make_tuple(color,
															CGAL::PLY_property<unsigned char>("red"),
															CGAL::PLY_property<unsigned char>("green"),
															CGAL::PLY_property<unsigned char>("blue")));
		
	return 0;
}
