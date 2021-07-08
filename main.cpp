#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/IO/Color.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include <iostream>
#include <fstream>
#include <array>

namespace OpenMesh {
	typedef TriMesh_ArrayKernelT<>	MyMesh;
}

namespace CGAL {
	typedef CGAL::Exact_predicates_exact_constructions_kernel	Kernel;
	typedef Arr_segment_traits_2<Kernel>						Traits;
	typedef Arrangement_2<Traits> 								Arrangement;
	typedef Traits::X_monotone_curve_2							Segment;
	typedef Point_set_3<Kernel::Point_3>						Point_set;
	typedef Point_set::Property_map<Color>						Color_map;
}

int main(int argc, char *argv[]) {
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " Mesh_File PC_File" << std::endl;
		return 1;
	}

	std::cout << "Convert " << argv[1] << " into " << argv[2] << std::endl;

	// Read Mesh

	OpenMesh::MyMesh mesh;
	
	mesh.request_halfedge_texcoords2D();
	mesh.request_face_texture_index();

	OpenMesh::IO::Options ropt;
	ropt += OpenMesh::IO::Options::FaceTexCoord;
	ropt += OpenMesh::IO::Options::Custom;

	if ( ! OpenMesh::IO::read_mesh(mesh, argv[1], ropt)) {
		std::cerr << "Error loading mesh from file " << argv[1] << std::endl;
		return 1;
	}

	OpenMesh::MPropHandleT< std::map< int, std::string > > textures_files;
	mesh.get_property_handle(textures_files, "TextureMapping");

	auto textures_index = mesh.face_texture_index_pph();

	OpenMesh::FPropHandleT<int> label;
	mesh.get_property_handle(label, "label");

	int num_textures = mesh.property(textures_files).size();

	// Arrangement on textures
	CGAL::Arrangement arrangements[num_textures];

	for (int i = 0; i < num_textures; i++) {
		arrangements[i] = CGAL::Arrangement();
	}

	for (auto edge: mesh.edges()) {
		auto tex0 = mesh.property(textures_index, edge.h0().face());
		auto from0 = mesh.texcoord2D(mesh.prev_halfedge_handle(edge.h0()));
		auto to0 = mesh.texcoord2D(edge.h0());
		auto tex1 = mesh.property(textures_index, edge.h1().face());
		auto from1 = mesh.texcoord2D(mesh.prev_halfedge_handle(edge.h1()));
		auto to1 = mesh.texcoord2D(edge.h1());

		bool one_edge = true;
		if (tex0 != tex1 || from0[0] != to1[0] || from0[1] != to1[1] || from1[0] != to0[0] || from1[1] != to0[1]) {
			one_edge = false;
		}

		if (one_edge) {
			insert(arrangements[tex0], CGAL::Segment(CGAL::Traits::Point_2(from0[0], from0[1]), CGAL::Traits::Point_2(to0[0], to0[1])));
		} else {
			insert(arrangements[tex0], CGAL::Segment(CGAL::Traits::Point_2(from0[0], from0[1]), CGAL::Traits::Point_2(to0[0], to0[1])));
			insert(arrangements[tex1], CGAL::Segment(CGAL::Traits::Point_2(from1[0], from1[1]), CGAL::Traits::Point_2(to1[0], to1[1])));
		}
	}

	// Write Point Cloud

	CGAL::Point_set point_set;
	CGAL::Color_map color;

	boost::tie (color, boost::tuples::ignore) = point_set.add_property_map<CGAL::Color>("color", CGAL::Color(0, 0, 0));

	for (int i = 0; i < num_textures; i++) {
		int width, height, n;
		unsigned char *data = stbi_load((std::string(argv[1]).substr(0, std::string(argv[1]).find_last_of('/')+1)+mesh.property(textures_files)[i]).c_str(), &width, &height, &n, 3);
		
		for (int v = 0; v < height; v++) {
			for (int u = 0; u < width; u++) {
				auto point = point_set.insert (CGAL::Kernel::Point_3(0., 0., 0.));
				color[*point] = CGAL::Color(data[3*(u + v*width) + 0], data[3*(u + v*width) + 1], data[3*(u + v*width) + 2]);
			}
		}
	}

	// Output file

	std::ofstream out_file (argv[2]);
	CGAL::write_ply_points_with_properties(out_file, point_set,
											CGAL::make_ply_point_writer (point_set.point_map()),
											std::make_tuple(color,
																CGAL::PLY_property<unsigned char>("red"),
																CGAL::PLY_property<unsigned char>("green"),
																CGAL::PLY_property<unsigned char>("blue")));
		
	return 0;
}
