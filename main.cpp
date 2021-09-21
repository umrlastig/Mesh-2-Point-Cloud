#define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 1

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/property_map.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

namespace OpenMesh {
	typedef TriMesh_ArrayKernelT<>	MyMesh;
}

namespace CGAL {
	typedef Simple_cartesian<double>							Point_set_kernel;
	typedef Point_set_3<Point_set_kernel::Point_3>				Point_set;

	typedef Exact_predicates_exact_constructions_kernel			Polygon_kernel;
	typedef CGAL::Polygon_2<Polygon_kernel> 					Polygon;
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
	mesh.request_vertex_normals();

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

	std::cout << "Mesh loaded" << std::endl;

	// Compute normal
	for (auto vertex: mesh.vertices()) {
		mesh.set_normal(vertex, OpenMesh::MyMesh::Normal(0));
	}

	std::map<OpenMesh::FaceHandle, bool> flat_faces;
	for (auto face: mesh.faces()) {
		std::vector<OpenMesh::MyMesh::Point> points;
		points.reserve(3);

		for (auto vertex: face.vertices()) {
			points.push_back(mesh.point(vertex));
		}
		
		auto normal = OpenMesh::cross(points[1] - points[0], points[2] - points[0]);

		if (normal.length() != 0) {
			flat_faces[face] = false;
			for (auto vertex: face.vertices()) {
				mesh.set_normal(vertex, mesh.normal(vertex) + normal);
			}
		} else {
			flat_faces[face] = true;
		}
	}

	for (auto vertex: mesh.vertices()) {
		auto normal = mesh.normal(vertex);
		if (normal.length() == 0) {
			for (auto v: vertex.vertices_ccw()) {
				auto edge = mesh.point(vertex) - mesh.point(v);
				if (edge.length() != 0) {
					normal += edge / edge.length();
				}
			}
		}
		mesh.set_normal(vertex, normal / normal.length());
	}

	std::cout << "Normal computed" << std::endl;

	//Open textures
	int num_textures = mesh.property(textures_files).size();

	int width[num_textures];
	int height[num_textures];
	int n[num_textures];
	unsigned char *data[num_textures];

	for (int i = 0; i < num_textures; i++) {
		data[i] = stbi_load((std::string(argv[1]).substr(0, std::string(argv[1]).find_last_of('/')+1)+mesh.property(textures_files)[i]).c_str(), &width[i], &height[i], &n[i], 3);
	}

	// Write Point Cloud
	CGAL::Point_set point_set;
	point_set.add_normal_map();

	CGAL::Point_set::Property_map<unsigned char> ps_red;
	CGAL::Point_set::Property_map<unsigned char> ps_green;
	CGAL::Point_set::Property_map<unsigned char> ps_blue;
	boost::tie (ps_red, boost::tuples::ignore) = point_set.add_property_map<unsigned char>("red", 0);
	boost::tie (ps_green, boost::tuples::ignore) = point_set.add_property_map<unsigned char>("green", 0);
	boost::tie (ps_blue, boost::tuples::ignore) = point_set.add_property_map<unsigned char>("blue", 0);
	
	CGAL::Point_set::Property_map<int> ps_label;
	boost::tie (ps_label, boost::tuples::ignore) = point_set.add_property_map<int>("label", 0);

	CGAL::Point_set::Property_map<int> ps_u;
	CGAL::Point_set::Property_map<int> ps_v;
	boost::tie (ps_u, boost::tuples::ignore) = point_set.add_property_map<int>("u", 0);
	boost::tie (ps_v, boost::tuples::ignore) = point_set.add_property_map<int>("v", 0);

	for (auto face: mesh.faces()) {
		if (!flat_faces[face]) {
			CGAL::Polygon polygon;
			CGAL::Point_set_kernel::Vector_3 points[3];
			CGAL::Point_set_kernel::Vector_3 normals[3];
			auto arrangement_id = mesh.property(textures_index, face);

			std::size_t i = 0;
			for (auto hh: face.halfedges()) {
				auto texCoord = mesh.texcoord2D(hh);
				polygon.push_back(CGAL::Polygon_kernel::Point_2(texCoord[0]*width[arrangement_id], texCoord[1]*height[arrangement_id]));
				auto point = mesh.point(hh.to());
				points[i] = CGAL::Point_set_kernel::Vector_3(point[0], point[1], point[2]);
				auto normal = mesh.normal(hh.to());
				normals[i] = CGAL::Point_set_kernel::Vector_3(normal[0], normal[1], normal[2]);
				i++;
			}

			if (polygon.is_simple()) {

				auto box = polygon.bbox();

				for (int u = ((int) box.xmin()); u < ((int) box.xmax()) + 1; u++) {
					for (int v = ((int) box.ymin()); v < ((int) box.ymax()) + 1; v++) {

						auto uv_point = CGAL::Polygon_kernel::Point_2(0.5 + u, 0.5 + v);

						if (!polygon.has_on_unbounded_side(uv_point)) {

							auto lambda = CGAL::Barycentric_coordinates::compute_triangle_coordinates_2(polygon[0], polygon[1], polygon[2], uv_point, CGAL::Polygon_kernel());

							CGAL::Point_set_kernel::Point_3 point = CGAL::Point_set_kernel::Point_3(0,0,0) + CGAL::to_double(lambda[0]) * points[0] + CGAL::to_double(lambda[1]) * points[1] + CGAL::to_double(lambda[2]) * points[2];
							CGAL::Point_set_kernel::Vector_3 normal = CGAL::to_double(lambda[0]) * normals[0] + CGAL::to_double(lambda[1]) * normals[1] + CGAL::to_double(lambda[2]) * normals[2];

							auto point_xyz = point_set.insert(point, normal);
							ps_label[*point_xyz] = mesh.property(label, face);
							ps_u[*point_xyz] = u;
							ps_v[*point_xyz] = v;

							ps_red[*point_xyz] = data[arrangement_id][3*(u + (height[arrangement_id] - v - 1)*width[arrangement_id]) + 0];
							ps_green[*point_xyz] = data[arrangement_id][3*(u + (height[arrangement_id] - v - 1)*width[arrangement_id]) + 1];
							ps_blue[*point_xyz] = data[arrangement_id][3*(u + (height[arrangement_id] - v - 1)*width[arrangement_id]) + 2];

						}
					}
				}
			}
		}
	}

	std::cout << "Point cloud computed" << std::endl;

	// Output file
	std::ofstream out_file (argv[2]);
	out_file << point_set;
		
	return 0;
}
