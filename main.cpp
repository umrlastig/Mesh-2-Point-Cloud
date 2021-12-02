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

#include <CGAL/grid_simplify_point_set.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include <getopt.h>
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

	int opt;
	const struct option options[] = {
		{"help", no_argument, NULL, 'h'},
		{"gridSubSampling", required_argument, NULL, 'g'},
		{"poissonDiskSampling", required_argument, NULL, 'p'},
		{NULL, 0, 0, '\0'}
	};

	double gridSubSampling = 0; //size of the grid use for grid sub-sampling
	double poissonDiskSampling = 0; //min distance for poisson disk sampling

	while ((opt = getopt_long(argc, argv, "hg:p:", options, NULL)) != -1) {
		switch(opt) {
			case 'h':
				std::cout << "Usage: " << argv[0] << " [OPTION] MESH_FILE PC_FILE" << std::endl;
				std::cout << "Convert the mesh from MESH_FILE into a point cloud in PC_FILE (PLY file)." << std::endl;
				std::cout << "By default, one point is sampled per texel. This behavior can be changed by using the -p option." << std::endl << std::endl;
				std::cout << "OPTIONS:" << std::endl;
				std::cout << " -h, --help                   Print this help anq quit." << std::endl;
				std::cout << " -g, --gridSubSampling=G      Subsample the final point cloud, keeping one point per voxel of size G." << std::endl;
				std::cout << " -p, --poissonDiskSampling=P  Use poisson isk sampling with parameter P." << std::endl << std::endl;
				return EXIT_SUCCESS;
				break;
			case 'g':
				gridSubSampling = strtod(optarg, NULL);
				break;
			case 'p':
				poissonDiskSampling = strtod(optarg, NULL);
				break;
		}
	}

	if (argc != optind + 2) {
		std::cerr << "MESH_FILE and PC_FILE are mandatory" << std::endl;
		std::cerr << "Usage: " << argv[0] << " [OPTION] MESH_FILE PC_FILE" << std::endl;
		std::cerr << "Use -h/--help to obtain more informations" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Convert " << argv[optind] << " into " << argv[optind+1] << std::endl;

	// Read Mesh
	OpenMesh::MyMesh mesh;
	
	mesh.request_halfedge_texcoords2D();
	mesh.request_face_texture_index();
	mesh.request_vertex_normals();
	mesh.request_face_normals();

	OpenMesh::IO::Options ropt;
	ropt += OpenMesh::IO::Options::FaceTexCoord;
	ropt += OpenMesh::IO::Options::Custom;

	if ( ! OpenMesh::IO::read_mesh(mesh, argv[optind], ropt)) {
		std::cerr << "Error loading mesh from file " << argv[optind] << std::endl;
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

	for (auto face: mesh.faces()) {
		std::vector<OpenMesh::MyMesh::Point> points;
		points.reserve(3);

		for (auto vertex: face.vertices()) {
			points.push_back(mesh.point(vertex));
		}
		
		auto normal = OpenMesh::cross(points[1] - points[0], points[2] - points[0]);

		if (normal.length() != 0) {
			mesh.set_normal(face, normal / normal.length());
			for (auto vertex: face.vertices()) {
				mesh.set_normal(vertex, mesh.normal(vertex) + normal);
			}
		} else {
			mesh.set_normal(face, normal);
		}
	}

	for (auto vertex: mesh.vertices()) {
		auto normal = mesh.normal(vertex);
		if (normal.length() == 0 and vertex.valence() > 0) {
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
		data[i] = stbi_load((std::string(argv[optind]).substr(0, std::string(argv[optind]).find_last_of('/')+1)+mesh.property(textures_files)[i]).c_str(), &width[i], &height[i], &n[i], 3);
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

	CGAL::Point_set::Property_map<int> face_id;
	boost::tie (face_id, boost::tuples::ignore) = point_set.add_property_map<int>("face_id", -1);

	CGAL::Point_set::Property_map<double> nx_0;
	CGAL::Point_set::Property_map<double> ny_0;
	CGAL::Point_set::Property_map<double> nz_0;
	boost::tie (nx_0, boost::tuples::ignore) = point_set.add_property_map<double>("fnx", 0);
	boost::tie (ny_0, boost::tuples::ignore) = point_set.add_property_map<double>("fny", 0);
	boost::tie (nz_0, boost::tuples::ignore) = point_set.add_property_map<double>("fnz", 0);

	for (auto face: mesh.faces()) {
		if (mesh.normal(face).length() > 0) {
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
							auto normal_0 = mesh.normal(face);

							auto point_xyz = point_set.insert(point, normal);
							ps_label[*point_xyz] = mesh.property(label, face);
							ps_u[*point_xyz] = u;
							ps_v[*point_xyz] = v;
							face_id[*point_xyz] = face.idx();

							nx_0[*point_xyz] = normal_0[0];
							ny_0[*point_xyz] = normal_0[1];
							nz_0[*point_xyz] = normal_0[2];

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

	// Simplify point set
	if (gridSubSampling > 0) {
		point_set.remove(CGAL::grid_simplify_point_set(point_set, gridSubSampling), point_set.end());
	}

	// Output file
	std::ofstream out_file (argv[optind+1]);
	out_file << point_set;
		
	return 0;
}
