#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>

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

namespace OpenMesh {
	typedef TriMesh_ArrayKernelT<>	MyMesh;
}

namespace CGAL {
	typedef Simple_cartesian<double>							Point_set_kernel;
	typedef Point_set_3<Point_set_kernel::Point_3>				Point_set;

	typedef Exact_predicates_exact_constructions_kernel			Arrangement_kernel;
	typedef Arr_segment_traits_2<Arrangement_kernel>			Traits;
	typedef Traits::X_monotone_curve_2							Segment;
	struct VertexData {
		Point_set_kernel::Point_3 point;
		Point_set_kernel::Vector_3 normal;
		VertexData () {};
		VertexData (Point_set_kernel::Point_3 point_, Point_set_kernel::Vector_3 normal_): point(point_), normal(normal_) {};
		VertexData (OpenMesh::DefaultTraits::Point point_, OpenMesh::DefaultTraits::Normal normal_){
			point = Point_set_kernel::Point_3(point_[0], point_[1], point_[2]);
			normal = Point_set_kernel::Vector_3(normal_[0], normal_[1], normal_[2]);
		};
	};
	struct HalfedgeData {
		int label;
		HalfedgeData (): label(-1) {};
		HalfedgeData (int label_): label(label_) {};
	};
	typedef Arr_dcel_base< Arr_extended_vertex< Arr_vertex_base< Traits::Point_2 >, VertexData>, Arr_extended_halfedge< Arr_halfedge_base< Segment>, HalfedgeData>, Arr_face_base >	Dcel;
	typedef Arrangement_2<Traits, Dcel>							Arrangement;
	typedef Arr_trapezoid_ric_point_location<Arrangement>		Trap_pl;
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

	for (auto face: mesh.faces()) {
		std::vector<OpenMesh::MyMesh::Point> points;
		points.reserve(3);

		for (auto vertex: face.vertices()) {
			points.push_back(mesh.point(vertex));
		}
		
		auto normal = OpenMesh::cross(points[1] - points[0], points[2] - points[0]);
		
		for (auto vertex: face.vertices()) {
			mesh.set_normal(vertex, mesh.normal(vertex) + normal);
		}
	}

	for (auto vertex: mesh.vertices()) {
		auto normal = mesh.normal(vertex);
		mesh.set_normal(vertex, normal / normal.length());
	}

	std::cout << "Normal computed" << std::endl;

	int num_textures = mesh.property(textures_files).size();

	// Arrangement on textures
	CGAL::Arrangement arrangements[num_textures];
	CGAL::Trap_pl* trap_pl[num_textures];

	for (int i = 0; i < num_textures; i++) {
		trap_pl[i] = new CGAL::Trap_pl(arrangements[i]);
	}

	for (auto edge: mesh.edges()) {
		auto p0 = mesh.point(edge.v0());
		auto p1 = mesh.point(edge.v1());
		auto n0 = mesh.normal(edge.v0());
		auto n1 = mesh.normal(edge.v1());

		if (!edge.is_boundary()) {
			auto tex0 = mesh.property(textures_index, edge.h0().face());
			auto from0 = mesh.texcoord2D(edge.h0().prev());
			auto to0 = mesh.texcoord2D(edge.h0());
			auto tex1 = mesh.property(textures_index, edge.h1().face());
			auto from1 = mesh.texcoord2D(edge.h1().prev());
			auto to1 = mesh.texcoord2D(edge.h1());

			bool one_edge = true;
			if (tex0 != tex1 || from0[0] != to1[0] || from0[1] != to1[1] || from1[0] != to0[0] || from1[1] != to0[1]) {
				one_edge = false;
			}

			if (one_edge) {
				if (from0[0] != to0[0] || from0[1] != to0[1]) {
					auto c = insert_non_intersecting_curve(arrangements[tex0], CGAL::Segment(CGAL::Traits::Point_2(from0[0], from0[1]), CGAL::Traits::Point_2(to0[0], to0[1])), *trap_pl[tex0]);
					c->source()->set_data(CGAL::VertexData(p0, n0));
					c->target()->set_data(CGAL::VertexData(p1, n1));
					c->set_data(CGAL::HalfedgeData(mesh.property(label, edge.h0().face())));
					c->twin()->set_data(CGAL::HalfedgeData(mesh.property(label, edge.h1().face())));
				}
			} else {
				if (from0[0] != to0[0] || from0[1] != to0[1]) {
					auto c = insert_non_intersecting_curve(arrangements[tex0], CGAL::Segment(CGAL::Traits::Point_2(from0[0], from0[1]), CGAL::Traits::Point_2(to0[0], to0[1])), *trap_pl[tex0]);
					c->source()->set_data(CGAL::VertexData(p0, n0));
					c->target()->set_data(CGAL::VertexData(p1, n1));
					c->set_data(CGAL::HalfedgeData(mesh.property(label, edge.h0().face())));
				}
				if (from1[0] != to1[0] || from1[1] != to1[1]) {
					auto c = insert_non_intersecting_curve(arrangements[tex1], CGAL::Segment(CGAL::Traits::Point_2(from1[0], from1[1]), CGAL::Traits::Point_2(to1[0], to1[1])), *trap_pl[tex1]);
					c->source()->set_data(CGAL::VertexData(p1, n1));
					c->target()->set_data(CGAL::VertexData(p0, n0));
					c->set_data(CGAL::HalfedgeData(mesh.property(label, edge.h1().face())));
				}
			}
		} else {
			if (!edge.h0().is_boundary()) {
				auto tex0 = mesh.property(textures_index, edge.h0().face());
				auto from0 = mesh.texcoord2D(edge.h0().prev());
				auto to0 = mesh.texcoord2D(edge.h0());
				if (from0[0] != to0[0] || from0[1] != to0[1]) {
					auto c = insert_non_intersecting_curve(arrangements[tex0], CGAL::Segment(CGAL::Traits::Point_2(from0[0], from0[1]), CGAL::Traits::Point_2(to0[0], to0[1])), *trap_pl[tex0]);
					c->source()->set_data(CGAL::VertexData(p0, n0));
					c->target()->set_data(CGAL::VertexData(p1, n1));
					c->set_data(CGAL::HalfedgeData(mesh.property(label, edge.h0().face())));
				}
			} else {
				auto tex1 = mesh.property(textures_index, edge.h1().face());
				auto from1 = mesh.texcoord2D(edge.h1().prev());
				auto to1 = mesh.texcoord2D(edge.h1());
				if (from1[0] != to1[0] || from1[1] != to1[1]) {
					auto c = insert_non_intersecting_curve(arrangements[tex1], CGAL::Segment(CGAL::Traits::Point_2(from1[0], from1[1]), CGAL::Traits::Point_2(to1[0], to1[1])), *trap_pl[tex1]);
					c->source()->set_data(CGAL::VertexData(p1, n1));
					c->target()->set_data(CGAL::VertexData(p0, n0));
					c->set_data(CGAL::HalfedgeData(mesh.property(label, edge.h1().face())));
				}
			}
		}
	}

	std::cout << "Arrangement computed" << std::endl;

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

	for (int i = 0; i < num_textures; i++) {
		int width, height, n;
		unsigned char *data = stbi_load((std::string(argv[1]).substr(0, std::string(argv[1]).find_last_of('/')+1)+mesh.property(textures_files)[i]).c_str(), &width, &height, &n, 3);
		
		for (int v = 0; v < height; v++) {
			for (int u = 0; u < width; u++) {
				auto uv_point = CGAL::Traits::Point_2((0.5+u)/width, (0.5+v)/height);
				auto obj = trap_pl[i]->locate(uv_point);
				const CGAL::Arrangement::Vertex_const_handle*   vh;
				const CGAL::Arrangement::Halfedge_const_handle* eh;
				const CGAL::Arrangement::Face_const_handle*     fh;

				CGAL::Point_set_kernel::Point_3 point;
				CGAL::Point_set_kernel::Vector_3 normal;
				int label;
				bool point_localised = true;

				if (fh = boost::get<CGAL::Arrangement::Face_const_handle>(&obj)) {// located inside a face
					if ((*fh)->has_outer_ccb()) {
						auto halfedge = (*fh)->outer_ccb();
						auto v0 = halfedge++->source();
						auto v1 = halfedge++->source();
						auto v2 = halfedge++->source();
						auto lambda = CGAL::Barycentric_coordinates::compute_triangle_coordinates_2<CGAL::Traits>(v0->point(), v1->point(), v2->point(), uv_point);
						point = CGAL::Point_set_kernel::Point_3(0,0,0) + CGAL::to_double(lambda[0]) * CGAL::Point_set_kernel::Vector_3(CGAL::Point_set_kernel::Point_3(0,0,0), v0->data().point) + CGAL::to_double(lambda[1]) * CGAL::Point_set_kernel::Vector_3(CGAL::Point_set_kernel::Point_3(0,0,0), v1->data().point) + CGAL::to_double(lambda[2]) * CGAL::Point_set_kernel::Vector_3(CGAL::Point_set_kernel::Point_3(0,0,0), v2->data().point);
						normal = CGAL::to_double(lambda[0]) * v0->data().normal + CGAL::to_double(lambda[1]) * v1->data().normal + CGAL::to_double(lambda[2]) * v2->data().normal;
						label = (*fh)->outer_ccb()->data().label;
					} else {
						point_localised = false;
					}
				} else if (eh = boost::get<CGAL::Arrangement::Halfedge_const_handle>(&obj)) {// located on an edge
					auto lambda = CGAL::Barycentric_coordinates::compute_segment_coordinates_2<CGAL::Traits>((*eh)->source()->point(), (*eh)->target()->point(), uv_point);
					point = CGAL::Point_set_kernel::Point_3(0,0,0) + CGAL::to_double(lambda[0]) * CGAL::Point_set_kernel::Vector_3(CGAL::Point_set_kernel::Point_3(0,0,0), (*eh)->source()->data().point) + CGAL::to_double(lambda[1]) * CGAL::Point_set_kernel::Vector_3(CGAL::Point_set_kernel::Point_3(0,0,0), (*eh)->target()->data().point);
					normal = CGAL::to_double(lambda[0]) * (*eh)->source()->data().normal + CGAL::to_double(lambda[1]) * (*eh)->target()->data().normal;
					label = (*eh)->data().label;
					if (label == -1) {
						label = (*eh)->twin()->data().label;
					}
				} else if (vh = boost::get<CGAL::Arrangement::Vertex_const_handle>(&obj)) {// located on a vertex
					point = (*vh)->data().point;
					normal = (*vh)->data().normal;
					label = (*vh)->incident_halfedges()->data().label;
					if (label == -1) {
						label = (*vh)->incident_halfedges()->twin()->data().label;
					}
				} else {CGAL_error_msg("Invalid object.");}

				if (point_localised) {
					auto point_xyz = point_set.insert(point, normal);
					ps_label[*point_xyz] = label;
					ps_u[*point_xyz] = u;
					ps_v[*point_xyz] = v;

					ps_red[*point_xyz] = data[3*(u + v*width) + 0];
					ps_green[*point_xyz] = data[3*(u + v*width) + 1];
					ps_blue[*point_xyz] = data[3*(u + v*width) + 2];
				}
			}
		}
	}

	std::cout << "Point cloud computed" << std::endl;

	// Output file
	std::ofstream out_file (argv[2]);
	CGAL::write_ply_point_set(out_file, point_set);
		
	return 0;
}
