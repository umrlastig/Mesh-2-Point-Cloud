#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

namespace OpenMesh {
	typedef TriMesh_ArrayKernelT<>	MyMesh;
}

namespace CGAL {
    typedef Simple_cartesian<double>		Kernel;
    typedef Surface_mesh<Kernel::Point_3>	Mesh;
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

	std::cout << mesh.n_vertices() << std::endl;
	
	return 0;
}
