# Mesh-2-Point-Cloud
Transform Mesh and Atlas into Point Cloud (X Y Z Intensity R G B)

## Requirements
- [CGAL](https://www.cgal.org/)
- [OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/)
 
## Compilation
```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

## Usage
```bash
$ ./mesh-2-pc [OPTION] MESH_FILE PC_FILE
```
OPTIONS:
- `-h`, `--help`                   Print this help anq quit.
- `-s S`, `--texelScale=S`           Adapt texel size with parameter S (default 1.0).
- `-p P`, `--poissonDiskSampling=P`  Use poisson isk sampling with parameter P.
- `-g G`, `--gridSubSampling=G`      Subsample the final point cloud, keeping one point per voxel of size G.
