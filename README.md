# Mesh-2-Point-Cloud
Transform Mesh and Atlas into Point Cloud (X Y Z Intensity R G B)

## Usage
```bash
$ mesh-2-pc [OPTION] MESH_FILE PC_FILE
```
OPTIONS:
 -h, --help                   Print this help anq quit.
 -s, --texelScale=S           Adapt texel size with parameter S (default 1.0).
 -p, --poissonDiskSampling=P  Use poisson isk sampling with parameter P.
 -g, --gridSubSampling=G      Subsample the final point cloud, keeping one point per voxel of size G.
 
