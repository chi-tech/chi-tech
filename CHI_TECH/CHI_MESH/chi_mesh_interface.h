#ifndef _chi_mesh_interface_h
#define _chi_mesh_interface_h

/**Structure for holding an interface. An interface can be surface-surface,
 * edge-edge, surface-edge or edge-surface.*/
struct chi_mesh::Interface
{
  chi_mesh::SurfaceMesh* master_surface;
  chi_mesh::LineMesh*    master_edge;

  chi_mesh::SurfaceMesh* slave_surface;
  chi_mesh::LineMesh*    slave_edge;

  Interface()
  {
    master_surface = nullptr;
    master_edge    = nullptr;
    slave_surface  = nullptr;
    slave_edge     = nullptr;
  }
};

#endif