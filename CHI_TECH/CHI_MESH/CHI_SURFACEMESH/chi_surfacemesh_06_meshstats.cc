#include "chi_surfacemesh.h"

#include <CHI_GRAPH/chi_graph.h>

#include <chi_log.h>

extern ChiLog chi_log;

//###################################################################
/**Checks for cyclic dependencies in this mesh.
 * Transport type sweeps have a step where the inter-cell dependency
 * is acyclically sorted. This step is repeated here.*/
void chi_mesh::SurfaceMesh::CheckCyclicDependencies(int num_angles)
{
  double tolerance = 1.0e-8;
  double dvarphi = 2.0*M_PI/num_angles;
  chi_mesh::Vector khat(0.0,0.0,1.0);


  //================================================== Loop over angles
  for (int a=0; a<num_angles; a++)
  {
    double varphi = 0.5*dvarphi + a*dvarphi;

    chi_mesh::Vector omega;
    omega.x = cos(varphi);
    omega.y = sin(varphi);
    omega.z = 0.0;

    //================================= Add all polyfaces to graph
    CHI_D_GRAPH G;
    size_t num_loc_cells = poly_faces.size();
    for (size_t c=0; c<num_loc_cells; c++)
      boost::add_vertex(G);

    //================================= Now construct dependencies
    for (size_t c=0; c<num_loc_cells; c++)
    {
      chi_mesh::PolyFace* face = poly_faces[c];

      size_t num_edges = face->edges.size();
      for (size_t e=0; e<num_edges; e++)
      {
        int v0i = face->edges[e][0];
        int v1i = face->edges[e][1];

        chi_mesh::Vector v01 = vertices[v1i] - vertices[v0i];
        chi_mesh::Vector n = v01.Cross(khat); n=n/n.Norm();

        double mu = omega.Dot(n);
        int neighbor = face->edges[e][2];
        if ( (mu > (0.0 + tolerance)) and (neighbor >= 0))
        {
          boost::add_edge(c,neighbor,G);
        }//if outgoing
      }//for edge
    }//for cell

    //================================================== Generic topological
    //                                                   sorting
    typedef boost::graph_traits<CHI_D_GRAPH>::vertex_descriptor gVertex;

    boost::property_map<CHI_D_GRAPH, boost::vertex_index_t>::type
      index_map = get(boost::vertex_index, G);

    std::vector<gVertex> sorted_list;
    try{
      boost::topological_sort(G,std::back_inserter(sorted_list));
    }
    catch (boost::bad_graph exc)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Function CheckCyclicDependencies. Detected cyclic depency.";
      exit(EXIT_FAILURE);
    }

    //================================= Cleanup
    G.clearing_graph();
  }//for angles

  GetMeshStats();

  chi_log.Log(LOG_0)
    << "Cyclic dependency check complete. No cycles or "
    << "bad mesh elements were detected";
}

//###################################################################
/**Gets simple mesh statistics.
 */
void chi_mesh::SurfaceMesh::GetMeshStats()
{
  std::vector<double> areas;
  std::vector<double> histo_bins;
  std::vector<int>    histo;

  int num_negative_areas = 0;


  //============================================= Compute areas for
  //                                              each polyface
  size_t num_loc_cells = poly_faces.size();
  areas.resize(num_loc_cells);
  double max_area = 0.0;
  for (size_t c=0; c<num_loc_cells; c++)
  {
    chi_mesh::PolyFace* face = poly_faces[c];

    size_t num_edges = face->edges.size();
    double area = 0.0;
    for (size_t e=0; e<num_edges; e++)
    {
      int v0i = face->edges[e][0];
      int v1i = face->edges[e][1];

      chi_mesh::Vector v01 = vertices[v1i] - vertices[v0i];
      chi_mesh::Vector v02 = face->face_centroid - vertices[v0i];

      //This is essentially the combine of the triangle for each side

      area += 0.5*(v01.x*v02.y - v01.y*v02.x);
    }//for edge

    areas[c] = area;
    if (area > max_area)
      max_area = area;

    if (area <= 0.0)
      num_negative_areas += 1;
  }//for cell

  //============================================= Sort the areas
  std::sort(areas.begin(),areas.end(),std::greater<double>());

  //============================================= Compute histogram bins
  histo_bins.resize(10);
  histo.resize(10,0);
  histo_bins[0] = max_area*1.05;
  for (int i=1;i<10; i++)
    histo_bins[i] = histo_bins[i-1]/2.0;

  //============================================= Polulate histogram
  for (auto area : areas)
  {
    int home_bin = 9;
    for (int i=0;i<10; i++)
    {

      if (area <= histo_bins[i])
        home_bin = i;

    }//check bins

    histo[home_bin] += 1;
  }//for areas


  std::stringstream output;
  for (int i=0;i<10; i++)
  {
    char buff[100];
    sprintf(buff,"%11.3e",histo_bins[i]);

    output << "Areas < " << buff << " = " << histo[i] << "\n";
  }
  output << "Number of negative or zero faces = " << num_negative_areas;

  chi_log.Log(LOG_0) << output.str();

}