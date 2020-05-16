#include "../chi_mesh.h"

#include "ChiMesh/SweepUtilities/SPDS/SPDS.h"
#include "ChiGraph/chi_directed_graph.h"

#include <chi_log.h>
#include <chi_mpi.h>
extern ChiMPI& chi_mpi;
extern ChiLog& chi_log;

#include <algorithm>

//###################################################################
/**Removes local cyclic dependencies.*/
void chi_mesh::sweep_management::
  RemoveLocalCyclicDependencies(chi_mesh::sweep_management::SPDS *sweep_order,
                                chi_graph::DirectedGraph &local_DG)
{
  //============================================= Utility lambdas
  auto IsInList = [](std::vector<int>& list, int val)
  {
    return std::find(list.begin(),list.end(),val) != list.end();
  };

  //============================================= Find initial SCCs
  auto SCCs = local_DG.FindStronglyConnectedComponents();

  //============================================= Remove bi-connected then
  //                                              tri-connected SCCs then
  //                                              n-connected
  for (auto& subDG : SCCs)
  {
    if (subDG.size()==2)
    {
      local_DG.RemoveEdge(subDG.front(), subDG.back());
      sweep_order->local_cyclic_dependencies.emplace_back(
        subDG.front(),
        subDG.back());
      if (chi_log.GetVerbosity() >= LOG_0VERBOSE_2)
        chi_log.Log(LOG_ALL)
          << "Bi-connected component removed: "
          << subDG.front() << "->"
          << subDG.back();
    }//bi-connected
    else if (subDG.size()==3)
    {
      bool found=false;
      for (int u : subDG)
      {
        for (int v : local_DG.vertices[u].ds_edge)
          if (IsInList(subDG,v))
          {
            found=true;
            local_DG.RemoveEdge(u, v);
            sweep_order->local_cyclic_dependencies.emplace_back(u, v);
            if (chi_log.GetVerbosity() >= LOG_0VERBOSE_2)
              chi_log.Log(LOG_ALL)
                << "Tri-connected component removed: "
                << u << "->" << v;
            break;
          }
        if (found) break;
      }//for u
    }//tri-connected
    else
    {
      //==================================== Add vertices to temporary graph
      chi_graph::DirectedGraph TG; //Temp Graph
      for (auto u : subDG)
        TG.AddVertex();

      //==================================== Add local connectivity
      int mapping_u = 0;
      for (auto u : subDG)
      {
        for (auto v : local_DG.vertices[u].ds_edge)
        {
          auto mapv = std::find(subDG.begin(),subDG.end(),v);
          if (mapv != subDG.end())
          {
            int mapping_v = mapv - subDG.begin();
            TG.AddEdge(mapping_u,mapping_v,local_DG.vertices[u].ds_weights[v]);
          }
        }//for v

        ++mapping_u;
      }//for u

      //==================================== Make a copy of the graph verts
      std::vector<chi_graph::GraphVertex> verts_copy;
      verts_copy.reserve(TG.vertices.size());
      for (auto& v : TG.vertices)
        verts_copy.push_back(v);

      //==================================== Solve the minimum Feedback
      //                                     Arc Set (FAS) problem
      auto s = TG.FindApproxMinimumFAS();

      //========================== Build a sequence map
      // This maps original sequence
      // to minFAS sequence. i.e. originally
      // we had v=0,1,2,3... and afterwards we
      // something like s=7,3,1,5,0,....
      // smap[v] then gives the position of v in s
      std::vector<int> smap(s.size(),-1);
      int count=0;
      for (int u: s)
        smap[u] = count++;

      //========================== Build edges to remove
      std::vector<std::pair<int,int>> edges_to_rem;
      for (auto& u : verts_copy)
      {
        int cur_map = smap[u.id];
        for (int v : u.ds_edge)
        {
          int adj_map = smap[v];
          if (adj_map < cur_map)
            edges_to_rem.emplace_back(u.id,v);
        }
      }

      for (auto& edge : edges_to_rem)
      {
        int u = subDG[edge.first];
        int v = subDG[edge.second];
        local_DG.RemoveEdge(u, v);
        sweep_order->local_cyclic_dependencies.emplace_back(u, v);
        if (chi_log.GetVerbosity() >= LOG_0VERBOSE_2)
          chi_log.Log(LOG_ALL)
            << "Removing edge " << u << " " << v;
      }

    }//n-connected
  }//for sub-DG


  //============================================= Find SSCs again
  // This step is like an insurance policy for if
  // something came through. There should be no SSCs
  // after the minFAS process, however, we just look
  // again in-case.
  SCCs = local_DG.FindStronglyConnectedComponents();

  //============================================= Rinse remove edges
  std::vector<std::pair<int,int>> edges_to_remove;
  int iter=0;
  while (not SCCs.empty())
  {
    if (chi_log.GetVerbosity() >= LOG_0VERBOSE_2)
      chi_log.Log(LOG_ALL)
        << "Inter cell cyclic dependency removal. Iteration " << ++iter;

    //=================================== Loop over sub-graphs
    edges_to_remove.clear();
    for (auto& subDG : SCCs)
    {
      //Identify edges to remove
      bool found_edge_to_remove = false;
      for (int v : subDG)
      {
        for (int w : local_DG.vertices[v].ds_edge)
          if (std::find(subDG.begin(), subDG.end(), w) != subDG.end())
          {
            found_edge_to_remove = true;
            edges_to_remove.emplace_back(v,w);
            break;
          }

        if (found_edge_to_remove) break;
      }//for v in subDG
    }//for subDG

    //Remove the edges
    for (auto& edge_to_remove : edges_to_remove)
    {
      local_DG.RemoveEdge(edge_to_remove.first, edge_to_remove.second);
      sweep_order->local_cyclic_dependencies.emplace_back(
        edge_to_remove.first,
        edge_to_remove.second);
    }

    // Refind SCCs
    SCCs = local_DG.FindStronglyConnectedComponents();
  }
}