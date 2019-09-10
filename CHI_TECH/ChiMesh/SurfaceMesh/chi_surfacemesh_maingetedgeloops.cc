#include <stdexcept>
#include"chi_surfacemesh.h"
#include <iostream>

/**Obtains a list of edges forming loops.*/
chi_mesh::EdgeLoopCollection* chi_mesh::SurfaceMesh::GetEdgeLoops()
{
  //================================================== Create new collection
  chi_mesh::EdgeLoopCollection* edge_loops = new chi_mesh::EdgeLoopCollection;

  //================================================== Build master list
  chi_mesh::EdgeList unused_edge_list;
  //============================================= Loop over faces
  std::vector<chi_mesh::Face>::iterator cur_face;
  for (cur_face = faces.begin();
          cur_face != faces.end();
          cur_face++)
  {
    //====================================== Loop over edges
    for (int e=0;e<3;e++)
    {
      if (cur_face->e_index[e][2]<0)
      {
        chi_mesh::Edge new_edge;
        new_edge.v_index[0] = cur_face->e_index[e][0];
        new_edge.v_index[1] = cur_face->e_index[e][1];

        new_edge.f_index[0] = cur_face - faces.begin();
        new_edge.f_index[2] = std::distance(faces.begin(),cur_face);

        try{
          new_edge.vertices[0] = vertices.at(new_edge.v_index[0]);
          new_edge.vertices[1] = vertices.at(new_edge.v_index[1]);

          unused_edge_list.push_back(new_edge);
        }
        catch(std::out_of_range o)
        {
          std::cerr << "EXCEPTION: Invalid vertex index\n";
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  //================================================== Check if open edges
//  if (unused_edge_list.size()==0)
//  {
//    std::cerr << "ERROR: No open edges found in function GetEdgeLoops.\n";
//    exit(EXIT_FAILURE);
//  }
//  else
//  {
//    std::cout << "Number of open edges found = " << unused_edge_list.size();
//    std::cout << std::endl;
//  }

  //============================================= Process lines as edgeloops
  for (unsigned ell=0; ell<lines.size(); ell++)
  {
    //chi_mesh::EdgeLoop* new_edge_loop = new EdgeLoop;
    chi_mesh::Edge new_edge;
    new_edge.v_index[0] = lines[ell].v_index[0];
    new_edge.v_index[1] = lines[ell].v_index[1];

    new_edge.vertices[0] = lines[ell].vertices[0];
    new_edge.vertices[1] = lines[ell].vertices[1];

    unused_edge_list.push_back(new_edge);

    //new_edge_loop->edges.push_back(new_edge);

    //edge_loops->push_back(new_edge_loop);
  }

  if (unused_edge_list.size()>0)
  {
    //================================================== Recursively find loops
    //============================================= Initial loop
    chi_mesh::EdgeLoop*  firstLoop = new chi_mesh::EdgeLoop;
    firstLoop->edges.push_back(unused_edge_list.back());
    unused_edge_list.pop_back();
    edge_loops->push_back(firstLoop);

    int iter=0;
    while (unused_edge_list.size()>0)
    {
      iter++;
      bool match_found = false;
      //=========================================== Reverse over each unused edge
      chi_mesh::EdgeList::reverse_iterator cur_edge,used_edge;
      for (cur_edge = unused_edge_list.rbegin();
           cur_edge != unused_edge_list.rend();
           ++cur_edge)
      {
        //==================================== Loop over each edge loop
        chi_mesh::EdgeLoopCollection::iterator cur_loop;
        for (cur_loop = edge_loops->begin();
             cur_loop != edge_loops->end();
             cur_loop++)
        {
          //============================= Get first and last edge
          chi_mesh::Edge* frst_edge = &(*cur_loop)->edges.front();
          chi_mesh::Edge* last_edge = &(*cur_loop)->edges.back();

          //============================= Compare to first edge
          if (cur_edge->v_index[1]==frst_edge->v_index[0])
          {

            match_found = true;
            (*cur_loop)->edges.insert((*cur_loop)->edges.begin(),*cur_edge);

            std::advance(cur_edge,1);
            unused_edge_list.erase(cur_edge.base());
            break;
          }

          //============================= Compare to last edge
          if (cur_edge->v_index[1]==last_edge->v_index[0])
          {

            match_found = true;
            (*cur_loop)->edges.push_back(*cur_edge);

            std::advance(cur_edge,1);
            unused_edge_list.erase(cur_edge.base());
            break;
          }

          if (match_found) {break;}
        }
        if (match_found) {break;}
      }

      if (!match_found)
      {
        chi_mesh::EdgeLoop*  new_loop = new chi_mesh::EdgeLoop;
        new_loop->edges.push_back(unused_edge_list.back());
        unused_edge_list.pop_back();
        edge_loops->push_back(new_loop);
      }

    }
  }

//  //============================================= Process lines as edgeloops
//  for (unsigned ell=0; ell<lines.size(); ell++)
//  {
//    chi_mesh::EdgeLoop* new_edge_loop = new EdgeLoop;
//    chi_mesh::Edge new_edge;
//    new_edge.v_index[0] = lines[ell].v_index[0];
//    new_edge.v_index[1] = lines[ell].v_index[1];
//
//    new_edge.vertices[0] = lines[ell].vertices[0];
//    new_edge.vertices[1] = lines[ell].vertices[1];
//
//    new_edge_loop->edges.push_back(new_edge);
//
//    edge_loops->push_back(new_edge_loop);
//  }

  //============================================= Verbose output

//  chi_mesh::EdgeLoopCollection::iterator cur_loop;
//  for (cur_loop = edge_loops->begin();
//       cur_loop != edge_loops->end();
//       cur_loop++)
//  {
//    printf("Edge Loop %d\n",cur_loop-edge_loops->begin());
//
//    chi_mesh::EdgeList::iterator used_edge;
//    for (used_edge = (*cur_loop)->edges.begin();
//         used_edge != (*cur_loop)->edges.end();
//         used_edge++)
//    {
//      printf("Used edge %d->%d\n",used_edge->v_index[0],used_edge->v_index[1]);
//    }
//
//  }



  return edge_loops;
}



