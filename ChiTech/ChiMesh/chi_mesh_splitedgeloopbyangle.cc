#include<math.h>
#include<iostream>
#include "chi_mesh.h"

chi_mesh::EdgeLoopCollection* chi_mesh::SplitEdgeLoopByAngle(EdgeLoop *input,
                                                   double angle)
{
  //================================================== Define cosine of angle
  double cosineangle = cos(angle*M_PI/180.0);
  //printf("costheta=%f\n",cosineangle);

  //================================================== Init the collection
  chi_mesh::EdgeLoopCollection* new_coll = new chi_mesh::EdgeLoopCollection;


  //================================================== Make unused list
  chi_mesh::EdgeList unused_edges;
  chi_mesh::EdgeList::iterator cur_edge;
  for (cur_edge = input->edges.begin();
          cur_edge != input->edges.end();
          cur_edge++)
  {
    unused_edges.push_back(*cur_edge);
  }

  //================================================== Init the first loop
  //Insert the last unused edge into the new loop
  chi_mesh::EdgeLoop* new_loop = new chi_mesh::EdgeLoop;
  new_loop->edges.insert(new_loop->edges.begin(),unused_edges.back());
  unused_edges.pop_back();
  new_coll->push_back(new_loop);

  //================================================== Iterate over all edges
  while (unused_edges.size()>0)
  {
    bool match_found = false;

    //=========================================== Loop over unused edges (revrs)
    chi_mesh::EdgeList::reverse_iterator rcur_edge;
    for (rcur_edge = unused_edges.rbegin();
            rcur_edge != unused_edges.rend();
            rcur_edge++)
    {
      //==================================== Loop over existing loops
      chi_mesh::EdgeLoopCollection::iterator cur_loop;
      for (cur_loop = new_coll->begin();
              cur_loop != new_coll->end();
              cur_loop++)
      {
        //printf("In loop %d\n",cur_loop-new_coll->begin());
        //============================= Loop over edges
        chi_mesh::EdgeList::iterator ref_edge;
        for (ref_edge = (*cur_loop)->edges.begin();
                ref_edge != (*cur_loop)->edges.end();
                ref_edge++)
        {
          chi_mesh::Vector3 vu = rcur_edge->vertices[1] -
                                 rcur_edge->vertices[0];
          chi_mesh::Vector3 vr = ref_edge->vertices[1] -
                                 ref_edge->vertices[0];
          vu = vu/vu.Norm();
          vr = vr/vr.Norm();

//          printf("Comparing edge %d->%d",rcur_edge->v_index[0],rcur_edge->v_index[1]);
//          printf(" to edge %d->%d",ref_edge->v_index[0],ref_edge->v_index[1]);
//          printf(" with cosineangle=%f",fabs(vr.Dot(vu)));
//
//          printf("    vu=%.3f %.3f  %.3f",vu.x,vu.y,vu.z);
//          printf("    vr=%.3f %.3f  %.3f\n",vr.x,vr.y,vr.z);
          if ((fabs(vr.Dot(vu)))>=cosineangle)
          {
            if (rcur_edge->v_index[1]==ref_edge->v_index[0])
            {
              //insert ahead
              (*cur_loop)->edges.insert(ref_edge,*rcur_edge);
              match_found=true;
              std::advance(rcur_edge,1);
              unused_edges.erase(rcur_edge.base());
              break;
            }
            else if (rcur_edge->v_index[0]==ref_edge->v_index[1])
            {
              //add behind
              (*cur_loop)->edges.insert(ref_edge+1,*rcur_edge);
              match_found=true;
              std::advance(rcur_edge,1);
              unused_edges.erase(rcur_edge.base());
              break;
            }
          }

        }
        if (match_found) {break;}
      }
      if (match_found) {break;}
    }

    if (!match_found)
    {
      new_loop = new chi_mesh::EdgeLoop;
      new_loop->edges.insert(new_loop->edges.begin(),unused_edges.back());
      unused_edges.pop_back();
      new_coll->push_back(new_loop);
      //printf("New loop created\n");
    }


  }


  //============================================= Verbose output
  /*chi_mesh::EdgeLoopCollection::iterator cur_loop;
  for (cur_loop = new_coll->begin();
       cur_loop != new_coll->end();
       cur_loop++)
  {
    printf("Edge Loop %d\n",(int)(cur_loop-new_coll->begin()));

    chi_mesh::EdgeList::iterator used_edge;
    for (used_edge = (*cur_loop)->edges.begin();
         used_edge != (*cur_loop)->edges.end();
         used_edge++)
    {
      printf("Used edge %d->%d\n",used_edge->v_index[0],used_edge->v_index[1]);
    }

  }*/



  return new_coll;
}