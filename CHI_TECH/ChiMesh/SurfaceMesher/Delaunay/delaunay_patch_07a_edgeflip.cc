#include"delaunay_mesher.h"


//###################################################################
/**Restores the Delaunay conditioning of the mesh using edge-flipping.
 **/
void chi_mesh::SurfaceMesherDelaunay::DelaunayPatch::
  EdgeFlip(Edge& edge_to_flip)
{
  printf("Flipping Edge %d->%d\n",edge_to_flip.v_index[0],edge_to_flip.v_index[1]);

  //=============================================================== Get edges triangles
  Tri* tau_1 = &triangles.at(edge_to_flip.f_index[0]);
  Tri* tau_2 = &triangles.at(edge_to_flip.f_index[2]);

//  printf("Tau1[%d] %d->%d->%d\n",edge_to_flip.f_index[0],tau_1->v_index[0],
//                                                         tau_1->v_index[1],
//                                                         tau_1->v_index[2]);
//  printf("Tau2[%d] %d->%d->%d\n",edge_to_flip.f_index[1],tau_2->v_index[0],
//                                                         tau_2->v_index[1],
//                                                         tau_2->v_index[2]);

  //=============================================================== Find neighbour A and B
  int neighborA = -1;
  int neighborB = -1;
  if (edge_to_flip.f_index[1]==0)
  {
    neighborA = tau_1->e_index[1][2];
    neighborB = tau_1->e_index[2][2];
  }
  else if (edge_to_flip.f_index[1]==1)
  {
    neighborA = tau_1->e_index[2][2];
    neighborB = tau_1->e_index[0][2];
  }
  else
  {
    neighborA = tau_1->e_index[0][2];
    neighborB = tau_1->e_index[1][2];
  }
//  printf("NeighborA = %d\n",neighborA);
//  printf("NeighborB = %d\n",neighborB);

  //=============================================================== Find neighbour C and D
  int neighborC = -1;
  int neighborD = -1;
  if (edge_to_flip.f_index[3]==0)
  {
    neighborC = tau_2->e_index[1][2];
    neighborD = tau_2->e_index[2][2];
  }
  else if (edge_to_flip.f_index[3]==1)
  {
    neighborC = tau_2->e_index[2][2];
    neighborD = tau_2->e_index[0][2];
  }
  else
  {
    neighborC = tau_2->e_index[0][2];
    neighborD = tau_2->e_index[1][2];
  }
//  printf("NeighborC = %d\n",neighborC);
//  printf("NeighborD = %d\n",neighborD);

  //=============================================================== Find vertex A and B
  int vertA = -1;
  for (int v=0; v<3; v++)
  {
    if ( (tau_1->v_index[v]!=edge_to_flip.v_index[0]) &&
         (tau_1->v_index[v]!=edge_to_flip.v_index[1]) )
    {
      vertA = tau_1->v_index[v];
      break;
    }
  }
  int vertB = -1;
  for (int v=0; v<3; v++)
  {
    if ( (tau_2->v_index[v]!=edge_to_flip.v_index[0]) &&
         (tau_2->v_index[v]!=edge_to_flip.v_index[1]) )
    {
      vertB = tau_2->v_index[v];
      break;
    }
  }

  //=============================================================== Modify tau_1 and tau_2
  tau_1->v_index[0] = vertA;
  tau_1->v_index[1] = vertB;
  tau_1->v_index[2] = edge_to_flip.v_index[1];

  tau_1->e_index[0][2]=edge_to_flip.f_index[2];
  tau_1->e_index[1][2]=neighborD;
  tau_1->e_index[2][2]=neighborA;

  tau_1->e_index[0][3]=-1;
  tau_1->e_index[1][3]=-1;
  tau_1->e_index[2][3]=-1;

  tau_1->e_index[0][0] = tau_1->v_index[0]; tau_1->e_index[0][1] = tau_1->v_index[1];
  tau_1->e_index[1][0] = tau_1->v_index[1]; tau_1->e_index[1][1] = tau_1->v_index[2];
  tau_1->e_index[2][0] = tau_1->v_index[2]; tau_1->e_index[2][1] = tau_1->v_index[0];


  tau_2->v_index[0] = vertB;
  tau_2->v_index[1] = vertA;
  tau_2->v_index[2] = edge_to_flip.v_index[0];

  tau_2->e_index[0][2]=edge_to_flip.f_index[0];
  tau_2->e_index[1][2]=neighborB;
  tau_2->e_index[2][2]=neighborC;

  tau_2->e_index[0][3]=-1;
  tau_2->e_index[1][3]=-1;
  tau_2->e_index[2][3]=-1;

  tau_2->e_index[0][0] = tau_2->v_index[0]; tau_2->e_index[0][1] = tau_2->v_index[1];
  tau_2->e_index[1][0] = tau_2->v_index[1]; tau_2->e_index[1][1] = tau_2->v_index[2];
  tau_2->e_index[2][0] = tau_2->v_index[2]; tau_2->e_index[2][1] = tau_2->v_index[0];

//  printf("Tau1 %3d,%3d,%3d -- %3d,%3d,%3d   Tau2 %3d,%3d,%3d -- %3d,%3d,%3d\n",
//         tau_1->v_index[0],tau_1->v_index[1],tau_1->v_index[2],tau_1->e_index[0][2],tau_1->e_index[1][2],tau_1->e_index[2][2],
//         tau_2->v_index[0],tau_2->v_index[1],tau_2->v_index[2],tau_2->e_index[0][2],tau_2->e_index[1][2],tau_2->e_index[2][2]);

  //=============================================================== Update neighbors B and D
  Tri* tau_neighborA= new Tri();
  Tri* tau_neighborB= new Tri();
  Tri* tau_neighborC= new Tri();
  Tri* tau_neighborD= new Tri();
  if (neighborA>=0)
  {
    tau_neighborA = &triangles.at(neighborA);
  }
  if (neighborB>=0)
  {
    tau_neighborB = &triangles.at(neighborB);
  }
  if (neighborC>=0)
  {
    tau_neighborC = &triangles.at(neighborC);
  }
  if (neighborD>=0)
  {
    tau_neighborD = &triangles.at(neighborD);
  }

  for (int e=0; e<3; e++)
  {
    if (neighborB>=0)
    {
      if (tau_neighborB->e_index[e][2]==edge_to_flip.f_index[0])
      {
        tau_neighborB->e_index[e][2]=edge_to_flip.f_index[2];
        //modded_triangles.push_back(*tau_neighborB);
        break;
      }
    }
  }
  for (int e=0; e<3; e++)
  {
    if (neighborD>=0)
    {
      if (tau_neighborD->e_index[e][2]==edge_to_flip.f_index[2])
      {
        tau_neighborD->e_index[e][2]=edge_to_flip.f_index[0];
        //modded_triangles.push_back(*tau_neighborD);
        break;
      }
    }
  }

  //=============================================================== Update edge connections
  for (int e1=0; e1<3; e1++)
  {
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborA->e_index[e1][0]==tau_1->e_index[e2][1]) &&
           (tau_neighborA->e_index[e1][1]==tau_1->e_index[e2][0])  )
      {
        tau_neighborA->e_index[e1][3] = e2;
        tau_1->e_index[e2][3]         = e1;
      }
    }
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborA->e_index[e1][0]==tau_2->e_index[e2][1]) &&
           (tau_neighborA->e_index[e1][1]==tau_2->e_index[e2][0])  )
      {
        tau_neighborA->e_index[e1][3] = e2;
        tau_2->e_index[e2][3]         = e1;
      }
    }
  }

  for (int e1=0; e1<3; e1++)
  {
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborB->e_index[e1][0]==tau_1->e_index[e2][1]) &&
           (tau_neighborB->e_index[e1][1]==tau_1->e_index[e2][0])  )
      {
        tau_neighborB->e_index[e1][3] = e2;
        tau_1->e_index[e2][3]         = e1;
      }
    }
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborB->e_index[e1][0]==tau_2->e_index[e2][1]) &&
           (tau_neighborB->e_index[e1][1]==tau_2->e_index[e2][0])  )
      {
        tau_neighborB->e_index[e1][3] = e2;
        tau_2->e_index[e2][3]         = e1;
      }
    }
  }

  for (int e1=0; e1<3; e1++)
  {
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborC->e_index[e1][0]==tau_1->e_index[e2][1]) &&
           (tau_neighborC->e_index[e1][1]==tau_1->e_index[e2][0])  )
      {
        tau_neighborC->e_index[e1][3] = e2;
        tau_1->e_index[e2][3]         = e1;
      }
    }
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborC->e_index[e1][0]==tau_2->e_index[e2][1]) &&
           (tau_neighborC->e_index[e1][1]==tau_2->e_index[e2][0])  )
      {
        tau_neighborC->e_index[e1][3] = e2;
        tau_2->e_index[e2][3]         = e1;
      }
    }
  }

  for (int e1=0; e1<3; e1++)
  {
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborD->e_index[e1][0]==tau_1->e_index[e2][1]) &&
           (tau_neighborD->e_index[e1][1]==tau_1->e_index[e2][0])  )
      {
        tau_neighborD->e_index[e1][3] = e2;
        tau_1->e_index[e2][3]         = e1;
      }
    }
    for (int e2=0; e2<3; e2++)
    {
      if ( (tau_neighborD->e_index[e1][0]==tau_2->e_index[e2][1]) &&
           (tau_neighborD->e_index[e1][1]==tau_2->e_index[e2][0])  )
      {
        tau_neighborD->e_index[e1][3] = e2;
        tau_2->e_index[e2][3]         = e1;
      }
    }
  }


}