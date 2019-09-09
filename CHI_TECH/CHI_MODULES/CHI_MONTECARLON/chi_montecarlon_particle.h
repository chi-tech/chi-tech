#ifndef _chi_montecarlon_particle_h
#define _chi_montecarlon_particle_h

#include <ChiMesh/chi_mesh.h>

/**Structure for storing particle information.*/
struct chi_montecarlon::Particle
{
  //Pos
  chi_mesh::Vector pos;

  //Dir
  chi_mesh::Vector dir;

  double w; //Weight
  int egrp; //Energy group
  int cur_cell_ind;

  bool alive;

  //=================================== Constructor
  Particle()
  {
    pos.x=0.0;   dir.x=0.0;
    pos.y=0.0;   dir.y=0.0;
    pos.z=0.0;   dir.z=1.0;

    w=1.0;
    egrp = 0;
    cur_cell_ind = -1;

    alive = true;
  }

  //=================================== Copy operator
  void operator=(const Particle& that)
  {
    this->pos.x = that.pos.x;  this->dir.x = that.dir.x;
    this->pos.y = that.pos.y;  this->dir.y = that.dir.y;
    this->pos.z = that.pos.z;  this->dir.z = that.dir.z;

    this->w = that.w;
    this->egrp = that.egrp;
    this->cur_cell_ind = that.cur_cell_ind;

    this->alive = that.alive;
  }
};


#endif