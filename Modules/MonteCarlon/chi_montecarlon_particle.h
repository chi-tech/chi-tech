#ifndef _chi_montecarlon_particle_h
#define _chi_montecarlon_particle_h

#include <ChiMesh/chi_mesh.h>

/**Structure for storing particle information.*/
struct chi_montecarlon::Particle
{
  //Pos
  chi_mesh::Vector pos;

  //Dir
  chi_mesh::Vector dir = chi_mesh::Vector(0.0,0.0,1.0);

  double w = 1.0; //Weight
  int egrp = 0; //Energy group
  int cur_cell_ind = -1;

  bool alive = true;
  bool banked = false;

  //=================================== Copy operator
  Particle& operator=(const Particle& that)
  {
    this->pos.x = that.pos.x;  this->dir.x = that.dir.x;
    this->pos.y = that.pos.y;  this->dir.y = that.dir.y;
    this->pos.z = that.pos.z;  this->dir.z = that.dir.z;

    this->w = that.w;
    this->egrp = that.egrp;
    this->cur_cell_ind = that.cur_cell_ind;

    this->alive = that.alive;

    return *this;
  }
};


#endif
