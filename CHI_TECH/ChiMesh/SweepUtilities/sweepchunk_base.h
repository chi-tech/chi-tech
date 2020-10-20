#ifndef _chi_sweepchunk_base_h
#define _chi_sweepchunk_base_h

//###################################################################
/**Sweep work function*/
class chi_mesh::sweep_management::SweepChunk
{
public:
  int                         angle_num;
  std::vector<double>*        x;
  bool                        suppress_surface_src;

public:
  /**Sets the location where flux moments are to be written.*/
  void SetDestinationPhi(std::vector<double>* destination_phi)
  {
    x = destination_phi;
  }

  virtual ~SweepChunk(){};

public:
  /**Sweep chunks should override this.*/
  virtual void Sweep(AngleSet* angle_set)
  {

  }
};

#endif
