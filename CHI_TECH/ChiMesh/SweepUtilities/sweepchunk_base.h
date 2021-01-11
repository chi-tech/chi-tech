#ifndef _chi_sweepchunk_base_h
#define _chi_sweepchunk_base_h

//###################################################################
/**Sweep work function*/
class chi_mesh::sweep_management::SweepChunk
{
public:
  std::vector<double>*        x;
  bool                        suppress_surface_src;

  SweepChunk(std::vector<double>* destination_phi, bool suppress_src)
    : x(destination_phi), suppress_surface_src(suppress_src)
  {}

  /**Sets the location where flux moments are to be written.*/
  void SetDestinationPhi(std::vector<double>* destination_phi)
  {
    x = destination_phi;
  }

  virtual ~SweepChunk()
  {};

  /**Sweep chunks should override this.*/
  virtual void Sweep(AngleSet* angle_set)
  {}
};

#endif
