#ifndef CHI_SWEEPCHUNK_BASE_H
#define CHI_SWEEPCHUNK_BASE_H

#include <functional>

//###################################################################
/**Sweep work function*/
class chi_mesh::sweep_management::SweepChunk
{
private:
  std::vector<double>* x;
  bool surface_source_active;

public:
  /**
   * Convenient typdef for the moment call back function. See moment_callbacks.
   *  Arguments are:
   *  chi_mesh::sweep_management::SweepChunk *
   *  chi_mesh::sweep_management::AngleSet *
   */
  typedef std::function<void(chi_mesh::sweep_management::SweepChunk * sweeper,
                             chi_mesh::sweep_management::AngleSet* angle_set)>
                             MomentCallbackF;

  /**
   * Functions of type MomentCallbackF can be added to the moment_callbacks
   * vector and these can be called from within functions taking a
   * LBSGroupset instance. The intention is that this function can
   * be used as a general interface to retrieve angular flux values
   */
  std::vector<MomentCallbackF> moment_callbacks;

  SweepChunk(std::vector<double>& destination_phi, bool suppress_src)
    : x(&destination_phi), surface_source_active(suppress_src)
  {}

  /**Sets the location where flux moments are to be written.*/
  void SetDestinationPhi(std::vector<double>& destination_phi)
  {
    x = (&destination_phi);
  }

  /**Returns a reference to the output vector.*/
  std::vector<double>& GetDestinationPhi()
  {
    return *x;
  }

  /**Activates or deactives the surface src flag.*/
  void SetSurfaceSourceActiveFlag(bool flag_value)
  {
    surface_source_active = flag_value;
  }

  /**Returns the surface src-active flag.*/
  bool IsSurfaceSourceActive() const
  {return surface_source_active;}

  virtual ~SweepChunk() = default;

  /**Sweep chunks should override this.*/
  virtual void Sweep(AngleSet* angle_set)
  {}
};

#endif //CHI_SWEEPCHUNK_BASE_H
