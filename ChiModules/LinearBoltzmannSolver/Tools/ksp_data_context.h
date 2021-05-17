

//###################################################################
/**This is a simple data structure of basically pointers to
 * objects needed in matrix free operations.*/
struct KSPDataContext
{
  LinearBoltzmann::Solver& solver;
  SweepChunk&      sweep_chunk;
  LBSGroupset&    groupset;
  Vec&            operating_vector;
  chi_mesh::sweep_management::SweepScheduler& sweepScheduler;
  int last_iteration = -1;

  KSPDataContext(LinearBoltzmann::Solver& in_solver,
                 SweepChunk& in_sweep_chunk,
                 LBSGroupset& in_groupset,
                 Vec& in_operating_vector,
                 chi_mesh::sweep_management::SweepScheduler& in_sweep_scheduler) :
                 solver(in_solver),
                 sweep_chunk(in_sweep_chunk),
                 groupset(in_groupset),
                 operating_vector(in_operating_vector),
                 sweepScheduler(in_sweep_scheduler) {}
};