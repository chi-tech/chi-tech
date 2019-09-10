

//###################################################################
/**This is a simple data structure of basically pointers to
 * objects needed in matrix free operations.*/
struct KSPDataContext
{
  LinearBoltzmanSolver* solver;
  SweepChunk*      sweep_chunk;
  int              group_set_num;
  LBSGroupset*    groupset;
  KSP              krylov_solver;
  Vec              x_temp;
  chi_mesh::SweepManagement::SweepScheduler* sweepScheduler;
};