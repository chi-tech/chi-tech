

//###################################################################
/**This is a simple data structure of basically pointers to
 * objects needed in matrix free operations.*/
struct KSP_DATA_CONTEXT
{
  CHI_NPTRANSPORT* solver;
  SweepChunk*      sweep_chunk;
  int              group_set_num;
  NPT_GROUPSET*    groupset;
  KSP              krylov_solver;
  Vec              x_temp;
  chi_mesh::SweepManagement::SweepScheduler* sweepScheduler;
};