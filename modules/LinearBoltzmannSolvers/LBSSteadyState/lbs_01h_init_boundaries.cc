#include "lbs_linear_boltzmann_solver.h"

#include "chi_log.h"

//###################################################################
/**Initializes transport related boundaries. */
void lbs::SteadyStateSolver::InitializeBoundaries()
{
  //================================================== Determine boundary-ids
  //                                                   involved in the problem
  std::set<uint64_t> globl_unique_bids_set;
  {
    std::set<uint64_t> local_unique_bids_set;
    for (const auto& cell : grid_ptr_->local_cells)
      for (const auto& face : cell.faces)
        if (not face.has_neighbor)
          local_unique_bids_set.insert(face.neighbor_id);

    std::vector<uint64_t> local_unique_bids(local_unique_bids_set.begin(),
                                            local_unique_bids_set.end());
    const int local_num_unique_bids = static_cast<int>(local_unique_bids.size());
    std::vector<int> recvcounts(chi::mpi.process_count, 0);

    MPI_Allgather(&local_num_unique_bids,      //sendbuf
                  1, MPI_INT,                  //sendcount+type
                  recvcounts.data(),           //recvbuf
                  1, MPI_INT,              //recvcount+type
                  MPI_COMM_WORLD);             //comm

    std::vector<int> recvdispls(chi::mpi.process_count, 0);

    int running_displacement = 0;
    for (int locI=0; locI<chi::mpi.process_count; ++locI)
    {
      recvdispls[locI] = running_displacement;
      running_displacement += recvcounts[locI];
    }

    std::vector<uint64_t> recvbuf(running_displacement, 0);

    MPI_Allgatherv(local_unique_bids.data(), //sendbuf
                   local_num_unique_bids,    //sendcount
                   MPI_UINT64_T,             //sendtype
                   recvbuf.data(),           //recvbuf
                   recvcounts.data(),        //recvcounts
                   recvdispls.data(),        //recvdispls
                   MPI_UINT64_T,             //recvtype
                   MPI_COMM_WORLD);          //comm

    globl_unique_bids_set = local_unique_bids_set; //give it a head start

    for (uint64_t bid : recvbuf)
      globl_unique_bids_set.insert(bid);
  }


  //================================================== Initialize default
  //                                                   incident boundary
  typedef chi_mesh::sweep_management::BoundaryVacuum SweepVacuumBndry;
  typedef chi_mesh::sweep_management::BoundaryIncidentHomogenous SweepIncHomoBndry;
  typedef chi_mesh::sweep_management::BoundaryReflecting SweepReflectingBndry;

  const std::vector<double> zero_Gvec(num_groups_, 0.0);

  const chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  const chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  sweep_boundaries_.clear();
  for (uint64_t bid : globl_unique_bids_set)
    if (boundary_preferences_.count(bid) == 0 and
        sweep_boundaries_.count(bid) == 0)
    {
      sweep_boundaries_[bid] = std::make_shared<SweepVacuumBndry>(zero_Gvec);
    }//defaulted
    else if (sweep_boundaries_.count(bid) == 0)
    {
      const auto& bndry_pref = boundary_preferences_.at(bid);
      if (bndry_pref.type == lbs::BoundaryType::VACUUM)
        sweep_boundaries_[bid] =
          std::make_shared<SweepVacuumBndry>(zero_Gvec);
      else if (bndry_pref.type == lbs::BoundaryType::INCIDENT_ISOTROPIC)
        sweep_boundaries_[bid] =
          std::make_shared<SweepIncHomoBndry>(bndry_pref.isotropic_mg_source);
      else if (bndry_pref.type == lbs::BoundaryType::REFLECTING)
      {
        chi_mesh::Normal normal;
        if (bid == 0) normal = ihat;
        if (bid == 1) normal = ihat * -1.0;
        if (bid == 2) normal = jhat;
        if (bid == 3) normal = jhat * -1.0;
        if (bid == 4) normal = khat;
        if (bid == 5) normal = khat * -1.0;

        sweep_boundaries_[bid] =
          std::make_shared<SweepReflectingBndry>(zero_Gvec, normal);
      }
    }//non-defaulted

}