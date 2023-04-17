#include "chi_object.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

/**Base constructor.*/
ChiObject::ChiObject() : log_(chi::log), mpi_(chi::mpi) {}

/**Returns a reference to the logger instance. The logger is a singleton.*/
chi_objects::ChiLog& ChiObject::Log() { return log_; }

/**Returns a reference to the MPI_Info instance. The MPI_Info object is a
 * singleton.*/
chi_objects::MPI_Info& ChiObject::MPI() { return mpi_; }
