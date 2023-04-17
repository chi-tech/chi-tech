#ifndef CHITECH_CHI_OBJECT_H
#define CHITECH_CHI_OBJECT_H



namespace chi_objects
{
  class ChiLog;
  class MPI_Info;

}

class ChiObject
{
private:
  chi_objects::ChiLog& log_;   ///< Reference to the logger
  chi_objects::MPI_Info& mpi_; ///< Reference to mpi information
public:
  ChiObject();

  /**Returns a reference to the logger instance. The logger is a singleton.*/
  chi_objects::ChiLog& Log();

  /**Returns a reference to the MPI_Info instance. The MPI_Info object is a
 * singleton.*/
  chi_objects::MPI_Info& MPI();

  virtual ~ChiObject() = default;
};

#endif // CHITECH_CHI_OBJECT_H
