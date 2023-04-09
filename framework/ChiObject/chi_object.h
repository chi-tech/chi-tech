#ifndef CHITECH_CHI_OBJECT_H
#define CHITECH_CHI_OBJECT_H



namespace chi_objects
{
  class ChiLog;
  class MPI_Info;

  class ChiObject
  {
  private:
    ChiLog& log_;
    MPI_Info& mpi_;
  public:
    ChiObject();

    ChiLog& Log();
    MPI_Info& MPI();

    virtual ~ChiObject() = default;
  };
}

#endif // CHITECH_CHI_OBJECT_H
