#ifndef chi_runtime_h
#define chi_runtime_h

/**General utilities in ChiTech*/
class ChiTech
{
private:
  static void ParseArguments(int argc, char** argv);

public:
  static int  RunInteractive(int argc, char** argv);
  static int  RunBatch(int argc, char** argv);
  static int  Initialize(int argc, char** argv);
  static void Finalize();
};

#endif