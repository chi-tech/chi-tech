#ifndef chi_runtime_h
#define chi_runtime_h

int  ChiTechRunInteractive(int argc, char** argv);
int  ChiTechRunBatch(int argc, char** argv);
int  ChiTechInitialize(int argc, char** argv);
void ChiTechFinalize();

#endif