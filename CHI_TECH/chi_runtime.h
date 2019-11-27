#ifndef chi_runtime_h
#define chi_runtime_h

void ChiTechParseArguments(int argc, char** argv);
void ChiTechRunInteractive(int argc, char** argv);
void ChiTechRunBatch(int argc, char** argv);
void ChiTechInitialize(int argc, char** argv);
void ChiTechFinalize();

#endif