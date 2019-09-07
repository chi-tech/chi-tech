/** \mainpage

ChiTech is an engineering development platform under development to support
various scientific simulations. At the core of the platform is the **Lua**
scripting engine which allows the user to interact with their C++ code using
a script-environment.\n
\n
Please be patient as we add more capabilities and tutorials to the code. For
now we have only a couple of tutorials.\n
\n

 # Index

 \subpage Tutorial00 \n
 \subpage MeshTutorials \n
 \subpage ProgrammerManual \n
 <a href="https://www.lua.org/manual/5.3/"><B>Lua 5.3 Reference</B></a>\n
 <a href="https://www.mcs.anl.gov/petsc/documentation/index.html"><B>PETSc Documentation</B></a>\n
 <a href="https://lorensen.github.io/VTKExamples/site/"><B>Visualization Tool Kit (VTK) documentation</B></a>\n


\n


 ## Input Functions

### Math Utilities
 
<table>
<tr><td width="33%">chiLegendre()</td><td width="33%">chiLegendreDerivative()</td><td width="33%">chiYlm()</td></tr>
<tr><td width="33%">chiCreateQuadrature()</td><td width="33%">chiCreateProductQuadrature()</td></tr>
</table>
### Mesh Utilities
 
<table>
<tr><td width="33%">chiEdgeLoopSplitByAngle()</td><td width="33%">chiLineMeshCreateFromLoop()</td><td width="33%">chiLineMeshCreateFromArray()</td></tr>
<tr><td width="33%">chiLogicalVolumeCreate()</td><td width="33%">chiMeshHandlerCreate()</td><td width="33%">chiMeshHandlerGetSurfaceFromCollection()</td></tr>
<tr><td width="33%">chiMeshHandlerSetCurrent()</td><td width="33%">chiRegionCreate()</td><td width="33%">chiRegionAddSurfaceBoundary()</td></tr>
<tr><td width="33%">chiRegionAddLineBoundary()</td><td width="33%">chiRegionGetBoundarySurfaceMesh()</td><td width="33%">chiRegionExportMeshToPython()</td></tr>
<tr><td width="33%">chiRegionExportMeshToObj()</td><td width="33%">chiRegionExportMeshToVTK()</td><td width="33%">chiSurfaceMeshCreate()</td></tr>
<tr><td width="33%">chiSurfaceMeshImportFromOBJFile()</td><td width="33%">chiSurfaceMeshImportFromTriangleFiles()</td><td width="33%">chiSurfaceMeshExportToObj()</td></tr>
<tr><td width="33%">chiSurfaceMeshExportPolyFile()</td><td width="33%">chiSurfaceMeshGetEdgeLoops()</td><td width="33%">chiSurfaceMeshGetEdgeLoopsPoly()</td></tr>
<tr><td width="33%">chiSurfaceMeshSplitByPatch()</td><td width="33%">chiSurfaceMeshExtractOpenEdgesToObj()</td><td width="33%">chiSurfaceMeshCheckCycles()</td></tr>
<tr><td width="33%">chiSurfaceMesherCreate()</td><td width="33%">chiSurfaceMesherExecute()</td><td width="33%">chiSurfaceMesherSetProperty()</td></tr>
<tr><td width="33%">chiSurfaceMesherExportToObj()</td><td width="33%">chiVolumeMesherCreate()</td><td width="33%">chiVolumeMesherExecute()</td></tr>
<tr><td width="33%">chiVolumeMesherSetProperty()</td><td width="33%">chiDomDecompose2D()</td><td width="33%">chiDecomposeSurfaceMeshPxPy()</td></tr>
</table>
### Field-function Manipulation
 
<table>
<tr><td width="33%">chiFFInterpolationCreate()</td><td width="33%">chiFFInterpolationSetProperty()</td><td width="33%">chiFFInterpolationInitialize()</td></tr>
<tr><td width="33%">chiFFInterpolationExecute()</td><td width="33%">chiFFInterpolationExportPython()</td><td width="33%">chiFFInterpolationGetValue()</td></tr>
</table>
### MPI Utilities
 
<table>
<tr><td width="33%">chiMPIBroadcastCellsets()</td><td width="33%">chiMPIReceiveCellsets()</td><td width="33%">chiMPIBarrier()</td></tr>
</table>
### Logging Utilities
 
<table>
<tr><td width="33%">chiLogSetVerbosity()</td><td width="33%">chiLog()</td></tr>
</table>
### Physics Utilities
 
<table>
<tr><td width="33%">chiSolverAddRegion()</td><td width="33%">chiSolverExecute()</td><td width="33%">chiPhysicsAddMaterial()</td></tr>
<tr><td width="33%">chiPhysicsMaterialAddProperty()</td><td width="33%">chiPhysicsMaterialSetProperty()</td><td width="33%">chiSolverAddFieldFunction()</td></tr>
<tr><td width="33%">chiGetFieldFunctionList()</td><td width="33%">chiExportFieldFunctionToVTK()</td><td width="33%">chiExportFieldFunctionToVTKG()</td></tr>
<tr><td width="33%">chiExportMultiFieldFunctionToVTKG()</td></tr>
</table>
### Monte-carlo Utilities
 
<table>
<tr><td width="33%">chiMonteCarlonCreateSolver()</td><td width="33%">chiMonteCarlonCreateSource()</td><td width="33%">chiMonteCarlonInitialize()</td></tr>
<tr><td width="33%">chiMonteCarlonExecute()</td><td width="33%">chiMonteCarlonSetProperty()</td></tr>
</table>
### Diffusion
 
<table>
<tr><td width="33%">chiDiffusionCreateSolver()</td><td width="33%">chiDiffusionInitialize()</td><td width="33%">chiDiffusionExecute()</td></tr>
<tr><td width="33%">chiDiffusionSetProperty()</td></tr>
</table>
### Linear Boltzman Solver
 
<table>
<tr><td width="33%">chiLBSransportCreateSolver()</td><td width="33%">chiLBSSetProperty()</td><td width="33%">chiLBSInitialize()</td></tr>
<tr><td width="33%">chiLBSExecute()</td><td width="33%">chiLBSGetFieldFunctionList()</td><td width="33%">chiLBSGetScalarFieldFunctionList()</td></tr>
</table>
### Linear Boltzman Solver - Groupset manipulation
\ref LuaLBSGroupsets Main page
<table>
<tr><td width="33%">chiLBSCreateGroup()</td><td width="33%">chiLBSCreateGroupset()</td><td width="33%">chiLBSGroupsetAddGroups()</td></tr>
<tr><td width="33%">chiLBSGroupsetSetQuadrature()</td><td width="33%">chiLBSGroupsetSetAngleAggDiv()</td><td width="33%">chiLBSGroupsetSetGroupSubsets()</td></tr>
<tr><td width="33%">chiLBSGroupsetSetIterativeMethod()</td><td width="33%">chiLBSGroupsetSetResidualTolerance()</td><td width="33%">chiLBSGroupsetSetMaxIterations()</td></tr>
<tr><td width="33%">chiLBSGroupsetSetGMRESRestartIntvl()</td><td width="33%">chiLBSGroupsetSetWGDSA()</td><td width="33%">chiLBSGroupsetSetTGDSA()</td></tr>
</table>


*/
