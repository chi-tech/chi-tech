/** \page Tutorial01 Tutorial 1: Basic Diffusion Simulation (No output)

 Let us tackle the basic diffusion equation as applied to steady state heat
 transfer.\n
 \f[
 \nabla k \nabla T = q
 \f]

 The purpose of this tutorial is to understand the basics of the
 input paradigm which is completely different from the traditional
 input-card based paradigm.
 ChiTech uses an imbedded Lua console. For a history of Lua please see its
 wikipedia page
 <a href="https://en.wikipedia.org/wiki/Lua_(programming_language)">here</a>.
 This paradigm allows the user to control many aspects of the program flow
 as well to be more versatile when creating inputs, i.e. reading files, for
  loops to create multiple entities etc. So let us get started.

 ## Step 1 - Create your input file, open it in a text editor

 Go to the directory of your
 choosing and make an input file of your choosing. The input file does not
 require a specific .lua or .inp extension although if you use .lua then
 there are many
 text editors that will be able to provide syntax highlighting.

## Step 2 - Create a mesh handler

\code
chiMeshHandlerCreate()
\endcode

This function call creates a new mesh handler and pushes it to a global
stack. It also make the newly created handler the "current" meaning all
 mesh operations will operate on this handler.
 The concept of creating, reading or manipulating computational meshes is
centered on operating on a specific chi_mesh::MeshHandler object. Physics/math
 objects can then freely operate on meshes by simply specifying which handler
 to use.

## Step 3 - Load a 2D surface mesh as template for extrusion

 The extrusion mesh generator requires a 2D mesh to generate the actual
 extrusion and this is how to load such a surface. For now we will use
 a test object included in ChiTech, which is an orthogonal grid of 32x32 cells
 (1024 total).


 \image html "Meshing/SquareMesh2x2Quads.png"  "Figure 1 - Mesh spanning from [-1,-1] to [1,1]" width=350px

\code
newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj")
\endcode

The first function call, CHI_LUA::chiSurfaceMeshCreate, returns a handle to an
 unpopulated surface mesh. This handle is essentially a stack index on the
 current mesh handler.
The second function call, CHI_LUA::chiSurfaceMeshImportFromOBJFile, receives
 as arguments the handle to the empty surface mesh and a filename pointing to
 a
 <a href="https://en.wikipedia.org/wiki/Wavefront_.obj_file">wavefront</a>
 ".obj" file.

<B>NOTE:</B> File names and paths are relative to where the executable is
 executed. i.e. If you are in your desktop folder then the path will be relative
 to it. Absolute paths can also be specified.

## Step 3 - Create a region and add the surface mesh as a boundary

A region is a weak concept for now but will later be used for transformation
 of meshes and boundaries. Boundaries get added to regions. A surface mesh
 boundary does not have a physical sense when using an extrusion.
 i.e It doesn't define the boundary on that surface. It gets assimilated
 during meshing process.

\code
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
\endcode

The first function call, CHI_LUA::chiRegionCreate, returns a handle to an
 empty Region similar to what we used with the surface mesh.
 The second function call, CHI_LUA::chiRegionAddSurfaceBoundary, receives
 as arguments a region handle and a surface mesh handle.

## Step 4 - Creating the meshers

Generating meshes in ChiTech follows the same paradigm as Star-CCM+. We first
 execute a surface meshing routing then a volume meshing routine. This is true
 whether you are doing a 1D, 2D or 3D simulation.

\code
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);
\endcode

A surface mesher is created using CHI_LUA::chiSurfaceMesherCreate and for our
purpose we are using a SURFACEMESHER_PREDEFINED which means we will be using
 predefined 2D cells and therefore will not be creating/modifying the surface
 mesh. The extruder is created with the call to CHI_LUA::chiVolumeMesherCreate
 with a type VOLUMEMESHER_EXTRUDER.

## Step 5 - Setup extrusion layers

Extrusion layers are needed by the extruder to define the 3D meshes. By default
a single layer of height 1.0 will be used provided that the user does not
 provide layer information.

\code
subdivisions=3
chiVolumeMesherSetProperty(EXTRUSION_LAYER,3.0,subdivisions,"My Layer");
\endcode

Extrusion layers are handled as a property of the extruder mesh and hence they
 are specified with a call to CHI_LUA::chiVolumeMesherSetProperty. It requires
 as arguments firstly the property index, which for this case is,
 EXTRUSION_LAYER, followed by the height of this layer (3.0), the number
 of sub-divisions, and a name for the layer.
 At this moment the name is of no
 physical meaning other than to help the user organize.

## Step 6 - Execute the meshers

By executing the meshers we create the geometry of the mesh. All material id's
 will be empty (-1) but the mesh can be visualized for error checking.

\code
chiSurfaceMesherExecute();
chiVolumeMesherExecute();
\endcode


## Step 7 - Assign material ID's

Materials ID's can conveniently be specified using logical volumes. There are
many options for logical volumes ranging from primitive parametric surfaces
 to non-convex user generated surfaces using surface meshes. For this tutorial
 we will use a rectangular paralellipiped (brick) as follows.

\code
material = chiPhysicsAddMaterial("Test Material");

vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)
\endcode

We first create a material using the CHI_LUA::chiPhysicsAddMaterial function.
 The handle of this material will essentially be zero but it is not technically
 required by the chiVolumeMesherSetProperty() function since this function
 operates on the integer supplied.
This material is added to the physics environment and therefore has scope over
all mesh handlers and all physics entities. We then create a logical volume
 using the function CHI_LUA::chiLogicalVolumeCreate with arguments RPP,
 specifying that we will be using a Rectangular Parallelipiped and then a
 series of dimensions specifying xmin-xmax-ymin-ymax-zmin-zmax. Every cell
 centroid within these dimensions will be flagged as being "within" the logical
 volume.

Actually setting the mesh's material id's is a utility facilitated by a
 volume mesher property and hence we call CHI_LUA::chiVolumeMesherSetProperty
 with a property index MATID_FROMLOGICAL. Next we provided a handle to the
 logical volume (vol0) and the desired material id. Logical volumes are very
 diverse and their uses are discussed elsewhere.

The culmination of this step is all done within a physics agnostic framework.
The user can even export the mesh for visualization using the function
 chiRegionExportMeshToObj(). The extruded mesh is shown below:


\image html "Meshing/SquareMeshExtruded.png"  "Figure 2 - Extruded mesh" width=350px

## Step 8 - Adding material properties

Now that the cells have been assigned a material id we need to add
 properties to the material conducive to a diffusion simulation. For a heat
 transfer diffusion simulation we will need to know the thermal conductivity
 "k" and the volumetric source strength "q". Both of these can be simple
 scalar values.

\code
chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
chiPhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)
\endcode

In this code we created the material properties using the function
 CHI_LUA::chiPhysicsMaterialAddProperty which requires a handle to the reference
 material, the property type (SCALAR_VALUE), and a name for the property.
 Unlike extrusion layers the property name can be used in further calls to refer
 to the specific property.

Material property values are set using the function
 CHI_LUA::chiPhysicsMaterialSetProperty which again expects a handle to the
 reference material, then either a material property id or name (in this case
 name), then an operation index and value(s). For this case we used an operation
 index SINGLE_VALUE which is the only operation supported by SCALAR_VALUE. In
 future the user can specify, as an example, temperature dependent values which
 will support the operation FROM_TABLE, but that is a topic for a different time.

## Step 9 - Setup the diffusion physics

The following sequence of function calls completely define the diffusion solver.

\code
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLC);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-4)
\endcode

We first create the diffusion solver with a call to
 CHI_LUA::chiDiffusionCreateSolver. This creates the solver and pushes it onto
 the physics handler. The function returns the handle.
 We then add our mesh region to this solver using the function
 CHI_LUA::chiSolverAddRegion which expects a handle to the reference solver
 followed by the handle to the relevant region.

Next we can set numerous diffusion solver properties which can comprehensively
 viewed in its specific documentation (CHI_LUA::chiDiffusionSetProperty).

## Step 10 - Initialize and Solve

The final step of this process is to initialize and execute the diffusion solver.

\code
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)
\endcode

## Complete input file

Here is the complete input file with comments

\code
--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj")

-- Setup Region
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);

-- Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_EXTRUDER);

-- Setup extrusion layers
subdivisions=3
chiVolumeMesherSetProperty(EXTRUSION_LAYER,3.0,subdivisions,"My Layer");

--  Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

material = chiPhysicsAddMaterial("Test Material");

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)

-- Export visualization
chiRegionExportMeshToObj(region1,"Tutorial01Mesh.obj",true)

--############################################### Add material properties


-- Set material properties
chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
chiPhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)


--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverAddRegion(phys1,region1)
chiDiffusionSetProperty(phys1,DISCRETIZATION_METHOD,PWLC);
chiDiffusionSetProperty(phys1,RESIDUAL_TOL,1.0e-4)

-- Initialize and Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)
\endcode


## Step 11 - Execute the code
Assuming you added the executable to your PATH environment variable, the code
can be executed by typing the executable name followed by the input file path
 (relative or absolute).

\verbatim
ChiTech Tutorial01.lua
\endverbatim

The output produced will look as follows:

\verbatim
[0]  2020-01-08 13:30:30 Running ChiTech in batch-mode with 1 processes.
[0]  ChiTech number of arguments supplied: 1
[0]  Surface mesh loaded with 0 triangle faces and 1024 polygon faces.
[0]  00:00:00 VolumeMesherExtruder executed. Memory in use = 12.3516 MB
[0]  VolumeMesherExtruder: Total number of cell layers is 3
[0]  VolumeMesherExtruder: Extruding cells
[0]  VolumeMesherExtruder: Cells extruded = 3072
[0]  VolumeMesherExtruder: Number of cells in region = 3072
[0]  VolumeMesherExtruder: Number of nodes in region = 4356
[0]  00:00:00 chiVolumeMesherExecute: Volume meshing completed. Memory used = 3.49 MB
[0]  Total process memory used after meshing 15.8 MB
[0]  00:00:00 Setting material id from logical volume.
[0]  00:00:00 Done setting material id from logical volume. Number of cells modified = 3072.
[0]  Exported Material Volume mesh to Tutorial01Mesh_m0.obj
[0]
[0]  00:00:00 Diffusion Solver: Initializing Diffusion solver PETSc
[0]  Computing nodal reorderings for CFEM
[0]  Time taken during nodal reordering 0.00318
[0]  Determining nodal connections
[0]  Building sparsity pattern.
[0]  Setting matrix preallocation.
[0]  Computing cell matrices
[0]  00:00:01 Diffusion Solver: Diffusion Solver initialization time 1.47667
[0]  Diffusion Solver: Assembling A and b
[0]  Diffusion Solver: Local matrix instructions
[0]  Diffusion Solver: Communicating matrix assembly
[0]  Diffusion Solver: Solving system
[0]  Diffusion iteration    0 - Residual 34.8524639
[0]  Diffusion iteration    1 - Residual 2.6205789
[0]  Diffusion iteration    2 - Residual 0.2875545
[0]  Diffusion iteration    3 - Residual 0.0372920
[0]  Diffusion iteration    4 - Residual 3.783e-03
[0]  Diffusion iteration    5 - Residual 3.774e-04
[0]  Diffusion iteration    6 - Residual 4.446e-05
[0]  Convergence reason: 2
[0]  Diffusion Solver: Number of iterations =6
[0]  Timing:
[0]  Assembling the matrix: 0.017682
[0]  Solving the system   : 0.01209
[0]  Diffusion Solver execution completed!
[0]  Final program time 00:00:01
[0]  2020-01-08 13:30:31 ChiTech finished execution of Tutorials/Tutorial01.lua
\endverbatim


 */