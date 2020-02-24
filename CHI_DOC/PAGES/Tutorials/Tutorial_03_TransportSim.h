/** \page Tutorial03 Tutorial 3: Basic transport simulation

Let us tackle a basic transport simulation on a 1D mesh\n

\f[
 \hat{\Omega} \nabla \Psi(\hat{\Omega}) + \sigma_t \Psi (\hat{\Omega}) =
 \int_{\hat{\Omega}} \sigma_s (\hat{\Omega}' \to \hat{\Omega})
 \Psi (\hat{\Omega}').d\hat{\Omega}' + q(\hat{\Omega})
\f]

## Step 1 - Create your input file, open it in a text editor

 As before, go to the directory of your
 choosing and make an input file of your choosing.

## Step 2 - Create a mesh

A 1D simulation is created from a Line-mesh. The line mesh expects to be
 populated by a number of points supplied by a lua table.

\code
mesh={}
N=5000
L=30.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
line_mesh = chiLineMeshCreateFromArray(mesh)
\endcode

As soon as we have a line mesh we can now proceed as before.

\code
region1 = chiRegionCreate()
chiRegionAddLineBoundary(region1,line_mesh);

-- Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_LINEMESH1D);

-- Execute meshing
chiSurfaceMesherExecute();
chiVolumeMesherExecute();

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)
\endcode

Some notably new items is the different volume mesher (VOLUMEMESHER_LINEMESH1D).
 Note also that even though there is no notion of a surfacemesher we still
 need to use a pass-through SurfaceMesher. This helps to keep the coding
 pipeline simple.

We also set the material id's the same we we do any mesh format.

## Step 3 - Add a material with a transport cross-section

Similar to the diffusion tutorial we have to create a material, then add a
 property to it, and then to set the property.

\code
material0 = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)

chiPhysicsMaterialSetProperty(material0,
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,
                              1,     --Num grps
                              1.0,   --Sigma_t
                              1.0)   --Scattering ratio
\endcode

The general material property TRANSPORT_XSECTIONS is used for
 non-fission transport materials. Simple cross-sections can be specified by
 using the operation SIMPLEXS1 which expects the number of groups,
 the total cross-section, and then the scattering ratio.

## Step 4 - Add Transport physics

\code
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

chiLBSCreateGroup(phys1)

pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,32)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
chiLBSGroupsetAddGroups(phys1,gs0,0,0)
chiLBSGroupsetSetQuadrature(phys1,gs0,pquad)

--========== Boundary conditions
bsrc = {0.5}
chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
                  ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD1D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)
\endcode

 A transport solver is invoked by using a call to chiLBSCreateSolver().
 This creates a derived physics solver so the mesh region gets added to the
 solver using the generic call chiSolverAddRegion(). Past this point we need
 to create the single required group with chiLBSCreateGroup() and then a
 quadrature rule for integration of the angular fluxes. Since we are dealing
 with a 1D simulation we will be integration over theta translating to a
 cosine from -1 to 1 and therefore the appropriate quadrature rule would be a
 Gauss-Legendre quadrature rule. For this we use a call to
 chiCreateProductQuadrature() and specifying GAUSS_LEGENDRE as the rule and
 32 number of discrete angles.

 The next step in the process is to assign a group-set. Group-sets are very
 useful aggregation features in higher dimension simulations but here we
 only have a single groupset. The group-set is created with a call to
 chiLBSCreateGroupset(). Next we add groups to the group-set using a range,
 however, since we only have one group here the range will be 0 to 0. The
 final piece of a groupset is to add a quadrature which is achieved with a
 call to chiLBSGroupsetSetQuadrature().


## Step 5 - Initialize and Execute

\code
chiLBSInitialize(phys1)
chiLBSExecute(phys1)
\endcode

This should be intuitive.

## Step 6 - Add output

\code
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, L+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])


chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)

if (chi_location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end
\endcode

Instead of using a SLICE interpolator we instead opt to use a LINE interpolator.
 A line-interpolator needs two points as the initial and final point (i.e.
 LINE_FIRSTPOINT and LINE_SECONDPOINT) as well as the total number of
 points LINE_NUMBEROFPOINTS. Finally add the first available scalar field
 field function to the interpolator.

 The output should look as follows:

  \image html "Physics/Tut3Output.png"  "Figure 1 - Output of the 1D simulation" width=350px




*/