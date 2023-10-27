--############################################### Setup mesh
-- This mesh is 2.0 x 2.0 offset to the right by 1.0
-- In cartesian coordinates,   V=4.0
-- In cylindrical coordinates, V=pi * (3.0^2 - 1.0^2) * 2.0 = 50.265482
nodes={}
N=10
L=2.0
xmin = 1.0
dx = L/N
for i=1,(N+1) do
  k=i-1
  nodes[i] = xmin + k*dx
end

meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes, nodes} })
chi_mesh.MeshGenerator.Execute(meshgen1)

chi_physics.FieldFunctionGridBased.Create
({
  name = "unit",
  sdm_type = "FV",
  coordinate_system = "cylindrical",
  initial_value = 1.0
})

chi.CellVolumeIntegralPostProcessor.Create
({
  name = "unit_pp_volume",
  field_function = "unit"
})
chi.ExecutePostProcessors({"unit_pp_volume"})
--trueV = math.pi * ((L+1.0)^2 - 1.0^2) * L
--compV = chi.PostProcessorGetValue("unit_pp_volume")
--print(trueV, compV, compV/trueV)