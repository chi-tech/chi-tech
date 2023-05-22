phys0 = prk.TransientSolver.Create(
{
  initial_rho = 0.0,
  initial_source = 0.0,
  dt = 0.01,
})

chiSolverInitialize(phys0)

for t=1,200 do
  chiSolverStep(phys0)
  time = prk.GetParam(phys0, "time_next")
  print(time,
    prk.GetParam(phys0, "population_next"),
    prk.GetParam(phys0, "period"),
    prk.GetParam(phys0, "period"))

  chiSolverAdvance(phys0)
  if (time > 0.1) then
    prk.SetParam(phys0, "rho", 0.21)
  end
end
