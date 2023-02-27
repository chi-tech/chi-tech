-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.51187 and 1.42458e-03
num_procs = 4
--Unstructured mesh




--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromVTU("ZPhi3D_g0_0.vtu")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh);

chiSurfaceMesherExecute();
chiVolumeMesherExecute();

--############################################### Exports
if master_export == nil then
    chiMeshHandlerExportMeshToVTK("ZVTU")
end
