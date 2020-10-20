--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

-- chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE, "one_grp.csx")

-- chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--         SIMPLEXS0,num_groups,0.1)

-- src={}
-- for g=1,num_groups do
--     src[g] = 0.0
-- end
-- src[1] = 1.0
-- chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

