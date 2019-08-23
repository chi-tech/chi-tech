PropertyClass = {}
PropertyClass.__index = PropertyClass;
PropertyCount = 0;

dofile(chinPropertyDir.."chin_property_00_constrdestr.lua")
dofile(chinPropertyDir.."chin_property_01_functions.lua")

--################ Design thoughts
--Properties are elements one level above features and 3D elements. They are 
--copy-paste-able, drag-and-drop-able and can be instanced. They can also
--be processed inside a node editor. They can be of different data formats:
--*Scalars
--*Vectors (size 2,3,4)
--*Matrices (size n)
--*Text values
--
--They can have different units (C,F, deg, rad, etc.), preferred units can
--be set.
--Bi-directional links can be assigned to them allowing multiple features
--to change them or receive changes from them.