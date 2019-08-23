TA_PRESSURE               = 1
TA_BCPRESSURE             = 2
TA_TEMPERATURE            = 3
TA_BCTEMPERATURE          = 4
TA_VOIDFRACTION           = 5
TA_BCVOIDFRACTION         = 6
TA_VELOCITY               = 7
TA_JVELOCITY              = 8
TA_PRESSURETEMPERATURE    = 9
TA_BCPRESSURETEMPERATURE  = 10
TA_JLOSSCOEFFICIENTS      = 11

--========================================================= Create system
systemA = chiThermoCreateSystem();

--========================================================= Create control volumes
cv000 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv001 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv002 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv003 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv004 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv005 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv006 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv007 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv008 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});
cv009 = chiThermoCreateVolumeFromCoordinates(systemA,{x=4,y=5,z=6},{x=7,y=8,z=9});

chiThermoSetComponentProperty(systemA,cv000,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv001,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv002,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv003,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv004,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv005,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv006,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv007,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv008,TA_PRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,cv009,TA_PRESSURETEMPERATURE,100.0e3,293.15);

--========================================================= Create BCs
bc000 = chiThermoCreateBC(systemA);
bc001 = chiThermoCreateBC(systemA);
chiThermoSetComponentProperty(systemA,bc000,TA_BCPRESSURETEMPERATURE,100.0e3,293.15);
chiThermoSetComponentProperty(systemA,bc001,TA_BCPRESSURETEMPERATURE,100.0e3,293.15);

--========================================================= Create Junctions
sj900 = chiThermoCreateSJunction(systemA);
sj901 = chiThermoCreateSJunction(systemA);

sj000 = chiThermoCreateSJunction(systemA);
sj001 = chiThermoCreateSJunction(systemA);
sj002 = chiThermoCreateSJunction(systemA);
sj003 = chiThermoCreateSJunction(systemA);
sj004 = chiThermoCreateSJunction(systemA);
sj005 = chiThermoCreateSJunction(systemA);
sj006 = chiThermoCreateSJunction(systemA);
sj007 = chiThermoCreateSJunction(systemA);
sj008 = chiThermoCreateSJunction(systemA);

chiThermoSetComponentProperty(systemA,sj900,TA_JLOSSCOEFFICIENTS,1000.0,1000.0);

--========================================================= Nodalize
chiThermoConnectTwoComponents(systemA,bc000,sj900,cv000,0);
chiThermoConnectTwoComponents(systemA,cv000,sj000,cv001,0);
chiThermoConnectTwoComponents(systemA,cv001,sj001,cv002,0);
chiThermoConnectTwoComponents(systemA,cv002,sj002,cv003,0);
chiThermoConnectTwoComponents(systemA,cv003,sj003,cv004,0);
chiThermoConnectTwoComponents(systemA,cv004,sj004,cv005,0);
chiThermoConnectTwoComponents(systemA,cv005,sj005,cv006,0);
chiThermoConnectTwoComponents(systemA,cv006,sj006,cv007,0);
chiThermoConnectTwoComponents(systemA,cv007,sj007,cv008,0);
chiThermoConnectTwoComponents(systemA,cv008,sj008,cv009,0);
chiThermoConnectTwoComponents(systemA,cv009,sj901,bc001,0);

--========================================================= Initialize
chiThermoInitialize(systemA);
