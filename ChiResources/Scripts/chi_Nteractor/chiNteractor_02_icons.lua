-------------------------------------------------------------------------------- Icon Global Variables
--========================================================= Begin GlobVar
chinGlobal.Icons = {};


--================================= Type ____ GlobalVar

--============== Icon Information
chinGlobal.Icons.iconTextureSize            = 1024;                                                                      -- Icon Texture Size
chinGlobal.Icons.iconsCutSize               = 32;                                                                       -- The cutout size
chinGlobal.Icons.iconPadding                = 1;                                                                        -- Padding inbetween icons
chinGlobal.Icons.iconOffset                 = chinGlobal.Icons.iconPadding*1 + chinGlobal.Icons.iconsCutSize;           -- Offset of icons starting from 0x0
chinGlobal.Icons.iconScale                  = chinGlobal.Icons.iconsCutSize/chinGlobal.Icons.iconTextureSize            -- Scales the icon accordingly
chinGlobal.Icons.iconShift                  = chinGlobal.Icons.iconOffset/chinGlobal.Icons.iconTextureSize;             -- Shift in texture location
chinGlobal.Icons.iconTranslation            = chinGlobal.Icons.iconPadding/chinGlobal.Icons.iconTextureSize;             -- Shift in texture location


--============== Icon Sizing
chinGlobal.Icons.iconSizeSmall              = 16;                                                                       -- Icon Size Small
chinGlobal.Icons.iconSizeMedium             = 25;                                                                       -- Icon Size Medium
chinGlobal.Icons.iconSizeLarge              = 35;                                                                       -- Icon Size Large
chinGlobal.Icons.iconSizeXlarge             = 45;                                                                       -- Icon Size X Large
chinGlobal.Icons.iconSizeCustom             = 79;                                                                       -- Icon Size Custom




--================================= Type ____ GlobalVar

chinIconExpander                            = {0,0};
chinIconReducer                             = {1,0};
chinIconFolder                              = {2,0};
chinIconTexture                             = {3,0};
chinIconObject                              = {4,0};
chinIconCheckboxU                           = {5,0};
chinIconCheckboxC                           = {6,0};
chinIconMaterial                            = {7,0};
chinIconTransform                           = {8,0};
chinIconArrowUp                             = {9,0};
chinIconArrowDown                           = {10,0};
chinIconArrowLeft                           = {11,0};
chinIconArrowRight                          = {12,0};
chinIconRedBall                             = {13,0};
chinIconFolderBlue                          = {14,0};
chinIconAtom                                = {15,0};
chinIconNothing                             = {16,0};
chinIconGrnBall                             = {16,0};
chinIconNothing                             = {17,0};
chinIconBlueBall                            = {18,0};
chinIconStep                                = {19,0};
chinIconRun                                 = {20,0};
chinIconStop                                = {21,0};
chinIconIni                                 = {22,0};


