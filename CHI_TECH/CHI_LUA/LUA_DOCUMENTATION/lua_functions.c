/** Sets the property of the text
\param handle int Handle to the text object in question.
\param property int Property number.
\param value varying Value of the property to set
Property numbers\n
 TEXT_VALUE    = text\n
 TEXT_POSITION = position requires float x, float y, float z\n
 TEXT_COLOR    = color requires float r, float g, float b, float a\n
 TEXT_FONT     = font int font index\n
 TEXT_SCALE    = scale requires float xscale, float yscale, float zscale\n
\n
Note if only xscale is provided then it is uniformly applied to yscale and zscale.
\ingroup LuaText
\author Jan */
int CHI_LUA::chiTextSetProperty(int handle, int property, varying value)
{return;} 
/** Creates a text object.
\param name char* Generic name for the object. Doesn't really do anything.
\return handle int Handle to the text object.
\ingroup LuaText
\author Jan */
int CHI_LUA::chiTextCreate(char* name)
{return;} 
/** \brief Assigns a transform to the lines.
\param index  int Index of the line.
\param transformIndex  int Index of the transform.
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineSetTransform(int index, int transformIndex)
{return;} 
/** \brief Sets the viewport for the lines.
\param index  int Index of the line.
\param viewportIndex  int Index of the viewport.
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineSetviewport(int index, int viewportIndex)
{return;} 
/** \brief Adjusts line color.
\param index  int Index of the line.
\param r float R Color component.
\param g float G Color component.
\param b float B Color component.
\param a float A Color component.
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineChangeColor(int index, float r, float g, float b, float a)
{return;} 
/** \defgroup LuaLine3D 3D Lines
 * Hallo*/
/** \brief Adds a vertex to the line vertex stack.
\param index  int Index of the line.
\param x float X-location.
\param y float Y-location.
\param z float Z-location.
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineAddVertex(int index, float x, float y, float z)
{return;} 
/** \brief Adjusts a specific vertex position.
\param index  int Index of the line.
\param vertNum int Vertex number.
\param x float X-location.
\param y float Y-location.
\param z float Z-location.
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineChangeVertex(int index, int vertNum, float x, float y, float z)
{return;} 
/** \brief Adjusts stipple options including line width.
\param index  int Index of the line.
\param stippleFlag bool true=Enabled, false=Disabled
\param stippleFactor int Stipple factor.
\param lineWidth float Line width.
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineSetStipple(int index, bool stippleFlag, int stippleFactor, float lineWidth)
{return;} 
/** \brief Creates a 3D line.
\param name   char*  Line name
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineCreate(char* name)
{return;} 
/** \brief Creates a 3D line from a chi_mesh::LineMesh.
\param name   char*  Line name
\ingroup LuaLine3D
\author Jan*/
int CHI_LUA::chi3DLineCreateFromLineMesh(char* name)
{return;} 
/** chiLightSetProperty(char* lightName, char* propertyName, variant value).
- lightName. Name of the object.
- propertyName. Name of the property.
- value. Value the property should take.
\ingroup LuaLights
\author Jan*/
int CHI_LUA::chiLightSetProperty()
{return;} 
/**\ingroup LuaLights Lights*/
/** chiLightCreate(char* name). Creates a new light.
 \ingroup LuaLights
\author Jan*/
int CHI_LUA::chiLightCreate()
{return;} 
/**\defgroup LuaShadows Shadows*/
/** Simple function to turn on shadows.
 *
\ingroup LuaShadows
\author Jan*/
int CHI_LUA::chiShadowsEnable()
{return;} 
/** Simple function to turn on shadows.
\ingroup LuaShadows
\author Jan*/
int CHI_LUA::chiShadowsDisable()
{return;} 
/** \brief Adds a vertex to the line vertex stack.
\param index  int Index of the line.
\param x float X-location.
\param y float Y-location.
\param z float Z-location.
\ingroup LuaLine
\author Jan*/
int CHI_LUA::chiLineAddVertex(int index, float x, float y, float z)
{return;} 
/** \brief Adjusts line color.
\param index  int Index of the line.
\param r float R Color component.
\param g float G Color component.
\param b float B Color component.
\param a float A Color component.
\ingroup LuaLine
\author Jan*/
int CHI_LUA::chiLineChangeColor(int index, float r, float g, float b, float a)
{return;} 
/** \brief Adjusts stipple options including line width.
\param index  int Index of the line.
\param stippleFlag bool true=Enabled, false=Disabled
\param stippleFactor int Stipple factor.
\param lineWidth float Line width.
\ingroup LuaLine
\author Jan*/
int CHI_LUA::chiLineSetStipple(int index, bool stippleFlag, int stippleFactor, float lineWidth)
{return;} 
/** \brief Sets the viewport for the lines.
\param index  int Index of the line.
\param viewportIndex  int Index of the viewport.
\ingroup LuaLine
\author Jan*/
int CHI_LUA::chiLineSetviewport(int index, int viewportIndex)
{return;} 
/** \defgroup LuaLine 2D Lines
 * Hallo*/
/** \brief Creates a 2D line.
\param name   char*  Line name
\ingroup LuaLine
\author Jan*/
int CHI_LUA::chiLineCreate(char* name)
{return;} 
/** \brief Adjusts a specific vertex position.
\param index  int Index of the line.
\param vertNum int Vertex number.
\param x float X-location.
\param y float Y-location.
\param z float Z-location.
\ingroup LuaLine
\author Jan*/
int CHI_LUA::chiLineChangeVertex(int index, int vertNum, float x, float y, float z)
{return;} 
/** chiCreateObject(char* name). Loads a surface mesh from a file.
- name. Name of the object (null by default).
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectCreate()
{return;} 
/**\defgroup LuaObjects Objects*/
/** Exports a graphical surface to a wavefront .obj file.
\param SurfaceHandle int Handle to the surface to be exported.
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectExportSurface(int SurfaceHandle)
{return;} 
/** chiLoadSurface(char* fileName). Loads a surface mesh from a file.
- fileName. Name of the .obj file that contains the surface mesh.
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectLoadSurface()
{return;} 
/** Loads a surface mesh from a chi_mesh::SurfaceMesh object.
\param Handle int Handle to the surface mesh.
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectLoadSurfaceFromSurfaceMesh(int Handle)
{return;} 
/**\defgroup LuaObjects Objects*/
/** chiObjectAddSurface(int objectIndex, int surfaceIndex). Loads a surface mesh from a file.
- objectIndex. Numerical index of the object in the object stack.
- surfaceIndex. Numerical index of the surface in the surfaceMesh stack.
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectAddSurface()
{return;} 
/**chiObjectQuery(int query, var param1, ..., var param N). Queries the
object stack for miscellaneous properties.
 \ingroup LuaObjects
\author G-mo*/
int CHI_LUA::chiObjectQuery()
{return;} 
/** chiObjectSetProperty(char* objectName, char* propertyName, variant value). Creates a material with default properties.
- objectName. Name of the object.
- propertyName. Name of the property.
- value. Value the property should take.
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectSetProperty()
{return;} 
/** chiObjectGetCentroid(char* surfaceName). Determines the centroid of a surface.
Returns x,y,z.
\ingroup LuaObjects
\author Jan */
int CHI_LUA::chiObjectGetCentroid()
{return;} 
/**\defgroup LuaTextures Textures */
/** chiTextureLoad(char* filePath). Loads a texture to the graphics
 * environment
- filePath. 	Path to the file.
\ingroup LuaTextures
\author Jan*/
int CHI_LUA::chiTextureLoad()
{return;} 
/** \defgroup Lua_General General
*/
/** Implements the sleep function in lua
\param time int Time to sleep in milliseconds.
\ingroup Lua_General
\author Jan*/
int CHI_LUA::chiSleep(int time)
{return;} 
/** \defgroup LuaDisplayerControl Displayer Control
 * */
/**Creates a new displayer and adds it to the currently selected scene.
\ingroup LuaDisplayerControl
\return int Displayer handle.
\author Jan*/
int CHI_LUA::chiDisplayerCreate()
{return;} 
/**Sets the viewport coordinates of the currently active displayer.
\param xmin double
\param ymin double
\param xmax double
\param ymax double
\ingroup LuaDisplayerControl
\author Jan*/
int CHI_LUA::chiDisplayerSetViewport(double xmin, double ymin, double xmax, double ymax)
{return;} 
/**Gets the viewport coordinates of the currently active displayer.
 *
\return xmin,ymin,xmax,ymax float Viewport dimensions.
\ingroup LuaDisplayerControl
\author Jan*/
int CHI_LUA::chiDisplayerGetViewport()
{return;} 
/** \defgroup LuaSceneControl Scene control
 * Hallo*/
/** \brief Simple function to enable the 2D scene.
\param flag bool Flag choice for 2D option. true=scene is 2D.
\ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiSet2D(bool flag)
{return;} 
/** \brief Simple function to enable a scene's 3D content.
 *
 \param flag bool Flag choice for 3D option. true=scene is 3D.
 *
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiSet3D(bool flag)
{return;} 
/** \brief Binds the scene for new object creation.
 *
\param number int Scene number to bind.
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiBindScene(int number)
{return;} 
/** \brief Gets the currently bound scene.
 *
 * \return sceneNumber int The number of the scene currently bound.
 * \return displayer   int Displayer number of the displayer currently selected.
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiGetScene()
{return;} 
/** \brief Gets the number of scene currently available.
 *
 * \return sceneCount int the number of scenes available.
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiGetSceneCount()
{return;} 
/** \brief Creates a new scene.
 *
 * \return sceneNumber int The number index of the newly created scene.
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiCreateScene()
{return;} 
/** \brief Requests a refresh of the scene.
 *
\param sceneNumber int Index of the scene for which the request is made.
 *
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiRequestSceneRefresh(int sceneNumber)
{return;} 
/** \brief Set the refresh mode of a scene.
 *
\param sceneMode int Update mode of the scene.
 * \ingroup LuaSceneControl
\author Jan*/
int CHI_LUA::chiSetSceneUpdateMode(int sceneMode)
{return;} 
/** bool chiViewportSetProperty() Set viewport dimensions.
 * Creates a viewport.
\param xmin   double  Viewport x min
\param ymin   double  Viewport y min
\param xmax   double  Viewport x max
\param ymax   double  Viewport y max
\ingroup LuaViewports
\author Gmo */
int CHI_LUA::chiViewportSetProperty(double xmin, double ymin, double xmax, double ymax)
{return;} 
/** \defgroup LuaViewports Viewports*/
/** int chiViewportCreate() Creates a viewport.
 *
\param  name    str Assigns a name to the viewport
 
\return viewportIndex int Index of the newly created viewport.
\ingroup LuaViewports
\author Gmo */
int CHI_LUA::chiViewportCreate(str name)
{return;} 
/** chiGraphicsTranslateCamera(char* name). Allows the camera to change its orientation.
\ingroup LuaCameras
\author Jan, GMO */
int CHI_LUA::chiGraphicsOrientCamera()
{return;} 
/** chiGraphicsCameraType(float value). Sets camera type (1=Perspective,!1=Orthographic).
\ingroup LuaCameras
\author Jan */
int CHI_LUA::chiGraphicsCameraType()
{return;} 
/** chiGraphicsTranslateCamera(char* name). Manual changes the camera's rotation.
\ingroup LuaCameras
\author Jan, GMO */
int CHI_LUA::chiGraphicsRotateCamera()
{return;} 
/** chiGraphicsGetCameraPosition(char* name). Retrieves the camera's rotation.
\author GMO */
int CHI_LUA::chiGraphicsGetCameraRotation()
{return;} 
/** chiGraphicsTranslateCamera(char* name). Allows the translation of the camera according to the normal.
\ingroup LuaCameras
\author Jan, GMO */
int CHI_LUA::chiGraphicsTranslateCamera()
{return;} 
/** chiGraphicsCameraOrthoWidth(float value). Sets camera orthographic width.
\param orthoWidth float       Width of the orthographic window.
\param orthoRound=false bool  Flag indicating if the orthoWidth should be rounded.
\note If the camera is part of a setup for displaying text then orthoRound should be set to true.
\ingroup LuaCameras
\author Jan */
int CHI_LUA::chiGraphicsCameraOrthoWidth(float orthoWidth, bool orthoRound=false)
{return;} 
/**\defgroup LuaCameras Cameras*/
/** chiGraphicsCameraFOV(char* name). Sets camera FOV.
\ingroup LuaCameras
\author Jan */
int CHI_LUA::chiGraphicsCameraFOV()
{return;} 
/** chiGraphicsTranslateCamera(char* name). Allows for a camera to be set.
\ingroup LuaCameras
\author Jan, GMO */
int CHI_LUA::chiGraphicsPositionCamera()
{return;} 
/** chiGraphicsGetCameraPosition(char* name). Retrieves the camera's position.
\author GMO */
int CHI_LUA::chiGraphicsGetCameraPosition()
{return;} 
/** \defgroup LuaWindowControl Window control
 * Hallo*/
/** \brief <B>chiWindowCreate()</B>
 * Creates a new window with a new scene.
\param name char Name of the form.
\param defer bool Flag=true Creates a window in a deferred sense.
\return newSceneNumber int Handle to the new scene
\ingroup LuaWindowControl
\author Jan*/
int CHI_LUA::chiWindowCreate(char name, bool defer)
{return;} 
/** chiSetWindowProperties(int xSize, int ySize, [int xPos], [int yPos]). Sets the window size for the
 * current scene. This does not change any displayers.
- xSize, the width of the window in pixels.
- ySize, the height of the window in pixels.
- xPos,  the x position of the window in pixels. (If omitted, current x location is used)
- yPos,  the y position of the window in pixels. (If omitted, current y location is used)
\ingroup LuaWindowControl
\author Jan */
int CHI_LUA::chiSetWindowProperties()
{return;} 
/** Sets the scene for the selected window.
\author Jan */
int CHI_LUA::chiSetWindowScene()
{return;} 
/** chiGetWindowProperties(). Gets the window size.
\author Jan */
int CHI_LUA::chiGetWindowProperties()
{return;} 
/** Maximizes a window.
\author Jan.*/
int CHI_LUA::chiWindowMaximize()
{return;} 
/** Normal displays a window (i.e. Not-maximized or Not-minimized).
\author Jan.*/
int CHI_LUA::chiWindowNormalize()
{return;} 
int CHI_LUA::chiWindowSetCursor()
{return;} 
/** \defgroup LuaText Text control
 * Hallo*/
/** chiSetLabel(char* name, char* string,int x,int y, [color_r,color_g,color_b,scale,font]).
Sets a label to be displayed.
*/
int CHI_LUA::chiSetLabel()
{return;} 
/** chiSetLabel(char* name, char* string,int x,int y, [color_r,color_g,color_b,scale,font]).
Sets a label to be displayed.
*/
int CHI_LUA::chiSetLabelProperty()
{return;} 
/**chiObjectQuery(int query, var param1, ..., var param N). Queries the
object stack for miscellaneous properties.
\param queryNumber int Number of the query.
\ingroup LuaTransforms
\author G-mo*/
int CHI_LUA::chiTransformQuery(int queryNumber)
{return;} 
/** chiTransformGet(). Get Transform parameters.
\param transformIndex int Index number of the transform.
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformGet(int transformIndex)
{return;} 
/** \defgroup LuaTransforms Transforms*/
/** CreateTransform(char* name). Creates a transform and places it on the tool transform-stack.
\param name char* Name of the transform.
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformCreate(char* name)
{return;} 
/** chiTransformSetScale(). Sets Transform scale
\param transformHandle int Handle to the transform
\param dx double Scale factor x
\param dy double Scale factor y
\param dz double Scale factor z
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformSetScale(int transformHandle, double dx, double dy, double dz)
{return;} 
/** chiTransformSetRotationPoint(char*\int transformName, double dx, double dy, double dz). Sets Transform rotation point.
\param transformHandle int Handle to the transform
\param dx double Rotation point x
\param dy double Rotation point y
\param dz double Rotation point z
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformSetRotationPoint(int transformHandle, double dx, double dy, double dz)
{return;} 
/** chiTransformSetRotation(char*\int transformName, double dx, double dy, double dz). Sets Transform rotation
\param transformHandle int Handle to the transform
\param dx double Rotation in degrees around axis
\param dy double Rotation in degrees around axis
\param dz double Rotation in degrees around axis
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformSetRotation(int transformHandle, double dx, double dy, double dz)
{return;} 
/** chiTransformSetParent() Sets the parent of an equation.
\param childIndex int Index of the child.
\param parentIndex int Index of the parent.
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformSetParent(int childIndex, int parentIndex)
{return;} 
/** chiTransformSetTranslation(char*\int transformName, double dx, double dy, double dz). Sets Transform translation
\param transformHandle int Handle to the transform
\param dx double Translation distance in x
\param dy double Translation distance in y
\param dz double Translation distance in z
\ingroup LuaTransforms
\author Jan*/
int CHI_LUA::chiTransformSetTranslation(int transformHandle, double dx, double dy, double dz)
{return;} 
/**chiMaterialQuery(int query, var param1, ..., var param N). Queries the
material stack for miscellaneous properties.
\ingroup LuaMaterials
\author Jan*/
int CHI_LUA::chiMaterialQuery()
{return;} 
/** chiMaterialGeneratePreview().
 *
\ingroup LuaMaterials
\author G-mo*/
int CHI_LUA::chiMaterialGeneratePreview()
{return;} 
/**\defgroup LuaMaterials Materials*/
/** chiCreateMaterial(char* name). Creates a material with default properties.
- name. Name of the material (null by default).
\ingroup LuaMaterials
\author Jan*/
int CHI_LUA::chiMaterialCreate()
{return;} 
/** chiMaterialSetProperty(char* materialName, char* propertyName, variant value). Creates a material with default properties.
- materialName. Name of the material.
- propertyName. Name of the property.
- value. Value the property should take.
### Properties available
|Property                  |Description [ValueType]                    |
|--------------------------|--------------------------------|
|DiffuseTexture			   |Name of the texture (normally the actual filepath) for the Diffuse texture. [char*]    	|
|AmbientTexture			   |Name of the texture (normally the actual filepath) for the Ambient texture.  [char*]  	|
|SpecularTexture		   |Name of the texture (normally the actual filepath) for the Specular texture. [char*]    |
|NormalTexture			   |Name of the texture (normally the actual filepath) for the Normal texture. [char*]   	|
|EnvironmentTexture		   |Name of the texture (normally the actual filepath) for the Environment texture. [char*]|
|EmissiveTexture		   |Name of the texture (normally the actual filepath) for the Emissive texture.   [char*] |
\ingroup LuaMaterials
\author Jan & Reptation*/
int CHI_LUA::chiMaterialSetProperty()
{return;} 
/**chiMaterialUpdate(int matlIndex, luatable materialTable). Assigns the
indicated material properties to the given table.
\ingroup LuaMaterials
\author Jan*/
int CHI_LUA::chiMaterialUpdate()
{return;} 
/** Converts a 3D point to window space.
x,y,z
returns x,y;
 0  1  2  3
 4  5  6  7
 8  9 10 11
12 13 14 15
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chi3DPointToScreen()
{return;} 
/** \defgroup LuaPoints Manipulating points */
/** Creates a point collection.
name    Line name
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chiPointCreate()
{return;} 
/** Adds a vertex to the point collection vertex stack.
name     Line name
float    x
float    y
float    z
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chiPointAddVertex()
{return;} 
/** Adjusts a specific vertex position.
name     Line name
int	     vertex number
float    x
float    y
float    z
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chiPointChangeVertex()
{return;} 
/** Converts a 3D point to window space.
x,y,z
returns x,y;
 0  1  2  3
 4  5  6  7
 8  9 10 11
12 13 14 15
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chiScreenTo3DPoint()
{return;} 
/** Creates a point collection.
name    Line name
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chi3DPointCreate()
{return;} 
/** Adds a vertex to the point collection vertex stack.
name     Line name
float    x
float    y
float    z
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chi3DPointAddVertex()
{return;} 
/** Adjusts a specific vertex position.
name     Line name
int	     vertex number
float    x
float    y
float    z
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chi3DPointChangeVertex()
{return;} 
/** Assigns a transform to the line
name     Line name
name     Transform name
\ingroup LuaPoints
\author Jan*/
int CHI_LUA::chi3DPointSetTransform()
{return;} 
/** Reloads the shaders
 *
 \ingroup LuaShaders
\author Gmo*/
int CHI_LUA::chiReloadShader()
{return;} 
/**\defgroup LuaShaders Shaders*/
/**chiObjectQuery(int query, var param1, ..., var param N). Queries the
object stack for miscellaneous properties.
 \ingroup LuaShaders
\author G-mo*/
int CHI_LUA::chiShaderQuery()
{return;} 
/** \defgroup LuaPie Raspberry Pie
 * Hallo*/
/** \brief Exports a GPIO pin.
\param pinNumber int GPIO pin number to be exported.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieExportPin(int pinNumber)
{return;} 
/** \brief Sets a pin to READ or WRITE.
\param pinNumber int GPIO pin number.
\param mode int 0=Read, 1=Write.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieSetPinMode(int pinNumber, int mode)
{return;} 
/** \brief Sets a pin to HIGH or LOW.
\param pinNumber int GPIO pin number.
\param mode int 0=Low, 1=High
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieSetPinValue(int pinNumber, int mode)
{return;} 
/** \brief Reads a pin.
\param pinNumber int GPIO pin number.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieGetPinValue(int pinNumber)
{return;} 
/** \brief Initializes SPI.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieInitSPI()
{return;} 
/** \brief Reads an SPI MCP3008 chip unbuffered.
\param channelNumber int Channel 0-7 to be read off the MCP3008 chip.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieReadSPIChannel(int channelNumber)
{return;} 
/** \brief Sets an SPI MCP3008 chip to either be buffered or not.
\param channelNumber int Channel 0-7 to be read off the MCP3008 chip.
\param bufferFlag bool true=Is buffered, false= Not buffered.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieSetSPIBuffer(int channelNumber, bool bufferFlag)
{return;} 
/** \brief Reads an SPI MCP3008 chip buffer.
\param channelNumber int Channel 0-7 to be read off the MCP3008 chip.
\param bufferPos int Buffer position (0 to 999).
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieGetSPIBuffer(int channelNumber, int bufferPos)
{return;} 
/** \brief Initializes the Pi's UART serial communication.
\param baudrate   int Can be any of the <I>Serial Options</I>.
### Serial options
PI3_BAUD_1200  0
PI3_BAUD_2400  1
PI3_BAUD_4800  2
PI3_BAUD_9600  3
PI3_BAUD_19200 4
PI3_BAUD_38400 5
\example
CHI_PI3 chipie;
chipie.InitializeSerial(PI3_BAUD_9600);
\endexample
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieInitializeSerial(int baudrate)
{return;} 
/** \brief Writes a message to the serial port.
\param message char Message to be sent.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieSerialWrite(char message)
{return;} 
/** \brief Reads a message from the serial port.
\ingroup LuaPie
\author Jan*/
int CHI_LUA::chiPieSerialRead()
{return;} 
/** bool chiThermoSetComponentProperty() Sets the property of a component.
 *
\param sysHndle  Handle to the system being referenced.
\param compHndle Handle to the component being referenced.
\param propCode  Property code.
 
 \return Success=true
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoSetComponentProperty(Handle sysHndle, Handle compHndle, Property propCode)
{return;} 
/** bool chiThermoGetComponentProperty() Sets the property of a component.
 *
\param sysHndle  Handle to the system being referenced.
\param compHndle Handle to the component being referenced.
\param propCode  Property code.
 
 \return Success=true
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoGetComponentProperty(Handle sysHndle, Handle compHndle, Property propCode)
{return;} 
/** \defgroup LuaThermoalpha Thermoalpha
 * Hallo*/
/** void chiThermoCreateSystem() Generate an empty Thermal-hydraulic system.
\return Returns a unique handle for the created system.
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoCreateSystem()
{return;} 
/** void chiThermoConnectTwoComponents() Connect two hydrodynamic components using a single junction.
\param systemHandle  Handle to the system to which all the components belong.
\param leftComponent Component to be connected to the left of the junction.
\param sjunc Single junction to be used for the connection.
\param rigtComponent Component to be connected to the right of the junction.
\param mode 0=end-begin, 1=begin-end, 2=end-end,
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoConnectTwoComponents(Handle systemHandle, Component leftComponent, Single sjunc, Component rigtComponent, 0=end-begin, mode)
{return;} 
/** bool chiThermoInitialize() Initializes system.
\param systemHandle int Handle to the system which should be initialized.
\return Returns true if successfully initialized and false otherwise.
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoInitialize(int systemHandle)
{return;} 
/** int chiThermoCreateVolumeFromCoordinates() Creates a hydrodynamic volume using a start and end coordinate.
\params systemHandle int Handle to the system to which the volume belongs.
\params point1		 Table Table with fields {x,y,z} containing the start location.
\params point2		 Table Table with fields {x,y,z} containing the end location.
\return Returns a unique handle for the created volume. (int)
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoCreateVolumeFromCoordinates(int systemHandle, Table point1, Table point2)
{return;} 
/** int chiThermoCreateSJunction() Creates a hydrodynamic single junction.
\params systemHandle Handle to the system to which the volume belongs.
\return Returns a unique handle for the created single junction.
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoCreateSJunction(Handle systemHandle)
{return;} 
/** int chiThermoCreateBC() Creates a hydrodynamic boundary condition.
\params systemHandle int Handle to the system to which the volume belongs.
\return Returns a unique handle for the created boundary condition.
\ingroup LuaThermoalpha
\author Jan*/
int CHI_LUA::chiThermoCreateBC(int systemHandle)
{return;} 
/** Splits an edge loop into edges if they differ by a certain angle.
\param LoopCollectionHandle int Handle to the Loop collection.
\param LoopHandle int Handle to the loop inside the collection.
\param Angle double (Optional) Value of the angle by which to split. Default 1 deg.
\return Handle int. Handle to the newly created LoopCollection.
\ingroup LuaMesh
\author Jan*/
int CHI_LUA::chiEdgeLoopSplitByAngle(int LoopCollectionHandle, int LoopHandle, double Angle)
{return;} 
/** \defgroup LuaMesh Meshing
 * \ingroup LuaPhysics
 *
 * ## 2D Mesh setup
 * Below is an example of setting up a 2D mesh.
 *
\code
chiMeshHandlerCreate()
--
newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,"CHI_RESOURCES/TestObjects/TestSurface5_simplices.obj")
loops,loop_count = chiSurfaceMeshGetEdgeLoops(newSurfMesh)
--
--
line_mesh = {};
line_mesh_count = 0;
--
for k=1,loop_count do
  split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
  for m=1,split_count do
    line_mesh_count = line_mesh_count + 1;
    line_mesh[line_mesh_count] = chiLineMeshCreateFromLoop(split_loops,m-1);
  end
--
end
\endcode
*/
/** Loads mesh data from a wavefront object.
\return success bool Return true if file was successfully loaded and false
 otherwise.
\ingroup LuaSurfaceMesh
\author Jan*/
int CHI_LUA::chiSurfaceMeshImportFromOBJFile()
{return;} 
/** Gets a list of edge loops for the given surface mesh.
\return Handle int Handle to the edge loops.
\return Count int Number of edge loops found.
\ingroup LuaSurfaceMesh
\author Jan*/
int CHI_LUA::chiSurfaceMeshGetEdgeLoops()
{return;} 
/** \defgroup LuaSurfaceMesh Surface Meshes
 * \ingroup LuaMesh
*/
/** Creates a new surface mesh.
\return Handle int Handle to the created surface mesh.
\ingroup LuaSurfaceMesh
\author Jan*/
int CHI_LUA::chiSurfaceMeshCreate()
{return;} 
/** Adds a surface mesh boundary to the region
\param RegionHandle int Handle to the region for which boundary is to be added.
\param SurfaceHandle int Handle to the surface mesh.
\ingroup LuaRegion
\author Jan*/
int CHI_LUA::chiRegionAddSurfaceBoundary(int RegionHandle, int SurfaceHandle)
{return;} 
/** Adds a line mesh boundary to the region
\param RegionHandle int Handle to the region for which boundary is to be added.
\param LineMeshHandle int Handle to the line mesh.
\ingroup LuaRegion
\author Jan*/
int CHI_LUA::chiRegionAddLineBoundary(int RegionHandle, int LineMeshHandle)
{return;} 
/** \defgroup LuaRegion Regions
 * \ingroup LuaMesh
 */
/** Creates a new region mesh.
\return Handle int Handle to the created region.
\ingroup LuaRegion
\author Jan*/
int CHI_LUA::chiRegionCreate()
{return;} 
/** \defgroup LuaMeshHandler Mesh Handler
 * \ingroup LuaMesh
*/
/** Creates a mesh handler.
\return Handle int Handle to the created mesh handler.
\ingroup LuaMeshHandler
\author Jan*/
int CHI_LUA::chiMeshHandlerCreate()
{return;} 
/** \defgroup LuaLineMesh Line Meshes
 * \ingroup LuaMesh
*/
/** Creates a new line mesh from a loop.
 *
\param LoopCollectionHandle int Handle to the Loop collection.
\param LoopHandle int Handle to the loop inside the collection.
\return Handle int Handle to the created line mesh.
\ingroup LuaLineMesh
\author Jan*/
int CHI_LUA::chiLineMeshCreateFromLoop(int LoopCollectionHandle, int LoopHandle)
{return;} 
/** \defgroup LuaSurfaceMesher Surface Re-meshers
 * \ingroup LuaMesh
*/
/** Creates a new surface mesher remeshing.
 *
\param Type int Surface Remesher type.
Remesher types:\n
 SURFACEMESHER_PREDEFINED = No remeshing is performed.\n
 SURFACEMESHER_DELAUNAY   = Delaunay surface remesher.
\ingroup LuaSurfaceMesher
\author Jan*/
int CHI_LUA::chiSurfaceMesherCreate(int Type)
{return;} 
/** Executes the surface meshing pipeline.
\ingroup LuaSurfaceMesher
\author Jan*/
int CHI_LUA::chiSurfaceMesherExecute()
{return;} 
/** \defgroup LuaVolumeMesher Volume Meshers
 * \ingroup LuaMesh
*/
/** Creates a new volume mesher.
 *
\param Type int Volume Remesher type.
Remesher types:\n
 VOLUMEMESHER_PREDEFINED2D = No remeshing is performed.\n
\ingroup LuaVolumeMesher
\author Jan*/
int CHI_LUA::chiVolumeMesherCreate(int Type)
{return;} 
/** Executes the volume meshing pipeline.
\ingroup LuaVolumeMesher
\author Jan*/
int CHI_LUA::chiVolumeMesherExecute()
{return;} 
/** Adds a region to a solver.
\ingroup LuaSolver
\author Jan*/
int CHI_LUA::chiSolverExecute()
{return;} 
/**\defgroup LuaPhysics Physics */
/** \defgroup LuaSolver Solvers
 * \ingroup LuaPhysics*/
/** Adds a region to a solver.
 *
\param SolverHandle int Handle to the solver.
\param RegionHandle int Handle to the region.
\ingroup LuaSolver
\author Jan*/
int CHI_LUA::chiSolverAddRegion(int SolverHandle, int RegionHandle)
{return;} 
/** Creates a simple point source at [0 0 -].
 *
\param SolverHandle int Handle to an existing montecarlo solver.
\return Handle int Handle to the created source.
\ingroup LuaMonteCarlon
\author Jan*/
int CHI_LUA::chiMonteCarlonCreateSource(int SolverHandle)
{return;} 
/**\defgroup LuaMonteCarlon Monte Carlo N-particle
 * \ingroup LuaSolver*/
/** Creates a MonteCarlon solver.
\return Handle int Handle to the created solver.
\ingroup LuaMonteCarlon
\author Jan*/
int CHI_LUA::chiMonteCarlonCreateSolver()
{return;} 
/** \defgroup LuaDiffusion Diffusion
 * \ingroup LuaSolver
 *
 *
 * Solves the general diffusion equation of the form
 *
 * \f{equation}{
 * -\nabla D \nabla \phi = q.
 * \f}
 *
 * Given a discretization method this solver will assemble the matrix $A$ using
 * a system of linear equations.
 *
 *
 * */
/** Initialize the Diffusion solver.
 *
\param SolverHandle int Handle to an existing diffusion solver.
\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
int CHI_LUA::chiDiffusionExecute(int SolverHandle)
{return;} 
/** Creates a Diffusion solver.
\return Handle int Handle to the created solver.
\ingroup LuaDiffusion
\author Jan*/
int CHI_LUA::chiDiffusionCreateSolver()
{return;} 
/** Initialize the Diffusion solver.
 *
\param SolverHandle int Handle to an existing diffusion solver.
\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
int CHI_LUA::chiDiffusionInitialize(int SolverHandle)
{return;} 
/** Sets a property of a Diffusion solver.
\param SolverHandle int Handle to an existing diffusion solver.
\param PropertyIndex int Code for a specific property.
\param Values varying Number of inputs associated with the index.<br>
##_
##PropertyIndex\n
DISCRETIZATION_METHOD =  Discretization method.\n
##Discretization methods
 CFEM2D = Linear Continues Finite Element 2D.
\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
int CHI_LUA::chiDiffusionSetProperty(int SolverHandle, int PropertyIndex, varying Values)
{return;} 
