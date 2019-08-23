#version 330
//#########################################################
//#                                                       #
//#                  PHONG SHADER (FRAGMENT)              #
//#                                                       #
//#########################################################






//=============================================== Interpolations
in      vec3        texCoord0;
in      vec3        normal;
in      vec3        bi_norm;
in      vec3        tangent;

in      vec3        pixelWorldSpace;     //Absolute world space
in      vec3        pixelViewSpace;      //Pixel to eye

//=============================================== Uniforms
uniform mat4        projectionMatrix;   // 0Projection matrix
uniform mat4        viewMatrix;         // 1View matrix
uniform mat4        modelMatrix;        // 2Model matrix
uniform mat4        normalMatrix;       // 3Normal matrix
uniform mat4        textureMatrix;      // 4Texture matrix

uniform int         lightCount;         // 5Number of lights
uniform int         lightEnabled[10];   // 6 0=false, 1=true
uniform int         lightType[10];      // 7 0=Sun,1=point,2=spot
uniform vec3        lightPosition[10];  // 8Doubles as Sun-direction
uniform vec4        lightParams[10];    // 9Misc. light params
uniform ivec4       lightBParams[10];   //10Light bool-like params
uniform vec4        lightColor[10];     //11Light color

uniform vec4        materialDiffuse;    //12Material diffuse color
uniform vec4        materialAmbient;    //13Material ambient color
uniform vec4        materialSpecular;   //14Material specular color
uniform float       materialShininess;  //15Material shininess

uniform sampler2D   textureIndex[6];    //16Texture indices
uniform int         textureFlags[12];   //17Texture enabled flags

uniform ivec4       vobjectID;          //18Object ID

//=============================================== Color output definitions
layout(location = 0) out vec4 color0;
layout(location = 1) out vec4 color1;
layout(location = 2) out vec4 color2;

//=============================================== Utility function prototypes
vec4 ApplyPhong(vec4 inputDiffuseColor, vec4 inputSpecularColor,
                float inputShininess,int lightIndex,
                vec3 N, vec3 E);








//############################################### Main function
void main (void)
{
    vec4  ambientColor       = materialAmbient;
    vec4  diffuseColor       = materialDiffuse;
    vec4  specularColor      = materialSpecular;
    float shininess          = materialShininess;

    vec3 N = normalize(normal);          //Normal vector
    vec3 E = normalize(-pixelViewSpace); //Pixel to eye vector

    mat3 tangentMatrix = mat3(normalize(tangent), normalize(bi_norm), N);

    //================================= Setting up textures
    if (textureFlags[0] == 1) {ambientColor    = texture2D(textureIndex[0],texCoord0.st);}
    if (textureFlags[1] == 1) {diffuseColor    = texture2D(textureIndex[1],texCoord0.st);}
    if (textureFlags[2] == 1) {shininess       = texture2D(textureIndex[2],texCoord0.st).r;}
    if (textureFlags[3] == 1)
    {
        N    = normalize(texture2D(textureIndex[3],texCoord0.st).xyz*2.0 -1.0);
        N.y  = N.y*-1.0;
        N    = normalize(tangentMatrix*N);
    }

    //================================= if multiplications are enabled
    if (textureFlags[6] == 1) {ambientColor   *= materialAmbient;}
    if (textureFlags[7] == 1) {diffuseColor   *= materialDiffuse;}

    //================================= Lights
    vec4 shadedColor = vec4(0.0,0.0,0.0,0.0);
    for (int k=0;k<lightCount;k++)
    {
        if (lightEnabled[k]==1)
        {
            shadedColor += ApplyPhong(diffuseColor,specularColor,shininess,k,N,E);
        }
    }
    if (lightCount>0) {diffuseColor = shadedColor;}


    //================================= Output
    color0 = diffuseColor + ambientColor;
    //color0 = materialAmbient;
}
//###############################################








//=============================================== Apply Phong Lighting
vec4 ApplyPhong(vec4 inputDiffuseColor, vec4 inputSpecularColor, float inputShininess,int lightIndex, vec3 N, vec3 E)
{
    vec4  outputColor = vec4(0.0,0.0,0.0,0.0);;
    float diffuseIntensity  = 0.0;
    float specularIntensity = 0.0;
    float specularFactor    = 0.0;

    vec3 L = vec3(0.0,0.0,0.0);          //Light direction
    vec3 R = vec3(0.0,0.0,0.0);          //Reflection vector

    L = normalize(lightPosition[lightIndex]-pixelViewSpace);
    diffuseIntensity = clamp(max(0.0,dot(N,L)),0.0,1.0);
    outputColor+=
    diffuseIntensity*lightColor[lightIndex]*inputDiffuseColor;
    outputColor.a = inputDiffuseColor.a;

    R = reflect(-L,N);
    specularIntensity = clamp(max(0.0,dot(E,R)),0.0,1.0);
    specularFactor = pow(specularIntensity,inputShininess);

    outputColor+= specularFactor*lightColor[lightIndex]*inputSpecularColor;

    return outputColor;
}