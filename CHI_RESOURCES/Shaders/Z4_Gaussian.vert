#version 130
//#########################################################
//#                                                       #
//#           GAUSSIAN BLUR SHADER (VERTEX)               #
//#                                                       #
//#########################################################

in          vec4    vVertex;

void main() {
	gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_Position = gl_ModelViewProjectionMatrix*vVertex;
}
