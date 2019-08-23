#version 130
//#########################################################
//#                                                       #
//#         POST-PROCESS SHADER (FRAGMENT)                #
//#                                                       #
//#########################################################
uniform sampler2D texIndex0;

void main (void)
{
    
      
      gl_FragData[0] = texture2D(texIndex0,gl_TexCoord[0].st);
}