#version 130
//#########################################################
//#                                                       #
//#         SHADOWMAP GENERATION SHADER (FRAGMENT)        #
//#                                                       #
//#########################################################
in      vec3 shadowCoord;
in      vec4 fragPos;
in      vec3 lightPos;


void main (void)
{
      //float squareDistance;
      //squareDistance = clamp(length(shadowCoord)/888.0,0.0,1.0); 
      //
      //
      //if (squareDistance < 1.0) 
      //{
      //  gl_FragColor = vec4(squareDistance,squareDistance,squareDistance,1.0);
      //}
      
      float squareDistance;
      squareDistance=length(shadowCoord);
      //squareDistance=length(fragPos.xyz-lightPos);
      
      float factor1;
      float factor2;
      float factor3;
      float factor4;
      factor1=0.0;
      factor2=0.0;
      factor3=0.0;
      
      
      
      if (squareDistance<=10.0) 
      {
        factor1=clamp(squareDistance/10.0,0.0,1.0);
        
      }
      if ((squareDistance>10.0) && (squareDistance<=20.0))
      {
        factor1=1.0;
        factor2=clamp((squareDistance-10.0)/10.0,0.0,1.0);
      }
      if ((squareDistance>20.0) && (squareDistance<=30.0))
      {
        factor1=1.0;
        factor2=1.0;
        factor3=clamp((squareDistance-20.0)/10.0,0.0,1.0);
      }
      if ((squareDistance>30.0))
      {
        factor1=1.0;
        factor2=1.0;
        factor3=1.0;
      }

      if (squareDistance<0.0) {gl_FragColor = vec4(1.0,1.0,1.0,1.0);}
      else{gl_FragColor = vec4(factor1,factor2,factor3,1.0);}
      gl_FragColor = vec4(factor1,factor2,factor3,1.0);
      //gl_FragColor = vec4(factor1,factor1,factor1,1.0);  
}
