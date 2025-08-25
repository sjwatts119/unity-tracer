Shader "RayTracer/AccumulationShader"
{
    Properties
    {
        _MainTex ("Current Texture", 2D) = "white" {} // Set automatically by Graphics.Blit
    }
    SubShader
    {
        Cull Off ZWrite Off ZTest Always
        Pass
        {
            CGPROGRAM
            #pragma vertex RayVertexShader
            #pragma fragment AccumulationFragmentShader
            #include "UnityCG.cginc"

            /*
             * Structs
             */
            
            struct VertexToFragment
            {
                float4 screenPosition : SV_POSITION;
                float2 pixelCoordinates : TEXCOORD0;
            };

            /*
             * Passed Parameters
             */
            
            sampler2D _MainTex;
            sampler2D PreviousTexture;
            int FrameNumber;

            /*
             * Shader entry points
             */

            // Runs per vertex
            VertexToFragment RayVertexShader(appdata_base meshVertexData)
            {
                VertexToFragment vertexOutput;
                vertexOutput.screenPosition = UnityObjectToClipPos(meshVertexData.vertex);
                vertexOutput.pixelCoordinates = meshVertexData.texcoord;
                return vertexOutput;
            }
            
            // Runs per pixel
            float4 AccumulationFragmentShader(VertexToFragment pixelData) : SV_Target
            {
                float4 current = tex2D(_MainTex, pixelData.pixelCoordinates);
                float4 previous = tex2D(PreviousTexture, pixelData.pixelCoordinates);
                
                // Calculate weight for blending
                float weight = 1.0 / (FrameNumber + 1);
                
                // Blend current frame with accumulated previous frames
                float4 accumulated = lerp(previous, current, weight);
                
                return accumulated;
            }
            ENDCG
        }
    }
}