Shader "RayTracer/RayShader"
{
    SubShader
    {
        Pass
        {
            CGPROGRAM
            #pragma vertex RayTracerVertexShader
            #pragma fragment RayTracerFragmentShader
            #include "UnityCG.cginc"

            struct VertexToFragment
            {
                float4 screenPosition : SV_POSITION;
                float2 pixelCoordinates : TEXCOORD0;
            };

            // Runs per vertex
            VertexToFragment RayTracerVertexShader(appdata_base meshVertexData)
            {
                VertexToFragment vertexOutput;
                vertexOutput.screenPosition = UnityObjectToClipPos(meshVertexData.vertex);
                vertexOutput.pixelCoordinates = meshVertexData.texcoord;
                return vertexOutput;
            }

            // Runs per pixel
            fixed4 RayTracerFragmentShader(VertexToFragment pixelData) : SV_Target
            {
                return fixed4(pixelData.pixelCoordinates.x, pixelData.pixelCoordinates.y, 0.0, 1.0);
            }
            ENDCG
        }
    }
}