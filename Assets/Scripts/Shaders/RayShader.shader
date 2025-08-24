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

            /*
             * Structs
             */
            
            struct Ray
            {
                float3 origin;
                float3 direction;
            };

            struct Sphere
            {
                float3 centre;
                float radius;
            };

            struct VertexToFragment
            {
                float4 screenPosition : SV_POSITION;
                float2 pixelCoordinates : TEXCOORD0;
            };

            /*
             * Struct "Constructors"
             */
            
            Ray CreateRay(float3 origin, float3 direction)
            {
                Ray ray;
                ray.origin = origin;
                ray.direction = direction;
                return ray;
            }

            Sphere CreateSphere(float3 centre, float radius)
            {
                Sphere sphere;
                sphere.centre = centre;
                sphere.radius = radius;
                return sphere;
            }

            /*
             * Buffers
             */

            StructuredBuffer<Sphere> SphereBuffer;

            /*
             * Methods
             */

            // Runs per vertex
            VertexToFragment RayTracerVertexShader(appdata_base meshVertexData)
            {
                VertexToFragment vertexOutput;
                vertexOutput.screenPosition = UnityObjectToClipPos(meshVertexData.vertex);
                vertexOutput.pixelCoordinates = meshVertexData.texcoord;
                return vertexOutput;
            }
            
            bool RayHitsSphere(Ray ray, Sphere sphere)
            {
                // Get the vector from the ray's origin to the sphere's centre
                float3 originToCentre = ray.origin - sphere.centre;

                // Calculate quadratic coefficients
                float a = dot(ray.direction, ray.direction);
                float b = 2.0 * dot(ray.direction, originToCentre);
                float c = dot(originToCentre, originToCentre) - sphere.radius * sphere.radius;

                // Now we can solve for three scenarios:
                // Negative: Ray misses the sphere entirely
                // Zero: Ray is tangent to the sphere (touches at exactly one point)
                // Positive: Ray intersects the sphere at two points (entry and exit)
                float discriminant = b * b - 4.0 * a * c;

                // Return whether or not the ray hits the sphere
                return (discriminant >= 0);
            }

            fixed4 GetRayColour(Ray ray)
            {
                // For now, we will just use the first sphere we have passed in via our buffer
                Sphere sphere = SphereBuffer[0];
                
                if (RayHitsSphere(ray, sphere))
                {
                    // Draw the sphere in red for now
                    return float4(1.0, 0.0, 0.0, 1.0);
                }
                
                // Get the unit vector of the ray direction
                float3 unitDirection = normalize(ray.direction);

                // Draw a background for now.
                float a = 0.5 * (unitDirection.y + 1.0);
                return (1.0 - a) * float4(1.0, 1.0, 1.0, 1.0) + a * float4(0.5, 0.7, 1.0, 1.0);
            }

            // Runs per pixel, requires return of RGBA colour with each channel in range 0-1
            fixed4 RayTracerFragmentShader(VertexToFragment pixelData) : SV_Target
            {
                // Convert pixel coordinates to device coordinates
                float4 deviceCoordinates = float4(pixelData.pixelCoordinates * 2.0 - 1.0, 1.0, 1.0);

                // Get ray direction in camera space by multiplying against the current position of the camera
                float3 cameraSpaceDirection = mul(unity_CameraInvProjection, deviceCoordinates).xyz;

                // Transform the ray's direction into world space
                float3 rayDirection = normalize(mul(unity_CameraToWorld, float4(cameraSpaceDirection, 0.0)).xyz);

                return GetRayColour(CreateRay(_WorldSpaceCameraPos,rayDirection));
            }
            ENDCG
        }
    }
}
