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
             * Constants
             */
            #define INFINITY (1.0 / 0.0) // Positive infinity
            #define EPSILON 0.001       // A small value to offset rays to avoid self-intersection

            /*
             * Structs
             */

            struct Interval
            {
                float min;
                float max;
            };
            
            struct Ray
            {
                float3 origin;
                float3 direction;
            };

            struct RayHit
            {
                bool hit;
                float3 position;
                float3 normal;
                float t;
                bool frontFace;
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

            Interval CreateInterval(float min, float max)
            {
                Interval interval;
                interval.min = min;
                interval.max = max;
                return interval;
            }

            /*
             * Passed in buffer data
             */

            // Sphere
            StructuredBuffer<Sphere> SphereBuffer;
            int SphereCount;

            // Camera
            float CameraFocalDistance;
            float CameraPlaneWidth;
            float CameraPlaneHeight;
            
            // Anti-aliasing
            int SamplesPerPixel;

            // Ray Tracing
            int RayMaxDepth;

            /*
             * Helper Methods
             */

            float DegreesToRadians(float degrees)
            {
                return degrees * UNITY_PI / 180.0;
            }

            // Generate next random uint32
            uint PCGNext(inout uint state)
            {
                state = state * 747796405 + 2891336453;
                uint result = ((state >> ((state >> 28) + 4)) ^ state) * 277803737;
                result = (result >> 22) ^ result;
                return result;
            }

            float PCGRandomFloat(inout uint state)
            {
                return PCGNext(state) / 4294967295.0;
            }

            float PCGRandomFloatInRange(inout uint state, float min, float max)
            {
                return min + (max - min) * PCGRandomFloat(state);
            }

            // Random offset in range [-0.5, 0.5] using PCG algorithm
            float2 PCGSampleSquare(inout uint seed)
            {
                // Generate two random floats in range [-0.5, 0.5]
                float offsetX = PCGRandomFloatInRange(seed, -0.5, 0.5);
                float offsetY = PCGRandomFloatInRange(seed, -0.5, 0.5);
                
                return float2(offsetX, offsetY);
            }

            float3 PCGRandomVector(inout uint state)
            {
                return float3(PCGRandomFloat(state), PCGRandomFloat(state), PCGRandomFloat(state));
            }

            float3 PCGRandomVectorInRange(inout uint state, float min, float max)
            {
                return float3(PCGRandomFloatInRange(state, min, max), PCGRandomFloatInRange(state, min, max), PCGRandomFloatInRange(state, min, max));
            }

            // Generate a random unit vector using the Marsaglia method
            float3 PCGRandomUnitVector(inout uint state)
            {
                float u1, u2, s;

                do {
                    u1 = PCGRandomFloatInRange(state, -1.0, 1.0);
                    u2 = PCGRandomFloatInRange(state, -1.0, 1.0);
                    s = u1 * u1 + u2 * u2;
                } while (s >= 1.0 || s == 0.0);

                float factor = 2.0 * sqrt(1.0 - s);
                return float3(u1 * factor, u2 * factor, 1.0 - 2.0 * s);
            }

            float3 PCGRandomUnitVectorOnHemisphere(float3 normal, inout uint state)
            {
                float3 inUnitSphere = PCGRandomUnitVector(state);

                // Check if the random vector is in the same hemisphere as the normal, if not, invert it
                if (dot(inUnitSphere, normal) > 0.0) 
                    return inUnitSphere;
                else
                    return -inUnitSphere;
            }

            // Get a point along the ray at distance t from the origin
            Ray RayAt(Ray ray, float t)
            {
                return CreateRay(ray.origin + t * ray.direction, ray.direction);
            }

            float IntervalSize(Interval interval)
            {
                return interval.max - interval.min;
            }

            bool IntervalContains(Interval interval, float value)
            {
                return value >= interval.min && value <= interval.max;
            }

            bool IntervalSurrounds(Interval interval, float value)
            {
                return interval.min < value && value < interval.max;
            }

            RayHit RayHitSetFaceNormal(RayHit hitRecord, Ray ray, float3 outwardNormal)
            {
                hitRecord.frontFace = dot(ray.direction, outwardNormal) < 0;
                hitRecord.normal = hitRecord.frontFace ? outwardNormal : -outwardNormal;
                return hitRecord;
            }

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
            
            RayHit RayHitsSphere(Ray ray, Interval rayInterval, Sphere sphere)
            {
                // Get the vector from the ray's origin to the sphere's centre
                float3 originToCentre = ray.origin - sphere.centre;

                // Calculate quadratic coefficients
                float a = dot(ray.direction, ray.direction);
                float h = dot(ray.direction, originToCentre);
                float c = dot(originToCentre, originToCentre) - sphere.radius * sphere.radius;

                // Now we can solve for three scenarios:
                // Negative: Ray misses the sphere entirely
                // Zero: Ray is tangent to the sphere (touches at exactly one point)
                // Positive: Ray intersects the sphere at two points (entry and exit)
                float discriminant = h * h - a * c;

                // Initialize hit info with hit set to false for now
                RayHit hit = (RayHit)0;

                // If the discriminant is negative, the ray misses the sphere, so return our miss
                if (discriminant < 0)
                {
                    return hit; // No intersection
                }

                // We have an intersection, so return the nearest intersection distance
                float sqrtDiscriminant = sqrt(discriminant);
                
                // Find the nearest root that lies in the acceptable range.
                float root = (-h - sqrtDiscriminant) / a;
                if (!IntervalContains(rayInterval, root))
                {
                    root = (-h + sqrtDiscriminant) / a;
                    if (!IntervalContains(rayInterval, root))
                    {
                        return hit;
                    }
                }

                // We have a valid intersection within the ray interval, so populate the hit info
                hit.hit = true;
                hit.t = root;
                hit.position = RayAt(ray, root).origin;
                float3 outwardNormal = (hit.position - sphere.centre) / sphere.radius;
                hit = RayHitSetFaceNormal(hit, ray, outwardNormal);
                return hit;
            }

            // Get a ray for anti-aliasing sampling
            Ray GetRay(float2 uv, float2 offset)
            {
                // Convert UV coordinates (0-1) to image plane coordinates with offset
                float2 offsetUV = uv + offset / float2(_ScreenParams.x, _ScreenParams.y);
                float2 imagePlanePosition = (offsetUV - 0.5) * float2(CameraPlaneWidth, CameraPlaneHeight);

                // Ray direction in camera space
                float3 cameraRayDirection = normalize(float3(imagePlanePosition.x, imagePlanePosition.y, CameraFocalDistance));

                // Transform ray direction to world space
                float3 rayDirection = normalize(mul((float3x3)unity_CameraToWorld, cameraRayDirection));

                // Ray origin is from the camera position in world space for now (subject to change for depth of field)
                float3 rayOrigin = _WorldSpaceCameraPos;

                return CreateRay(rayOrigin, rayDirection);
            }

            // Get the closest hit for a ray
            RayHit GetHit(Ray ray)
            {
                // Start with no hit
                RayHit closestHit = (RayHit)0;
                closestHit.t = INFINITY;

                // Iterate over all spheres to find the closest hit
                for (int i = 0; i < SphereCount; i++)
                {
                    Sphere sphere = SphereBuffer[i];

                    // Check for intersection with the sphere
                    RayHit hit = RayHitsSphere(ray, CreateInterval(EPSILON, closestHit.t), sphere);

                    // If we have a hit and it's closer than our current closest hit, update closest hit
                    if (hit.hit && hit.t < closestHit.t)
                    {
                        closestHit = hit;
                    }
                }

                // In the future, we can add more geometry type collision checks here

                return closestHit;
            }

            // Get the colour for a ray by tracing it through the scene
            fixed3 GetRayColour(Ray ray, inout uint seed)
            {
                float3 rayColour = float3(1, 1, 1);
                float3 backgroundColour = float3(0.5, 0.7, 1.0); // Light blue background
                
                for (int depth = 0; depth < RayMaxDepth; depth++)
                {
                    RayHit hit = GetHit(ray);
                    
                    if (!hit.hit)
                    {
                        // If we hit nothing, return our current colour multiplied by the background gradient
                        float t = 0.5 * (normalize(ray.direction).y + 1.0);
                        rayColour *= backgroundColour;
                        break;
                    }
                    
                    ray.origin = hit.position;
                    ray.direction = PCGRandomUnitVectorOnHemisphere(hit.normal, seed);
                    rayColour *= 0.5;
                }
                
                return rayColour;
            }

            // Runs per pixel
            fixed4 RayTracerFragmentShader(VertexToFragment pixelData) : SV_Target
            {
                // Get pixel coordinates
                uint2 pixelCoord = uint2(pixelData.pixelCoordinates * _ScreenParams.xy);

                // Generate an RNG seed for this pixel
                uint pixelSeed = pixelCoord.x * 73856093u ^ pixelCoord.y * 19349663u;

                // Start at black as if we have no intersections, light wouldn't be reflected to the camera
                float3 pixelColor = float3(0, 0, 0);
                
                // Sample multiple rays per pixel for anti-aliasing
                for (int sample = 0; sample < SamplesPerPixel; sample++)
                {
                    // Make a seed for this sample
                    uint sampleSeed = pixelSeed + sample * 12345u;
                    
                    // Generate PCG-based random offset for this sample
                    float2 offset = PCGSampleSquare(sampleSeed);
                    
                    // Get ray with random offset
                    Ray ray = GetRay(pixelData.pixelCoordinates, offset);
                    
                    // Accumulate color using the same seed state
                    pixelColor += GetRayColour(ray, sampleSeed);
                }
                
                // Average the samples
                float pixelSampleScale = 1.0 / SamplesPerPixel;
                pixelColor *= pixelSampleScale;

                return fixed4(pixelColor, 1.0);
            }
            ENDCG
        }
    }
}