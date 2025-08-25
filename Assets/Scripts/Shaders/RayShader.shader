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
             * Material Types (instead of enum)
             */
            #define MATERIAL_LAMBERTIAN 0
            #define MATERIAL_METAL 1
            #define MATERIAL_GLASS 2
            #define MATERIAL_LIGHT 3

            /*
             * Structs
             */

            struct Material
            {
                int type;
                float3 albedo;
                float fuzz;
                float refractiveIndex;
                float3 emission;
            };

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
                bool didHit;
                float3 position;
                float3 normal;
                float t;
                bool frontFace;
                Material material;
            };

            struct RayScatter
            {
                Ray scatteredRay;
                float3 attenuation;
                float3 emission;
                bool didScatter;
            };

            struct Sphere
            {
                float3 centre;
                float radius;
                Material material;
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

            Sphere CreateSphere(float3 centre, float radius, Material material)
            {
                Sphere sphere;
                sphere.centre = centre;
                sphere.radius = radius;
                sphere.material = material;
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

            /*
             * Material scattering
             */

            // Scatter the ray for a Lambertian (diffuse) material
            RayScatter ScatterLambertian(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;

                // Scatter the ray in a random direction against the normal
                float3 scatterDirection = hit.normal + PCGRandomUnitVector(seed);

                // Handle the degenerate case where the scatter direction is near zero
                if (any(abs(scatterDirection) < EPSILON))
                {
                    scatterDirection = hit.normal; // Fallback to normal direction
                }

                // Build up the scatter info
                scatter.scatteredRay = CreateRay(hit.position, normalize(scatterDirection));
                scatter.attenuation = material.albedo;
                scatter.emission = material.emission;
                scatter.didScatter = true; // Lambertian always scatters

                return scatter;
            }

            // Scatter the ray for a Metal material
            RayScatter ScatterMetal(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;

                // Reflect the ray direction around the normal
                float3 reflectDirection = reflect(ray.direction, hit.normal);

                if (material.fuzz > 0.0)
                {
                    reflectDirection = normalize(reflectDirection) + (PCGRandomUnitVector(seed) * material.fuzz);
                }

                scatter.scatteredRay = CreateRay(hit.position, normalize(reflectDirection));
                scatter.attenuation = material.albedo;
                scatter.emission = material.emission;

                // Our ray was scattered if the scattered direction is in the same hemisphere as the normal (it bounced off instead of being absorbed)
                scatter.didScatter = dot(scatter.scatteredRay.direction, hit.normal) > 0;

                return scatter;
            }

            // Get a scatter object based on the material type
            RayScatter ScatterByMaterial(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                if (material.type == MATERIAL_METAL)
                {
                    return ScatterMetal(ray, hit, material, seed);
                } else if (material.type == MATERIAL_LAMBERTIAN)
                {
                    return ScatterLambertian(ray, hit, material, seed);
                }

                return ScatterLambertian(ray, hit, material, seed); // Default to Lambertian if for some reason we get an unknown material
            }

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

                // Initialise hit info with hit set to false for now
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
                hit.didHit = true;
                hit.t = root;
                hit.position = RayAt(ray, root).origin;
                hit.material = sphere.material;
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
                    if (hit.didHit && hit.t < closestHit.t)
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
                
                for (int depth = 0; depth < RayMaxDepth; depth++)
                {
                    RayHit hit = GetHit(ray);

                    // Did this ray hit anything?
                    if (!hit.didHit)
                    {
                        float3 backgroundGradient = lerp(float3(1.0, 1.0, 1.0), float3(0.8, 0.1, 0.5), 0.5 * (normalize(ray.direction).y + 1.0));
                        
                        // If we hit nothing, so simply return our calculated colour so far attenuated by our background gradient
                        rayColour *= backgroundGradient;
                        break;
                    }

                    // We hit something, so scatter the ray based on the material
                    RayScatter scatter = ScatterByMaterial(ray, hit, hit.material, seed);

                    // If the material didn't scatter, it was absorbed so simply return black
                    if (!scatter.didScatter) 
                    {
                        rayColour = float3(0, 0, 0);
                        break;
                    }

                    // Now we can alter the ray based on the effects of the material
                    ray = scatter.scatteredRay;
                    rayColour *= scatter.attenuation;
                    
                    // TODO this is wrong, need to add a scale factor for emission and it should be additive, not multiplicative
                    // rayColour += scatter.emission * rayColour;

                    // Russian roulette termination to prevent spending resources on hardly contributing rays
                    if (depth < 3) continue; // Don't terminate the first few bounces

                    // Don't terminate if the ray is still bright (over 10% in any channel)
                    float maxComponent = max(rayColour.r, max(rayColour.g, rayColour.b));
                    if (maxComponent > 0.1) continue;

                    // Survival probability is directly proportional to brightness
                    // Dimmer rays are more likely to be terminated
                    if (PCGRandomFloat(seed) > maxComponent) break;

                    // Survived, so scale up the colour to maintain energy
                    rayColour /= maxComponent;
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
                float3 pixelColour = float3(0, 0, 0);
                
                // Sample multiple rays per pixel for anti-aliasing
                for (int sample = 0; sample < SamplesPerPixel; sample++)
                {
                    // Make a seed for this sample
                    uint sampleSeed = pixelSeed + sample * 12345u;
                    
                    // Generate random offset for this sample
                    float2 offset = PCGSampleSquare(sampleSeed);
                    
                    // Get ray with random offset
                    Ray ray = GetRay(pixelData.pixelCoordinates, offset);
                    
                    // Accumulate colour
                    pixelColour += GetRayColour(ray, sampleSeed);
                }
                
                // Average the samples
                float pixelSampleScale = 1.0 / SamplesPerPixel;
                pixelColour *= pixelSampleScale;

                return fixed4(pixelColour, 1.0);
            }
            ENDCG
        }
    }
}