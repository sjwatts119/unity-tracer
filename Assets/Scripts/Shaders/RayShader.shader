Shader "RayTracer/RayShader"
{
    SubShader
    {
        Pass
        {
            CGPROGRAM
            #pragma vertex RayVertexShader
            #pragma fragment RayFragmentShader
            #include "UnityCG.cginc"

            /*
             * Constants
             */
            
            #define INFINITY (1.0 / 0.0) // Positive infinity
            #define EPSILON 0.001       // A small value to offset rays to avoid self-intersection

            /*
             * Material Types
             */
            
            #define MATERIAL_LAMBERTIAN 0
            #define MATERIAL_METAL 1
            #define MATERIAL_DIELECTRIC 2
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
                float emissionStrength;
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
                bool didScatter;
                Ray scatteredRay;
                float3 attenuation;
                float3 emission;
            };

            struct Sphere
            {
                float3 centre;
                float radius;
                Material material;
            };

            struct Quad
            {
                float3 q;
                float3 u;
                float3 v;
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
             * Passed Parameters
             */

            // Sphere
            StructuredBuffer<Sphere> SphereBuffer;
            int SphereCount;

            // Quad
            StructuredBuffer<Quad> QuadBuffer;
            int QuadCount;

            // Camera
            float CameraFocalDistance;
            float CameraDefocusAngle;
            float CameraPlaneWidth;
            float CameraPlaneHeight;
            float4x4 CameraLocalToWorld;
            
            // Anti-aliasing
            int SamplesPerPixel;

            // Ray Tracing
            int RayMaxDepth;

            // Other
            int FrameNumber;

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

            // Generate a random unit vector inside a unit sphere with the Marsaglia method
            float2 PCGMarsagliaUnitCircle(inout uint state) {
                float u1, u2, s;
                
                do {
                    u1 = PCGRandomFloatInRange(state, -1.0, 1.0);
                    u2 = PCGRandomFloatInRange(state, -1.0, 1.0);
                    s = u1 * u1 + u2 * u2;
                } while (s >= 1.0 || s == 0.0);
                
                return float2(u1, u2);
            }

            float3 PCGRandomUnitVector(inout uint state) {
                float2 p = PCGMarsagliaUnitCircle(state);
                float s = p.x * p.x + p.y * p.y;
                
                float factor = 2.0 * sqrt(1.0 - s);
                return float3(p.x * factor, p.y * factor, 1.0 - 2.0 * s);
            }

            // Generate a random unit vector using the Marsaglia method
            float3 PCGRandomUnitVectorInUnitDisk(inout uint state) {
                float2 p = PCGMarsagliaUnitCircle(state);
                return float3(p.x, p.y, 0.0);
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

            // Thanks https://github.com/SebLague for this
			// Random value in normal distribution (with mean=0 and sd=1)
			float PCGRandomValueNormalDistribution(inout uint state)
			{
				float theta = 2 * 3.1415926 * PCGRandomFloat(state);
				float rho = sqrt(-2 * log(PCGRandomFloat(state)));
				return rho * cos(theta);
			}

            // Thanks https://github.com/SebLague for this
			// Calculate a random direction
			float3 PCGRandomDirection(inout uint state)
			{
				float x = PCGRandomValueNormalDistribution(state);
				float y = PCGRandomValueNormalDistribution(state);
				float z = PCGRandomValueNormalDistribution(state);
				return normalize(float3(x, y, z));
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

            // Set the normal and front face of a ray hit based on the ray direction and outward normal
            RayHit RayHitSetFaceNormal(RayHit hitRecord, Ray ray, float3 outwardNormal)
            {
                hitRecord.frontFace = dot(ray.direction, outwardNormal) < 0;
                hitRecord.normal = hitRecord.frontFace ? outwardNormal : -outwardNormal;
                return hitRecord;
            }

            // Use the schlick approximation to estimate the reflectance at a dielectric surface
            float SchlickReflectance(float cosine, float refractiveIndex)
            {
                float r0 = (1.0 - refractiveIndex) / (1.0 + refractiveIndex);
                r0 = r0 * r0;
                return r0 + (1.0 - r0) * pow((1.0 - cosine), 5.0);
            }

            // Check if a point defined by (a, b) is inside the quad defined by the origin, u, and v vectors
            bool IsInteriorToQuad(Quad quad, float a, float b)
            {
                Interval unitInterval = CreateInterval(0.0, 1.0);

                return IntervalContains(unitInterval, a) && IntervalContains(unitInterval, b);
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

                // Scatter the ray in a random direction against the normal (using cosine-weighted hemisphere sampling)
                float3 scatterDirection = hit.normal + PCGRandomDirection(seed);

                // Handle the degenerate case where the scatter direction is near zero
                if (any(abs(scatterDirection) < EPSILON))
                {
                    scatterDirection = hit.normal; // Fallback to normal direction
                }

                // Build up the scatter info
                scatter.scatteredRay = CreateRay(hit.position, normalize(scatterDirection));
                scatter.attenuation = material.albedo;
                scatter.emission = float3 (0, 0, 0); // Lambertian surfaces
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
                scatter.emission = float3 (0, 0, 0);

                // Our ray was scattered if the scattered direction is in the same hemisphere as the normal (it bounced off instead of being absorbed)
                scatter.didScatter = dot(scatter.scatteredRay.direction, hit.normal) > 0;

                return scatter;
            }

            // Scatter the ray for a Dielectric material
            RayScatter ScatterDielectric(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;
                
                float3 attenuation = float3(1, 1, 1); // For now, dielectrics don't attenuate the ray (no colour change)
                float refractionRatio = hit.frontFace ? (1.0 / material.refractiveIndex) : material.refractiveIndex;

                // Calculate reflection and refraction directions
                float3 unitDirection = normalize(ray.direction);
                float cosTheta = min(dot(-unitDirection, hit.normal), 1.0);
                float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

                float3 direction;

                // Check for total internal reflection
                bool cannotRefract = refractionRatio * sinTheta > 1.0;
                
                // Use Schlick's approximation to determine reflection probability
                float reflectance = SchlickReflectance(cosTheta, refractionRatio);
                bool shouldReflect = cannotRefract || (reflectance > PCGRandomFloat(seed));

                if (shouldReflect){
                    direction = reflect(unitDirection, hit.normal);
                } else {
                    direction = refract(unitDirection, hit.normal, refractionRatio);
                }

                // Populate the scatter info
                scatter.scatteredRay = CreateRay(hit.position, normalize(direction));
                scatter.attenuation = attenuation;
                scatter.emission = float3 (0, 0, 0);
                scatter.didScatter = true; // Dielectric always scatters
                return scatter;
            }

            RayScatter ScatterLight(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;
                scatter.scatteredRay = ray; // Purely because this property needs to be set, it gets ignored when checking didScatter
                scatter.didScatter = false; // We can consider this ray absorbed, its tracing terminates at the light source
                scatter.attenuation = float3(0, 0, 0);
                scatter.emission = material.emission * material.emissionStrength;
                return scatter;
            }

            // Get a scatter object based on the material type
            RayScatter ScatterByMaterial(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                switch (material.type)
                {
                    case MATERIAL_LAMBERTIAN:
                        return ScatterLambertian(ray, hit, material, seed);
                    case MATERIAL_METAL:
                        return ScatterMetal(ray, hit, material, seed);
                    case MATERIAL_DIELECTRIC:
                        return ScatterDielectric(ray, hit, material, seed);
                    case MATERIAL_LIGHT:
                        return ScatterLight(ray, hit, material, seed);
                    default:
                        RayScatter noScatter;
                        noScatter.scatteredRay = ray; // Purely because this property needs to be set, it gets ignored when checking didScatter
                        noScatter.didScatter = false; // Unknown material, so we consider the ray absorbed
                        noScatter.attenuation = float3(0, 0, 0);
                        noScatter.emission = float3(0, 0, 0);
                        return noScatter;
                }
            }

            /*
             * Ray-Geometry intersection
             */

            RayHit RayHitsQuad(Ray ray, Interval rayInterval, Quad quad)
            {
                // Calculate the plane normal and D constant for the quad plane
                float3 normal = normalize(cross(quad.u, quad.v));
                float D = dot(normal, quad.q);
                float denom = dot(normal, ray.direction);

                // Initialsze hit info with hit set to false for now
                RayHit hit = (RayHit)0;

                if (abs(denom) < EPSILON)
                {
                    return hit; // Ray is parallel to the quad plane, so no intersection
                }

                // Calculate the intersection distance t along the ray
                float t = (D - dot(normal, ray.origin)) / denom;
                
                // Check if the intersection is within the ray interval
                if (!IntervalContains(rayInterval, t))
                {
                    return hit; // Intersection is outside the ray interval
                }

                // Calculate the hit point
                float3 hitPoint = ray.origin + t * ray.direction;
                float3 planarHitpoint = hitPoint - quad.q;

                // Calculate barycentric coordinates (alpha, beta) for the hit point
                float barycentricDenom = dot(cross(quad.u, quad.v), cross(quad.u, quad.v));
                float alpha = dot(cross(planarHitpoint, quad.v), cross(quad.u, quad.v)) / barycentricDenom;
                float beta = dot(cross(quad.u, planarHitpoint), cross(quad.u, quad.v)) / barycentricDenom;

                // Check if the hit point is inside the quad bounds
                if (!IsInteriorToQuad(quad, alpha, beta))
                {
                    return hit; // Intersection is outside the quad bounds
                }

                // We have a valid intersection, so populate the hit info
                hit.didHit = true;
                hit.t = t;
                hit.position = hitPoint;
                hit.material = quad.material;
                hit = RayHitSetFaceNormal(hit, ray, normal);

                return hit;
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
                hit.position = ray.origin + hit.t * ray.direction;
                hit.material = sphere.material;
                float3 outwardNormal = (hit.position - sphere.centre) / sphere.radius;
                hit = RayHitSetFaceNormal(hit, ray, outwardNormal);
                return hit;
            }

            /*
             * Ray Tracing
             */

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

                // Iterate over all quads to find the closest hit
                for (int i = 0; i < QuadCount; i++)
                {
                    Quad quad = QuadBuffer[i];

                    // Check for intersection with the quad
                    RayHit hit = RayHitsQuad(ray, CreateInterval(EPSILON, closestHit.t), quad);

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
                float3 accumulatedLight = float3(0, 0, 0); // Track accumulated light
                
                for (int depth = 0; depth < RayMaxDepth; depth++)
                {
                    RayHit hit = GetHit(ray);

                   // Did this ray hit anything?
                    if (!hit.didHit)
                    {
                        // Sky gradient from white to light blue based on ray direction
                        // float t = 0.5 * (ray.direction.y + 1.0); // Convert Y from [-1,1] to [0,1]
                        // float3 skyColor = lerp(float3(1.0, 1.0, 1.0), float3(0.5, 0.7, 1.0), t);
                        // return accumulatedLight + (skyColor * rayColour);
                        return accumulatedLight * float3(0, 0, 0); // Black background
                    }

                    // We hit something, so scatter the ray based on the material
                    RayScatter scatter = ScatterByMaterial(ray, hit, hit.material, seed);

                    // Add emission from this surface, weighted by current ray color
                    accumulatedLight += scatter.emission * rayColour;

                    // If the material didn't scatter (like a light source), we're done
                    if (!scatter.didScatter) 
                    {
                        return accumulatedLight;
                    }

                    // Move to the scattered ray for the next bounce
                    ray = scatter.scatteredRay;

                    // Apply attenuation from our material
                    rayColour *= scatter.attenuation;

                    // Russian roulette termination to prevent spending resources on hardly contributing rays
                    if (depth < 3) continue; // Don't terminate the first few bounces

                    // Don't terminate if the ray is still bright (over 1% in any channel)
                    float maxComponent = max(rayColour.r, max(rayColour.g, rayColour.b));
                    if (maxComponent > 0.01) continue;

                    // Survival probability is directly proportional to brightness
                    if (PCGRandomFloat(seed) > maxComponent) break;

                    // Survived, so scale up the colour to maintain energy
                    rayColour /= maxComponent;
                }
                
                return accumulatedLight;
            }

            // Generate a ray from either the camera or a defocus disk and apply anti-aliasing jitter
            Ray GetRay(inout uint sampleSeed, float3 camRight, float3 camUp, VertexToFragment pixelData)
            {
                // Our ray starts at the camera position
                float3 rayOrigin = _WorldSpaceCameraPos;

                // If we have defocus enabled, jitter the ray origin within a disk perpendicular to the view direction
                if (CameraDefocusAngle > 0.0)
                {
                    float2 defocusJitter = PCGRandomUnitVectorInUnitDisk(sampleSeed) * CameraDefocusAngle * (CameraFocalDistance / _ScreenParams.xy);
                    rayOrigin += (camRight * defocusJitter.x) + (camUp * defocusJitter.y);
                }

                // Apply anti-aliasing jitter within the pixel
                float2 antialiasingJitter = PCGSampleSquare(sampleSeed) / _ScreenParams.xy;
                float2 jitteredCoord = pixelData.pixelCoordinates + antialiasingJitter;

                // Calculate the jittered focus point in local camera space, then transform to world space
                float3 jitteredFocusPointLocal = float3(jitteredCoord - 0.5, 1) * float3(CameraPlaneWidth, CameraPlaneHeight, CameraFocalDistance);
                float3 jitteredFocusPoint = mul(CameraLocalToWorld, float4(jitteredFocusPointLocal, 1));

                // Calculate ray direction from origin to jittered focus point
                float3 rayDirection = normalize(jitteredFocusPoint - rayOrigin);

                return CreateRay(rayOrigin, rayDirection);
            }

            /*
             * Shader entry points
             */

            // Runs per pixel
            fixed4 RayFragmentShader(VertexToFragment pixelData) : SV_Target
            {
                // Get pixel coordinates
                uint2 pixelCoord = uint2(pixelData.pixelCoordinates * _ScreenParams.xy);

                // Get camera basis vectors from the local to world matrix (Thanks to https://github.com/SebLague for this strategy)
                float3 camRight = CameraLocalToWorld._m00_m10_m20;
                float3 camUp = CameraLocalToWorld._m01_m11_m21;

                // Generate an RNG seed for this pixel
                uint pixelSeed = pixelCoord.x * 73856093u ^ pixelCoord.y * 19349663u + FrameNumber * 834123;

                // Start at black as if we have no intersections, light wouldn't be reflected to the camera
                float3 pixelColour = float3(0, 0, 0);
                
                // Sample multiple rays per pixel for anti-aliasing
                for (int sample = 0; sample < SamplesPerPixel; sample++)
                {
                    // Make a seed for this sample
                    uint sampleSeed = pixelSeed + sample * 12345u + FrameNumber * 834123;

                    // Create the ray
                    Ray ray = GetRay(sampleSeed, camRight, camUp, pixelData);
                    
                    // Accumulate colour
                    pixelColour += GetRayColour(ray, sampleSeed);
                }
                
                // Average the samples
                float pixelSampleScale = 1.0 / SamplesPerPixel;
                pixelColour *= pixelSampleScale;

                return fixed4(pixelColour, 1.0);
            }

            // Runs per vertex
            VertexToFragment RayVertexShader(appdata_base meshVertexData)
            {
                VertexToFragment vertexOutput;
                vertexOutput.screenPosition = UnityObjectToClipPos(meshVertexData.vertex);
                vertexOutput.pixelCoordinates = meshVertexData.texcoord;
                return vertexOutput;
            }
            ENDCG
        }
    }
}