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
            #define INTERSECTION_EPSILON 1e-5 // Epsilon for intersection tests
            #define OFFSET_EPSILON 1e-3 // Epsilon to avoid self-intersection when offsetting ray origins

            /*
             * Material Types
             */
            
            #define MATERIAL_SOLID 0
            #define MATERIAL_DIELECTRIC 1
            #define MATERIAL_EMISSIVE 2

            /*
             * Structs
             */

            struct VertexToFragment
            {
                float4 screenPosition : SV_POSITION;
                float2 pixelCoordinates : TEXCOORD0;
            };

            struct Interval
            {
                float min;
                float max;
            };
            
            struct Material
            {
                int type;
                float3 colour;
                
                float reflectivity;
                float roughness;
                
                float refractiveIndex;
                float absorptionStrength;
                
                float emissionStrength;
            };
            
            struct Ray
            {
                float3 origin;
                float3 direction;
                float3 inverseDirection;
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

            struct AABB
            {
                float3 min;
                float3 max;
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

            struct Cuboid
            {
                float3 centre;
                float3 size;
                float4x4 worldToLocal;
                float4x4 localToWorld;
                Material material;
            };

            struct Triangle
            {
                float3 v0;
                float3 v1;
                float3 v2;
                Material material;
            };

            struct BvhNode
            {
                float3 aabbMin;
                float3 aabbMax;
                int leftFirst; // If leaf node, index of first primitive. If internal node, index of first child node.
                int primitiveCount; // If leaf node, number of primitives. If internal node, 0.
            };

            /*
             * Struct "Constructors"
             */
            
            Ray CreateRay(float3 origin, float3 direction)
            {
                Ray ray;
                ray.origin = origin;
                ray.direction = direction;
                ray.inverseDirection = rcp(direction);
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

            AABB CreateAABB(float3 min, float3 max)
            {
                AABB aabb;
                aabb.min = min;
                aabb.max = max;
                return aabb;
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

            // Cuboid
            StructuredBuffer<Cuboid> CuboidBuffer;
            int CuboidCount;

            // Triangle
            StructuredBuffer<Triangle> TriangleBuffer;
            int TriangleCount;

            // Bvh Nodes
            StructuredBuffer<BvhNode> BvhNodeBuffer;
            int BvhNodeCount;

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
				float theta = 2 * UNITY_PI * PCGRandomFloat(state);
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

            // Scatter the ray for a Metal material
            RayScatter ScatterMetal(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;

                // Reflect the ray direction around the normal
                float3 reflectDirection = reflect(ray.direction, hit.normal);

                if (material.roughness > 0.0)
                {
                    reflectDirection = normalize(reflectDirection) + (PCGRandomUnitVector(seed) * material.roughness);
                }

                scatter.scatteredRay = CreateRay(hit.position + hit.normal * OFFSET_EPSILON, normalize(reflectDirection));
                scatter.attenuation = material.colour;
                scatter.emission = float3 (0, 0, 0);

                // Our ray was scattered if the scattered direction is in the same hemisphere as the normal (it bounced off instead of being absorbed)
                scatter.didScatter = dot(scatter.scatteredRay.direction, hit.normal) > 0;

                return scatter;
            }

            // Scatter the ray for a Lambertian (diffuse) material
            RayScatter ScatterSolid(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                // If the material has reflectivity, we need to decide whether to reflect or scatter diffusely
                if (material.reflectivity > PCGRandomFloat(seed))
                {
                    return ScatterMetal(ray, hit, material, seed);
                }

                RayScatter scatter;

                // Scatter the ray in a random direction against the normal (using cosine-weighted hemisphere sampling)
                float3 scatterDirection = hit.normal + PCGRandomDirection(seed);

                // Handle the degenerate case where the scatter direction is near zero
                if (any(abs(scatterDirection) < INTERSECTION_EPSILON))
                {
                    scatterDirection = hit.normal; // Fallback to normal direction
                }

                // Build up the scatter info
                // scatter.scatteredRay = CreateRay(hit.position, normalize(scatterDirection));
                scatter.scatteredRay = CreateRay(hit.position + hit.normal * OFFSET_EPSILON, normalize(scatterDirection));
                scatter.attenuation = material.colour;
                scatter.emission = float3 (0, 0, 0); // Lambertian surfaces
                scatter.didScatter = true; // Lambertian always scatters

                return scatter;
            }

            // Scatter the ray for a Dielectric material
            RayScatter ScatterDielectric(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;
                
                float3 attenuation = float3(1, 1, 1);

                // Handle absorption for rays inside the material
                if (!hit.frontFace && material.absorptionStrength > 0.0)
                {
                    // Invert the colour to get the absorption coefficient
                    float3 absorptionCoeff = (1.0 - material.colour) * material.absorptionStrength;
                    attenuation = exp(-absorptionCoeff * hit.t);
                }
                
                float refractionRatio = hit.frontFace ? (rcp(material.refractiveIndex)) : material.refractiveIndex;

                // Calculate reflection and refraction directions
                float3 unitDirection = normalize(ray.direction);
                float cosTheta = min(dot(-unitDirection, hit.normal), 1.0);
                float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
                
                // Check for total internal reflection
                bool cannotRefract = refractionRatio * sinTheta > 1.0;
                
                // Use Schlick's approximation to determine reflection probability
                float reflectance = SchlickReflectance(cosTheta, refractionRatio);
                bool shouldReflect = cannotRefract || (reflectance > PCGRandomFloat(seed));

                float3 direction;
                float3 offsetDirection;

                if (shouldReflect){
                    direction = reflect(unitDirection, hit.normal);
                    offsetDirection = hit.normal;
                } else {
                    direction = refract(unitDirection, hit.normal, refractionRatio);
                    offsetDirection = direction;
                }

                // Populate the scatter info
                scatter.scatteredRay = CreateRay(hit.position + offsetDirection * OFFSET_EPSILON, direction);
                scatter.attenuation = attenuation;
                scatter.emission = float3 (0, 0, 0);
                scatter.didScatter = true; // Dielectric always scatters
                return scatter;
            }

            RayScatter ScatterEmissive(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                RayScatter scatter;
                scatter.scatteredRay = ray; // Purely because this property needs to be set, it gets ignored when checking didScatter
                scatter.didScatter = false; // We can consider this ray absorbed, its tracing terminates at the light source
                scatter.attenuation = float3(0, 0, 0);
                scatter.emission = material.colour * material.emissionStrength;
                return scatter;
            }

            // Get a scatter object based on the material type
            RayScatter ScatterByMaterial(Ray ray, RayHit hit, Material material, inout uint seed)
            {
                switch (material.type)
                {
                    case MATERIAL_SOLID:
                        return ScatterSolid(ray, hit, material, seed);
                    case MATERIAL_DIELECTRIC:
                        return ScatterDielectric(ray, hit, material, seed);
                    case MATERIAL_EMISSIVE:
                        return ScatterEmissive(ray, hit, material, seed);
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

            // HLSL Optimised slab intersection method, thanks https://medium.com/@bromanz
            bool SlabIntersection(float3 origin, float3 invDir, float3 boxMin, float3 boxMax, Interval rayInterval, out float tNear, out float tFar)
            {
                float3 t0 = (boxMin - origin) * invDir;
                float3 t1 = (boxMax - origin) * invDir;
                
                float3 tMin = min(t0, t1);
                float3 tMax = max(t0, t1);
                
                tNear = max(rayInterval.min, max(tMin.x, max(tMin.y, tMin.z)));
                tFar = min(rayInterval.max, min(tMax.x, min(tMax.y, tMax.z)));
                
                return tNear < tFar;
            }

            // Check if a ray hits an AABB
            bool RayHitsAABB(Ray ray, Interval rayInterval, AABB aabb)
            {
                float tNear, tFar; // Unused out parameters
                return SlabIntersection(ray.origin, ray.inverseDirection, aabb.min, aabb.max, rayInterval, tNear, tFar);
            }

            // Check if a ray hits a cuboid
            RayHit RayHitsCuboid(Ray ray, Interval rayInterval, Cuboid cuboid)
            {
                RayHit hit = (RayHit)0;

                // Transform the ray into the cuboid's local space
                float3 rayToCentre = ray.origin - cuboid.centre;
                float3 localOrigin = mul((float3x3)cuboid.worldToLocal, rayToCentre);
                float3 localDirection = mul((float3x3)cuboid.worldToLocal, ray.direction);

                // Define the half-size of the cuboid in local space
                float3 halfSize = cuboid.size * 0.5;
                float3 inverseDirection = rcp(localDirection);
                
                float tNear, tFar;
                if (!SlabIntersection(localOrigin, inverseDirection, -halfSize, halfSize, rayInterval, tNear, tFar))
                {
                    return hit;
                }
                
                // Check if ray origin is inside the cube
                bool rayStartsInside = all(abs(localOrigin) <= halfSize);

                // If the ray starts inside the cuboid, position this interval to start at the far face
                float t = rayStartsInside ? tFar : (tNear >= rayInterval.min) ? tNear : tFar;

                // Check if t is within the ray interval
                if (t < rayInterval.min || t > rayInterval.max)
                {
                    return hit;
                }

                // Get hit points in local and world space
                float3 localHitPoint = localOrigin + t * localDirection;
                float3 worldHitPoint = ray.origin + t * ray.direction;
                
                // Calculate normal based on which face was hit
                float3 faceDistance = abs(localHitPoint / halfSize);
                
                float3 localNormal;

                // Determine which face was hit based on the largest component of the normalised hit point
                if (faceDistance.x > faceDistance.y - INTERSECTION_EPSILON && faceDistance.x > faceDistance.z - INTERSECTION_EPSILON)
                    localNormal = float3(sign(localHitPoint.x), 0, 0);
                else if (faceDistance.y > faceDistance.z - INTERSECTION_EPSILON)  
                    localNormal = float3(0, sign(localHitPoint.y), 0);
                else
                    localNormal = float3(0, 0, sign(localHitPoint.z));

                // Transform the local normal back to world space
                float3 worldNormal = normalize(mul(cuboid.localToWorld, float4(localNormal, 0.0)).xyz);

                // We have a valid intersection, so populate the hit info
                hit.didHit = true;
                hit.t = t;
                hit.position = worldHitPoint;
                hit.material = cuboid.material;
                hit = RayHitSetFaceNormal(hit, ray, worldNormal);
                
                return hit;
            }

            RayHit RayHitsQuad(Ray ray, Interval rayInterval, Quad quad)
            {
                // Calculate the plane normal and D constant for the quad plane
                float3 normal = normalize(cross(quad.u, quad.v));
                float D = dot(normal, quad.q);
                float denom = dot(normal, ray.direction);

                // Initialise hit info with hit set to false for now
                RayHit hit = (RayHit)0;

                if (abs(denom) < INTERSECTION_EPSILON)
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
                float3 crossUV = cross(quad.u, quad.v);
                float barycentricDenom = dot(crossUV, crossUV);
                float alpha = dot(cross(planarHitpoint, quad.v), crossUV) / barycentricDenom;
                float beta = dot(cross(quad.u, planarHitpoint), crossUV) / barycentricDenom;

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

            // Moller Trumbore intersection algorithm
            RayHit RayHitsTriangle(Ray ray, Interval rayInterval, Triangle tri)
            {
                // Initialise hit info with hit set to false for now
                RayHit hit = (RayHit)0;

                // Calculate the two edge vectors of the triangle
                float3 edge1 = tri.v1 - tri.v0;
                float3 edge2 = tri.v2 - tri.v0;
                
                // Begin calculating the determinant - also used to calculate u parameter
                float3 h = cross(ray.direction, edge2);
                float a = dot(edge1, h);

                // If the determinant is near zero, the ray is parallel to the triangle plane
                if (abs(a) < INTERSECTION_EPSILON)
                {
                    return hit; // Ray is parallel to the triangle, so no intersection
                }

                float f = 1.0 / a;
                float3 s = ray.origin - tri.v0;
                
                // Calculate the u parameter and test bounds
                float u = f * dot(s, h);

                // Check if the intersection is outside the triangle (u coordinate)
                if (u < 0.0 || u > 1.0)
                {
                    return hit; // Intersection is outside the triangle bounds
                }

                // Prepare to test the v parameter
                float3 q = cross(s, edge1);
                
                // Calculate the v parameter and test bounds
                float v = f * dot(ray.direction, q);

                // Check if the intersection is outside the triangle (v coordinate and u+v constraint)
                if (v < 0.0 || u + v > 1.0)
                {
                    return hit; // Intersection is outside the triangle bounds
                }

                // Calculate the intersection distance t along the ray
                float t = f * dot(edge2, q);

                // Check if the intersection is within the ray interval
                if (!IntervalContains(rayInterval, t))
                {
                    return hit; // Intersection is outside the ray interval
                }

                // We have a valid intersection, so populate the hit info
                hit.didHit = true;
                hit.t = t;
                hit.position = ray.origin + t * ray.direction;
                hit.material = tri.material;
                hit = RayHitSetFaceNormal(hit, ray, normalize(cross(edge1, edge2)));

                return hit;
            }
            
            /*
             * Ray Tracing
             */

            // Get the closest hit for a ray
            RayHit GetHit(Ray ray)
            {
                RayHit closestHit = (RayHit)0;
                closestHit.t = INFINITY;
                Interval rayInterval = CreateInterval(INTERSECTION_EPSILON, closestHit.t);

                // Test individual primitives first
                for (int i = 0; i < SphereCount; i++)
                {
                    Sphere sphere = SphereBuffer[i];
                    rayInterval.max = closestHit.t;
                    RayHit hit = RayHitsSphere(ray, rayInterval, sphere);
                    if (hit.didHit && hit.t < closestHit.t)
                    {
                        closestHit = hit;
                    }
                }

                for (int i = 0; i < CuboidCount; i++)
                {
                    Cuboid cuboid = CuboidBuffer[i];
                    rayInterval.max = closestHit.t;
                    RayHit hit = RayHitsCuboid(ray, rayInterval, cuboid);
                    if (hit.didHit && hit.t < closestHit.t)
                    {
                        closestHit = hit;
                    }
                }

                for (int i = 0; i < QuadCount; i++)
                {
                    Quad quad = QuadBuffer[i];
                    rayInterval.max = closestHit.t;
                    RayHit hit = RayHitsQuad(ray, rayInterval, quad);
                    if (hit.didHit && hit.t < closestHit.t)
                    {
                        closestHit = hit;
                    }
                }
                
                // Traverse BVH with distance-ordered traversal
                int stack[32];
                int stackPtr = 0;
                stack[stackPtr++] = 0;

                while (stackPtr > 0)
                {
                    int nodeIndex = stack[--stackPtr];
                    BvhNode node = BvhNodeBuffer[nodeIndex];
                    
                    if (node.primitiveCount > 0) // Leaf node
                    {
                        for (int i = 0; i < node.primitiveCount; i++)
                        {
                            Triangle tri = TriangleBuffer[node.leftFirst + i];
                            Interval triInterval = CreateInterval(INTERSECTION_EPSILON, closestHit.t);
                            RayHit hit = RayHitsTriangle(ray, triInterval, tri);
                            if (hit.didHit && hit.t < closestHit.t)
                            {
                                closestHit = hit;
                            }
                        }
                    }
                    else // Internal node
                    {
                        int leftIdx = node.leftFirst;
                        int rightIdx = node.leftFirst + 1;
                        
                        // Get distances to both children
                        float tNearL, tFarL, tNearR, tFarR;
                        BvhNode leftChild = BvhNodeBuffer[leftIdx];
                        BvhNode rightChild = BvhNodeBuffer[rightIdx];
                        
                        bool hitLeft = SlabIntersection(ray.origin, ray.inverseDirection, 
                            leftChild.aabbMin, leftChild.aabbMax, 
                            CreateInterval(INTERSECTION_EPSILON, closestHit.t), tNearL, tFarL);
                        bool hitRight = SlabIntersection(ray.origin, ray.inverseDirection, 
                            rightChild.aabbMin, rightChild.aabbMax, 
                            CreateInterval(INTERSECTION_EPSILON, closestHit.t), tNearR, tFarR);
                        
                        // Push children in distance order (furthest first, nearest last)
                        if (hitLeft && hitRight)
                        {
                            if (tNearL < tNearR)
                            {
                                if (tNearR < closestHit.t) stack[stackPtr++] = rightIdx;
                                if (tNearL < closestHit.t) stack[stackPtr++] = leftIdx;
                            }
                            else
                            {
                                if (tNearL < closestHit.t) stack[stackPtr++] = leftIdx;
                                if (tNearR < closestHit.t) stack[stackPtr++] = rightIdx;
                            }
                        }
                        else if (hitLeft && tNearL < closestHit.t)
                        {
                            stack[stackPtr++] = leftIdx;
                        }
                        else if (hitRight && tNearR < closestHit.t)
                        {
                            stack[stackPtr++] = rightIdx;
                        }
                    }
                }
                
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
                        return accumulatedLight; // Black background
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
                    if (depth < 5) continue;

                    // Calculate the maximum component of the current ray colour
                    float maxComponent = max(rayColour.r, max(rayColour.g, rayColour.b));
                    
                    // Clamp our minimum survival probability to avoid creating excessively bright or dark pixels
                    float survivalProbability = max(maxComponent, 0.1);  
                    
                    if (PCGRandomFloat(seed) > survivalProbability) {
                        break;
                    }
                    
                    // Survived, so scale up the colour to keep the same amount of energy in the scene on average
                    rayColour /= survivalProbability;
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
                    float2 defocusJitter = PCGRandomUnitVectorInUnitDisk(sampleSeed).xy * CameraDefocusAngle * (CameraFocalDistance / _ScreenParams.xy);
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
                uint pixelSeed = pixelCoord.x * 7919u + pixelCoord.y * 7927u + FrameNumber * 1009u;

                // Start at black as if we have no intersections, light wouldn't be reflected to the camera
                float3 pixelColour = float3(0, 0, 0);
                
                // Sample multiple rays per pixel for anti-aliasing
                for (int sample = 0; sample < SamplesPerPixel; sample++)
                {
                    // Create the ray
                    Ray ray = GetRay(pixelSeed, camRight, camUp, pixelData);
                    
                    // Accumulate colour
                    pixelColour += GetRayColour(ray, pixelSeed);
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