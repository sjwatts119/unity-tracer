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
            #define EPSILON 0.0001       // A small value to offset rays to avoid self-intersection

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

            /*
             * Helper Methods
             */

            float DegreesToRadians(float degrees)
            {
                return degrees * UNITY_PI / 180.0;
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

            fixed4 GetRayColour(Ray ray)
            {
                // Set up variables to track the closest hit, starting at infinity distance
                RayHit closestHit = (RayHit)0;
                float closestT = INFINITY;

                // Check for intersection with all spheres in the scene
                for (int i = 0; i < SphereCount; i++)
                {
                    // Check for intersection with this sphere
                    RayHit hitRecord = RayHitsSphere(ray, CreateInterval(EPSILON, closestT), SphereBuffer[i]);

                    // If we hit something and it's closer than our previous closest hit, update our closest hit
                    if (hitRecord.hit && hitRecord.t < closestT)
                    {
                        closestT = hitRecord.t;
                        closestHit = hitRecord;
                    }
                }

                // If we hit something, shade based on the normal at the hit point
                if (closestHit.hit)
                {
                    return 0.5f * (fixed4(closestHit.normal, 1.0) + fixed4(1.0, 1.0, 1.0, 1.0));
                }

                // If we didn't hit anything, return a gradient background colour
                float3 unitDirection = normalize(ray.direction);
                float a = 0.5 * (unitDirection.y + 1.0);
                return (1.0 - a) * float4(1.0, 1.0, 1.0, 1.0) + a * float4(0.5, 0.7, 1.0, 1.0);
            }

            // Runs per pixel, requires return of RGBA colour with each channel in range 0-1
            fixed4 RayTracerFragmentShader(VertexToFragment pixelData) : SV_Target
            {
                // Convert pixel coordinates (0-1) to image plane coordinates
                float2 imagePlanePosition = (pixelData.pixelCoordinates - 0.5) * float2(CameraPlaneWidth, CameraPlaneHeight);

                // Ray direction in camera space
                float3 cameraRayDirection = normalize(float3(imagePlanePosition.x, imagePlanePosition.y, CameraFocalDistance));

                // Transform ray directionto world space
                float3 rayDirection = normalize(mul((float3x3)unity_CameraToWorld, cameraRayDirection));

                // Ray origin is from the camera position in world space for now (subject to change for depth of field)
                float3 rayOrigin = _WorldSpaceCameraPos;

                return GetRayColour(CreateRay(rayOrigin, rayDirection));
            }
            ENDCG
        }
    }
}
