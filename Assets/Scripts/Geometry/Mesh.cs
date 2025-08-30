using Geometry.Abstract;
using Geometry.Structs;
using UnityEngine;

namespace Geometry
{
    public class Mesh : TraceableGroup
    {
        public Triangle[] Triangles => GetTriangles();
        
        // Probably add other primitive types later (spheres, quads, cuboids)
        
        private Triangle[] GetTriangles()
        {
            // Get the mesh from the MeshFilter component
            var mesh = GetComponent<MeshFilter>().sharedMesh;
            
            var triangles = mesh.triangles;
            var vertices = mesh.vertices;
            var primitives = new Triangle[triangles.Length / 3];
            
            // Convert each triangle to world space and create a Triangle primitive
            for (var i = 0; i < triangles.Length; i += 3)
            {
                var v0 = transform.TransformPoint(vertices[triangles[i]]);
                var v1 = transform.TransformPoint(vertices[triangles[i + 1]]);
                var v2 = transform.TransformPoint(vertices[triangles[i + 2]]);
                primitives[i / 3] = new Triangle
                {
                    v0 = v0,
                    v1 = v1,
                    v2 = v2,
                    material = GetMaterial()
                };
            }
            
            return primitives;
        }

        public override AABB ToAABB()
        {
            // for now we are just looking at triangles, this is really scuffed, but it works for testing
            var triangles = Triangles;
            
            if (triangles.Length == 0)
            {
                return new AABB
                {
                    min = Vector3.zero,
                    max = Vector3.zero
                };
            }
            var min = triangles[0].v0;
            var max = triangles[0].v0;
            
            foreach (var triangle in triangles)
            {
                min = Vector3.Min(min, Vector3.Min(triangle.v0, Vector3.Min(triangle.v1, triangle.v2)));
                max = Vector3.Max(max, Vector3.Max(triangle.v0, Vector3.Max(triangle.v1, triangle.v2)));
            }
            
            return new AABB
            {
                min = min,
                max = max
            };
        }

        public override PrimitiveGroup ToPrimitiveGroup(int sphereCount, int quadCount, int cuboidCount, int triangleCount)
        {
            return new PrimitiveGroup
            {
                sphereStart = sphereCount,
                sphereCount = 0,
                quadStart = quadCount,
                quadCount = 0,
                cuboidStart = cuboidCount,
                cuboidCount = 0,
                triangleStart = triangleCount,
                triangleCount = Triangles.Length,
                boundingBox = ToAABB()
            };
        }
        
        void OnDrawGizmosSelected()
        {
            var aabb = ToAABB();
    
            Vector3 center = (aabb.min + aabb.max) * 0.5f;
            Vector3 size = aabb.max - aabb.min;
    
            Gizmos.color = Color.green;
            Gizmos.DrawWireCube(center, size);
        }
    }
    
    
}