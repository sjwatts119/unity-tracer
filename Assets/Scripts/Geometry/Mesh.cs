using System;
using Geometry.Abstract;
using Geometry.Structs;
using UnityEngine;

namespace Geometry
{
    public class Mesh : Traceable
    {
        public Triangle[] GetPrimitives()
        {
            var mesh = GetComponent<MeshFilter>().sharedMesh;
            var triangles = mesh.triangles;
            var vertices = mesh.vertices;
            var primitives = new Triangle[triangles.Length / 3];
            
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

        public AABB ToAABB()
        {
            var triangles = GetPrimitives();
            
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