using System;
using Geometry.Structs;

namespace Core
{
    [System.Serializable]
    public class SceneData
    {
        public Sphere[] spheres = Array.Empty<Sphere>();
        public Quad[] quads = Array.Empty<Quad>();
        public Cuboid[] cuboids = Array.Empty<Cuboid>();
        public Triangle[] triangles = Array.Empty<Triangle>();

        // Add BVH nodes TODO

        public bool HasData => spheres.Length > 0 || quads.Length > 0 || cuboids.Length > 0 || triangles.Length > 0;

        public int TotalPrimitiveCount =>  spheres.Length + quads.Length + cuboids.Length + triangles.Length;
    }
}
