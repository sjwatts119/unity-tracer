// Assets/Scripts/RayTracing/Data/SceneData.cs

using System;
using Geometry.Structs;

namespace Core
{
    [System.Serializable]
    public class SceneData
    {
        // Individual primitives (not in groups)
        public Sphere[] individualSpheres = Array.Empty<Sphere>();
        public Quad[] individualQuads = Array.Empty<Quad>();
        public Cuboid[] individualCuboids = Array.Empty<Cuboid>();
        public Triangle[] individualTriangles = Array.Empty<Triangle>();

        // Grouped primitives (referenced by groups)
        public Sphere[] groupedSpheres = Array.Empty<Sphere>();
        public Quad[] groupedQuads = Array.Empty<Quad>();
        public Cuboid[] groupedCuboids = Array.Empty<Cuboid>();
        public Triangle[] groupedTriangles = Array.Empty<Triangle>();

        // Groups that reference the grouped primitives
        public PrimitiveGroup[] primitiveGroups = Array.Empty<PrimitiveGroup>();
        
        // Add BVH nodes TODO

        public bool HasData =>
            individualSpheres.Length > 0 || individualQuads.Length > 0 ||
            individualCuboids.Length > 0 || individualTriangles.Length > 0 ||
            primitiveGroups.Length > 0;

        public int TotalPrimitiveCount =>
            individualSpheres.Length + individualQuads.Length +
            individualCuboids.Length + individualTriangles.Length +
            groupedSpheres.Length + groupedQuads.Length +
            groupedCuboids.Length + groupedTriangles.Length;
    }
}
