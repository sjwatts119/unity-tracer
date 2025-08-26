using Geometry.Interfaces;
using UnityEngine;
using Vector3 = UnityEngine.Vector3;

namespace Geometry.Structs
{
    [System.Serializable]
    public struct Sphere : IPrimitive
    {
        public Vector3 centre;
        public float radius;
        public Materials.Structs.Material material;
    }

    [System.Serializable]
    public struct Quad : IPrimitive
    {
        public Vector3 q; // Bottom-left corner
        public Vector3 u; // Edge vector along the width ( q ---u--- )
        public Vector3 v; // Edge vector along the height (rotate above example -90 degrees)
        public Materials.Structs.Material material;
    }
    
    [System.Serializable]
    public struct Cuboid : IPrimitive
    {
        public Vector3 centre;
        public Vector3 size; // Scale in each dimension
        public Matrix4x4 worldToLocal; // Inverse of localToWorld
        public Matrix4x4 localToWorld; // Transformation matrix from local to world space
        public Materials.Structs.Material material;
    }
    
    // Implement triangle primitive
    
    // An abstract group of primitives that can be traced against as one object
    // Bounding box is calculated from the extremities of all primitives in the group
    [System.Serializable]
    public struct PrimitiveGroup
    {
        public int sphereStart; public int sphereCount; // Number of spheres and starting index in the buffer
        public int quadStart; public int quadCount; // Number of quads and starting index in the buffer
        public int cuboidStart; public int cuboidCount; // Number of cuboids and starting index in the buffer
        public int triangleStart; public int triangleCount; // Number of triangles and starting index in the buffer
        // public AABB boundingBox;
    }
}