using System.Numerics;
using Geometry;
using Vector3 = UnityEngine.Vector3;

namespace Core
{
    public static class ShaderStructs
    {
        [System.Serializable]
        public struct Material
        {
            public MaterialType type;
            public Vector3 albedo;
            public float fuzz;
            public float refractiveIndex;
            public Vector3 emission;
        }
        
        [System.Serializable]
        public struct Sphere : HasMaterial
        {
            public Vector3 centre;
            public float radius;
            public Material material { get; set; }
        }

        [System.Serializable]
        public struct Quad : HasMaterial
        {
            public Vector3 q; // Bottom-left corner
            public Vector3 u; // Edge vector along the width ( q ---u--- )
            public Vector3 v; // Edge vector along the height (rotate above example 90 degrees)
            public Material material { get; set; }
        }
    }
}