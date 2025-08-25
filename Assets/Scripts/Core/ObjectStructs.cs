using Geometry;
using UnityEngine;

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
    }
}