using UnityEngine;
using Vector3 = UnityEngine.Vector3;

namespace Materials.Structs
{
    [System.Serializable]
    public enum MaterialType
    {
        Lambertian = 0,
        Metal = 1,
        Dielectric = 2,
        Light = 3,
    }

    [System.Serializable]
    public struct Material
    {
        public MaterialType type;
        public Vector3 albedo;
        public float fuzz;
        public float refractiveIndex;
        public Vector3 emission;
        public float emissionStrength;
    }
}
