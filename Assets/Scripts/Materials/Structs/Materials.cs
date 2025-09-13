using UnityEngine;
using Vector3 = UnityEngine.Vector3;

namespace Materials.Structs
{
    [System.Serializable]
    public enum MaterialType
    {
        Solid = 0,
        Dielectric = 1,
        Emissive = 2,
    }
    
    [System.Serializable]
    public struct Material
    {
        public MaterialType type;
        public Vector3 colour;
        
        public float reflectivity; // Likelihood of diffuse reflection vs specular reflection
        public float roughness; // Affects the spread of specular reflections, does not affect diffuse reflections
        
        public float refractiveIndex; // The refractive index of the material
        public float absorptionStrength; // The amount the colour attenuates the ray
        
        public float emissionStrength; // The brightness of the emission
    }
}
