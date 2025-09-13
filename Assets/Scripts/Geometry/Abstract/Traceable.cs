using UnityEngine;
using Vector3 = UnityEngine.Vector3;
using Materials.Structs;

namespace Geometry.Abstract
{
    public abstract class Traceable : UnityEngine.MonoBehaviour
    {
        [Header("Material Properties")]
        [SerializeField]
        private MaterialType materialType = MaterialType.Solid;
        
        [Header("Solid Properties")]
        [SerializeField, ColorUsage(false)]
        private Color colour = Color.white;
        
        [Header("Solid Properties")]
        [SerializeField, Range(0f, 1f)]
        private float roughness = 0.0f;
        
        [Header("Solid Properties")]
        [SerializeField, Range(0f, 1f)]
        private float reflectivity = 0.5f;
        
        [Header("Dielectric Properties")]
        [SerializeField, Range(1.0f, 3.0f)]
        private float refractiveIndex = 1.5f;

        [Header("Dielectric Properties")] 
        [SerializeField, Range(0f, 3f)]
        private float absorptionStrength = 0f;
        
        [Header("Emission Properties")]
        [SerializeField, Range(0f, 30f)]
        private float emissionStrength = 1f;

        protected Materials.Structs.Material GetMaterial()
        {
            return new Materials.Structs.Material
            {
                type = materialType,
                colour = new Vector3(colour.r, colour.g, colour.b),
                reflectivity = reflectivity,
                roughness = roughness,
                refractiveIndex = refractiveIndex,
                absorptionStrength = absorptionStrength,
                emissionStrength = emissionStrength
            };
        }

        protected virtual void OnValidate()
        {
            ValidateMaterial();
        }

        private void ValidateMaterial()
        {
            switch (materialType) 
            {
                case MaterialType.Solid:
                    ValidateSolid();
                    break;
                case MaterialType.Dielectric:
                    ValidateDielectric();
                    break;
                case MaterialType.Emissive:
                    ValidateEmissive();
                    break;
                default:
                    Debug.LogError($"Unknown material type on {gameObject.name}");
                    break;
            }
        }

        private void ValidateSolid()
        {
            absorptionStrength = 0f;
            emissionStrength = 0f;
            refractiveIndex = 1f;
            roughness = Mathf.Clamp01(roughness);
            
            colour = new Color(
                Mathf.Clamp01(colour.r),
                Mathf.Clamp01(colour.g),
                Mathf.Clamp01(colour.b),
                colour.a
            );
        }

        private void ValidateDielectric()
        {
            emissionStrength = 0f;
            roughness = 0f;
            reflectivity = 0f;
            
            colour = new Color(
                Mathf.Clamp01(colour.r),
                Mathf.Clamp01(colour.g),
                Mathf.Clamp01(colour.b),
                colour.a
            );
        }
        
        private void ValidateEmissive()
        {
            refractiveIndex = 1f;
            roughness = 0f;
            reflectivity = 0f;
            absorptionStrength = 0f;
            
            colour = new Color(
                Mathf.Clamp01(colour.r),
                Mathf.Clamp01(colour.g),
                Mathf.Clamp01(colour.b)
            );
        }
    }
}