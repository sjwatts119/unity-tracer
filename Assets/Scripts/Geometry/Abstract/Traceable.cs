using UnityEngine;
using Vector3 = UnityEngine.Vector3;
using Materials.Structs;

namespace Geometry.Abstract
{
    public abstract class Traceable : UnityEngine.MonoBehaviour
    {
        [Header("Material Properties")]
        [SerializeField]
        private MaterialType materialType = MaterialType.Lambertian;
        
        [Header("Base Properties")]
        [SerializeField, ColorUsage(false)]
        private Color albedo = Color.white;
        
        [Header("Metal Properties")]
        [SerializeField, Range(0f, 1f)]
        private float fuzz = 0.0f;
        
        [Header("Dielectric Properties")]
        [SerializeField, Range(1.0f, 3.0f)]
        private float refractiveIndex = 1.5f;
        
        [Header("Emission Properties")]
        [SerializeField, ColorUsage(false)]
        private Color emission = Color.black;
        
        [Header("Emission Properties")]
        [SerializeField, Range(0f, 30f)]
        private float emissionStrength = 1f;

        protected Materials.Structs.Material GetMaterial()
        {
            return new Materials.Structs.Material
            {
                type = materialType,
                albedo = new Vector3(albedo.r, albedo.g, albedo.b),
                fuzz = fuzz,
                refractiveIndex = refractiveIndex,
                emission = new Vector3(emission.r, emission.g, emission.b),
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
                case MaterialType.Lambertian:
                    ValidateLambertian();
                    break;
                case MaterialType.Metal:
                    ValidateMetal();
                    break;
                case MaterialType.Dielectric:
                    ValidateDielectric();
                    break;
                case MaterialType.Light:
                    ValidateLight();
                    break;
                default:
                    Debug.LogError($"Unknown material type on {gameObject.name}");
                    break;
            }
        }

        private void ValidateLambertian()
        {
            emission = Color.black;
            emissionStrength = 0f;
            fuzz = 0f;
            refractiveIndex = 1f;
            
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b),
                albedo.a
            );
        }

        private void ValidateMetal()
        {
            emission = Color.black;
            emissionStrength = 0f;
            refractiveIndex = 1f;
            
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b),
                albedo.a
            );
            
            fuzz = Mathf.Clamp01(fuzz);
        }

        private void ValidateDielectric()
        {
            emission = Color.black;
            emissionStrength = 0f;
            fuzz = 0f;
            
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b),
                albedo.a
            );
        }
        
        private void ValidateLight()
        {
            fuzz = 0f;
            refractiveIndex = 1f;
            
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b)
            );
            
            emission = new Color(
                Mathf.Max(0f, emission.r),
                Mathf.Max(0f, emission.g),
                Mathf.Max(0f, emission.b),
                emission.a
            );
        }
    }
}