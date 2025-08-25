using Core;
using UnityEngine;

namespace Geometry
{
    public abstract class Object<T> : MonoBehaviour where T : struct, HasMaterial
    {
        [Header("Material Properties")]
        [SerializeField] private MaterialType materialType = MaterialType.Lambertian;
        
        [Header("Base Properties")]
        [SerializeField, ColorUsage(false)]
        private Color albedo = Color.white;
        
        [Header("Metal Properties")]
        [SerializeField, Range(0f, 1f)]
        private float fuzz = 0.0f;
        
        [Header("Glass Properties")]
        [SerializeField, Range(0.1f, 3.0f)]
        private float refractiveIndex = 1.5f;
        
        [Header("Emission Properties")]
        [SerializeField, ColorUsage(false)]
        private Color emission = Color.black;

        public MaterialType MaterialType => materialType;
        public Color Albedo => albedo;
        public float Fuzz => fuzz;
        public float RefractiveIndex => refractiveIndex;
        public Color Emission => emission;
        
        public abstract T ToShaderData();
        
        protected ShaderStructs.Material GetMaterial()
        {
            return new ShaderStructs.Material
            {
                type = materialType,
                albedo = new Vector3(albedo.r, albedo.g, albedo.b),
                fuzz = fuzz,
                refractiveIndex = refractiveIndex,
                emission = new Vector3(emission.r, emission.g, emission.b)
            };
        }

        protected virtual void OnValidate()
        {
            ValidateMaterial();
        }

        protected void ValidateMaterial()
        {
            switch (materialType) 
            {
                case MaterialType.Lambertian:
                    ValidateLambertian();
                    break;
                case MaterialType.Metal:
                    ValidateMetal();
                    break;
                case MaterialType.Glass:
                    ValidateGlass();
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
            // Lambertian shouldn't have fuzz or emission or a non-default refractive index
            emission = Color.black;
            fuzz = 0f;
            refractiveIndex = 1f;
            
            // Clamp albedo components to valid range [0,1]
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b),
                albedo.a
            );
        }

        private void ValidateMetal()
        {
            // Metal shouldn't have emission or a non-default refractive index
            emission = Color.black;
            refractiveIndex = 1f;
            
            // Clamp albedo components to valid range [0,1]
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b),
                albedo.a
            );
            
            fuzz = Mathf.Clamp01(fuzz);
        }

        private void ValidateGlass()
        {
            // Glass shouldn't have fuzz or emission
            emission = Color.black;
            fuzz = 0f;
            
            // Clamp albedo components to valid range [0,1]
            albedo = new Color(
                Mathf.Clamp01(albedo.r),
                Mathf.Clamp01(albedo.g),
                Mathf.Clamp01(albedo.b),
                albedo.a
            );
        }
        
        private void ValidateLight()
        {
            // Light shouldn't have fuzz or a non-default refractive index
            fuzz = 0f;
            refractiveIndex = 1f;
            
            // Clamp albedo components to valid range [0,1]
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