using UnityEngine;
using Core;

namespace Geometry
{
    public class Sphere : Object<ShaderStructs.Sphere>
    {
        [Header("Sphere Geometry")]
        [SerializeField, Min(0.01f)]
        private float radius = 0.5f;

        public float Radius => radius;

        public override ShaderStructs.Sphere ToShaderData()
        {
            return new ShaderStructs.Sphere
            {
                centre = transform.position,
                radius = radius,
                material = GetMaterial()
            };
        }

        protected override void OnValidate()
        {
            // Call the parent's validation first
            base.OnValidate();
            
            // Validate sphere-specific properties
            radius = Mathf.Max(0.01f, radius);
            
            if (radius > 100f)
            {
                Debug.LogWarning($"Sphere {gameObject.name} has very large radius ({radius}) - may cause performance issues");
            }
        }

        protected void OnDrawGizmos()
        {
            // Set gizmo color based on material properties
            Color gizmoColor = MaterialType == MaterialType.Light ? Emission : Albedo;
            gizmoColor.a = 1f;
            Gizmos.color = gizmoColor;
            
            // Draw sphere wireframe
            Gizmos.DrawWireSphere(transform.position, radius);
            
            // Draw solid sphere for light materials
            if (MaterialType == MaterialType.Light && Emission.maxColorComponent > 0.1f)
            {
                var emissionColor = Emission;
                emissionColor.a = 0.2f;
                Gizmos.color = emissionColor;
                Gizmos.DrawSphere(transform.position, radius);
            }
        }
    }
}