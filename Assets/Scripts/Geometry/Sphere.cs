using UnityEngine;
using Core;

namespace Geometry
{
    public class Sphere : Object<ShaderStructs.Sphere>
    {
        private float CalculateRadius()
        {
            Vector3 scale = transform.localScale;
            return Mathf.Max(scale.x, scale.y, scale.z) * 0.5f;
        }

        public override ShaderStructs.Sphere ToShaderData()
        {
            return new ShaderStructs.Sphere
            {
                centre = transform.position,
                radius = CalculateRadius(),
                material = GetMaterial()
            };
        }
    }
}