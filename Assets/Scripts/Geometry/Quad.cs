using UnityEngine;
using Core;

namespace Geometry
{
    public class Quad : Object<ShaderStructs.Quad>
    {
        public override ShaderStructs.Quad[] ToShaderData()
        {
            // Get the scaled right and up vectors
            Vector3 scaledRight = transform.right * transform.localScale.x;
            Vector3 scaledUp = transform.up * transform.localScale.y;
            
            // Calculate corner as bottom-left of the quad
            Vector3 corner = transform.position - (scaledRight * 0.5f) - (scaledUp * 0.5f);
            
            return new ShaderStructs.Quad[]
            {
                new ShaderStructs.Quad
                {
                    q = corner,
                    u = scaledRight,
                    v = scaledUp,
                    material = GetMaterial()
                }
            };
        }
    }
}
