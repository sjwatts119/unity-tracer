using UnityEngine;
using Core;

namespace Geometry
{
    public class Plane : Object<ShaderStructs.Quad>
    {
        public override ShaderStructs.Quad[] ToShaderData()
        {
            // Unity's default plane is 10x10 versus a quad being 1x1, just hack in a scale factor lol
            const float planeSize = 10f;
            
            Vector3 scaledRight = transform.right * transform.localScale.y * planeSize;
            Vector3 scaledForward = transform.forward * transform.localScale.z * planeSize;
            
            // Calculate corner as bottom-left of the plane
            Vector3 corner = transform.position - (scaledRight * 0.5f) - (scaledForward * 0.5f);
            
            return new ShaderStructs.Quad[]
            {
                new ShaderStructs.Quad
                {
                    q = corner,
                    u = scaledRight,
                    v = scaledForward,
                    material = GetMaterial()
                }
            };
        }
    }
}
