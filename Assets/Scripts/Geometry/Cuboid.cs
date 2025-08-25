using UnityEngine;
using Core;

namespace Geometry
{
    public class Cuboid : Object<ShaderStructs.Quad>
    {
        public override ShaderStructs.Quad[] ToShaderData()
        {
            var material = GetMaterial();
            Vector3 pos = transform.position;
            Vector3 scale = transform.localScale;
            
            var quads = new ShaderStructs.Quad[6];
            
            Vector3 right = transform.right * scale.x;
            Vector3 up = transform.up * scale.y;
            Vector3 forward = transform.forward * scale.z;
            
            Vector3 halfRight = right * 0.5f;
            Vector3 halfUp = up * 0.5f;
            Vector3 halfForward = forward * 0.5f;
            
            // Front face
            quads[0] = new ShaderStructs.Quad
            {
                q = pos + halfForward - halfRight - halfUp,
                u = right,
                v = up,
                material = material
            };
            
            // Back face
            quads[1] = new ShaderStructs.Quad
            {
                q = pos - halfForward + halfRight - halfUp,
                u = -right,
                v = up,
                material = material
            };
            
            // Right face
            quads[2] = new ShaderStructs.Quad
            {
                q = pos + halfRight - halfForward - halfUp,
                u = forward,
                v = up,
                material = material
            };
            
            // Left face
            quads[3] = new ShaderStructs.Quad
            {
                q = pos - halfRight + halfForward - halfUp,
                u = -forward,
                v = up,
                material = material
            };
            
            // Top face
            quads[4] = new ShaderStructs.Quad
            {
                q = pos + halfUp - halfRight - halfForward,
                u = right,
                v = forward,
                material = material
            };
            
            // Bottom face
            quads[5] = new ShaderStructs.Quad
            {
                q = pos - halfUp - halfRight + halfForward,
                u = right,
                v = -forward,
                material = material
            };
            
            return quads;
        }
    }
}