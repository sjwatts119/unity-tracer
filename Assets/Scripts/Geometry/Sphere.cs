using UnityEngine;

namespace Geometry
{
    [System.Serializable]
    public class Sphere : MonoBehaviour
    {
        public float radius = 0.5f;
        public Color color = Color.magenta;
    
        // Convert to struct for shader use
        public Core.ShaderStructs.Sphere ToShaderData()
        {
            return new Core.ShaderStructs.Sphere
            {
                centre = transform.position,
                radius = radius
            };
        }
    }
}

