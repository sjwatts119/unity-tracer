using System.IO;
using UnityEngine;

namespace Shaders
{
    public class RayShader
    {
        private const string Name = "RayTracer/RayShader";
        private static Material _material;

        // Get a cached material instance from the shader
        public static Material Material
        {
            get
            {
                if (_material)
                {
                    return _material;
                }

                var shader = Shader.Find(Name);
                
                if (!shader) 
                {
                    throw new FileNotFoundException("Failed to load shader " + Name);
                }
                
                return new Material(shader);
            }
        }
    }
}