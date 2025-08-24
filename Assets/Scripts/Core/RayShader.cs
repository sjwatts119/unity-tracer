using System.IO;
using UnityEngine;

namespace Core
{
    public class RayShader
    {
        private const string Name = "RayTracer/RayShader";
        private static Material _material;

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
                    throw new FileNotFoundException("Failed to load shader at " + Name);
                }
                
                return new Material(shader);
            }
        }
    }
}