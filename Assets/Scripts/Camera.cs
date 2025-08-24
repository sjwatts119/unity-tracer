using System.IO;
using Core;
using UnityEngine;

[ExecuteAlways, ImageEffectAllowedInSceneView]
public class CameraOverride : MonoBehaviour
{
    public void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        var mat = RayShader.GetMaterial();

        if (!mat)
        {
            throw new FileNotFoundException("Failed to load shader at " + RayShader.Name);
        }
        
        Graphics.Blit(source, destination, mat);
    }
}