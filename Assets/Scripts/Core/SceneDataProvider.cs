using System;
using System.Collections.Generic;
using System.Linq;
using Geometry.Abstract;
using Geometry.Structs;
using UnityEngine;
using Object = UnityEngine.Object;

namespace Core
{
public class SceneDataProvider
{
    public SceneData GatherSceneData()
    {
        var data = new SceneData();
        
        // Get all individual primitives
        var spheres = Object.FindObjectsByType<Geometry.Sphere>(FindObjectsSortMode.None);
        var quads = Object.FindObjectsByType<Geometry.Quad>(FindObjectsSortMode.None);
        var planes = Object.FindObjectsByType<Geometry.Plane>(FindObjectsSortMode.None);
        var cuboids = Object.FindObjectsByType<Geometry.Cuboid>(FindObjectsSortMode.None);
        
        // Convert to struct arrays
        var individualSpheres = spheres.Select(sphere => (Sphere)sphere.ToPrimitive()).ToArray();
        var individualQuads = quads.Select(quad => (Quad)quad.ToPrimitive())
                                 .Concat(planes.Select(plane => (Quad)plane.ToPrimitive()))
                                 .ToArray();
        var individualCuboids = cuboids.Select(cuboid => (Cuboid)cuboid.ToPrimitive()).ToArray();
        
        // Get ALL TraceableGroups
        var groups = Object.FindObjectsByType<TraceableGroup>(FindObjectsSortMode.None);
        
        // Track grouped primitives
        var groupedSpheresList = new List<Sphere>();
        var groupedQuadsList = new List<Quad>();
        var groupedCuboidsList = new List<Cuboid>();
        var groupedTrianglesList = new List<Triangle>();
        var primitiveGroups = new List<PrimitiveGroup>();
        
        // Process each group
        foreach (var group in groups)
        {
            // Get primitives from this group
            var primitiveArrays = group.GetPrimitives();
            
            // Only create a group if it has primitives
            if (primitiveArrays.spheres.Length == 0 && 
                primitiveArrays.quads.Length == 0 && 
                primitiveArrays.cuboids.Length == 0 && 
                primitiveArrays.triangles.Length == 0)
                continue;
            
            // Create the group with current indices
            var primitiveGroup = new PrimitiveGroup
            {
                sphereStart = groupedSpheresList.Count,
                sphereCount = primitiveArrays.spheres.Length,
                quadStart = groupedQuadsList.Count,
                quadCount = primitiveArrays.quads.Length,
                cuboidStart = groupedCuboidsList.Count,
                cuboidCount = primitiveArrays.cuboids.Length,
                triangleStart = groupedTrianglesList.Count,
                triangleCount = primitiveArrays.triangles.Length,
                boundingBox = group.ToAABB()
            };
            
            // Add primitives to the grouped lists
            groupedSpheresList.AddRange(primitiveArrays.spheres);
            groupedQuadsList.AddRange(primitiveArrays.quads);
            groupedCuboidsList.AddRange(primitiveArrays.cuboids);
            groupedTrianglesList.AddRange(primitiveArrays.triangles);
            
            primitiveGroups.Add(primitiveGroup);
        }
        
        // Assign to data structure
        data.individualSpheres = individualSpheres;
        data.individualQuads = individualQuads;
        data.individualCuboids = individualCuboids;
        data.individualTriangles = Array.Empty<Triangle>(); // No individual triangles for now
        
        data.groupedSpheres = groupedSpheresList.ToArray();
        data.groupedQuads = groupedQuadsList.ToArray();
        data.groupedCuboids = groupedCuboidsList.ToArray();
        data.groupedTriangles = groupedTrianglesList.ToArray();
        
        data.primitiveGroups = primitiveGroups.ToArray();
        
        return data;
    }
}
}