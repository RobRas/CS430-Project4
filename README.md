This program uses a raytracer to create 3D images from a json file of objects. The image is of PPM P6 format. This version includes spot lights and point lights, with diffuse and specular reflection. Objects can both reflect and refract colors.

To run: raytrace width height input.json output.ppm

The input file should have one camera object. It supports up to 128 additional spheres and planes, as well as 128 additional light sources.

This program was written by Robert Rasmussen - rsr47
