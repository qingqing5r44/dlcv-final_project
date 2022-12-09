import sys
import trimesh
import pyrender
import numpy as np
file = str(sys.argv[1])
fuze_trimesh = trimesh.load(file)
material = pyrender.MetallicRoughnessMaterial(
    metallicFactor=0.0,
    alphaMode='OPAQUE',
    baseColorFactor=(1.0, 1.0, 0.9, 1.0))
mesh = pyrender.Mesh.from_trimesh(fuze_trimesh, material=material)
scene = pyrender.Scene(bg_color=[0.0, 0.0, 0.0, 0.0],
                        ambient_light=(0.3, 0.3, 0.3))
scene.add(mesh,'mesh')

camera = pyrender.camera.IntrinsicsCamera(
            fx=2600, fy=2600,
            cx=310, cy=330)
camera_pose = [[ 1,  0,  0, -1.04828188e-02],
               [ 0,  1,  0,  6.06489088e-03],
               [ 0,  0,  1,  1.05632458e+01],
               [ 0,  0,  0,  1]]        
scene.add(camera, pose=camera_pose)
pyrender.Viewer(scene, use_raymond_lighting=True, viewport_size = [650,600])