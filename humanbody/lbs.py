'''
Copyright 2015 Matthew Loper, Naureen Mahmood and the Max Planck Gesellschaft.  All rights reserved.
This software is provided for research purposes only.
By using this software you agree to the terms of the SMPL Model license here http://smpl.is.tue.mpg.de/license

More information about SMPL is available here http://smpl.is.tue.mpg.
For comments or questions, please email us at: smpl@tuebingen.mpg.de


About this file:
================
This file defines linear blend skinning for the SMPL loader which 
defines the effect of bones and blendshapes on the vertices of the template mesh.

Modules included:
- global_rigid_transformation: 
  computes global rotation & translation of the model
- verts_core: [overloaded function inherited from verts.verts_core]
  computes the blending of joint-influences for each vertex based on type of skinning

'''


import numpy as np

def global_rigid_transformation(pose, J, kintree_table, xp):
    results = {}
    pose = pose.reshape((-1,3))
    id_to_col = {kintree_table[1,i] : i for i in range(kintree_table.shape[1])}
    parent = {i : id_to_col[kintree_table[0,i]] for i in range(1, kintree_table.shape[1])}


    import cv2
    rodrigues = lambda x : cv2.Rodrigues(x)[0]

    with_zeros = lambda x : xp.vstack((x, xp.array([[0.0, 0.0, 0.0, 1.0]])))
    results[0] = with_zeros(xp.hstack((rodrigues(pose[0,:]), J[0,:].reshape((3,1)))))        
        
    for i in range(1, kintree_table.shape[1]):
        results[i] = results[parent[i]].dot(with_zeros(xp.hstack((
            rodrigues(pose[i,:]),
            ((J[i,:] - J[parent[i],:]).reshape((3,1)))
            ))))

    pack = lambda x : xp.hstack([np.zeros((4, 3)), x.reshape((4,1))])
    
    results = [results[i] for i in sorted(results.keys())]
    results_global = results

    if True:
        results2 = [results[i] - (pack(
            results[i].dot(xp.concatenate( ( (J[i,:]), [0] ) )))
            ) for i in range(len(results))]
        results = results2
    result = xp.dstack(results)
    return result, results_global


    
