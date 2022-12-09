from geometry import perspective_projection

model_joints = 
rotation = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
camera_t = [-0.14640908,  0.21162955,  8.92637   ]
focal_length = 5000
camera_center = [800., 600.]
projected_joints = perspective_projection(model_joints, rotation, camera_t, focal_length, camera_center)

#