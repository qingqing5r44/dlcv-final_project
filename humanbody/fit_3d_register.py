# encoding: utf-8
"""
Copyright 2016 Max Planck Society, Federica Bogo, Angjoo Kanazawa. All rights reserved.
This software is provided for research purposes only.
By using this software you agree to the terms of the SMPLify license here:
     http://smplify.is.tue.mpg.de/license

About this Script:
============
This is a demo version of the algorithm implemented in the paper,
which fits the SMPL body model to the image given the joint detections.
The code is organized to be run on the LSP dataset.
See README to see how to download images and the detected joints.
"""

from os.path import join, exists, abspath, dirname
from os import makedirs
import logging
import cPickle as pickle
from time import time
from glob import glob
import argparse
import sys
sys.path.append('../../../smpl')
import cv2
import numpy as np
import chumpy as ch
import numpy.linalg as LA
from opendr.camera import ProjectPoints
from lib.robustifiers import GMOf
from smpl_webuser.serialization import load_model
from smpl_webuser.lbs import global_rigid_transformation
from smpl_webuser.verts import verts_decorated
from lib.sphere_collisions import SphereCollisions
from lib.max_mixture_prior import MaxMixtureCompletePrior
from render_model import render_model
import json
from scipy import sparse
from scipy import optimize 
data_path = '/media/xyz/本地磁盘/panoptic-toolbox/'
seq_name = '160224_haggling1'
hd_skel_json_path = data_path+seq_name+'/hdPose3d_stage1_coco19/hd/'
hd_face_json_path = data_path+seq_name+'/hdFace3d/'
hd_hand_json_path = data_path+seq_name+'/hdHand3d/'
from smpl_webuser.serialization import load_model

face_detect = [48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 7, 8, 9, 17, 18, 19, 20, 21, 36, 37, 38, 39, 40, 41, 22, 23, 24, 25, 26, 42, 43, 44, 45, 46, 47, 27, 28, 29, 30, 31, 32, 33, 34, 35]
fm_face_idx = [3555, 8992, 8993, 10781, 8559, 6907, 7242, 91, 7331, 10811, 10614, 3576, 9432, 10619, 10832, 294, 7349, 77, 10942, 10617, 10791, 8802, 6997, 9217, 9572, 3590, 9566, 9792, 3725, 3591, 9108, 9117, 10698, 10670, 7510, 338, 7003, 7493, 7134, 7032, 124, 7792, 8627, 2781, 2785, 408, 329, 8776, 331, 3518, 9495, 10772, 315, 35]

fm_left_idx = np.array([2264, 2210, 1985, 2226, 2113, 2250, 2743, 2707, 2764, 2718, 2122, 2135, 2124, 2161, 2037, 2027, 2307, 2329, 2271, 1999, 2372, 2368, 2360, 2396, 2439, 2431, 2188, 2198, 2580, 2480, 2514, 2477, 2574, 2519, 2770, 2193, 2604, 2591, 2631, 2589, 2690, 2652]).reshape(21, 2)
fm_right_idx =np.array( [5524, 5361, 5528, 5503, 5422, 5687, 6188, 6143, 6175, 6151, 5477, 5569, 5620, 5695, 5667, 5717, 5778, 5747, 5723, 5703, 5810, 5823, 5791, 5840, 5889, 5857, 5622, 5609, 5931, 5933, 5948, 5911, 6008, 5969, 6204, 5628, 6037, 6051, 6066, 6023, 6122, 6086]).reshape(21, 2)

m = load_model( '../../../smpl/models/basicModel_m_lbs_10_207_0_v1.0.0.pkl' )
mh = load_model( '../../../smpl/models/basicModel_m_lbs_10_207_0_v1.0.0.pkl' )
def rodrigues(x):
    return cv2.Rodrigues(x)[0]

def featureEnergyGivenshape(faceFeatureCoord, pi, lamb, t, tempFeatureCoord):
    R = rodrigues(pi)
    selectedFacePoints = faceFeatureCoord
    selectedFacePointxy = lamb * ((R.dot(selectedFacePoints.T))) + t.reshape([3,1])
    errorMatrix = selectedFacePointxy - tempFeatureCoord.T
    return LA.norm(errorMatrix, "fro")

def optRSTfromShape_t(faceFeatureCoord, R0, s0, t0, tempFeatureCoord):
    global h_6, h_7, h_8
    x0 = np.zeros([6])
    def energyFunc (x):
        pi = (x[0:3])
        lamb = 100.
        t = x[3:]
        return featureEnergyGivenshape(faceFeatureCoord, pi, lamb, t, tempFeatureCoord) 
    x0[0:3] = rodrigues(R0).reshape(3)
 
    x0[3:] = t0.reshape(3)
    x = optimize.fmin_bfgs(energyFunc, x0, maxiter = 25)
    R = rodrigues(x[0:3])
    s = 100.
    t = x[3:]
    return R, s, t   

_LOGGER = logging.getLogger(__name__)

# Mapping from LSP joints to SMPL joints.
# 0 Right ankle  8
# 1 Right knee   5
# 2 Right hip    2
# 3 Left hip     1
# 4 Left knee    4
# 5 Left ankle   7
# 6 Right wrist  21
# 7 Right elbow  19
# 8 Right shoulder 17
# 9 Left shoulder  16
# 10 Left elbow    18
# 11 Left wrist    20
# 12 Neck           -
# 13 Head top       added


# --------------------Camera estimation --------------------
def optimize_on_joints3d(j3d, model, n_betas = 10):
    prior = MaxMixtureCompletePrior(n_gaussians=8).get_gmm_prior()
    # get the mean pose as our initial pose
    init_pose = np.hstack((np.zeros(3), prior.weights.dot(prior.means)))
    # define the mapping LSP joints -> SMPL joints
    # cids are joints ids for LSP:
    cids = range(15) + [15, 16, 17, 18]
    #cids = range(15)
    # joint ids for SMPL
    # SMPL does not have a joint for head, instead we use a vertex for the head
    # and append it later.
    smpl_ids = [12, 15, 0, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8]
    # the vertex id for the joint corresponding to the head
    head_id = [122, 448, 3619, 3941]

    # weights assigned to each joint during optimization;
    # the definition of hips in SMPL and LSP is significantly different so set
    # their weights to zero
    base_weights = np.array(
        [1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1], dtype=np.float64)
    betas = ch.zeros(n_betas)


    # instantiate the model:
    # verts_decorated allows us to define how many
    # shape coefficients (directions) we want to consider (here, n_betas)
    sv = verts_decorated(
            trans=ch.zeros(3),
            pose=ch.array(init_pose),
            v_template=model.v_template,
            J=model.J_regressor,
            betas=betas,
            shapedirs=model.shapedirs[:, :, :n_betas],
            weights=model.weights,
            kintree_table=model.kintree_table,
            bs_style=model.bs_style,
            f=model.f,
            bs_type=model.bs_type,
            posedirs=model.posedirs)
    # make the SMPL joints depend on betas
    Jdirs = np.dstack([model.J_regressor.dot(model.shapedirs[:, :, i])
                           for i in range(len(betas))])
    J_onbetas = ch.array(Jdirs).dot(betas) + model.J_regressor.dot(
            model.v_template.r)

    # get joint positions as a function of model pose, betas and trans
    (_, A_global) = global_rigid_transformation(
            sv.pose, J_onbetas, model.kintree_table, xp=ch)
    Jtr = ch.vstack([g[:3, 3] for g in A_global]) + sv.trans

    # add the head joint, corresponding to a vertex...
    Jtr = ch.vstack((Jtr, sv[head_id[0]]))
    Jtr = ch.vstack((Jtr, sv[head_id[1]]))
    Jtr = ch.vstack((Jtr, sv[head_id[2]]))
    Jtr = ch.vstack((Jtr, sv[head_id[3]]))
    # ... and add the joint id to the list
    smpl_ids.append(len(Jtr) - 4)
    smpl_ids.append(len(Jtr) - 3)
    smpl_ids.append(len(Jtr) - 2)
    smpl_ids.append(len(Jtr) - 1)
    # update the weights using confidence values
    weights =  base_weights
    # data term: distance between observed and estimated joints in 2D
    #obj_j2d = lambda w: w * (weights.reshape(1, -1)).dot( ( ( j3d - Jtr[smpl_ids] ) ** 2).sum( axis = 1 )) 
    obj_j2d = lambda w, sigma:  (w * weights.reshape((-1, 1)) * GMOf((j3d - Jtr[smpl_ids]), sigma))
          # mixture of gaussians pose prior
    pprior = lambda w: w * prior(sv.pose)
        # joint angles pose prior, defined over a subset of pose parameters:
        # 55: left elbow,  90deg bend at -np.pi/2
        # 58: right elbow, 90deg bend at np.pi/2
        # 12: left knee,   90deg bend at np.pi/2
        # 15: right knee,  90deg bend at np.pi/2
    alpha = 10
    my_exp = lambda x: alpha * ch.exp(x)
    obj_angle = lambda w: w * ch.concatenate([my_exp(sv.pose[55]), my_exp(-sv.pose[
                                                 58]), my_exp(-sv.pose[12]), my_exp(-sv.pose[15])])



    opt_weights = zip([4.04 * 1e2, 4.04 * 1e2, 57.4, 4.78],
                          [1e2, 5 * 1e1, 1e1, .5 * 1e1])

    # run the optimization in 4 stages, progressively decreasing the
    # weights for the priors
    #for stage, (w, wbetas) in enumerate(opt_weights):
    objs = {}
    w = 4.04 * 1e2
    wbetas = 5 * 1e1
    objs['j2d'] = obj_j2d(1e4, 10)

    objs['pose'] = pprior(w)

    objs['pose_exp'] = obj_angle(0.317 * w)

    objs['betas'] = wbetas * betas

    ch.minimize(
                objs,
                x0=[sv.betas, sv.pose, sv.trans],
                method='dogleg',
                options={'maxiter': 100,
                         'e_3': .0001,
                         'disp': 0})
    return sv


def pointError (sets1, sets2, weights, ids):
    error = sets1[ids] - sets2
    errorT = error.dot(error.T)
    errorT = ch.diag(errorT).reshape(19, 1)
    return weights.reshape(1, 19).dot(errorT)
    

# --------------------Core optimization --------------------
def optimize_on_joints(j2d,
                       model,
                       cam,
                       img,
                       prior,
                       try_both_orient,
                       body_orient,
                       n_betas=10,
                       regs=None,
                       conf=None,
                       viz=False):
    """Fit the model to the given set of joints, given the estimated camera
    :param j2d: 14x2 array of CNN joints
    :param model: SMPL model
    :param cam: estimated camera
    :param img: h x w x 3 image 
    :param prior: mixture of gaussians pose prior
    :param try_both_orient: boolean, if True both body_orient and its flip are considered for the fit
    :param body_orient: 3D vector, initialization for the body orientation
    :param n_betas: number of shape coefficients considered during optimization
    :param regs: regressors for capsules' axis and radius, if not None enables the interpenetration error term
    :param conf: 14D vector storing the confidence values from the CNN
    :param viz: boolean, if True enables visualization during optimization
    :returns: a tuple containing the optimized model, its joints projected on image space, the camera translation
    """
    t0 = time()
    # define the mapping LSP joints -> SMPL joints
    # cids are joints ids for LSP:
    cids = range(12) + [13]
    # joint ids for SMPL
    # SMPL does not have a joint for head, instead we use a vertex for the head
    # and append it later.
    smpl_ids = [8, 5, 2, 1, 4, 7, 21, 19, 17, 16, 18, 20]

    # the vertex id for the joint corresponding to the head
    head_id = 411

    # weights assigned to each joint during optimization;
    # the definition of hips in SMPL and LSP is significantly different so set
    # their weights to zero
    base_weights = np.array(
        [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.float64)

    if try_both_orient:
        flipped_orient = cv2.Rodrigues(body_orient)[0].dot(
            cv2.Rodrigues(np.array([0., np.pi, 0]))[0])
        flipped_orient = cv2.Rodrigues(flipped_orient)[0].ravel()
        orientations = [body_orient, flipped_orient]
    else:
        orientations = [body_orient]

    if try_both_orient:
        # store here the final error for both orientations,
        # and pick the orientation resulting in the lowest error
        errors = []

    svs = []
    cams = []
    for o_id, orient in enumerate(orientations):
        # initialize the shape to the mean shape in the SMPL training set
        betas = ch.zeros(n_betas)

        # initialize the pose by using the optimized body orientation and the
        # pose prior
        init_pose = np.hstack((orient, prior.weights.dot(prior.means)))

        # instantiate the model:
        # verts_decorated allows us to define how many
        # shape coefficients (directions) we want to consider (here, n_betas)
        sv = verts_decorated(
            trans=ch.zeros(3),
            pose=ch.array(init_pose),
            v_template=model.v_template, 
            J=model.J_regressor,
            betas=betas,
            shapedirs=model.shapedirs[:, :, :n_betas],
            weights=model.weights,
            kintree_table=model.kintree_table,
            bs_style=model.bs_style,
            f=model.f,
            bs_type= model.bs_type,
            posedirs=model.posedirs)

        # make the SMPL joints depend on betas
        Jdirs = np.dstack([model.J_regressor.dot(model.shapedirs[:, :, i])
                           for i in range(len(betas))])
        J_onbetas = ch.array(Jdirs).dot(betas) + model.J_regressor.dot(
            model.v_template.r)

        # get joint positions as a function of model pose, betas and trans
        (_, A_global) = global_rigid_transformation(
            sv.pose, J_onbetas, model.kintree_table, xp=ch)
        Jtr = ch.vstack([g[:3, 3] for g in A_global]) + sv.trans

        # add the head joint, corresponding to a vertex...
        Jtr = ch.vstack((Jtr, sv[head_id]))

        # ... and add the joint id to the list
        if o_id == 0:
            smpl_ids.append(len(Jtr) - 1)

        # update the weights using confidence values
        weights = base_weights * conf[
            cids] if conf is not None else base_weights

        # project SMPL joints on the image plane using the estimated camera
        cam.v = Jtr

        # data term: distance between observed and estimated joints in 2D
        obj_j2d = lambda w, sigma: (
            w * weights.reshape((-1, 1)) * GMOf((j2d[cids] - cam[smpl_ids]), sigma))

        # mixture of gaussians pose prior
        pprior = lambda w: w * prior(sv.pose)
        # joint angles pose prior, defined over a subset of pose parameters:
        # 55: left elbow,  90deg bend at -np.pi/2
        # 58: right elbow, 90deg bend at np.pi/2
        # 12: left knee,   90deg bend at np.pi/2
        # 15: right knee,  90deg bend at np.pi/2
        alpha = 10
        my_exp = lambda x: alpha * ch.exp(x)
        obj_angle = lambda w: w * ch.concatenate([my_exp(sv.pose[55]), my_exp(-sv.pose[
                                                 58]), my_exp(-sv.pose[12]), my_exp(-sv.pose[15])])

        if viz:
            import matplotlib.pyplot as plt
            plt.ion()

            def on_step(_):
                """Create visualization."""
                plt.figure(1, figsize=(10, 10))
                plt.subplot(1, 2, 1)
                # show optimized joints in 2D
                tmp_img = img.copy()
                for coord, target_coord in zip(
                        np.around(cam.r[smpl_ids]).astype(int),
                        np.around(j2d[cids]).astype(int)):
                    if (coord[0] < tmp_img.shape[1] and coord[0] >= 0 and
                            coord[1] < tmp_img.shape[0] and coord[1] >= 0):
                        cv2.circle(tmp_img, tuple(coord), 3, [0, 0, 255])
                    if (target_coord[0] < tmp_img.shape[1] and
                            target_coord[0] >= 0 and
                            target_coord[1] < tmp_img.shape[0] and
                            target_coord[1] >= 0):
                        cv2.circle(tmp_img, tuple(target_coord), 3,
                                   [0, 255, 0])
                plt.imshow(tmp_img[:, :, ::-1])
                plt.draw()
                plt.show()
                plt.pause(1e-2)

            on_step(_)
        else:
            on_step = None

        if regs is not None:
            # interpenetration term
            sp = SphereCollisions(
                pose=sv.pose, betas=sv.betas, model=model, regs=regs)
            sp.no_hands = True
        # weight configuration used in the paper, with joints + confidence values from the CNN
        # (all the weights used in the code were obtained via grid search, see the paper for more details)
        # the first list contains the weights for the pose priors,
        # the second list contains the weights for the shape prior
        opt_weights = zip([4.04 * 1e2, 4.04 * 1e2, 57.4, 4.78],
                          [1e2, 5 * 1e1, 1e1, .5 * 1e1])

        # run the optimization in 4 stages, progressively decreasing the
        # weights for the priors
        for stage, (w, wbetas) in enumerate(opt_weights):
            _LOGGER.info('stage %01d', stage)
            objs = {}

            objs['j2d'] = obj_j2d(1e7, 100)

            objs['pose'] = pprior(w)

            objs['pose_exp'] = obj_angle(0.317 * w)

            objs['betas'] = wbetas * betas

            if regs is not None:
                objs['sph_coll'] = 1e3 * sp

            ch.minimize(
                objs,
                x0=[sv.betas, sv.pose],
                method='dogleg',
                callback=on_step,
                options={'maxiter': 100,
                         'e_3': .0001,
                         'disp': 0})

        t1 = time()
        _LOGGER.info('elapsed %.05f', (t1 - t0))
        if try_both_orient:
            errors.append((objs['j2d'].r**2).sum())
        svs.append(sv)
        cams.append(cam)

    if try_both_orient and errors[0] > errors[1]:
        choose_id = 1
    else:
        choose_id = 0
    if viz:
        plt.ioff()
    return (svs[choose_id], cams[choose_id].r, cams[choose_id].t.r)


# Select HD frame
hd_idx_set = np.arange(5757,5758)
bodyPointNum = [2, 1, 2, 2, 5, 2, 2, 5, 2, 2, 5, 2, 2, 6, 2, 1, 1, 1, 1]
row = [0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 13, 13, 13, 14, 14, 15, 16, 17, 18]
col = [3067, 1305, 331, 3510, 6541, 1291, 1891, 1650, 1622, 1666, 1661, 1725, 1961, 2242, 1512, 3108, 1055, 1012, 1005, 1523, 1148, 3322, 3202, 5354, 5287, 5125, 5133, 5130, 5127, 5091, 5669, 5565, 4399, 6538, 4532, 4495, 4491, 4548, 4505, 4500, 6735, 6730, 122, 448, 3619, 3941]
data = [0.5, 0.5, 1., 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 1.0/6, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0]
selectM = sparse.csr_matrix((np.array(data),(np.array(row), np.array(col))), shape = (19, 6890))
for hd_idx in hd_idx_set:

	# Load the json file with this frame's skeletons
	skel_json_fname = hd_skel_json_path+'body3DScene_{0:08d}.json'.format(hd_idx)
	with open(skel_json_fname) as dfile:
	    bframe = json.load(dfile)
        # Load the json file with this frame's face
        face_json_fname = hd_face_json_path+'faceRecon3D_hd{0:08d}.json'.format(hd_idx)
        with open(face_json_fname) as dfile:
	    fframe = json.load(dfile)

	hand_json_fname = hd_hand_json_path+'handRecon3D_hd{0:08d}.json'.format(hd_idx)
	with open(hand_json_fname) as dfile:
	    hframe = json.load(dfile)
        for bodyidx in xrange(len(bframe['bodies'])):
		body = bframe['bodies'][bodyidx]
		hand = hframe['people'][bodyidx]
		face = fframe['people'][bodyidx]
		righthand3d = np.array(hand['right_hand']['landmarks']).reshape((-1,3))
		skel = np.array(body['joints19']).reshape((-1,4))
		skeleton_data = skel[:, :3]
		leftHand3d = np.array(hand['left_hand']['landmarks']).reshape((-1,3))
		face3d = np.array(face['face70']['landmarks']).reshape((-1,3))[face_detect]
		firstselect = [0, 2, 3, 9, 12, 1] 
		featureNum = len(firstselect)
		smplFeaturePairs = selectM.dot(m.r)[firstselect, :]
		skeleton_data1 = skeleton_data[firstselect, :]
		A_square = np.zeros([12, 12]) 
		b_square = np.zeros([12, 1])
		for i in xrange(featureNum):
		    B = np.zeros([3, 12])
		    y = np.zeros([3,1])
		    B[0, 0:3] = smplFeaturePairs[i, :]
		    B[0, 9] = 1
		    y[:, 0] =  skeleton_data1[i, :]
		    B[1, 3:6] = smplFeaturePairs[i, :]
		    B[1, 10] = 1
		    B[2, 6:9] = smplFeaturePairs[i, :]
		    B[2, 11] = 1
		    B_square = np.dot(B.T, B)
		    y_square = np.dot(B.T, y)
		    A_square = A_square + B_square
		    b_square = b_square + y_square
		A_square_mat = np.mat(A_square)
		b_square_arr = np.array(b_square)
		x = np.linalg.solve(A_square_mat, b_square_arr)
		TT = np.array([x[0:3], x[3:6], x[6:9]]).reshape(3,3)
		b = np.array(x[9:12])

		U,S,V=LA.svd(TT)
		R = U.dot(V.T)
		s = (S[0]+S[1]+S[2])/3
		R, s, b = optRSTfromShape_t(smplFeaturePairs, R, s, b, skeleton_data1)
		TT = s * R    
		b= b.reshape([3,1])
		deformedsmplshape = (np.dot(TT, (m.r).T) +  b).T             
		f_3d = (1.0/s)*(LA.inv(R)).dot((skeleton_data.T - b))
		f_3d = f_3d.T
		f_3d = ch.array(f_3d)
		face3d = (1.0/s)*(LA.inv(R)).dot((face3d.T - b))
		face3d = face3d.T
		righthand3d = (1.0/s)*(LA.inv(R)).dot((righthand3d.T - b))
		righthand3d = righthand3d.T
		leftHand3d = (1.0/s)*(LA.inv(R)).dot((leftHand3d.T - b))
		leftHand3d = leftHand3d.T
		sv = optimize_on_joints3d(f_3d, m, 10)   
		b= b.reshape([3,1])
		deformedsmplshape = (np.dot(TT, (sv.r).T) +  b).T     
		outmesh_path = '%d/%d/deformed.obj'%(hd_idx, bodyidx)
		outmesh_path1 = '%d/%d/coef.txt'%(hd_idx, bodyidx)
		outmesh_path2 = '%d/%d/global.txt'%(hd_idx, bodyidx)
		outmesh_path3 = '%d/%d/anchorpoint.txt'%(hd_idx, bodyidx)
		outmesh_path4 = '%d/%d/faceconset.txt'%(hd_idx, bodyidx)
		outmesh_path5 = '%d/%d/handconset.txt'%(hd_idx, bodyidx)
		outmesh_path6 = '%d/%d/bodyconset.txt'%(hd_idx, bodyidx)
		globalRotate = rodrigues(np.array(sv.pose[:3])).reshape(9)
		globaltrans = np.array(sv.trans).reshape(3)
		fromPctosmpl = LA.inv(R).reshape(9)
		fromsmpltoPc = R.reshape(9)
		mh.betas[:] = np.array(sv.betas)
		mh.pose[:] = np.zeros(72)
		with open( outmesh_path, 'w') as fp:
		    for v in deformedsmplshape:
			fp.write( 'v %f %f %f\n' % ( v[0], v[1], v[2]) )

		    for f in m.f+1: # Faces are 1-based, not 0-based in obj files
			fp.write( 'f %d %d %d\n' %  (f[0], f[1], f[2]) )
		with open( outmesh_path1, 'w') as fp1:
		    for v in sv.betas:
			fp1.write( '%f \n' % ( v) )

		    for f in sv.pose[3:]: # Faces are 1-based, not 0-based in obj files
			fp1.write( '%f\n' %  (f) )
		with open( outmesh_path2, 'w') as fp2:
		    for v in globalRotate:
			fp2.write( '%f \n' % ( v) )
		    for f in fromPctosmpl: # Faces are 1-based, not 0-based in obj files
			fp2.write( '%f\n' %  (f) )
		    for f in fromsmpltoPc: # Faces are 1-based, not 0-based in obj files
			fp2.write( '%f\n' %  (f) )
		    fp2.write('%f %f %f %f %f %f %f \n'%(s, b[0,0], b[1, 0], b[2,0], globaltrans[0], globaltrans[1], globaltrans[2]))
		with open( outmesh_path3, 'w') as fp3:
		    fp3.write('1768 %f %f %f \n'%(mh.r[3500, 0], mh.r[3500, 1], mh.r[3500, 2]))
		    fp3.write('3113 %f %f %f \n'%(mh.r[3141, 0], mh.r[3141, 1], mh.r[3141, 2]))
		    fp3.write('-1')

		with open( outmesh_path4, 'w') as fp4:
		    for i in xrange(len(face_detect)):
			fp4.write( '1 %d %f %f %f\n' % ( fm_face_idx[i], face3d[i, 0], face3d[i,1], face3d[i, 2]) )
		    fp4.write('-1')	
		with open( outmesh_path5, 'w') as fp5:
		    for i in xrange(21):
			fp5.write( '2 %d %d %f %f %f\n' % ( fm_left_idx[i,0], fm_left_idx[i, 1], leftHand3d[i,0], leftHand3d[i, 1], leftHand3d[i, 2]) )
		    for i in xrange(21):
			fp5.write( '2 %d %d %f %f %f\n' % ( fm_right_idx[i,0], fm_right_idx[i, 1], righthand3d[i,0], righthand3d[i, 1], righthand3d[i, 2]) )
		    fp5.write('-1')	
		with open( outmesh_path6, 'w') as fp6:
		    for i in xrange(f_3d.shape[0]):
			fp6.write( ' %f %f %f\n' % ( f_3d[i, 0], f_3d[i,1], f_3d[i, 2]) )
		    fp6.write('-1')			
