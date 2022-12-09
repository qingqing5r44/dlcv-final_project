# encoding: utf-8
"""
Created on Wed Nov 14 08:26:32 2018

@author: asus
"""

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
import numpy.linalg as LA
import joblib
import tensorflow as tf
import sys
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"]="3"
from lbs import global_rigid_transformation
sys.path.append('E:\Graphics\bachelor_thesis\Code\splocs-master')
sys.path.append('E:\Graphics\bachelor_thesis\Code\humanbody\humanbody')
C_body = np.load('female_pose_training_data_all.npz')['C']
gmm = joblib.load('gmm.pkl')
motion2posemodel = joblib.load('woman_pose_regression.pkl')
pose2posmodel = joblib.load("Fromcof2posmodel.pickle")
fromposeparam2faceanchormodel = joblib.load("fromposeparam2faceanchor.pkl")
hhh = np.load('smpl_param.npz')
shape_dirs = hhh['shapedirs']
v_template = hhh['v_template']
kintree_table = hhh['kintree_table']
J_regressor = np.load('J_regressor2.npz')['J_regressor']
#hhh = np.load('male_body_training_data_all.npz')
#Input_ = hhh['coef_input']
#Out_error = hhh['DE_error']
#train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder = train_test_split (Input_, Out_error, train_size = 0.9, random_state = 33)
import numpy.linalg as LA
import cv2
from scipy import optimize 

hhh = np.load('pose2faceanchordataset.npz')
Input = hhh['Input']
Output = hhh['Output']
train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder = train_test_split (Input, Output, train_size = 0.9, random_state = 33)

ss_x = preprocessing.StandardScaler()
train_x_disorder = ss_x.fit_transform(train_x_disorder)
test_x_disorder = ss_x.transform(test_x_disorder)
#
ss_y = preprocessing.StandardScaler()
train_y_disorder = ss_y.fit_transform(train_y_disorder)
test_y_disorder = ss_y.transform(test_y_disorder)
J_dirs = np.dstack([J_regressor.dot(shape_dirs[:, :, i]) for i in range(10)])


#rg = joblib.load('face_anchor_predict.pkl')
#hhh = np.load('DE_training_data.npz')
#Input_ = hhh['Input_']
#Out_error = hhh['Out_error']
#C_body = hhh['C']
#model = joblib.load('pose_regression_model.m')
#print(4)
#train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder = train_test_split (Input_, Out_error, train_size = 0.9, random_state = 33)
#ss_x = preprocessing.StandardScaler()
#train_x_disorder = ss_x.fit_transform(train_x_disorder)
#test_x_disorder = ss_x.transform(test_x_disorder)
#
#ss_y = preprocessing.StandardScaler()
#train_y_disorder = ss_y.fit_transform(train_y_disorder)
#test_y_disorder = ss_y.transform(test_y_disorder)
#del train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder, Input_, Out_error

def getJointsPos (x, y):
    betas = np.zeros(10)
    pose = np.zeros(72)
    betas[:] = np.array(x)
    pose[3:] = np.array(y)
    J_onbetas = J_dirs.dot(betas) + J_regressor.dot(v_template)
    (_, A_global) = global_rigid_transformation(pose, J_onbetas, kintree_table, xp=np)
    Jtr = np.vstack([g[:3, 3] for g in A_global])
    Jtr = Jtr.reshape(72)
    return list(Jtr)
    

def optShapePoseFromFeaturePoints(FeaturePos):
    idxset = [12, 16, 18, 20, 1, 4, 7, 17, 19, 21, 2, 5, 8, 24, 25, 26, 27]
    FeaturePos = np.array(FeaturePos).reshape(17, 3)
    predict = lambda x:pose2posmodel.predict(x.reshape(1, 410)).reshape(28, 3)[idxset,:]
    rotate = lambda x: rodrigues(x[:3]).dot(predict(x[6:]).T) + x[3:6].reshape(3, 1)
    computeEnergy1 = lambda x: LA.norm(rotate(x) - FeaturePos.T, "fro")
    computeEnergy2 = lambda x: gmm.score(x.reshape(1, 400))
    computeEnergy =lambda x: 1e5 * computeEnergy1(x) + 1 * computeEnergy2(x[16:])
    x0 = np.zeros(416)
    t = optimize.fmin_bfgs(computeEnergy, x0, maxiter = 30)
    return t

def rodrigues(x):
    return cv2.Rodrigues(x)[0]
def featureEnergyGivenshape(faceFeatureCoord, pi, t, tempFeatureCoord):
    R = rodrigues(pi)
    selectedFacePoints = faceFeatureCoord
    selectedFacePointxy = ((R.dot(selectedFacePoints.T))) + t.reshape([3,1])
    errorMatrix = selectedFacePointxy - tempFeatureCoord.T
    return LA.norm(errorMatrix, "fro")

def findRigidtransform(faceFeatureCoord, R0, t0, tempFeatureCoord):
    x0 = np.zeros([6])
    def energyFunc (x):
        pi = (x[0:3])
        t = x[3:]
        return featureEnergyGivenshape(faceFeatureCoord, pi, t, tempFeatureCoord) 
    faceFeatureCoord = np.array(faceFeatureCoord).reshape(-1, 3)
    tempFeatureCoord = np.array(tempFeatureCoord).reshape(-1, 3)
    A_square = np.zeros([12, 12]) 
    b_square = np.zeros([12, 1])
    for i in range(faceFeatureCoord.shape[0]):
        B = np.zeros([3, 12])
        y = np.zeros([3,1])
        B[0, 0:3] = faceFeatureCoord[i, :]
        B[0, 9] = 1
        y[:, 0] =  tempFeatureCoord[i, :]
        B[1, 3:6] = faceFeatureCoord[i, :]
        B[1, 10] = 1
        B[2, 6:9] = faceFeatureCoord[i, :]
        B[2, 11] = 1
        B_square = np.dot(B.T, B)
        y_square = np.dot(B.T, y)
        A_square = A_square + B_square
        b_square = b_square + y_square
    A_square_mat = np.mat(A_square)
    b_square_arr = np.array(b_square)
    xl = np.linalg.solve(A_square_mat, b_square_arr)
    TT = np.array([xl[0:3], xl[3:6], xl[6:9]]).reshape(3,3)
    b = np.array(xl[9:12])
    U,S,V=LA.svd(TT)
    R = U.dot(V.T)
    R0 = np.array(R).reshape(3, 3)
    x0[0:3] = rodrigues(R0).reshape(3)
    x0[3:] = np.array(b).reshape(3)
    x = optimize.fmin_bfgs(energyFunc, x0, maxiter = 50)
    R = rodrigues(x[0:3])
    t = x[3:]
    return list(R.reshape(9)), list(t)   

def load_body_motion1():

    import tensorflow.compat.v1 as tf
    #ss_L_y = preprocessing.StandardScaler()
    #train_L_y_disorder = ss_L_y.fit_transform(train_L_y_disorder)
    #test_L_y_disorder = ss_L_y.transform(test_L_y_disorder)
    tf_x = tf.placeholder(tf.float32, [None,79], name = 'tf_x')     # input x
    tf_y = tf.placeholder(tf.float32, [None,65380], name = 'tf_y')     # input y
    
    # neural network layers
    l1 = tf.layers.dense(tf_x, 1000, tf.nn.relu, name = 'l1')          # hidden layer
    output = tf.layers.dense(l1, 65380, name = 'output')                     # output layer
    loss=tf.reduce_mean(tf.reduce_sum(tf.square(tf_y-output),
    								  reduction_indices=[1]), name = 'loss')
    sess = tf.Session()                                 # control training and others
    sess.run(tf.global_variables_initializer())         # initialize var in graph
    saver = tf.train.Saver()
    saver.restore(sess, './model1/deerror.ckpt') 
    return sess, tf_x,  output
	
def load_body_motion():
	n_input = 410
	n_output = 65380
	n_hidden_1 = 1500  # 1st layer num features

	tf_body_x = tf.placeholder(tf.float32, [None,n_input], name = 'tf_body_x')     # input x
	tf_body_y = tf.placeholder(tf.float32, [None,n_output], name = 'tf_body_y')     # input y

	# neural network layers
	b_l1 = tf.layers.dense(tf_body_x, n_hidden_1, tf.nn.relu, name = 'e_l1')          # hidden layer
	b_output = tf.layers.dense(b_l1, n_output, name = 'e_output')                     # output layer
	b_loss=tf.reduce_mean(tf.reduce_sum(tf.square(tf_body_y-b_output),reduction_indices=[1]), name = 'loss')
	sess = tf.Session()                                 # control training and others
	sess.run(tf.global_variables_initializer())         # initialize var in graph
	saver = tf.train.Saver()
	save_path = saver.restore(sess, 'femalemodel/deerror1.ckpt')
#    de = tf.matmul(tf_body_x[:,10:], np.array(C_body,dtype = 'float32'))
#    de_Input = tf.placeholder(tf.float32, [None,65380], name = 'deinput')
#    de_ref = tf.placeholder(tf.float32, [None,65380], name = 'de_ref')
#    de_output = tf.add(de,b_output )
	return sess, tf_body_x, b_output

def outputSkeletonPos(x, y):
    t = np.zeros([1, 410])
    t[0, :10] = np.array(x).reshape(10)
    t[0, 10:] = np.array(y).reshape(400)
    h = pose2posmodel.predict(t).reshape(84)
    return list(h)

def outputPoseParm(x, y):
    t = np.zeros([1, 79])
    t[0, :10] = np.array(x).reshape(10)
    t[0, 10:] = np.array(y).reshape(69)
    h = motion2posemodel.predict(t).reshape(400)
    return list(h)

def recons1(x, tf_x, output, sess):
    #x = Input_[frame, :].reshape([1, 79])
    x = np.array(x).reshape(1, 79)
    x_trans = ss_x.transform(x)
    y = model.predict(x)
    DE = y.dot(C_body)
    y_l = sess.run(output, {tf_x: x_trans})
    y_trans = 1 *ss_y.inverse_transform(y_l)
    DE_predict = DE + y_trans
    DE_predict = DE_predict.reshape(65380)
    return list(DE_predict)
	
def recons(x, tf_x, output, sess):
    #x = Input_[2, :].reshape(1, 410)
    x = np.array(x).reshape(1, 410)
    DE = x[[0], 10:].dot(C_body)
    y_l = sess.run(output, {tf_x: x})
    DE_predict = DE + y_l
    DE_predict = DE_predict.reshape(65380)
    return list(DE_predict)
    
def reconsFace1(x):
	x = np.array(x).reshape(1, 79)
	y_l = rg.predict(x)
	y = y_l.reshape(12, 3)
	y_L = np.zeros([12, 3, 3])
	for i in range(12):
		y_L[i, :, :] = rodrigues(y[i, :].reshape(3,1))
	y_L = y_L.reshape(36 * 3)
	return list(y_L)
	
def reconsFace(x):
	x = ss_x.transform(np.array(x).reshape(1, 410))
	y_l = fromposeparam2faceanchormodel.predict(x)
	y = (ss_y.inverse_transform(y_l)).reshape(12, 3)
	y_L = np.zeros([12, 3, 3])
	for i in range(12):
		y_L[i, :, :] = rodrigues(y[i, :].reshape(3,1))
	y_L = y_L.reshape(36 * 3)
	return list(y_L)

def optimizex(x0, DEopt, tf_x, output, sess):
    t0 = np.array(x0)[10:]
    DEopt = np.array(DEopt)
    def computEnergy(x):
        l = np.zeros([410])
        l[10:] = x
        l[:10] = x0[:10]
        error = np.array(recons(l, tf_x, output, sess)) - DEopt
        return error.T.dot(error)
    t = optimize.fmin_bfgs(computEnergy, t0, maxiter = 20)
    t1 = np.zeros(410)
    t1[:10] = x0[:10]
    t1[10:] = t
    y_l = recons(t1, tf_x, output, sess)
    return list(t1), y_l
        

