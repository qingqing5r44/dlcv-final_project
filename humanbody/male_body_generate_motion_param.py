# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 08:58:40 2018

@author: asus
"""

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.neural_network import MLPRegressor
from sklearn.externals import joblib
import numpy.linalg as LA
import cv2
from scipy import optimize 
def rodrigues(x):
    return cv2.Rodrigues(x)[0]

rg = joblib.load('face_anchor_predict.pkl')
hhh = np.load('DE_training_data.npz')
Input_ = hhh['Input_']
Out_error = hhh['Out_error']
C_body = hhh['C']
model = joblib.load('pose_regression_model.m')
print(4)
train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder = train_test_split (Input_, Out_error, train_size = 0.9, random_state = 33)
ss_x = preprocessing.StandardScaler()
train_x_disorder = ss_x.fit_transform(train_x_disorder)
test_x_disorder = ss_x.transform(test_x_disorder)

ss_y = preprocessing.StandardScaler()
train_y_disorder = ss_y.fit_transform(train_y_disorder)
test_y_disorder = ss_y.transform(test_y_disorder)
del train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder, Input_, Out_error

def load_body_motion():

    import tensorflow as tf
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
    saver.restore(sess, './model/deerror.ckpt') 
    return sess, tf_x,  output

#def load_body_motion():
#	sess = tf.Session()                                 # control training and others
#	sess.run(tf.global_variables_initializer())         # initialize var in graph
#	saver = tf.train.Saver()
#	saver.restore(sess, './model/deerror.ckpt')
#	return sess
        
#def load():
#    sess = tf.Session()                                 # control training and others
#    sess.run(tf.global_variables_initializer())         # initialize var in graph
#    saver = tf.train.Saver()
#    saver.restore(sess, './model/deerror.ckpt')
#    return sess

#model MLP
#rg = MLPRegressor(tol = 0.00001, solver = 'lbfgs', hidden_layer_sizes = (100, 200, 400,), activation='logistic',max_iter=7000)
#rg.fit(train_x_disorder, train_y_disorder)
#mlp_score = rg.score(test_x_disorder, test_y_disorder)
#print('mlp: ', mlp_score)
#mlp_score = rg.score(train_x_disorder, train_y_disorder)
#print('mlp: ', mlp_score)
    
def reconsFace(x):
	x = np.array(x).reshape(1, 79)
	y_l = rg.predict(x)
	y = y_l.reshape(12, 3)
	y_L = np.zeros([12, 3, 3])
	for i in range(12):
		y_L[i, :, :] = rodrigues(y[i, :].reshape(3,1))
	y_L = y_L.reshape(36 * 3)
	return list(y_L)

def recons(x, tf_x, output, sess):
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
#    f = open("F:\project\humanbody\humanbody\\de\\%04d.txt"%frame, 'w')
#        #f = open("F:\\data\\" + moduleType +"\\IICData\\" + "\\dC"+ str(nComp) +".txt", 'w')
#    for i in DE_predict:    
#        f.write(str(i) + "\n")
#    f.close()
    
def optimizeXY(x0, y0, de_target,  tf_x, output, sess):
    def compute(l):
        t = np.zeros([1, 79])
        t[0, :10]  = x0
        t[0, 10:] = l
        x_trans = ss_x.transform(t)
        y = model.predict(t)
        de = y.dot(C_body)
        y_l = sess.run(output, {tf_x:x_trans})
        y_trans = 1 * ss_y.inverse_transform(y_l)
        DE_predict = (de + y_trans).reshape(65380)
        return LA.norm(DE_predict - de_target, "fro")
    t0= np.zeros([69])
    t0[:] = y0
    t_final = optimize.fmin_bfgs(compute, t0, maxiter = 25)
    x_final = np.zeros([1, 79])
    x_final[0, :10] = x0
    x_final[0, 10:] = t_final
    return recons(x_final, tf_x, output, sess), t_final
        
    
    