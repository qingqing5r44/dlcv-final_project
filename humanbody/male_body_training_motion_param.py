# -*- coding: utf-8 -*-
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
#from sklearn.neural_network import MLPRegressor
from sklearn.externals import joblib
model = joblib.load("pose_regression_model.m")
#rg = joblib.load('face_anchor_model.m')
hhh = np.load('DE_training_data.npz')
Input_ = hhh['Input_']
Out_error = hhh['Out_error']
C = hhh['C']
train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder = train_test_split (Input_, Out_error, train_size = 0.9, random_state = 33)

#train_D_y_disorder = train_y_disorder[:, 2 * np.arange(32690)]
#train_L_y_disorder = train_y_disorder[:, 1 + 2 * np.arange(32690)]
#
#test_D_y_disorder = test_y_disorder[:, 2 * np.arange(32690)]
#test_L_y_disorder = test_y_disorder[:, 1 + 2 * np.arange(32690)]

#数据标准化
ss_x = preprocessing.StandardScaler()
train_x_disorder = ss_x.fit_transform(train_x_disorder)
test_x_disorder = ss_x.transform(test_x_disorder)

ss_y = preprocessing.StandardScaler()
train_y_disorder = ss_y.fit_transform(train_y_disorder)
test_y_disorder = ss_y.transform(test_y_disorder)

#ss_L_y = preprocessing.StandardScaler()
#train_L_y_disorder = ss_L_y.fit_transform(train_L_y_disorder)
#test_L_y_disorder = ss_L_y.transform(test_L_y_disorder)



import tensorflow as tf

        

tf_x = tf.placeholder(tf.float32, [None,79], name = 'tf_x')     # input x
tf_y = tf.placeholder(tf.float32, [None,65380], name = 'tf_y')     # input y

# neural network layers
l1 = tf.layers.dense(tf_x, 1000, tf.nn.relu, name = 'l1')          # hidden layer
output = tf.layers.dense(l1, 65380, name = 'output')                     # output layer
loss=tf.reduce_mean(tf.reduce_sum(tf.square(tf_y-output),
                                  reduction_indices=[1]), name = 'loss')

sample_size = train_x_disorder.shape[0]
batch = 1000
start = 0
end = batch
#loss = tf.losses.mean_squared_error(tf_y, output)   # compute cost
global_step = tf.Variable(0)
learning_rate = tf.train.exponential_decay(3e-3, global_step, decay_steps = 200, decay_rate = 0.99, staircase = True)
optimizer = tf.train.GradientDescentOptimizer(learning_rate)
train_op = optimizer.minimize(loss, global_step= global_step)

sess = tf.Session()                                 # control training and others
sess.run(tf.global_variables_initializer())         # initialize var in graph
saver = tf.train.Saver()
#saver.restore(sess, './model/deerror.ckpt') 
total = 8000

for step in range(total):
    
    # train and net output
    _, l, pred = sess.run([train_op, loss, output], {tf_x: train_x_disorder[start:end, :], tf_y: train_y_disorder[start: end, :]})
    if step % 50 == 0:
        print(str(step) + ', start:' + str(start)+ ', end: ' + str(end))
        print('loss is: ' + str(l))
        print('train loss:' + str(sess.run(loss, {tf_y:train_y_disorder, tf_x: train_x_disorder})))
        print('test loss:' + str(sess.run(loss, {tf_y:test_y_disorder, tf_x: test_x_disorder})))
        # print('prediction is:' + str(pred))
    start = end if end<sample_size else 0
    end = start + batch

saver = tf.train.Saver()
saver.save(sess, './model/deerror.ckpt')

#model MLP
#rg = MLPRegressor(tol = 0.00001, solver = 'lbfgs', hidden_layer_sizes = (100, 200, 400,), activation='logistic',max_iter=7000)
#rg.fit(train_x_disorder, train_y_disorder)
#mlp_score = rg.score(test_x_disorder, test_y_disorder)
#print('mlp: ', mlp_score)
#mlp_score = rg.score(train_x_disorder, train_y_disorder)
#print('mlp: ', mlp_score)

def recons():
    for frame in range(4312):
        x = Input_[frame, :].reshape([1, 79])
        x_trans = ss_x.transform(x)
        y = model.predict(x)
        DE = y.dot(C)
        y_l = sess.run(output, {tf_x: x_trans})
        y_trans = 1 *ss_y.inverse_transform(y_l)
        DE_predict = DE + y_trans
        DE_predict = DE_predict.reshape(65380)
        f = open("F:\project\humanbody\humanbody\\de\\%04d.txt"%frame, 'w')
            #f = open("F:\\data\\" + moduleType +"\\IICData\\" + "\\dC"+ str(nComp) +".txt", 'w')
        for i in DE_predict:    
            f.write(str(i) + "\n")
        f.close()
    
#model_gbr_GridSearch = GradientBoostingRegressor()
#param_grid = {'n_estimators':range(20, 81, 10),
#              'learning_rate': [0.2, 0.1, 0.05, 0.02, 0.01],
#              'max_depth':[4, 6, 8],
#              'min_samples_leaf': [3, 5, 9, 14],
#              'max_features':[0.8, 0.5, 0.3, 0.1]}
#from sklearn.model_selection import GridSearchCV
#estimator = GridSearchCV(model_gbr_GridSearch, param_grid)
#estimator.fit(train_x_disorder, train_y_disorder)
#print('zuiyou: ',estimator.best_params_)
#print('score ', estimator.score(test_x_disorder, test_y_disorder))


