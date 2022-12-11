# -*- coding: utf-8 -*-
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
#from sklearn.externals import joblib
import tensorflow as tf

hhh = np.load('male_body_training_data_all.npz')
Input_ = hhh['coef_input']
Out_error = hhh['DE_error']
train_x_disorder, test_x_disorder, train_y_disorder, test_y_disorder = train_test_split (Input_, Out_error, train_size = 0.9, random_state = 33)
#ss_x = preprocessing.StandardScaler()
#train_x_disorder = ss_x.fit_transform(train_x_disorder)
#test_x_disorder = ss_x.transform(test_x_disorder)



n_input = 410
n_output = 65380
n_hidden_1 = 1500  # 1st layer num features


tf_body_x = tf.placeholder(tf.float32, [None,n_input], name = 'tf_body_x')     # input x
tf_body_y = tf.placeholder(tf.float32, [None,n_output], name = 'tf_body_y')     # input y



# neural network layers
b_l1 = tf.layers.dense(tf_body_x, n_hidden_1, tf.nn.relu, name = 'e_l1')          # hidden layer
b_output = tf.layers.dense(b_l1, n_output, name = 'e_output')                     # output layer
b_loss=tf.reduce_mean(tf.reduce_sum(tf.square(tf_body_y-b_output),
                                  reduction_indices=[1]), name = 'loss')
sample_size = train_x_disorder.shape[0]
batch = 1000
start = 0
end = batch
#loss = tf.losses.mean_squared_error(tf_y, output)   # compute cost
global_step = tf.Variable(0)
learning_rate = tf.train.exponential_decay(3e-2, global_step, decay_steps = 200, decay_rate = 0.995, staircase = True)
#learning_rate = tf.placeholder(tf.float32)
optimizer = tf.train.GradientDescentOptimizer(learning_rate)
train_op = optimizer.minimize(b_loss, global_step= global_step)

sess = tf.Session()                                 # control training and others
sess.run(tf.global_variables_initializer())         # initialize var in graph

total = 15000
lr = 3e-7
for step in range(total):
    
    # train and net output

    _, l, pred = sess.run([train_op, b_loss, b_output], {tf_body_x: train_x_disorder[start:end, :], tf_body_y: train_y_disorder[start: end, :]})
    if step % 50 == 0:
        print(str(step) + ', start:' + str(start)+ ', end: ' + str(end))
        print('loss is: ' + str(l))
        print('train loss:' + str(sess.run(b_loss, {tf_body_y:train_y_disorder, tf_body_x: train_x_disorder})))
        print('test loss:' + str(sess.run(b_loss, {tf_body_y:test_y_disorder, tf_body_x: test_x_disorder})))
        # print('prediction is:' + str(pred))
    start = end if end<sample_size else 0
    end = start + batch

saver = tf.train.Saver()
save_path = saver.save(sess, 'male_deep_network/deerror1.ckpt')
# 输出保存路径
print ('Save to path: ', save_path)
