# -*- coding: utf-8 -*-
import tensorflow as tf
sess = tf.compat.v1.Session()
a = tf.constant(2)
b = tf.constant(2)
c = tf.add(a, b)
print(c)

def Hello():
    d=sess.run(c)
    return d
def Add(a, b):
    return a+b