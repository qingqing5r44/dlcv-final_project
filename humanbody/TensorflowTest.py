# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:57:03 2018

@author: asus
"""

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