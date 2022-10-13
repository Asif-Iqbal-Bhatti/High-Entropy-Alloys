#!/usr/bin/env python3

'''
#==================================================
# USAGE  :: Regression analysis using TENSORFLOW
#           predicting SRO and Energy for HEAs
# AUTHOR :: ASIF IQBAL
# DATED  :: 12/10/2022

# How can linear regression be done with TensorFlow?
# Linear Regression It is the relationship between the 
# dependent and independent variables, where the independent 
# variable is denoted by the letter "x" and the dependent 
# variable is the response variable, which is denoted by the 
# letter "y". y = mx + c.
#==================================================
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow.compat.v1 as tf 
from tensorflow import keras
from tensorflow.keras import layers

print("TENSORFLOW VERSION",tf.__version__)
tf.compat.v1.disable_eager_execution()


NUM_PARALLEL_EXEC_UNITS = 2
config = tf.compat.v1.ConfigProto(
	intra_op_parallelism_threads=NUM_PARALLEL_EXEC_UNITS, 
	inter_op_parallelism_threads=2, 
	allow_soft_placement=True,
	device_count = {'CPU': NUM_PARALLEL_EXEC_UNITS})
	
rate_of_learning = 0.008
epochs_for_training = 5000

# *.CSV FILENAME TO READ DATA
orig = pd.read_csv('VASPRun_x0.20_800_1st_2nd_3rd.csv')
orig_sor = orig.sort_values(by=['TE/at'], ascending=True)

x_data = orig_sor["TE/at"].astype("float")
y_data = orig_sor["Total_SRO"].astype("float")

length_x = len(x_data)
New_X_data = tf.placeholder(tf.float32, shape=(None, ), name = 'x_data') 
New_Y_data = tf.placeholder(tf.float32, shape=(None, ), name = 'y_data') 

var_W = tf.Variable(np.random.normal(), name = "W") 
var_b = tf.Variable(np.random.normal(), name = "b")

y_predictions = tf.add(tf.multiply(New_X_data, var_W), var_b) 
function_cost = tf.reduce_mean(tf.square(y_predictions-New_Y_data)) # OR LOSS
optimizer_gradient = tf.train.GradientDescentOptimizer(rate_of_learning).minimize(function_cost) 


with tf.Session(config=config) as sess: 
	sess.run(tf.global_variables_initializer())
		
	for epoch in range(epochs_for_training):
		sess.run(optimizer_gradient, feed_dict = {New_X_data: x_data, New_Y_data: y_data})
		if (epoch + 1) % 500 == 0:
			c = sess.run(function_cost, feed_dict = {New_X_data : x_data, New_Y_data : y_data}) 
			print("Epoch", (epoch + 1), ": cost =", c, "var_1 =", sess.run(var_W), "var_2 =", sess.run(var_b)) 
      
	cost_training = sess.run(function_cost, feed_dict = {New_X_data: x_data, New_Y_data: y_data}) 
	weight = sess.run(var_W) 
	bias = sess.run(var_b)

#print(type(weight), type(bias), type(x_data))
make_predictions = weight * x_data + bias 
print(f'Training cost = {cost_training}, Weight = {weight}, bias = {bias}')
plt.plot(x_data, y_data, 'bo',label ='Sample data taken') 
plt.plot(x_data, make_predictions, label ='line fitted') 
plt.title('Result for linear regression') 
plt.legend() 
plt.savefig("test.png", dpi=400,)
plt.show()


