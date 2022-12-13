#!/usr/bin/env python3

'''
#==================================================
# USAGE  :: DNN using TENSORFLOW
#           predicting SRO and Energy for HEAs
# AUTHOR :: ASIF IQBAL
# DATED  :: 12/10/2022
#==================================================
'''

import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import tensorflow.compat.v1 as tf 
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow.keras import Sequential
from tensorflow.keras import optimizers
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras.utils import plot_model
from tensorflow.keras.callbacks import TensorBoard

print("TENSORFLOW VERSION",tf.__version__)
tf.disable_eager_execution()
logdir="logs/fit/" + datetime.now().strftime("%Y%m%d-%H%M%S")
tensorboard_callback = TensorBoard(log_dir=logdir)

# *.CSV FILENAME TO READ DATA
df = pd.read_csv('VASPRun_x0.20_800_1st_2nd_3rd.csv')
labels=df['Total_SRO']
features = df.iloc[:,0:2]
X=features
#df = np.loadtxt('VASPRun_x0.20_800_1st_2nd_3rd.csv', dtype='float',delimiter=',',skiprows=1)
#X=df[:,:]
y=np.ravel(labels)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.80, random_state=42)

scaler = StandardScaler().fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)

#========= MODEL
model = Sequential()
model.add(Dense(64, activation='tanh', input_shape=(2,)))
model.add(Dropout(0.2))
model.add(Dense(48, activation='tanh'))
model.add(Dense(24, activation='tanh'))
model.add(Dense(1, activation='tanh'))

opt = optimizers.SGD(lr = 0.01)

model.compile(
loss='binary_crossentropy',
optimizer=opt,
metrics=['mean_squared_error'],)

model.fit(X_train, y_train, epochs=900, batch_size=100, verbose=1, shuffle=False,
callbacks=[tensorboard_callback], validation_data=(X_test, y_test))

y_pred = model.predict(X_test)
score = model.evaluate(X_test, y_test, verbose=1)
print(score)

#========= PLOTTING
plot_model(
    model,
    to_file="model.png",
    show_shapes=True,
    show_dtype=False,
    show_layer_names=True,
    rankdir="TB",
    expand_nested=True,
    dpi=300,
)

model.save("model.h5")
#tensorboard --logdir=logs/ --host localhost --port 8088
