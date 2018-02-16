#!/usr/bin/env python

from __future__ import print_function
import os, sys
import numpy as np

from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
#from keras.utils import np_utils
from keras import metrics
from keras import optimizers

USAGE="""
prediction.py weight_dir
"""

def create_model():
    model = Sequential()

    model.add(Dense(units=1500, input_dim=1129))
    model.add(Activation("relu"))
    model.add(Dropout(0.5))

    for i in range(0,3):
        model.add(Dense(units=1500, input_dim=1500))
        model.add(Activation("relu"))
    model.add(Dropout(0.5))

    model.add(Dense(units=4, input_dim=1500))
    model.add(Activation("sigmoid"))

    model.compile(loss='binary_crossentropy',
              optimizer='adadelta',
              metrics=['accuracy'])
    return model


def main ():

    if len(sys.argv)<4:
        print(USAGE)
        sys.exit()

    data_dir="."
    data_file=sys.argv[1]
    pred_file=sys.argv[2]
    weight_file=sys.argv[3]

    fea_npz=data_dir+"/data_feature.dat.npz"
    fea=np.load(fea_npz)    #1129 featrues

    X_pred=fea['feature']

#    N_pred=len(X_pred)

#    model = None 
    model = create_model()
    model.load_weights(weight_file)

    Y_pred= model.predict(X_pred)
    np.savez(pred_file,y_pred=Y_pred)

if __name__ == '__main__':

    main ()

