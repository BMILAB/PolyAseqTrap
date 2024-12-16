#!/usr/bin/env python

## usage ##
## DeepIP_train (~30m 8k seqs)
# (1) label 0/1 in >title of trainSeq
# python DeepIP_train.py -trainSeq Newtrain.fa -trainedModel Newtrain.train.model_I.hdf5 -epoch 1

# (2) provided labels half=1 / half=0,  and not use label 1/0 in trainSeq
# python DeepIP_train.py -trainSeq Newtrain.fa -trainedModel Newtrain.train.model_II.hdf5 -seqLabel 1/0 -epoch 1

# (3) provided labels half=0 / half=1, and not use label 01 in trainSeq
# python DeepIP_train.py -trainSeq Newtrain.fa -trainedModel Newtrain.train.model_III.hdf5 -seqLabel 01 -epoch 1

# (4) set saveFreq to save check point
# python DeepIP_train.py -trainSeq Newtrain.fa -trainedModel Newtrain.train.model_III.hdf5 -seqLabel 01 -epoch 10 -saveFreq 2
# will output files like: <trainmodel>.e001.hdf5, e.g., Newtrain.train.model_III.e001.hdf5

## (4) Run from R
# py_run_string("trainSeq='Newtrain.fa'; trainedModel='Newtrain.train.model_I.hdf5'; epoch=1; seqLabel='1/0'")
# py_run_file("DeepIP_train.py") 

## input ##
# trainSeq: should have negative_0/...e_1 at the end
# trainedModel (optional): default = DeepIP_trained_model.hdf5
# epoch (optional): default 100
# seqLabel (optional): default ''

## Output ##
# Newtrain.train.model.hdf5
# Newtrain.train.model.hdf5.history.csv

import pydoc
import sys
import numpy as np
import math
import random

from keras.preprocessing import sequence
from keras.models import Model
from keras.layers import Bidirectional, Input, concatenate, add
from keras.layers.core import Dense, Dropout, Flatten
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import MaxPooling1D, AveragePooling1D
#from keras.layers.recurrent import LSTM
from keras.layers import LSTM
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.layers import PReLU
#from keras.layers.advanced_activations import PReLU
from sklearn.metrics import average_precision_score, roc_auc_score
import pandas as pd
#from matplotlib import pyplot as plt

from keras.callbacks import ModelCheckpoint

from sklearn.metrics import average_precision_score, roc_auc_score,confusion_matrix,accuracy_score,precision_score,recall_score,f1_score
np.random.seed(1337)
random.seed(1337)


def oneHotEncodingForSeq(rawSeqList):
    if len(rawSeqList) != 0:
        encodedSeq = np.zeros((len(rawSeqList), len(rawSeqList[0]), 5))
        for i in range(len(rawSeqList)):
            sequence = rawSeqList[i]
            j = 0
            for s in sequence:
                if s == 'A' or s == 'a':
                    encodedSeq[i][j] = [1, 0, 0, 0, 0]
                elif s == 'T' or s == 't':
                    encodedSeq[i][j] = [0, 1, 0, 0, 0]
                elif s == 'C' or s == 'c':
                    encodedSeq[i][j] = [0, 0, 1, 0, 0]
                elif s == 'G' or s == 'g':
                    encodedSeq[i][j] = [0, 0, 0, 1, 0]
                elif s == 'N' or s == 'n':
                    encodedSeq[i][j] = [0, 0, 0, 0, 1]
                else:
                    print>> sys.stderr, 'ERROR: Unwanted nucleotide: ' + s
                j = j + 1
        return encodedSeq
    else:
        return 0

def sequenceModel(seqInput):
    seqCov = Conv1D(filters=512,
                    kernel_size=8,
                    padding="valid",
                    input_shape=(200, 5),
                    activation="relu",
                    strides=1)(seqInput)

    seqPool = MaxPooling1D(pool_size=3, strides=3)(seqCov)
    seqDout1 = Dropout(rate=0.7)(seqPool)
    seqBiLstm = Bidirectional(LSTM(units=128, return_sequences=True))(seqDout1)
    seqDout2 = Dropout(rate=0.7)(seqBiLstm)
    seqFlat = Flatten()(seqDout2)
    seqDen2 = Dense(256, kernel_initializer='glorot_uniform', activation='relu')(seqFlat)
    seqDout4 = Dropout(rate=0.7)(seqDen2)

    return seqDout4

## check last two chars in a seqfile, return 0 (not label), 1 (is label)  
def titleLabeled(seqfile):
  yes=0
  with open(seqfile, "r") as file:
    # Read the first line and strip any leading or trailing whitespace
    header = file.readline().strip()
    if header[0]=='>':
      if header[-2:]=='/1' or header[-2:]=='/0' or header[-2:]=='_1' or header[-2:]=='_0' or header[-2:]==':1' or header[-2:]==':0':
        yes=1
  return yes

## check the existence of a var 
def isVarDefined(var_name):
    vars_dict = globals() if var_name in globals() else locals()
    if var_name in vars_dict:
        return True
    else:
        return False

### main function ####

## arguments ##
# trainSeq
# trainedModel
# seqLabel
# epoch

print("Starting DeepIP_train ...")
for i in range(len(sys.argv)):
    if i < len(sys.argv) - 1:
        if sys.argv[i] == '-trainSeq' or sys.argv[i] == '-train':
            trainSeq = sys.argv[i + 1]
            
        #.hdf5
        if sys.argv[i] == '-trainedModel' or sys.argv[i] == '-model':
            trainedModel = sys.argv[i + 1]
            
        if sys.argv[i] == '-epoch' or sys.argv[i] == '-N':
            epoch = int(sys.argv[i + 1])
            
        if sys.argv[i] == '-seqLabel' or sys.argv[i] == '-label':
            seqLabel = sys.argv[i + 1]
            
        if sys.argv[i] == '-saveFreq' or sys.argv[i] == '-F':
            saveFreq = int(sys.argv[i + 1])

if not isVarDefined('trainedModel'):
  trainedModel = 'DeepIP_trained_model.hdf5'
  
if not isVarDefined('epoch'):
  epoch=100
  
if not isVarDefined('seqLabel'):
  seqLabel=''  

if not isVarDefined('saveFreq'):
  saveFreq=epoch
  
if saveFreq>epoch:
  saveFreq=epoch
  
if not isVarDefined('trainSeq'):
    print('Please provide trainSeq for model training!')
    sys.exit();
  
if trainedModel[-4:]!='hdf5':
  trainedModel=trainedModel+'.hdf5'

trainedRes=trainedModel+".history.csv"

print(">>> trainSeq (.fa): " + trainSeq)
print(">>> trainedModel (.hdf5): " + trainedModel)
print(">>> training history (.csv): " + trainedRes)
print(">>> epoch (iteration number): " + str(epoch))
print(">>> saveFreq (save check point): " + str(saveFreq))
print(">>> seqLabel: " + seqLabel)

## check seqLabel
if seqLabel == '':
  print('>>>>>> seqLabel not provided, will read true labels from last char of fa header.')
else:
  print('>>>>>> seqLabel provided, will use true labels as given.')
  if seqLabel != '0/1' and seqLabel != '1/0' and seqLabel != '0:1' and seqLabel != '1:0' and seqLabel != '01' and seqLabel != '10':
    print('Error: seqLabel should be like 0/1, 1/0, 0:1, 1:0, 01, 10!')
    sys.exit()

if trainSeq == '':
    print('Please provide trainSeq for model training!')
    sys.exit();
    
##1.1process trainSeq
if seqLabel != '':
  labeled=0  #use provided labels
else:
  labeled=titleLabeled(trainSeq)
  if labeled==0:
    print("seqLabel not provided, but trainSeq is not labeled (last char=0/1)!")
    sys.exit()

testingSequenceList = []
testingLabelList = []
for line in open(trainSeq):
  info = line[0:(len(line) - 1)]
  if '>' in info:
    if labeled == 1:
      testingLabelList.append(int(info[-1:]))
  else:
    testingSequenceList.append(info)

# user provided labels, half 0 half 1
if labeled==0:
  nseq=len(testingSequenceList)
  if nseq % 2 !=0:
    print(f"seqLabel is provided, but the sequence count ({nseq}) is not even (half true; half false)!")
    sys.exit()
  nhalf=int(nseq/2)
  if seqLabel == '0/1' or seqLabel == '0:1' or seqLabel == '01':
    testingLabelList=[0]*nhalf + [1]*nhalf
  else:
    testingLabelList=[1]*nhalf + [0]*nhalf

if len(testingLabelList) != len(testingSequenceList):
  print(f"sequence count ({nseq}) is not equal to label count ({len(testingLabelList)})!")
  sys.exit()

encodedTestingSeq = oneHotEncodingForSeq(testingSequenceList) #get the seq encoding to binary

testingData = []
testingData.append(encodedTestingSeq)

# Building deep learning model
training_net = []
# deep learning sub-model for sequence
seqInput = Input(shape=(200, 5))
seqModel = sequenceModel(seqInput)

den1 = Dense(256, kernel_initializer='glorot_uniform', activation='relu')(seqModel)
dout1 = Dropout(rate=0.7)(den1)
den2 = Dense(128, kernel_initializer='glorot_uniform', activation='relu')(dout1)
dout2 = Dropout(rate=0.7)(den2)
den3 = Dense(64, activation='relu')(dout2)
den4 = Dense(1, activation='sigmoid')(den3)
# model = Model(inputs=[seqInput, ssInput1, ssInput2, ssInput3], outputs=den4)
model = Model(inputs=seqInput, outputs=den4)
#  Graph disconnected: cannot obtain value for tensor Tensor
model.compile(loss='binary_crossentropy', optimizer='nadam', metrics=['accuracy'])

###save the best model,need chanage the fix name
# filepath = 'train_weights_improment_{epoch:02d}_{acc:.2f}.hdf5'
#checkpoint = ModelCheckpoint(filepath, monitor='acc', verbose=1,save_weights_only=True, mode='max',period=2)

#tensorflow:`period` argument is deprecated. Please use `save_freq` to specify the frequency in number of batches seen.

# save check point when saveFreq is set.
if saveFreq<epoch:
  trainedModel = trainedModel[:-4]+"e{epoch:03d}.hdf5"
  checkpoint = ModelCheckpoint(trainedModel, monitor='acc', verbose=1, save_weights_only=True, mode='max', save_freq=saveFreq)
else:
  checkpoint = ModelCheckpoint(trainedModel, monitor='acc', verbose=1, save_weights_only=True, mode='max', save_freq=epoch) 

# verbose:show progress; if monitor='val_acc'/acc,then mode select 'max',if monitor select loss,the mode select 'min';period indicate interval   savy_best_only=True
callbacks_list = [checkpoint] #need to put callback_list
#the first method : fit the model
testingLabelList=np.array(testingLabelList)
result=model.fit(testingData,testingLabelList,batch_size=2042,epochs=epoch,validation_split=0.33,
          shuffle=1,verbose=1,callbacks=callbacks_list)

# loss	accuracy	val_loss	val_accuracy
DF=pd.DataFrame(result.history)
DF.to_csv(trainedRes, index=False)
print(">>> Output training history  to ", trainedRes)


