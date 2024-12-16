#!/usr/bin/env python

## usage ##
# DeepIP_test
# (1) label 0/1 in >title of testSeq
# python DeepIP_test.py -testSeq Newtest.fa -trainedModel Newtrain.train.model_I.hdf5 -outputFile Newtest.predicted_I.csv
# (2) provided labels and not use label 1/0 in trainSeq; half=1, half=0
# python DeepIP_test.py -testSeq Newtest.fa -trainedModel Newtrain.train.model_I.hdf5 -outputFile Newtest.predicted_II.csv -seqLabel 1/0
# (3) half=0, half=1
# python DeepIP_test.py -testSeq Newtest.fa -trainedModel Newtrain.train.model_I.hdf5 -outputFile Newtest.predicted_III.csv -seqLabel 0/1
# (4) all=0
# python DeepIP_test.py -testSeq Newtest.fa -trainedModel Newtrain.train.model_I.hdf5 -outputFile Newtest.predicted_IV.csv -seqLabel 0
# (5) all=1
# python DeepIP_test.py -testSeq Newtest.fa -trainedModel Newtrain.train.model_I.hdf5 -outputFile Newtest.predicted_IV.csv -seqLabel 1

## Run from R
# py_run_string("testSeq='Newtest.fa'; trainedModel='Newtrain.train.model.hdf5'; outputFile='Newtest.predicted_III.csv'; seqLabel='1/0'")
# py_run_file("DeepIP_test.py") 

## input ##
# testSeq: if not have negative_0/...e_1 at the end, then set all label to 1
# trainedModel: .hdf5
# outputFile: .csv
# seqLabel (optional): default ''

## Output ##
# Newtest.predicted.csv: title, score, true_label(0/1), predict_label(0/1), res(TP/FP/TN/FN)


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
from keras.layers import LSTM
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.layers import PReLU
from sklearn.metrics import average_precision_score, roc_auc_score
import pandas as pd
#from matplotlib import pyplot as plt

from keras.callbacks import ModelCheckpoint

from sklearn.metrics import average_precision_score, roc_auc_score,confusion_matrix,accuracy_score,precision_score,recall_score,f1_score
np.random.seed(1337)
random.seed(1337)


## output each sequence: title, score, true_label(0/1), predict_label(0/1), res(TP/FP/TN/FN)
def addtag(DataFrame, file):
    result = DataFrame
    
    result['true_label'] = result['true_label'].astype(int)

    conditions = [((result['predict_label'] == 1) & (result['true_label'] == 1)),
                  ((result['predict_label'] == 1) & (result['true_label'] == 0)),
                  ((result['predict_label'] == 0) & (result['true_label'] == 1)),
                  ((result['predict_label'] == 0) & (result['true_label'] == 0))
                  ]
                  
    choice = ['TP', 'FP', 'FN', 'TN']
    result['res'] = np.select(conditions, choice, default='invalid')

    result.to_csv(file, index=False)

    return 0

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
# testSeq = ''
# trainedModel = ''
# outputFile = ''
# seqLabel = ''

print(">>> Start DeepIP_test...")

for i in range(len(sys.argv)):
    if i < len(sys.argv) - 1:
        if sys.argv[i] == '-testSeq' or sys.argv[i] == '-test':
            testSeq = sys.argv[i + 1]
        if sys.argv[i] == '-outputFile' or sys.argv[i] == '-output':
            outputFile = sys.argv[i + 1]
        if sys.argv[i] == '-trainedModel' or sys.argv[i] == '-model':
            trainedModel = sys.argv[i + 1]
        if sys.argv[i] == '-seqLabel' or sys.argv[i] == '-label':
            seqLabel = sys.argv[i + 1]

if not isVarDefined('trainedModel'):
    print('Please provide trainedModel (.hdf5)!')
    sys.exit();  
if not isVarDefined('seqLabel'):
  seqLabel=''  
if not isVarDefined('testSeq'):
    print('Please provide testSeq (.fa)!')
    sys.exit();  
if not isVarDefined('outputFile'):
    print('Please provide outputFile (.csv)!')
    sys.exit();      
if trainedModel == '':
    print('Please provide trainedModel (.hdf5)!')
    sys.exit();  
if outputFile == '':
    print('Please provide outputFile (.csv)!')
    sys.exit();  
if testSeq == '':
    print('Please provide testSeq (.fa)!')
    sys.exit();  

if outputFile[-3:]!='csv':
  outputFile=outputFile+'.csv'

seqLabel=str(seqLabel)

print(">>> testSeq (.fa): " + testSeq)
print(">>> trainedModel (.hdf5): " + trainedModel)
print(">>> outputFile (.csv): " + outputFile)
print(">>> seqLabel: " + seqLabel)

## check seqLabel
if seqLabel == '':
  print('>>>>>> seqLabel not provided, will read true labels from last char of fa header.')
else:
  print('>>>>>> seqLabel provided, will use true labels as given.')
  lbls=['0/1','0:1','01','1/0','1:0','10','00','0:0','0/0','11','1/1','1:1','1','0']
  if not any(lbl in seqLabel for lbl in lbls):
    print(f"Error seqLabel, should be like {lbls}!")
    sys.exit()

testingSequenceList = []
testingLabelList = []
testingTitleList = []

##1.1process testSeq
if seqLabel != '':
  labeled=0  #use provided labels
else:
  labeled=titleLabeled(testSeq)
  if labeled==0:
    print("seqLabel not provided, but trainSeq is not labeled (last char=0/1)!")
    sys.exit()
    
for line in open(testSeq):
    info = line[0:(len(line) - 1)]
    if '>' in info:
      if labeled == 1: #label from the seq title
        testingLabelList.append(int(info[-1:]))
      testingTitleList.append(info)
    else:
      testingSequenceList.append(info)
  
# user provided labels, half 0 + half 1, or all 0/0
nseq=len(testingSequenceList)
if labeled==0:
  nhalf=int(nseq/2)
  if seqLabel == '0' or seqLabel == '00' or seqLabel == '0:0' or seqLabel == '0/0':
    testingLabelList=[0]*nseq
  if seqLabel == '1' or seqLabel == '11' or seqLabel == '1/1' or seqLabel == '1:1':
    testingLabelList=[1]*nseq
  if seqLabel == '01' or seqLabel == '0/1' or seqLabel == '0:1':
    if nseq % 2 !=0:
      print(f"seqLabel is like 0/1, but the sequence count ({nseq}) is not even (half 0; half 1)!")
      sys.exit()
    testingLabelList=[0]*nhalf + [1]*nhalf
  if seqLabel == '10' or seqLabel == '1/0' or seqLabel == '1:0':
    if nseq % 2 !=0:
      print(f"seqLabel is like 1/0, but the sequence count ({nseq}) is not even (half 1; half 01)!")
      sys.exit()
    testingLabelList=[1]*nhalf + [0]*nhalf

if len(testingLabelList) != len(testingSequenceList):
  print(f"sequence count ({nseq}) is not equal to label count ({len(testingLabelList)})!")
  sys.exit()

## testSeq is all sequence, then add >s1/0 like titles
if len(testingTitleList)==0:
  prelist = []
  for i in range(1, nseq):
    prelist.append(">s"+str(i)+'/'+str(testingLabelList[i-1]))

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

model.load_weights(trainedModel)

print ('================= predicting =====================')
test_predict_result = model.predict(testingData, batch_size=2042, verbose=0)

## process the predict result
test_predict_result = test_predict_result.tolist()  #process the predict score
test_predict_result=[x[0] for x in test_predict_result] #strip the '[]'

DF=pd.DataFrame(data=testingTitleList,columns=['title']) #add title
DF['score']=test_predict_result   #add score
DF['true_label']=testingLabelList 
DF['predict_label']=DF['score'].apply(lambda x:1 if x>0.5 else 0) #add tag
addtag(DF, outputFile)

print ('confusion_matrix: [0/1; 0/1]' )
print ( str(confusion_matrix(testingLabelList, DF['predict_label'].values)) )
print ('Acc: ' + str(accuracy_score(testingLabelList, DF['predict_label'].values)))
print ('precison: ' + str(precision_score(testingLabelList, DF['predict_label'].values)))
print ('F1: ' + str(f1_score(testingLabelList, DF['predict_label'].values)))
print ('recall: ' + str(recall_score(testingLabelList, DF['predict_label'].values)))
print ('AUC: ' + str(roc_auc_score(testingLabelList, test_predict_result)))
print ('AUPRC: ' + str(average_precision_score(testingLabelList, test_predict_result)))

print("Finish DeepIP_test.")
