library(reticulate) 
use_condaenv("/path/anaconda3/envs/DeepPASTA/")
setwd('/path/polyAseqTrap/DeepIP_path/')
## --------- train_1w on test.all.fa --------
##Test data: test.all.fa, length 200, with PA located at position 101.
##seqLabel='' indicates that labels such as :1 (real PA) and :0 (not a real PA) are provided at the end of the sequence header.
##If the sequence header does not have labels, you can use something like seqLabel=1 to indicate that all provided sequences are real PA.
cmd=sprintf("testSeq='%s'; trainedModel='%s'; outputFile='%s'; seqLabel=''",
            'test.all.fa',
            'ath.train.model.hdf5',
            'ath.train.model_ON_test.all.csv')
py_run_string(cmd)
py_run_file(DeepIP_TEST) 








