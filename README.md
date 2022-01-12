This is the implementation and dataset "App Miscategorization Detection: A Case Study on Google Play". Didi Surian, Suranga Seneviratne, Aruna Seneviratne and Sanjay Chawla. IEEE Transactions on Knowledge and Data Engineering 29 (8), 1591-1604

------------------
Codes
------------------
- core.m
  Topic model based on the von Mises-Fisher (vMF) distribution
- utils (folder):
  All supporting Matlab scripts (data conversion, synthetic data generator, etc.)

-------------------------------------------------------------------------
NOTES
-------------------------------------------------------------------------
The experiments were done using:
- MATLAB ver. 8.2.0.701 (R2013b)
- LIBSVM 3.18: 
       http://www.csie.ntu.edu.tw/~cjlin/libsvm/    
  We have used matlab interface for LIBSVM
- Stanford Topic Modeling toolbox v0.4.0: 
       http://nlp.stanford.edu/software/tmt/tmt-0.4/
- k-means++ vers. 2013-02-08

The experiments were performed on a machine with Intel Core(TM) Duo CPU T6400 @2.00 GHz,1.75 GB or RAM, running on Microsoft Windows XP SP3.

Please refer to the paper about all parameter settings.
