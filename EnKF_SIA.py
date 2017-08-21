# -*- coding: utf-8 -*-
"""
Created on Tue July 28 11:04:46 2017

@author: Siavash

EnKF with damping factor
"""
def HEnKF(A,D,R,a):
    import numpy as np
    from numpy.linalg import inv
    NENS = len(A.T)
    NOBS = len(D)
    NP = len(A) - NOBS
    N1 = 1./NENS*np.ones((NENS,NENS))
    B = (np.identity(NENS) - N1)
    DA = np.dot(A,B)
    G = DA[NP:len(DA), :].T
    GT = DA[NP:len(DA), :]
    HPHp = np.dot(GT,G) /(NENS-1)
    rhop=HPHp
    K2 = rhop+R
    U, s, V = np.linalg.svd(K2, full_matrices=True)
    S = np.diag(s)
    K2inv = V*inv(S)*U.T
    PHt = np.dot(DA,G) /(NENS-1)
    K = np.dot(PHt,K2inv)
    At=a*K
    Af=(D-A[NP:len(A), :])
    AUPDATE = A + np.dot(At,Af)
    return AUPDATE
        

   
# Function to calculate the perturbed observations for EnKF
def std_maker(REFOBS):
    import numpy as np
    refobs_rng = max(REFOBS)-min(REFOBS);
    nn = len(REFOBS)
    QCPMAP_temp = (0.1 *refobs_rng)*np.ones((nn,1));
    QCPMAP = np.ones((len(QCPMAP_temp)))
    QCPMAP = QCPMAP_temp
    return QCPMAP

    
# Function to calculate the error and covariance of error for EnKF   
def OBSERVATION2(NENS,R,NOBS,OBSN):
    import numpy as np
    OBS = OBSN
    RSQR = R
    OBST = len(R.T)
    for OBSI in xrange(0,OBST):
        print "iter is", repr(OBSI)
        A = RSQR[:,OBSI:OBSI+1]*np.ones((1,NENS))
        B1 = np.random.normal(0, 1, NOBS*NENS)
        B = np.reshape(B1, (NOBS, NENS))
        EPSILON_temp = A*B
        EPSILON = EPSILON_temp        
        AD = OBS[:,OBSI:OBSI+1]*np.ones((1,NENS))
        D_temp = AD+EPSILON
        D = D_temp
        
    return[D,EPSILON]    