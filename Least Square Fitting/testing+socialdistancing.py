import numpy as np
from scipy.integrate import odeint

def model_testing_distancing(time,*p,changeATE=0.1,changeS=0,changeDA=1/8,changeAI=0,lag=None):

    E0 = 10**p[2]
    IS0 = 10**p[3]
    IA0 = r*IS0
    QS0 = df.active[0]
    RS0 = df.recovered[0]
    D0 = df.dead[0]
    RU0 = 0
    P0 = 0
    QA0 = 0
    RA0 = 0
    b = p[0]
    d = p[1]
    ss = 0
    aa = p[4]
    l = p[5]
    k = p[6]

    def dX(X,t):
        if (lag is None and t>df.shape[0]) or (lag is not None and t>df.shape[0]+lag):
            a_t=changeATE/(1-changeATE)
            a = 0
            s = 0.66
            d_a_ = changeDA
            b_ = changeAI*b
        else:
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
            b_ = b
        M = np.array([[ -a,       0,     0,              0,     0,     0,     0,     0,     s,     0,     0],
                      [  0,-(1+r)*g,     0,              0,     0,     0,     0,     0,     0,     0,     0],
                      [  0,       g,    -d,              0,     0,     0,     0,     0,     0,     0,     0],
                      [  0,     r*g,     0, -d_a_-d_a_*a_t,     0,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     d,              0,  -l-k,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,              0,     l,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,              0,     k,     0,     0,     0,     0,     0,     0],
                      [  0,       0,     0,           d_a_,     0,     0,     0,     0,     0,     0,     0],
                      [  a,       0,     0,              0,     0,     0,     0,     0,    -s,     0,     0],
                      [  0,       0,     0,       d_a_*a_t,     0,     0,     0,     0,     0,    -l,     0],
                      [  0,       0,     0,              0,     0,     0,     0,     0,     0,     l,     0]])
        C = X[0]*X[2]*np.array([[  -b_/N],
                                [   b_/N],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0]])
        C = C + X[0]*X[3]*np.array([[-a_i*b_/N],
                                    [ a_i*b_/N],
                                    [        0],
                                    [        0],
                                    [        0],
                                    [        0],
                                    [        0],
                                    [        0],
                                    [        0],
                                    [        0],
                                    [        0]])
        return (np.matmul(M,X.reshape(1,-1).T)+C).flatten()

    output = odeint(dX,[N-E0-IS0-QS0-RA0-RS0-D0-RU0-P0-QA0,E0,IS0,IA0,QS0,RA0,D0,RU0,P0,QA0,RS0],time)
    return output
