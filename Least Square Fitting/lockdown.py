import numpy as np
from scipy.integrate import odeint

def model_lockdown_fast(time,*p,changeS=None,lag=None,lag2=None):

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
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
            if changeS is not None:
                if lag2 is not None and t>df.shape[0]+lag+lag2:
                    a_t=a_test
                    s = 0
                    d_a_ = d_a
                    a = aa
                else:
                    a_t=a_test
                    s= 0.66
                    d_a_ = d_a
                    if changeS!=0:
                        a = 0
                    else:
                        a = aa
        else:
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
        a_t=a_test
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
        C = X[0]*X[2]*np.array([[  -b/N],
                                [   b/N],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0]])
        C = C + X[0]*X[3]*np.array([[-a_i*b/N],
                                    [ a_i*b/N],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0]])
        return (np.matmul(M,X.reshape(1,-1).T)+C).flatten()

    output = odeint(dX,[N-E0-IS0-QS0-RA0-RS0-D0-RU0-P0-QA0,E0,IS0,IA0,QS0,RA0,D0,RU0,P0,QA0,RS0],time)
    return output

def model_lockdown_slow(time,*p,changeS=None,lag=None,lag2=None):

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
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
            if changeS is not None:
                if lag2 is not None and t>df.shape[0]+lag+lag2:
                    a_t=a_test
                    s = 0
                    d_a_ = d_a
                    a = aa
                else:
                    a_t=a_test
                    s= changeS*aa
                    d_a_ = d_a
                    if changeS!=0:
                        a = 0
                    else:
                        a = aa
        else:
            a_t=a_test
            a = aa
            s = ss
            d_a_ = d_a
        a_t=a_test
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
        C = X[0]*X[2]*np.array([[  -b/N],
                                [   b/N],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0],
                                [     0]])
        C = C + X[0]*X[3]*np.array([[-a_i*b/N],
                                    [ a_i*b/N],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0],
                                    [       0]])
        return (np.matmul(M,X.reshape(1,-1).T)+C).flatten()

    output = odeint(dX,[N-E0-IS0-QS0-RA0-RS0-D0-RU0-P0-QA0,E0,IS0,A0,QS0,RA0,D0,RU0,P0,QA0,RS0],time)
    return output
