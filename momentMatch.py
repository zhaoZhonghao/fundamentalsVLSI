#! /usr/bin/python
'''for demo of delay calculation on model-reduced ladder(moment matching); the
original ladder code is in demo_ladder.py, which is imported for waveform
comparing. Numpy is used for matrix manipulation, Matplotlib for plloting
and Control for s-domain analysis'''


import numpy as np
import matplotlib.pyplot as plt
import control
import demo_ladder as demo

def mtp(target, a, b):
    return float(np.dot(np.dot(a.T, target), b))


def momentMatching(N=10, qParam=3, endTime=300, deltaTime=0.1):
    '''
    Model-reduced ladder waveform simulator(moment matching):

    '''
    ## Everything is the same as in the laderWave() for matrix initialization
    c_val = 1
    g_val = 1

    C = np.zeros([N, N])
    G = np.zeros([N, N])

    diagC = [c_val] * N
    diagG = [-2 * g_val] * N
    sub_diagG = [g_val] * (N - 1)

    C = np.diag(diagC)
    G = np.diag(diagG)
    G = G + np.diag(sub_diagG, -1) + np.diag(sub_diagG, 1)
    ## Last node has only one resistor connected, so change -2 to -1 on diagnal
    G[N - 1][N - 1] += g_val

    ## Input step stimulus. The input node is kept as 1V
    ## Per Norton's Theorem, the equivalent input current source is 1 amp
    i_source = np.zeros([N, 1])
    i_source[0] = 1

    ## state-space MIT slides' version, NOT common one!!
    ## interpret the EXACT system in form of state space
    ## which is elaborated below:
    ## Conventionally and initially,
    ##     C * [dv/dt] = G * [v] + inv(C) * [i_source]
    ##     Let
    ##        Amat = inv(C) * G
    ##        bmat = inv(C) * [i_source]
    ##        cmat = [0, 0, ..., 0, 1] (of the same dimension)
    ##     The last ladder node Voltage can be rewritten as
    ##        [dv/dt] = A * [v] + bmat
    ##        y = cmat * [v]  ## NOTE y is the ! scalar ! output
    Amat = np.dot(np.linalg.inv(C), G)
    bmat = np.dot(np.linalg.inv(C), i_source)
    cmat = np.zeros([N, 1])
    cmat[N - 1] = 1

    ## ->-> MOMENT MATCHING Begins <-<- ##
    midinvpwrA = np.linalg.inv(Amat)  # middle inversed powered Amatrix
    invA = midinvpwrA

    ## moment matching is based the assumption that Hr(s) = H(s)
        ## expand Hr(s)=H(s)
        ## with H(s) defined as:
    ##     H(s) = m0 + m1*s + m2*s^2 + ... + mk*s^ k + (inf)
    ##     where mk is called moment and decided as cmat * Amat^(-k) * bmat
    ## with Hr(s) defined as:
    ##     Hr(s) = (b0 + b1*s + ... + b(q-1)*s^(q-1))/(1 + a1*s + ... + aq*s^q)
    ## CROSS MULTiPLYING and match terms of same order:
        ## SECTION I : amMAT * [a] = amRVCT

        ## amMAT is the leftside Matrix, that relates [a] and [m] # thus amMAT
        ## amRVCT is the *R*ightside *vector* # thuse amRVCT
        ## both amMAT and amRVCT is calculable with the moment definition
        ## based on above definitions, we can get coefficient a's
        ## with [a] = inv(amMAT) * amRVCT
##--------
    ##     m0        m1    ...    m(q-1)    aq        mq
    ##     m1        m2    ...    mq        a(q-1)    m(q+1)
    ##    [    ....... .....             ]* [ .. ] = [ ..   ]
    ##     m(q-1)    mq ...       m(2q-2)   a1       m(2q-1)
##--------
        ## NOTICE ascending or descending
        ## NOTICE that amMAT all m(i) are in sub-diags(rightUP -> leftDOWN)
        ## that we first generate a MIRRORED amMAT(Alias) i.e.,
##--------
    ##     m(q-1)        m(q-2)    ...    m0
    ##     mq            m(q-1)    ...    m1
    ##    [    ....... .....                 ]
    ##     m(2q-2)       m(2q-3) ...      m(q-1)
##--------
        ## this way, we can use numpy.diag easily
        ## Once MIRRORED matrix's generation is achived, horizontally flip it
    ## coefficient a's are attainable under above guidance and instruction

        ## SECTION II : bmMAT * bmRVCT = [b]
        ## bmMAT is the matrix that relates [b] and [m]
        ## bmLVCT is the *L*eft *vector* in the equation
    ## bmMAT, bmLVCT and [b] shares a similar relationship
##--------
    ##     1       0    ...0     ..    0         m0        b
    ##     a1      1    ...0     ..    0         m1        b(q+1)
    ##     a2      a1      1           0
    ##    [    ....... ...a1           0 ] *   [ ..   ] = [ ..   ]
    ##     a(q-1)  a(q-2) ...    a1    1         m(q-1)    b(2q-1)
##--------
        ## this way, we can generate bmMAT, bmRVCT after [a] is obtained
        ## Once bmMAT orginal matrix's generation is achived,
        ## [b] = bmMAT * [m]
    ## coefficient b's are attainable under above guidance and instruction

    amMAT = np.diag([0.0] * qParam)
    amRVCT = np.zeros([qParam, 1])
    bmMAT = np.diag([1.0] * qParam)
    bmLVCT = np.zeros([qParam, 1])

    for i in range(qParam * 2):
        loadCell = mtp(midinvpwrA, cmat, bmat)
        midinvpwrA = np.dot(midinvpwrA, invA)
        if i < 2 * qParam - 1:
            celldiagPos = qParam - i - 1
            celldiagCtt = [loadCell] * (qParam - abs(qParam - i - 1))
            amMAT += np.diag(celldiagCtt, celldiagPos)
        if i >= qParam:
            amRVCT[i - qParam][0] = -loadCell  # notice the negative logo!!
        if i < qParam:
            bmLVCT[i][0] = loadCell

    amMAT = np.fliplr(amMAT)
    amRes = np.dot(np.linalg.inv(amMAT), amRVCT)
    den = []
    for each in amRes:
        den.append(float(each))
    den.append(1)
    # print(den)  # denominator a_q, a_q-1, ..., a_2, a_1, 1
    amRes = (np.fliplr(amRes.T)).T  # ascending order -> for transferFunction

    ## bmMAT SECTION
    bmMAT = np.zeros([qParam, qParam])
    loadCell = 1
    for i in range(qParam):
        celldiagPos = -i
        celldiagCtt = [loadCell] * (qParam - i)
        bmMAT += np.diag(celldiagCtt, celldiagPos)
        loadCell = float(amRes[i][0])  # notice int or float choosela   ()
    bmRes = np.dot(bmMAT, bmLVCT)
    nom = []
    for each in bmRes:
        nom.append(float(each) * -1.0)
    # print(nom)
    nom.reverse()
    
    ## with a's and b's coefficient vectors, i.e. [amRes] and [bmRes]
    ## stored in their coresponding list nom and den,
    ## elements in ascending form -->
        ## nom = -[b(q-1), b(q-2), ..., b2, b1, b0]
            ## nom = -[b(q-1)*s^(q-1), b(q-2)*s^(q-2), ..., b2*s^2, b1*s, b0]
        ## den =  [aq, a(q-1), ..., a2, a1, 1]
            ## den =  [aq*s^q, a(q-1)*s^(q-1), ..., a2*s^2, a1*s, 1]
    G = control.TransferFunction(nom, den) # transfer function
    U = control.TransferFunction([1], [1, 0]) # stimulus
    Y = G * U
    t = np.arange(0, endTime, deltaTime)
    t, ymm = control.impulse_response(Y, t) # inverse Laplace
    return (ymm, t)


def plot_wave(x, y, color='b', name="", pen=3):
    '''plot the waveform'''
    plt.plot(x, y, color, label=name, linewidth=pen)


if __name__ == '__main__':
    ## If this program is not 'imported' by other program,
    ## call momentMatching() some times to get voltage waves and time intervals
    ## also get corresponding exact curve of the same N(order)
    ## Here, for clarity and full covergency (set time=20000)
        ##       | momentMatching   |   ladderWave()
        ##__________________________________________
        ## N=10  |        Y         |       Y
        ## N=100 |        Y         |       Y
    (v, t) = momentMatching(endTime=20000)
    plot_wave(t, v, 'y:', name="N=10")
    (v, t) = demo.ladderWave(endTime=20000)
    plot_wave(t, v, name="N=10(ori)", pen=1)
    
    (v, t) = momentMatching(N=100, endTime=20000)
    plot_wave(t, v, 'c:', "N=100")
    (v, t) = demo.ladderWave(N=100, endTime=20000)
    plot_wave(t, v, name="N=100(ori)", pen=1)

    plt.title("moment match(dashline) | original(solidline)", fontsize=14)
    plt.ylabel("Output(V)")
    plt.xlabel("Time: in your deltaTime unit")
    plt.legend()
    plt.show()
