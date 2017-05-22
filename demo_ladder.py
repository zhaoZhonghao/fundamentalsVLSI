#!/usr/bin/env python
'''for R/C ladder delay calculation demo, using Backward Euler method.
NumPy is used for matrix operation, MatplotLib for plotting'''

import numpy as np
import matplotlib.pyplot as plt

def elmoreDelay(n):
    '''Elmore Delay calculator for uniform R/C ladder
    ---
    + Feedback ladder delay, given n as the ladder order;
    0.69 is a constant factor, which is emphasized in lecture slides
    '''
    d = n * (n + 1) / 2
    return d * 0.69


def ladderWave(N=10, endTime=1000, deltaTime=0.01, tpdMode=False):
    '''
    Ladder waveform simulator: returns transient waveform on ladder end
    ---
    + N -> order of uniform RC ladder, default is 10-order;
    + endTime determines how long the transient simulation takes;
    + deltaTime is the Euler step-length;
    + when tpdMode is True, simulator stops whenever above 0.5*Vin;
    Returns a (t, v) data pair <- time points and waveform data points
    '''
    c_val = 1
    g_val = 1

    C = np.zeros([N, N])
    G = np.zeros([N, N])

    diagC = [c_val] * N
    diagG = [-2 * g_val] * N
    sub_diagG = [g_val] * (N-1)
    print(diagC)
    print(diagG)
    print(sub_diagG, "\n")

    C = np.diag(diagC)
    G = np.diag(diagG)
    G = G + np.diag(sub_diagG, -1) + np.diag(sub_diagG, 1)
    # Last node has only 1 resistor connected, so change -2 to -1 on diagonal
    G[N - 1][N - 1] += g_val
    print(C)
    print(G, "\n")

    ## Input step stimulus. The input node is kept as 1V.
    ## Per Norton's Theorem, the equivalent input current source is 1A.
    i_source = np.zeros([N, 1])
    i_source[0] = 1

    ## Discretize time; 'time' is a list for discrete intervals
    start = 0
    time = [x * deltaTime for x in range(start, int(endTime / deltaTime))]

    ## Begin simulation.
    ## v_nodes is a list of voltages on all nodes, for just one moment
    ## vout is a list for voltage on the last stage, for all time points
    v_nodes = np.zeros([N, 1])
    vout = []
    i = 0
    ## Calculate all node voltages in time t
    ## Check matrices C and G
    ## per KCL, Iout=Iin, [C]*[dv/dt] = [G]*[v] + i_source
    ## Forward Euler => [C]*[(vt1-vt0)/deltaT] = [G]*[vt0] + i_source
    ## Backward Euler => [C]*[(vt1-vt0)/deltaT] = [G]*[vt1] + i_source
    ##     vt1-vt0 = deltaT * C-1 * G * vt1 + deltaT * C-1 * i_source
    ##     (I - deltaT * C-1 * G) * vt1 = vt0 + deltaT * C-1 * i_source
    ## where C-1 is a simple notation of inv(C)
    ## Now solve vt1 = ?
    ## Rewrite the equation into:
    ## A * [vt1] = [vt0] + B => [vt1] = A-1 * ([vt0] + B)
    Ainv = np.linalg.inv(np.eye(N) - deltaTime * np.dot(np.linalg.inv(C), G))

    ## If i_source is time-variant, B should be changed in each time point;
    ## in this demo, i_source is a constant
    B = deltaTime * np.dot(np.linalg.inv(C), i_source)

    ## Give an anlysis on the above formula, and try Forward Euler?

    for tt in time:
        i += 1
        v_nodes = np.dot(Ainv, v_nodes + B)

        ## Add the last stage voltage to vout, for this moment
        vout.append(v_nodes[N-1][0])

        if tpdMode:
            if v_nodes[N-1] > 0.5:
                print("As time=%.2f, Node%d reaches 0.5Vin (%.7f) V\n"
                      % (tt, N, v_nodes[N-1]))
                print("Elmore Delay= %.2f, of identical order N=%d"
                      % (elmoreDelay(N), N))
                break

    return (vout, time[0:i])


def plot_wave(x, y, color='b', name=""):
    '''plot the waveform'''
    plt.plot(x, y, color, label=name, linewidth=2)


if __name__ == '__main__':
    ## If this program is not 'imported' by other program,
    ## call ladderWave() some times to get voltage waves and time intervals
    (v, t) = ladderWave()
    plot_wave(t, v, name="N=10")

    (v, t) = ladderWave(N=5)
    plot_wave(t, v, 'c', "N=5")

    (v, t) = ladderWave(N=20)
    plot_wave(t, v, 'g', "N=20")

    (v, t) = ladderWave(N=50)
    plot_wave(t, v, 'y', "N=50")

    ## Try other commands for tpdMode?
    #(v, t) = ladderWave(tpdMode=True, N=1, deltaTime=0.01)
    #(v, t) = ladderWave(tpdMode=True, N=2, deltaTime=0.01)
    #(v, t) = ladderWave(deltaTime=0.01, tpdMode=True, N=15)
    #(v, t) = ladderWave(N=100, endTime=5000, tpdMode=True)

    ## Or try other commands for showing longer ladder?
    #(v, t) = ladderWave(N=30, endTime=5000)
    #plot_wave(t, v, 'y')
    #(v, t) = ladderWave(N=100, endTime=5000)
    #plot_wave(t, v, 'r')

    plt.title("R/C Ladder Waveform", fontsize=14)
    plt.ylabel("Output(V)")
    plt.xlabel("Time: in your deltaTime unit")
    plt.legend()
    plt.show()
