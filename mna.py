#! /usr/bin/python
"""Modified Nodal Analysis Demo"""
import re as re
import numpy as np


def charConvert(target):
    if re.match(r"[a-z]", target):
        return ord(target) - 96
    else:
        return target


netlistFile = input("Input Netlist File (Abs)Path:\n")

with open(netlistFile) as foo:
    N = 0
    M = 0
    for line in foo:
        node = re.match(r"^[rvi](\d+|\w)\s(\d+|\w)\s(\d+|\w)\s(\d+\.*\d*)", line)
        volt = re.match(r"^v(\d+|\w)\s(\d+|\w)\s(\d+|\w)\s(\d+\.*\d*)", line)
        if node:
            nodeA = charConvert(node.group(2))
            nodeB = charConvert(node.group(3))
            nodeBuf = max(int(nodeA), int(nodeB))
        N = max(N, nodeBuf)

        if volt:
            M = M + 1
    # End of Counting Nodes(N, excluding GND) & VolSources(M).

with open(netlistFile) as foo:
    D = np.zeros([M, M])
    G = np.zeros([N, N])
    B = np.zeros([N, M])
    vectorI = np.zeros([N, 1])
    vectorE = np.zeros([M, 1])

    for line in foo:
        valid = re.match(r"^r(\d+|\w)\s(\d+|\w)\s(\d+|\w)\s(\d+\.*\d*)", line)
        if valid:  # deal with Resistors
            nodeA = int(charConvert(valid.group(2))) - 1
            nodeB = int(charConvert(valid.group(3))) - 1
            value = 1 / float(valid.group(4))
            if nodeA > -1 and nodeB > -1:
                G[nodeA][nodeA] = G[nodeA][nodeA] + value
                G[nodeB][nodeB] = G[nodeB][nodeB] + value
                G[nodeA][nodeB] = G[nodeA][nodeB] - value
                G[nodeB][nodeA] = G[nodeB][nodeA] - value
            elif nodeA > -1:
                G[nodeA][nodeA] = G[nodeA][nodeA] + value
            elif nodeB > -1:
                G[nodeB][nodeB] = G[nodeB][nodeB] + value

        valid = re.match(r"^v(\d+|\w)\s(\d+|\w)\s(\d+|\w)\s(\d+\.*\d*)", line)
        if valid:  # deal with voltage Sources
            voltIndex = int(charConvert(valid.group(1))) - 1
            posNode = int(charConvert(valid.group(2))) - 1
            negNode = int(charConvert(valid.group(3))) - 1
            voltValue = float(valid.group(4))
            vectorE[voltIndex][0] = voltValue

            if posNode > -1:
                B[posNode][voltIndex] = 1
            if negNode > -1:
                B[negNode][voltIndex] = -1

        valid = re.match(r"^i(\d+|\w)\s(\d+|\w)\s(\d+|\w)\s(\d+\.*\d*)", line)
        if valid:  # deal with current Sources
            posNode = int(charConvert(valid.group(2))) - 1
            negNode = int(charConvert(valid.group(3))) - 1
            currentValue = float(valid.group(4))
            if posNode > -1:
                vectorI[posNode][0] = vectorI[posNode][0] + currentValue
            if negNode > -1:
                vectorI[negNode][0] = vectorI[negNode][0] - currentValue

    C = B.T
    A = np.vstack((np.hstack((G, B)), np.hstack((C, D))))
    vectorZ = np.vstack((vectorI, vectorE))
    vectorResult = np.dot(np.linalg.inv(A), vectorZ)
    for i in range(N):
        print(("Nodal Potential%d is %.2f V" % ((i + 1), vectorResult[i][0])))
    for i in range(M):
        print("Current thru Vs%d is %.2f A" % ((i + 1), vectorResult[N + i][0]))


    #print("A-matrix\n", A, "\n")
    #print("vector_Z:\n", vectorZ, "\n")
    #print("result:\n", res, "\n")
    #print("A inverse\n", np.linalg.inv(A),"\n")
    #print("vector_I:\n", vectorI, "\n")
    #print("vector_E:\n", vectorE, "\n")
    #print("B-matrix\n",B,"\n")
    #print("G-matrix\n",G,"\n")
    #print("D-matrix\n",D)
    ##print("C-matrix\n",C,"\n")
