import numpy as np
import math


def numbers_to_strings(argument):
    switcher = {
        0: "zero",
        1: "one",
        2: "two",
    }


class Component:
    def __init__(self, Type, Node1, Node2, Value, InitialValue):
        self.Type = Type
        self.Node1 = Node1
        self.Node2 = Node2
        self.Value = Value
        self.InitialValue = InitialValue


def ParsingFile(FileName):
    f = open(FileName, "r")
    if f.mode == "r":  # check if file is open
        FileContents = f.readlines()  # read file line by line
    # Removing \n at each line
   # print(FileContents)
   # print(np.shape(FileContents))
    for i in range(np.shape(FileContents)[0] - 1):  # loop for the lines before last
       # print(len(FileContents[i]))
       # print(FileContents[i])
        FileContents[i] = FileContents[i][0:len(FileContents[i]) - 1]
    print(np.shape(FileContents))
    print(FileContents)
    TimeStamp = FileContents[0]
    Dummy = FileContents[1]
    FileContents = FileContents[2:]
    ComponentList = []
    for i in range(np.shape(FileContents)[0]):
        Line = FileContents[i].split()
        ComponentList.append(Component(Line[0], Line[1], Line[2], Line[3], Line[4]))
    return (ComponentList, TimeStamp)


ComponentList = []
TimeStamp = 0
ComponentList, TimeStamp = ParsingFile("1.txt")
n = 0  # representing Number of Nodes
m = 0  # representing Number of ID voltage Source
for mComponent in ComponentList:
    n = max(int(mComponent.Node1[1]), int(mComponent.Node2[1]))+1  # To Get Number of Nodes
    if mComponent.Type == "Vsrc":
        m = m + 1
# INITIALIZING the Matrices
G = np.zeros((n, n))  # for A resistance
B = np.zeros((n, m))  # connection of the voltage sources
C = np.zeros((m, n))  # Transpose of B
D = np.zeros((m, m))  # is a zero matrix
I = np.zeros((n,1))
E = np.zeros(((m,1)))
A = np.zeros((m+n,m+n))
Z = np.zeros((n+m,1))

# print(G)
# print(B)
# print(C)
# print(D)

# Ctrl/ to comment multiple line selection

def initmatb(matrixb):
    #matrixb=np.array(matrixb)
    #print(type(matrixb))
    # Component_Type | Node1 | Node2 | Value | Initial_Value
    # print(ComponentList)
    for component in ComponentList:
        if component.Type == "Vsrc":
            for index, node in enumerate(matrixb):
                #print(index)
                pos=int(component.Node1[1])
                neg=int(component.Node2[1])
                if index == pos:
                    matrixb[index][0] = 1
                elif index == neg:
                    matrixb[index][0] = -1
                else:
                    matrixb[index][0] = 0


# initmatb(B)
# print(B)

def initmatg(matrixG):
    for component in ComponentList:
        if(component.Type == "R"):
            node1 = int(component.Node1[1])
            node2 = int(component.Node2[1])
            matrixG[node1][node1]+=1/float(component.Value)
            matrixG[node2][node2]+=1/float(component.Value)
            matrixG[node1][node2]=-1/float(component.Value)
            matrixG[node2][node1]=-1/float(component.Value)


# initmatg(G)
# print(G)

def initmatc():
    return  B.transpose()

def initmate(matrixe):
    for component in ComponentList:
        if(component.Type == "Vsrc"):
            volt=int(component.Value)
            for index,node in enumerate(matrixe):
                matrixe[index][0]=volt

def initmati(matrixi):
    for component in ComponentList:
        if component.Type == "Isrc":
            node1 = int(component.Node1[1])
            node2 = int(component.Node2[1])
            if(node1 !=0):
                matrixi[node1][0]=float(component.Value)
            if(node2 !=0):
                matrixi[node2][0]=float(component.Value)

def initmata(matrixa):
        for i in range(n):
           for j in range(n):
               matrixa[i][j]=G[i][j]

        for i in range(n):
            for j in range(m):
                matrixa[i][j+n]=B[i][j]

        for i in range(m):
            for j in range(n):
                matrixa[i+n][j]=C[i][j]

        for i in range(m):
            for j in range(m):
                matrixa[i+n][j+n]=D[i][j]
def initmatz(matrixz):
    for i in range(n+m):
        iteri=0
        itere=0
        if i<n:
            matrixz[i][0]=I[iteri][0]
            iteri+=1
        else:
            matrixz[i][0]=E[itere][0]
            itere+=1

def calculatex():
    a=np.linalg.inv(A)
    return np.matmul(a,Z)

initmatg(G)
initmatb(B)
C =initmatc()
initmate(E)
initmati(I)
initmata(A)
initmatz(Z)
X=calculatex()
# print(G)
# print(B)
# print(C)
# print(D)
# print(E)
# print(I)
# print(A)
#print(Z)
print(X)