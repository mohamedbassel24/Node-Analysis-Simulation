import numpy as np
class Component:
    def __init__(self, Type, Node1, Node2, Value, InitialValue):
        self.Type = Type
        self.Node1 = Node1
        self.Node2 = Node2
        self.Value = Value
        self.InitialValue = InitialValue


def ParsingFile(FileName):
    f = open(FileName + ".txt", "r")
    if f.mode == "r":  # check if file is open
        FileContents = f.readlines()  # read file line by line
    for i in range(np.shape(FileContents)[0] - 1):  # loop for the lines before last
        FileContents[i] = FileContents[i][0:len(FileContents[i]) - 1]
    mTimeStamp = FileContents[0]
    Dummy = FileContents[1]
    FileContents = FileContents[2:]
    nComponentList = []
    for i in range(np.shape(FileContents)[0]):
        Line = FileContents[i].split()
        if Line[0] != 1 or Line[0] != -1:
            nComponentList.append(Component(Line[0], Line[1], Line[2], Line[3], Line[4]))
    return nComponentList, mTimeStamp


def initmatg(matrixG, ComponentList):
    for component in ComponentList:
        if component.Type == "R":
            node1 = int(component.Node1[1]) - 1
            node2 = int(component.Node2[1]) - 1
            if node1 < 0 or node2 < 0:  # for Ground Logic
                mMax = max(node1, node2)
                node1 = node2 = mMax
                matrixG[node1][node1] += 1 / float(component.Value)
            else:
                matrixG[node1][node1] += 1 / float(component.Value)
                matrixG[node2][node2] += 1 / float(component.Value)
                matrixG[node1][node2] = -1 / float(component.Value)  # leh msh += heya kman
                matrixG[node2][node1] = -1 / float(component.Value)


def initmatb(matrixb, ComponentList):
    VolCounter = 0
    for component in ComponentList:
        if component.Type == "Vsrc":
            #       for index, node in enumerate(matrixb):
            # print(index)
            pos = int(component.Node1[1]) - 1
            neg = int(component.Node2[1]) - 1
            if pos >= 0:
                matrixb[pos][VolCounter] = 1
            if neg >= 0:
                matrixb[neg][VolCounter] = -1
            VolCounter = VolCounter + 1


def initmatc(B):
    return B.transpose()


def IniMatA(G, B, C, D):
    UpperA = np.hstack((G, B))
    DownA = np.hstack((C, D))
    return np.vstack((UpperA, DownA))


def initmate(matrixe, ComponentList):
    Index = 0
    for component in ComponentList:
        if component.Type == "Vsrc":
            volt = int(component.Value)

            # for index, node in enumerate(matrixe):
            matrixe[Index][0] = volt
            Index = Index + 1


def initmati(matrixi, ComponentList):
    for component in ComponentList:
        if component.Type == "Isrc":
            node1 = int(component.Node1[1]) - 1
            node2 = int(component.Node2[1]) - 1
            if node1 >= 0:
                matrixi[node1][0] = float(component.Value)
            if node2 >= 0:
                matrixi[node2][0] = float(component.Value)


def WriteToFile(FileName, Step, Values, n, m):
    f = open("RCircuit " + FileName + " output.txt", "w+")
    for i in range(n):
        mString = "V" + str(i + 1) + "\n"
        mString += Step + " "
        mString += str(Values[i]) + "\n"
        f.write(mString)
    for i in range(m):
        mString = "I_Vsrc" + str(i + 1) + "\n"
        mString += Step + " "
        mString += str(Values[i + n]) + "\n"
        f.write(mString)
    #  f.write("%f \n",  Values[i + n - 1])
    f.close()


ComponentList = []
TimeStamp = 0
FileNumber = 0
while 1:
    g = input("Enter File# : ")
    if 4 >= int(g) > 0:
        FileNumber = g
        break
    else:
        print("Enter a Valid File#")
ComponentList, TimeStamp = ParsingFile(FileNumber)
n = 0  # representing Number of Nodes
m = 0  # representing Number of ID voltage Source

for mComponent in ComponentList:
    n = max(n, int(mComponent.Node1[1]), int(mComponent.Node2[1]))  # To Get Number of Nodes
    if mComponent.Type == "Vsrc":  # Count The No. of Voltage src
        m = m + 1


#                                   Mat A
# INITIALIZING the Matrices
print(n)
G = np.zeros((n, n))  # for A resistance
B = np.zeros((n, m))  # connection of the voltage sources
C = np.zeros((m, n))  # Transpose of B
D = np.zeros((m, m))  # is a zero matrix
# Calculting The Matrices Values:

initmatg(G, ComponentList)
initmatb(B, ComponentList)
C = initmatc(B)
D = np.zeros((m, m))
print("MatG:\n", G)
print("MatB:\n", B)
print("MatC:\n", C)
print("MatD:\n", D)
A = IniMatA(G, B, C, D)
#                                 Mat X(Unknown)
V = np.zeros((n, 1))  # hold the unknown voltages at each node
J = np.zeros((m, 1))  # holds the unknown currents through the voltage sources.
X = np.vstack((V, J))
#                                  Mat Z
I = np.zeros((n, 1))
E = np.zeros((m, 1))
initmate(E, ComponentList)
initmati(I, ComponentList)
Z = np.vstack((I, E))
print("Z:", Z)
# Solving The AX=Z
print("MatA:", A, "\nMatZ", Z)
X = np.linalg.solve(A, Z)
WriteToFile(FileNumber, TimeStamp, X, n, m)
