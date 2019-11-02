import numpy as np


class Component:
    def __init__(self, Type, Node1, Node2, Value, InitialValue, NumVoltageSource=0):
        self.Type = Type
        self.Node1 = Node1
        self.Node2 = Node2
        self.Value = Value
        self.InitialValue = InitialValue
        self.mVoltageSrcNumber = NumVoltageSource


def ParsingFile(FileName):
    f = open(FileName + ".txt", "r")
    if f.mode == "r":  # check if file is open
        FileContents = f.readlines()  # read file line by line
    for i in range(np.shape(FileContents)[0] - 1):  # loop for the lines before last
        FileContents[i] = FileContents[i][0:len(FileContents[i]) - 1]
    mTimeStamp = FileContents[0]
    NumberOfIterations = FileContents[1]
    FileContents = FileContents[2:]
    nComponentList = []
    for i in range(np.shape(FileContents)[0]):
        Line = FileContents[i].split()
        if Line[0] != 1 or Line[0] != -1:
            nComponentList.append(Component(Line[0], Line[1], Line[2], Line[3], Line[4]))
    return nComponentList, mTimeStamp, NumberOfIterations


def initmatg(matrixG, ComponentList):
    for component in ComponentList:
        if component.Type == "R":
            node1 = int(component.Node1[1]) - 1
            node2 = int(component.Node2[1]) - 1
            if node1 < 0 or node2 < 0:  # for Ground Logic
                mMax = max(node1, node2)
                node1 = mMax
                matrixG[node1][node1] += 1 / float(component.Value)
            else:
                matrixG[node1][node1] += 1 / float(component.Value)
                matrixG[node2][node2] += 1 / float(component.Value)
                matrixG[node1][node2] += -1 / float(component.Value)
                matrixG[node2][node1] += -1 / float(component.Value)


def initmatb(matrixb, ComponentList):
    VolCounter = 0
    for component in ComponentList:
        if component.Type == "Vsrc" or component.Type == "I":
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


def initmate(matrixe, ComponentList, rStep):
    Index = 0
    for component in ComponentList:
        if component.Type == "Vsrc":
            volt = float(component.Value)
            matrixe[Index][0] = volt
            Index = Index + 1
        elif component.Type == "I":
            volt = float(component.Value)
            matrixe[Index][0] = (float(component.Value) / rStep) * (float(component.InitialValue))
            Index = Index + 1


def initmati(matrixi, ComponentList, rTimeStamp):
    for component in ComponentList:
        if component.Type == "Isrc":
            node1 = int(component.Node1[1]) - 1
            node2 = int(component.Node2[1]) - 1
            if node1 >= 0:
                matrixi[node1][0] += float(component.Value)
            if node2 >= 0:
                matrixi[node2][0] += -1 * float(component.Value)
        elif component.Type == "C":
            node1 = int(component.Node1[1]) - 1
            node2 = int(component.Node2[1]) - 1
            if node1 >= 0:
                matrixi[node1][0] += (float(component.Value) / rTimeStamp) * (float(component.InitialValue))
            if node2 >= 0:
                matrixi[node2][0] += -1 * ((float(component.Value) / rTimeStamp) * (float(component.InitialValue)))


def WriteToFile(FileName, Step, List, n, m):
    f = open("RCircuit " + FileName + " output.txt", "w+")
    for i in range(n):
        mString = "V" + str(i + 1) + "\n"
        for j in range(len(List)):
            mString += str((float(Step) * (j + 1)))[0:3] + " "
            mString += str(List[j][i]) + "\n"
            f.write(mString)
            mString = ""
    for i in range(m):
        mString = "I_Vsrc" + str(i + 1) + "\n"
        for j in range(len(List)):
            mString += str(float(Step) * (j + 1))[0:3] + " "
            mString += str(List[j][i + n]) + "\n"
            f.write(mString)
            mString = ""

    f.close()


def ConvertCap_Res(mList, rStep):
    """This Function Convert each capacitor into Resistance"""
    for mComp in mList:
        if mComp.Type == "C":
            mList.append(Component("R", mComp.Node1, mComp.Node2, str(float(rStep) / float(mComp.Value)), "0"))

    return mList


def ConvertInd_Res(mList, rStep):
    """This Function Convert each capacitor into Resistance"""
    for mComp in mList:
        if mComp.Type == "I":  # NOT G
            mList.append(Component("R", mComp.Node1, mComp.Node2, str(float(mComp.Value) / float(rStep)), "0"))

    return mList


def UpdateCInitalValue(mList, X):
    """This Function Convert each capacitor into Resistance"""
    for i in range(len(mList)):
        if mList[i].Type == "C":
            mLeft = int(mList[i].Node1[1]) - 1
            mRight = int(mList[i].Node2[1]) - 1
            pos = 0
            Neg = 0
            if mLeft < 0:
                pos = 0
            else:
                pos = X[int(mList[i].Node1[1]) - 1]
            if mRight < 0:
                Neg = 0
            else:
                Neg = X[int(mList[i].Node2[1]) - 1]
            mList[i].InitialValue = pos - Neg

    return


def UpdateIInitalValue(mList, X):
    """This Function Updates the Value of Inductor """
    for i in range(len(mList)):
        if mList[i].Type == "I":
            mList[i].InitialValue = float(X[int(mList[i].mVoltageSrcNumber)])

    return


ComponentList = []
TimeStamp = 0
FileNumber = 0
NumberOfIterations = 0
while 1:
    g = input("Enter File# : ")
    if 8 >= int(g) > 0:
        FileNumber = g
        break
    else:
        print("Enter a Valid File#")

ComponentList, TimeStamp, NumberOfIterations = ParsingFile(FileNumber)
n = 0  # representing Number of Nodes
m = 0  # representing Number of ID voltage Source
m = 0  # representing Number of ID voltage Source

for mComponent in ComponentList:
    n = max(n, int(mComponent.Node1[1]), int(mComponent.Node2[1]))  # To Get Number of Nodes
    if mComponent.Type == "Vsrc":  # Count The No. of Voltage src
        m = m + 1
    elif mComponent.Type == "I":
        m = m + 1
        mComponent.mVoltageSrcNumber = m - 1

#                                   Mat A
# INITIALIZING the Matrices
ComponentList = ConvertCap_Res(ComponentList, TimeStamp)
ComponentList = ConvertInd_Res(ComponentList, TimeStamp)
G = np.zeros((n, n))  # for A resistance
B = np.zeros((n, m))  # connection of the voltage sources
C = np.zeros((m, n))  # Transpose of B
D = np.zeros((m, m))  # is a zero matrix
# Calculting The Matrices Values:
D[0, 0] = -1 / float(TimeStamp)
D[1, 1] = -1 / float(TimeStamp)
initmatg(G, ComponentList)
initmatb(B, ComponentList)
C = initmatc(B)
D = np.zeros((m, m))
# print("MatG:\n", G)
# print("MatB:\n", B)
# print("MatC:\n", C)
# print("MatD:\n", D)
A = IniMatA(G, B, C, D)
#                                 Mat X(Unknown)
V = np.zeros((n, 1))  # hold the unknown voltages at each node
J = np.zeros((m, 1))  # holds the unknown currents through the voltage sources.
X = np.vstack((V, J))
#                                  Mat Z

mList = []
TimeStamp = float(TimeStamp)
NumberOfIterations = int(NumberOfIterations)
ActualTime = TimeStamp

while ActualTime <= (NumberOfIterations * TimeStamp):
    I = np.zeros((n, 1))
    E = np.zeros((m, 1))
    initmate(E, ComponentList, TimeStamp)
    initmati(I, ComponentList, TimeStamp)
    Z = np.vstack((I, E))
    X = np.linalg.solve(A, Z)
    mList.append(X)
    UpdateCInitalValue(ComponentList, X)
    UpdateIInitalValue(ComponentList, X[n:])
    ActualTime += TimeStamp

WriteToFile(FileNumber, TimeStamp, mList, n, m)
