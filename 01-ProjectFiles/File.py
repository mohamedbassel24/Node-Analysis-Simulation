import numpy as np
import math
def numbers_to_strings(argument):
    switcher = {
        0: "zero",
        1: "one",
        2: "two",
    }
class Component:
  def __init__(self,Type, Node1,Node2,Value,InitialValue):
    self.Type=Type
    self. Node1=Node1
    self.Node2=Node2
    self.Value=Value
    self.InitialValue=InitialValue


def ParsingFile(FileName):
    f = open(FileName, "r")
    if f.mode == "r":# check if file is open
        FileContents = f.readlines() # read file line by line
    #Removing \n at each line
    for i in range (np.shape(FileContents)[0]-1): # loop for the lines before last
        FileContents[i]=FileContents[i][0:len(FileContents[i])-1]
    print(np.shape(FileContents))
    print(FileContents)
    TimeStamp=FileContents[0]
    Dummy=FileContents[1]
    FileContents=FileContents[2:]
    ComponentList=[]
    for i in range (np.shape(FileContents)[0]):
        Line=FileContents[i].split()
        ComponentList.append(Component(Line[0],Line[1],Line[2],Line[3],Line[4]))
    return(ComponentList,TimeStamp)

ComponentList=[]
TimeStamp=0
ComponentList,TimeStamp=ParsingFile("1.txt")
