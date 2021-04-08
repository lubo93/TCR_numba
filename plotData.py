class PlotData:
    def __init__(self):
        self.xArgsName = ""
        self.clustersizesName = ""
        self.clustersizesAllErrName = ""
        self.diametersName = ""
        self.diametersAllErrName = ""
        self.eccentritiesName = ""
        self.degreeMeanName = ""
        self.degreeMeanAllErrName = ""
        
        self.xArgs = list()
        self.clustersizes = list()
        self.clustersizesAllErr = list()
        self.diameters = list()
        self.diametersAllErr = list()
        self.eccentrities = list()
        self.degreeMean = list()
        self.degreeMeanAllErr = list()
 
class PlotDataSingleNet(PlotData):
    def __init__(self):
        self.cPre = list()
        self.cPost = list()
        self.clustersizesPre = list()
        self.clustersizesPost = list()
        self.GnamesHist1Pre = list()
        self.GnamesHist1Post = list()
        
        self.N = -1
        self.min_ldVal = -1
        self.max_ldVal = -1
        self.extractNum = -1
    
class PlotDataManyNets():
    def __init__(self):
        self.clustersizesMaxPre = list()
        self.clustersizesMaxPost = list()
        self.degMaxsPre = list()
        self.degMaxsPost = list()
        self.Ns = list()
        self.isConnectedsPre = list()
        self.isConnectedsPost = list()
        self.clustersizesMaxPost = list()
       
        self.maxValIndices = []
 
        self.Ns = []
        #info: number of extracted letters from the DNA sequence 
        self.min_ldVal = -1
        self.max_ldVal = -1
        self.extractNum = -1

        # info: samples = 6 is usual
        self.samples = -1

        # info: maximum levenshtein distance
        self.maxValIndices = []

class PlotDataLevenshtein: 
    def __init__(self):
        self.xArgsName = ""
        self.clustersizesName = ""
        self.clustersizesAllErrName = ""
        self.diametersName = ""
        self.diametersAllErrName = ""
        self.eccentritiesName = ""
        self.degreeMeanName = ""
        self.degreeMeanAllErrName = ""

        self.xArgs = list()
        self.clustersizes = list()
        self.clustersizesAllErr = list()
        self.diameters = list()
        self.diametersAllErr = list()
        self.eccentrities = list()
        self.degreeMean = list()
        self.degreeMeanAllErr = list()
  
        self.clustersizesAllArray = list()
        self.diametersAllArray = list()
        self.centersAllArray = list()
        self.eccentritiesAllArray = list()
        self.degreeMeanAllArray = list()
        self.lengthsAllArray = list()

        # info: only needed for lenNet2.py: 
        self.xArgsAllArray = list()

class InnerDataLevenshtein:
     def __init__(self):
         self.clustersizes = list()
         self.diameters = list()
         self.centers = list()
         self.eccentrities = list()
         self.degreeMean = list()
         
         # info: only needded for lenNet2:
         self.clustersizesAllArray = list()
         self.diametersAllArray = list()
         self.centersAllArray = list()
         self.eccentritiesAllArray = list()
         self.degreeMeanAllArray = list()
         self.xArgsAllArray = list()

         self.lengths = list()

class InnerDataManyNets: # I think, we should find away to reduce the number of members here
     def __init__(self):
         self.eccs = list()
         self.dias = list()
         self.degMaxs = list()
         self.nodes = list()
         self.cons = list()
         self.clustersizes = list()
         self.clustersizesMax = list()
         self.isConnected = list()
         
         self.eccsPre = list()
         self.diasPre = list()
         self.degMaxsPre = list()
         self.nodesPre = list()
         self.consPre = list()
         self.clustersizesPre = list()
         self.clustersizesMaxPre = list()
         self.isConnectedsPre = list()

         self.eccsPost = list()
         self.diasPost = list()
         self.degMaxsPost = list()
         self.nodesPost = list()
         self.consPost = list()
         self.clustersizesPost = list()
         self.clustersizesMaxPost = list()
         self.isConnectedsPost = list()

class DataObject:
    def __init__(self):
        self.clustersizes = list()
        self.diameters = list()
        self.centers = list()
        self.eccentricities = list()

        self.clustersizesAllArray = list()
        self.diametersAllArray = list()
        self.centersAllArray = list()
        self.eccentricitiesAllArray = list()
        self.degreeMeanAllArray = list()
        self.xArgsAllArray = list()
