import math
import random
import numpy as np
from numpy.linalg import inv
import heapq
#from scipy.spatial import Voronoi, voronoi_plot_2d
#import matplotlib.pyplot as plt
from sys import exit
import pdb

class Macro:
    def __init__(self):
        self.ActualNum = 0
        self.density = 0
        self.Radius = 0
        self.Position_x = []
        self.Position_y = []
        self.numAntenna = 0
        self.power_dbm = 0
        self.power_W = 0

        self.UEPosition_x = []
        self.UEPosition_y = []

class Pico:
    def __init__(self):
        self.ActualNum = 0
        self.density = 0
        self.Position_x = []
        self.Position_y = []
        self.numAntenna = 0
        self.power_dbm = 0
        self.power_W = 0


class Channel:
    def __init__(self):
        self.Dis_User_MBSs = []
        self.Dis_User_PBSs = []
        self.User_BSs = []
        self.User_MBSs = []
        self.User_PBSs = []
        self.User_BSs_normalize = []
        self.G = []
        self.G_Hermitian = []
        self.G_psuedo = []
        self.identity_matrix = []
        self.ZF_beam_weight = []
        self.ZF_beam_weight_normalize = []

class Para:
    def __init__(self):
        self.SimulationRegion = 0
        self.SINRthreshold = 0

        self.PathLoss_exponent = 0
        self.PassLoss_reference_d = 0

        self.SNRGap_db = 0
        self.SNRGap_W = 0

        self.LoadingFactor = 0

        self.TotalCapacity = 0
        self.TotalInterference = 0
        self.TotalSignal = 0
        self.Capacity_PerBS = 0
        self.successnum = []

        self.selectedMBSnum = 0

class USER:
    def __init__(self):
        self.ActualNum = 0
        self.Position_x = []
        self.Position_y = []
        self.Alluser_selectedBS_index = []

class Cluster_father:
    def __init__(self):
        self.ue = UE()
        self.mbs = MBS()
        self.pbs = PBS()
        self.channel = Channel()

class UE:
    def __init__(self):
        self.Num = 0
        self.index = []
        self.ActualNum = 0
        self.Position_x = []
        self.Position_y = []
        self.received_signal_power = []
        self.received_interference_power = []
        self.received_SINR_power = []
        self.received_capacity_power = []
        self.SINR_W = []

        self.ActualPosition_x = []
        self.ActualPosition_y = []

class MBS:
    def __init__(self):
        self.index = []
        self.Position_inCluster_x = []
        self.Position_inCluster_y = []

class PBS:
    def __init__(self):
        self.Num = 0
        self.Position_inCluster_x = []
        self.Position_inCluster_y = []

class Channel:
    def __init__(self):
        self.Dis_User_MBSs = []
        self.Dis_User_PBSs = []
        self.User_BSs = []
        self.User_MBSs = []
        self.User_PBSs = []
        self.User_BSs_normalize = []
        self.G = []
        self.G_Hermitian = []
        self.G_psuedo = []
        self.identity_matrix = []
        self.ZF_beam_weight = []
        self.ZF_beam_weight_normalize = []

def GenerateMBSPosition(para, macro):
    Position_x = []
    Position_y = []

    macro.ActualNum = 3 #np.random.poisson(para.SimulationRegion * para.SimulationRegion * math.pi* macro.density)  # Total number of MBS
    #while macro.ActualNum < para.selectedMBSnum:
    #    macro.ActualNum = np.random.poisson(para.SimulationRegion * para.SimulationRegion * math.pi* macro.density)  

    Position_x, Position_y = GeneratePosition(macro.ActualNum, para.SimulationRegion)
    print("Number of BS=", macro.ActualNum)

    dis_MBS_center = []
    for i in range(macro.ActualNum):
        #dis_MBS_center.append(((Position_x[i])**2 + (Position_y[i])**2)**(1/2)) #--> The result is always 1.0
        dis_MBS_center.append( math.sqrt( (Position_x[i])**2 + (Position_y[i])**2 ) )

    for i in sorted(dis_MBS_center):
        macro.Position_x.append(Position_x[dis_MBS_center.index(i)])
        macro.Position_y.append(Position_y[dis_MBS_center.index(i)])

    return macro

# def GenerateUSERPosition(para, macro, user):
#     user.Position_x = []
#     user.Position_y = []
#
#     user.ActualNum = int(macro.ActualNum * para.LoadingFactor * macro.numAntenna)
#     user.Position_x, user.Position_y = GeneratePosition(user.ActualNum, para.SimulationRegion)

#    return user

def GenerateUserPosition(para, macro, user):
    # Divide simulation region in grids, each with size 50x50
    # Place one UE at the intersection point 
    # __________
    # |__|__|__|
    # |__|__|__|
    # |  |  |  |
    #

    Dis_BS_ClusterCenter = np.zeros(macro.ActualNum)
    user_index_x = np.arange((-para.SimulationRegion/2), (para.SimulationRegion/2 + 1), 50)
    user_index_y = np.arange((-para.SimulationRegion/2), (para.SimulationRegion/2 + 1), 50)

    # Create an empty array for BS to store the attached UEs
    macro.UEPosition_x = []
    macro.UEPosition_y = []
    for j in range(macro.ActualNum):
        macro.UEPosition_x.append([])
        macro.UEPosition_y.append([])


    for k in range(len(user_index_y)):
        for i in range(len(user_index_x)):
            for j in range(macro.ActualNum):
                # For MBS j, determine its distance to every UE in row k
                Dis_BS_ClusterCenter[j] = math.sqrt(math.pow((user_index_x[i]) - (macro.Position_x[j]), 2) + math.pow((user_index_y[k]) - (macro.Position_y[j]), 2))            
            minindex = Dis_BS_ClusterCenter.argmin() # Find the nearest UE in row k             
            macro.UEPosition_x[minindex].append(user_index_x[i])
            macro.UEPosition_y[minindex].append(user_index_y[k])
 
    usernum = int(macro.numAntenna * para.LoadingFactor) # Number of users to serve per BS
        
    user.Position_x = []
    user.Position_y = []
    for j in range(macro.ActualNum):
        selecteduser = []  
        #print("len(macro.UEPosition_x[j])=", len(macro.UEPosition_x[j]))      
        # For each BS, randomly choose (usernum) UE from the attached ones to serve
        if len(macro.UEPosition_x[j]) >= usernum:
            selecteduser = random.sample(range(len(macro.UEPosition_x[j])), usernum)
        else:
            selecteduser = random.sample(range(1000), usernum)
            
        
        # Determine the locations of served UEs
        for a in selecteduser:
            user.Position_x.append(macro.UEPosition_x[j][a])
            user.Position_y.append(macro.UEPosition_y[j][a])

    user.ActualNum = len(user.Position_x) # Total number of served users
    print("Number of users=", user.ActualNum)

    return macro, user

def GeneratePosition(ActualNum, SimulationRegion):
    Position_x = []
    Position_y = []
    for i in range(ActualNum):
        Position_x.append(np.random.uniform((-SimulationRegion / 2), (SimulationRegion / 2)))
        Position_y.append(np.random.uniform((-SimulationRegion / 2), (SimulationRegion / 2)))

    return Position_x, Position_y

def SelectCooperatedBS(para, macro, pico, user):

    user.Alluser_selectedBS_index = []
    Alluser_selectedBS_index_copy = []
    selectedBS_index = []

    # For each UE, determine the received power from every BS
    # for a in range(user.ActualNum):
        # avg_receivedpower_MBS = []
        # for i in range(macro.ActualNum):            
        #     #avg_receivedpower_MBS.append(macro.power_W * ((macro.Position_x[i] - user.Position_x[a])**2 + (macro.Position_y[i] - user.Position_y[a])**2)**(1/2 * (- para.PathLoss_exponent)))
        #     avg_receivedpower_MBS.append( macro.power_W * math.pow( math.sqrt( (macro.Position_x[i] - user.Position_x[a])**2 + (macro.Position_y[i] - user.Position_y[a])**2 ), - para.PathLoss_exponent) )

        # # Select the best (selectedMBSnum) BSs to from a cluster    
        # selectedMBS = heapq.nlargest(para.selectedMBSnum, avg_receivedpower_MBS)
        # selectedBS_index = [avg_receivedpower_MBS.index(x) for x in selectedMBS]   
        # user.Alluser_selectedBS_index.append(sorted(selectedBS_index))

    #Alluser_selectedBS_index_copy = user.Alluser_selectedBS_index[:]
    Alluser_selectedBS_index_copy = [ [0, 1], [0, 1], [1, 2]]
    user.Alluser_selectedBS_index = [ [0, 1], [0, 1], [1, 2]]

    # If two users select the same set of BSs, remove the redudant one
    #print("Alluser_selectedBS_index_copy=", Alluser_selectedBS_index_copy)
    # for a in range(user.ActualNum):
    #     for i in range(1,user.ActualNum - a):
    #         if user.Alluser_selectedBS_index[a] == user.Alluser_selectedBS_index[a + i]:
    #             Alluser_selectedBS_index_copy[a + i] = 0
    #print("Alluser_selectedBS_index_copy1=", Alluser_selectedBS_index_copy) 
              
    # Copy Alluser_selectedBS_index_copy to cluster
    cluster = [] 
    # a = 0
    # 
    # for j in Alluser_selectedBS_index_copy:
    #     if j != 0:
    #         cluster.append(Cluster_father())
    #         cluster[a].mbs.index = j
    #         for k in cluster[a].mbs.index:
    #             cluster[a].mbs.Position_inCiuster_x.append(macro.Position_x[k])
    #             cluster[a].mbs.Position_inCiuster_y.append(macro.Position_y[k])
    #             #print("a=", a, cluster[a].mbs.index)
    #         a = a + 1
	
    cluster.append(Cluster_father())
    cluster[0].mbs.index = [0, 1]
    cluster.append(Cluster_father())
    cluster[1].mbs.index = [1, 2]
    for j in range(len(cluster)):
        for m in cluster[j].mbs.index:
            cluster[j].mbs.Position_inCluster_x.append(macro.Position_x[m])
            cluster[j].mbs.Position_inCluster_y.append(macro.Position_y[m])
    
            
    #print(Alluser_selectedBS_index_copy.__len__())
    print(len(cluster))    
    
    
    for i in range(user.ActualNum):   #Attach users 
        for l in range(len(cluster)):
            if cluster[l].mbs.index == user.Alluser_selectedBS_index[i]: #%%
                cluster[l].ue.Position_x.append(user.Position_x[i])
                cluster[l].ue.Position_y.append(user.Position_y[i])
                cluster[l].ue.index.append(i)
    
    for l in range(len(cluster)):
        print("UEs in cluster=", l, cluster[l].ue.index)

    pdb.set_trace()
    # for i in range(len(cluster)):
    #     for a in range(len(cluster[i].ue.Position_x)):
    #         dis_user_pbs = []
    #         for k in range(pico.ActualNum):
    #             dis_user_pbs.append(((pico.Position_x[k] - cluster[i].ue.Position_x)**2 + (pico.Position_y[k] - cluster[i].ue.Position_y[a])**2)**(1/2))
    #         cluster[a].pbs.Position_inCiuster_x.append(pico.Position_x[dis_user_pbs.index(min(dis_user_pbs))])
    #         cluster[a].pbs.Position_inCiuster_y.append(pico.Position_y[dis_user_pbs.index(min(dis_user_pbs))])


    return para, cluster, macro, user

def Cluster_GenerateChannelZFbeamParameter(para, cluster, macro):
    for l in range(len(cluster)):  # Initialize list (size = # of users)
        cluster[l].channel.Dis_User_MBSs = []
        cluster[l].channel.Dis_User_PBSs = []
        cluster[l].channel.User_BSs = []
        cluster[l].channel.User_MBSs = []
        cluster[l].channel.User_PBSs = []
        cluster[l].channel.User_BSs_normalize = []
        cluster[l].channel.G = []
        cluster[l].channel.G_Hermitian = []
        cluster[l].channel.G_psuedo = []
        cluster[l].channel.identity_matrix = []
        cluster[l].channel.ZF_beam_weight = []
        cluster[l].channel.ZF_beam_weight_normalize = []
        for i in range(len(cluster[l].ue.Position_x)):
            cluster[l].channel.Dis_User_MBSs.append([])
            cluster[l].channel.Dis_User_PBSs.append([])
            cluster[l].channel.User_BSs.append([])
            cluster[l].channel.User_MBSs.append([])
            cluster[l].channel.User_PBSs.append([])
            cluster[l].channel.User_BSs_normalize.append([])
            cluster[l].channel.G.append([])
            cluster[l].channel.G_Hermitian.append([])
            cluster[l].channel.G_psuedo.append([])
            cluster[l].channel.identity_matrix.append([])
            cluster[l].channel.ZF_beam_weight.append([])
            cluster[l].channel.ZF_beam_weight_normalize.append([])        

        for i in range(len(cluster[l].ue.Position_x)):
            for j in range(len(cluster)):  
                cluster[l].channel.Dis_User_MBSs[i].append([])
                cluster[l].channel.User_BSs[i].append([])
                cluster[l].channel.User_MBSs[i].append([])
                cluster[l].channel.User_BSs_normalize[i].append([])
            for j in range(len(cluster)):  # Initialize user channel as an array
                cluster[l].channel.Dis_User_MBSs[i][j] = np.zeros(macro.numAntenna * len(cluster[j].mbs.Position_inCluster_x))
                cluster[l].channel.User_MBSs[i][j] = np.zeros(macro.numAntenna * len(cluster[j].mbs.Position_inCluster_x), 'complex')
            cluster[l].channel.G[i] = np.zeros(((macro.numAntenna * len(cluster[l].mbs.Position_inCluster_x)), 1), 'complex')
            cluster[l].channel.identity_matrix[i] = np.identity((macro.numAntenna * len(cluster[l].mbs.Position_inCluster_x)))
    
        

    return para, cluster

def Cluster_GenerateChannel(para, cluster, macro):
    # MBS (PBS) channel
    norm = 0
    # for i in range(len(cluster)):
    #     for k in range(len(cluster[i].ue.Position_x)):
    #         for j in range(len(cluster)):            	
    #             for n in range(len(cluster[j].mbs.Position_inCluster_x)):                	
    #                 cluster[i].channel.Dis_User_MBSs[k][j][ n*macro.numAntenna ] = math.sqrt(math.pow((cluster[i].ue.Position_x[k]) - (cluster[j].mbs.Position_inCluster_x[n]),2) + \
    #                                                                                        math.pow((cluster[i].ue.Position_y[k]) - (cluster[j].mbs.Position_inCluster_y[n]),2))
    #                 cluster[i].channel.User_MBSs[k][j][n * macro.numAntenna] = math.sqrt(math.pow((cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna]),-para.PathLoss_exponent)) * \
    #                                                                          (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))                                                                                                 
    #                 for q in range(1, macro.numAntenna):
    #                     cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna + q] = cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna]
    #                     cluster[i].channel.User_MBSs[k][j][n * macro.numAntenna + q] = math.sqrt(math.pow((cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna + q]),-para.PathLoss_exponent)) * \
    #                                                                                    (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))
  
    #             #cluster[i].channel.User_BSs[k][j] = np.append(cluster[i].channel.User_MBSs[k][j])
    #             norm = np.linalg.norm(cluster[i].channel.User_MBSs[k][j])
    #             cluster[i].channel.User_BSs_normalize[k][j] = cluster[i].channel.User_MBSs[k][j] / norm

    for l in range(len(cluster)):
        for i in range(len(cluster[l].ue.Position_x)):
            for j in range(len(cluster)):               
                for m in range(len(cluster[j].mbs.Position_inCluster_x)):                   
                    cluster[l].channel.Dis_User_MBSs[i][j][ m*macro.numAntenna ] = math.sqrt(math.pow((cluster[l].ue.Position_x[i]) - (cluster[j].mbs.Position_inCluster_x[m]),2) + \
                                                                                           math.pow((cluster[l].ue.Position_y[i]) - (cluster[j].mbs.Position_inCluster_y[m]),2))
                    cluster[l].channel.User_MBSs[i][j][m * macro.numAntenna] = math.sqrt(math.pow((cluster[l].channel.Dis_User_MBSs[i][j][m * macro.numAntenna]),-para.PathLoss_exponent)) * \
                                                                             (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))                                                                                                 
                    for q in range(1, macro.numAntenna):
                        cluster[l].channel.Dis_User_MBSs[i][j][m * macro.numAntenna + q] = cluster[l].channel.Dis_User_MBSs[i][j][m * macro.numAntenna]
                        cluster[l].channel.User_MBSs[i][j][m * macro.numAntenna + q] = math.sqrt(math.pow((cluster[l].channel.Dis_User_MBSs[i][j][m * macro.numAntenna + q]),-para.PathLoss_exponent)) * \
                                                                                       (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))
  
                #cluster[l].channel.User_BSs[i][j] = np.append(cluster[l].channel.User_MBSs[i][j])
                norm = np.linalg.norm(cluster[l].channel.User_MBSs[i][j])
                cluster[l].channel.User_BSs_normalize[i][j] = cluster[l].channel.User_MBSs[i][j] / norm
    
    return para, cluster

def Cluster_GenerateZFbeam(para, cluster, macro):    
    for l in range(len(cluster)):  # G_i
        for i in range(len(cluster[l].ue.Position_x)): # target user
            for k in range(len(cluster[l].ue.Position_x)):  # intra-cluster user
                if i != k:
                    cluster[l].channel.G[i] = np.c_[ cluster[l].channel.G[i], cluster[l].channel.User_BSs_normalize[k][l] ]
            
            cluster[l].channel.G[i] = np.delete(cluster[l].channel.G[i], 0, 1)
            cluster[l].channel.G_Hermitian[i] = np.conj(cluster[l].channel.G[i]).transpose()
            cluster[l].channel.G_psuedo[i] = np.dot(inv(np.dot(cluster[l].channel.G_Hermitian[i], cluster[l].channel.G[i])), cluster[l].channel.G_Hermitian[i])
            cluster[l].channel.ZF_beam_weight[i] = np.dot((cluster[l].channel.identity_matrix[i] - np.dot(cluster[l].channel.G[i], cluster[l].channel.G_psuedo[i])), \
                                                          cluster[l].channel.User_BSs_normalize[i][l])
            norm = np.linalg.norm(cluster[l].channel.ZF_beam_weight[i])
            cluster[l].channel.ZF_beam_weight[i] = cluster[l].channel.ZF_beam_weight[i] / norm

    return para, cluster

def Cluster_ReceivedPower(para, cluster, macro):    
    # Initial empty array
    for i in range(len(cluster[0].ue.Position_x)):        
        cluster[0].ue.received_signal_power.append([])
        cluster[0].ue.received_interference_power.append([])
        cluster[0].ue.received_SINR_power.append([])
        cluster[0].ue.received_capacity_power.append([])

    # Create array with size (Num of users in cluster 0)x1
    for i in range(len(cluster[0].ue.Position_x)):
        cluster[0].ue.received_signal_power[i] = 0
        cluster[0].ue.received_interference_power[i] = 0
    
    para.TotalSignal = 0
    for i in range(len(cluster[0].ue.Position_x)):  # signal in Watt
        # print(len(cluster[0].mbs.Position_inCluster_x))
        cluster[0].ue.received_signal_power[i] = ((len(cluster[0].mbs.Position_inCluster_x) * macro.power_W)/ len(cluster[0].ue.Position_x)) * \
                                                   math.pow(abs(np.dot(np.conj(cluster[0].channel.User_MBSs[i][0]).transpose(),cluster[0].channel.ZF_beam_weight[i])), 2)


    para.TotalInterference = 0
    # for k in range(len(cluster[0].ue.Position_x)):  # intereference in Watt
    #     for j in range(1, len(cluster)):
    #         # print(len(cluster[j].ue.Position_x))
    #         for i in range(len(cluster[j].ue.Position_x)):
    #             cluster[0].ue.received_interference_power[k] = ((len(cluster[j].mbs.Position_inCluster_x) * macro.power_W ) / len(cluster[j].ue.Position_x)) * \
    #                                                              math.pow(abs(np.dot(np.conj(cluster[0].channel.User_BSs[k][j]).transpose(),cluster[j].channel.ZF_beam_weight[i])), 2) \
    #                                                              + cluster[0].ue.received_interference_power[k]
    
    # 
    for i in range(len(cluster[0].ue.Position_x)):  # intereference in Watt
        for j in range(1, len(cluster)): # other clusters than cluster 0
            # print(len(cluster[j].ue.Position_x))
            for k in range(len(cluster[j].ue.Position_x)): # inter-cluster users in cluster j
                cluster[0].ue.received_interference_power[i] = ((len(cluster[j].mbs.Position_inCluster_x) * macro.power_W ) / len(cluster[j].ue.Position_x)) * \
                                                                 math.pow(abs(np.dot(np.conj(cluster[0].channel.User_MBSs[i][j]).transpose(),cluster[j].channel.ZF_beam_weight[k])), 2) \
                                                                 + cluster[0].ue.received_interference_power[i]


    para.TotalCapacity = 0
    #SINRthreshold = [10] #[-10, -5, 0, 5, 10, 15]
    para.successnum = []
    for i in range(len(para.SINRthreshold)):
        para.successnum.append([])

    for i in range(len(cluster[0].ue.Position_x)):  # SINR, ergodic sum rate, average ergodic sum rate
        cluster[0].ue.received_SINR_power[i] = cluster [0].ue.received_signal_power[i] / (cluster[0].ue.received_interference_power[i] + 10**(-40))
        for a in range(len(para.SINRthreshold)):
            if cluster[0].ue.received_SINR_power[i] > dbtonumber(para.SINRthreshold[a]):
                para.successnum[a].append(1)
            else:
                para.successnum[a].append(0)
        cluster[0].ue.received_capacity_power[i] = math.log((1 + (cluster[0].ue.received_SINR_power[i] / dbtonumber(para.SNRGap_db))), 2)
        para.TotalCapacity = cluster[0].ue.received_capacity_power[i] + para.TotalCapacity

    if len(cluster[0].mbs.Position_inCluster_x) == 0:
        para.TotalCapacity = 0
    else:
        para.TotalCapacity = para.TotalCapacity / len(cluster[0].ue.Position_x)

    # if cluster[0].ue.Num == 0:
    #     para.TotalInterference = 0
    #     para.TotalSignal = 0
    # else:
    #     para.TotalInterference = para.TotalInterference / cluster[0].ue.Num
    #     para.TotalSignal = para.TotalSignal / cluster[0].ue.Num

    

    return para, cluster

def dbmtoW(a):
    b = 0.001 * math.pow(10, (a / 10))
    return b

def converttodb(a):
    b = 10 * math.log10(a)
    return b

def dbtonumber(a):
    b = math.pow(10, (a / 10))
    return b

# ************************************************************************************************************************************#
# *************************************************BS Parameters***************************************************************************#
SimulationTime = 1#100

# Macro
macro = Macro()
macro.numAntenna = 2
macro.Radius = 500                                                                      # Cell radius
macro.density = (macro.Radius ** 2 * math.pi) ** (-1)                                   # BS density
macro.power_dbm = 43                                                                    # BS power (dbm)
macro.power_W = dbmtoW(macro.power_dbm)                                                 # BS power(W)

# Pico
pico = Pico()
pico.numAntenna = 1
pico.density = macro.density * 0                                                       # BS density
pico.power_dbm = 24                                                                    # BS power (dbm)
pico.power_W = dbmtoW(macro.power_dbm)                                                 # BS power (W)

# User
user = USER()

# Para
para = Para()
para.LoadingFactor = 0.5
para.SINRthreshold = [10]
para.selectedMBSnum = 2 # cluster size

# SimulationRegion
para.SimulationRegion = 100#6000
# pathloss
para.PathLoss_exponent = 3.76
para.PassLoss_reference_d = 0.3920

# ************************************************* Main ************************************************************************#
selectedMBSnum = para.selectedMBSnum # [2, 4, 6, 8, 10]
SINRthreshold = para.SINRthreshold  #[-10, -5, 0, 5, 10, 15]
# loading = [0.6, 0.2, 0.6, 0.8, 1]
simulationcapacity = []
simulationcoverage = []
meancapacity = []
meancoverage= []

# Initialize result holder
simulatemeancapacity = []
simulatemeancoverage = []
for i in range(len(SINRthreshold)):
    meancoverage.append([])
    simulatemeancoverage.append([])

for simulationtime in range(SimulationTime):
    macro = GenerateMBSPosition(para, macro)
    user = USER()
    macro, user = GenerateUserPosition(para, macro, user)
    para, cluster, macro, user = SelectCooperatedBS(para, macro, pico, user)
    para, cluster = Cluster_GenerateChannelZFbeamParameter(para, cluster, macro)
    para, cluster = Cluster_GenerateChannel(para, cluster, macro)
    para, cluster = Cluster_GenerateZFbeam(para, cluster,macro)
    para, cluster = Cluster_ReceivedPower(para, cluster, macro)
    simulatemeancapacity.append(para.TotalCapacity)
    for i in range(len(SINRthreshold)):
        simulatemeancoverage[i] = simulatemeancoverage[i] + para.successnum[i]
        print(simulationtime, simulatemeancoverage)

meancapacity.append(np.mean(simulatemeancapacity))
for i in range(len(SINRthreshold)):
    meancoverage[i].append(np.mean(simulatemeancoverage[i]))

print(meancapacity)
print(meancoverage)



# for i in range(len(macro.UEPosition_x[0])):
#      plt.plot(macro.UEPosition_x[0][i], macro.UEPosition_y[0][i], 'ro-')
# for i in range(len(macro.UEPosition_x[1])):
#      plt.plot(macro.UEPosition_x[1][i], macro.UEPosition_y[1][i], 'yo-')
# for i in range(len(user.Position_y)):
#      plt.plot(user.Position_x[i], user.Position_y[i], 'bo-')
# for i in range(len(macro.Position_x)):
#      plt.plot(macro.Position_x[i], macro.Position_y[i], 'go-')
# for i in range(len(cluster)):
#     for j in range(len(cluster[i].ue.Position_x)):
#         plt.plot(cluster[i].ue.Position_x[j], cluster[i].ue.Position_y[j], 'go-')
# plt.axis([-para.SimulationRegion/2,para.SimulationRegion/2,-para.SimulationRegion/2,para.SimulationRegion/2])
# plt.show()
