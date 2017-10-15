import math
import random
import numpy as np
from numpy.linalg import inv
import heapq
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from sys import exit

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
        self.density = 0
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
        self.Position_inCiuster_x = []
        self.Position_inCiuster_y = []

class PBS:
    def __init__(self):
        self.Num = 0
        self.Position_inCiuster_x = []
        self.Position_inCiuster_y = []

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

    macro.ActualNum = np.random.poisson(para.SimulationRegion * para.SimulationRegion * macro.density)  # Total number of MBS
    while macro.ActualNum < para.selectedMBSnum:
        macro.ActualNum = np.random.poisson(para.SimulationRegion * para.SimulationRegion * macro.density)  

    print("Num of MBS=", macro.ActualNum)

    Position_x, Position_y = GeneratePosition(macro.ActualNum, para.SimulationRegion)

    # Distance from each MBS to origin
    dis_MBS_center = []
    for i in range(macro.ActualNum):
        #dis_MBS_center.append(((Position_x[i])**2 + (Position_y[i])**2)**(1/2)) #--> The result is always 1.0
        dis_MBS_center.append( math.sqrt( (Position_x[i])**2 + (Position_y[i])**2 ) )

    # Arrange MBS according to their distance to origin from nearest to farest
    for i in sorted(dis_MBS_center):
        macro.Position_x.append(Position_x[dis_MBS_center.index(i)])
        macro.Position_y.append(Position_y[dis_MBS_center.index(i)])

    

    return macro

def GeneratePBSPosition(para, pico):
    Position_x = []
    Position_y = []

    pico.ActualNum = np.random.poisson(para.SimulationRegion * para.SimulationRegion * pico.density)  # Total number of PMS
    pico.Position_x, pico.Position_y = GeneratePosition(pico.ActualNum, para.SimulationRegion)

    return pico

def GenerateUserPosition(para, user):
    user.Position_x = []
    user.Position_y = []

    user.ActualNum = np.random.poisson(para.SimulationRegion * para.SimulationRegion * user.density)
    while user.ActualNum < 1:
        user.ActualNum = np.random.poisson(para.SimulationRegion * para.SimulationRegion * user.density)
        print(user.ActualNum)

    user.Position_x, user.Position_y = GeneratePosition(user.ActualNum, para.SimulationRegion)

    print("Num of users=", user.ActualNum)
    return user

# def GenerateUserPosition(para, macro, user):
#     # Divide simulation region in grids, each with size 50x50
#     # Place one UE at the intersection point 
#     # __________
#     # |__|__|__|
#     # |__|__|__|
#     # |  |  |  |
#     #

#     Dis_BS_ClusterCenter = np.zeros(macro.ActualNum)
#     user_index_x = np.arange((-para.SimulationRegion/2), (para.SimulationRegion/2 + 1), 50)
#     user_index_y = np.arange((-para.SimulationRegion/2), (para.SimulationRegion/2 + 1), 50)

#     # Create an empty array for BS to store the attached UEs
#     macro.UEPosition_x = []
#     macro.UEPosition_y = []
#     for j in range(macro.ActualNum):
#         macro.UEPosition_x.append([])
#         macro.UEPosition_y.append([])


#     for k in range(len(user_index_y)):
#         for i in range(len(user_index_x)):
#             for j in range(macro.ActualNum):
#                 # For MBS j, determine its distance to every UE in row k
#                 Dis_BS_ClusterCenter[j] = math.sqrt(math.pow((user_index_x[i]) - (macro.Position_x[j]), 2) + math.pow((user_index_y[k]) - (macro.Position_y[j]), 2))            
#             minindex = Dis_BS_ClusterCenter.argmin() # Find the nearest UE in row k             
#             macro.UEPosition_x[minindex].append(user_index_x[i])
#             macro.UEPosition_y[minindex].append(user_index_y[k])
 
#     usernum = int(macro.numAntenna * para.LoadingFactor) # Number of users to serve per BS
#     user.Position_x = []
#     user.Position_y = []
#     for j in range(macro.ActualNum):
#         selecteduser = []  
#         #print("len(macro.UEPosition_x[j])=", len(macro.UEPosition_x[j]))      
#         # For each BS, randomly choose (usernum) UE from the attached ones to serve
#         if len(macro.UEPosition_x[j]) >= usernum:
#             selecteduser = random.sample(range(len(macro.UEPosition_x[j])), usernum)
#         else:
#             selecteduser = random.sample(range(1000), usernum)
            
        
#         # Determine the locations of served UEs
#         for a in selecteduser:
#             user.Position_x.append(macro.UEPosition_x[j][a])
#             user.Position_y.append(macro.UEPosition_y[j][a])

#     user.ActualNum = len(user.Position_x) # Total number of served users

#     return macro, user

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
    for a in range(user.ActualNum):
        avg_receivedpower_MBS = []
        for i in range(macro.ActualNum):            
            #avg_receivedpower_MBS.append(macro.power_W * ((macro.Position_x[i] - user.Position_x[a])**2 + (macro.Position_y[i] - user.Position_y[a])**2)**(1/2 * (- para.PathLoss_exponent)))
            avg_receivedpower_MBS.append( macro.power_W * math.pow( math.sqrt( (macro.Position_x[i] - user.Position_x[a])**2 + (macro.Position_y[i] - user.Position_y[a])**2 ), - para.PathLoss_exponent) )

        # Select the best (selectedMBSnum) BSs to form a cluster    
        selectedMBS = heapq.nlargest(para.selectedMBSnum, avg_receivedpower_MBS)
        selectedBS_index = [avg_receivedpower_MBS.index(x) for x in selectedMBS]   
        user.Alluser_selectedBS_index.append(sorted(selectedBS_index))

    Alluser_selectedBS_index_copy = user.Alluser_selectedBS_index[:]

    # If two users select the same set of BSs, remove the redudant one
    #print("Alluser_selectedBS_index_copy=", Alluser_selectedBS_index_copy)
    for a in range(user.ActualNum):
        for i in range(1,user.ActualNum - a):
            if user.Alluser_selectedBS_index[a] == user.Alluser_selectedBS_index[a + i]:
                Alluser_selectedBS_index_copy[a + i] = 0
    #print("Alluser_selectedBS_index_copy1=", Alluser_selectedBS_index_copy) 
              

    # Copy Alluser_selectedBS_index_copy to cluster
    cluster = [] 
    a = 0
    for j in Alluser_selectedBS_index_copy:
        if j != 0:
            cluster.append(Cluster_father())
            cluster[a].mbs.index = j
            for k in cluster[a].mbs.index:
                cluster[a].mbs.Position_inCiuster_x.append(macro.Position_x[k])
                cluster[a].mbs.Position_inCiuster_y.append(macro.Position_y[k])
            a = a + 1
    #print(Alluser_selectedBS_index_copy.__len__())
    #print(len(cluster))
    #print(user.ActualNum)


    for j in range(user.ActualNum):   #Attach users 
        for i in range(len(cluster)):
            if cluster[i].mbs.index == user.Alluser_selectedBS_index[j]:
                cluster[i].ue.Position_x.append(user.Position_x[j])
                cluster[i].ue.Position_y.append(user.Position_y[j])

    # for i in range(len(cluster)):
    #     for a in range(len(cluster[i].ue.Position_x)):
    #         dis_user_pbs = []
    #         for k in range(pico.ActualNum):
    #             dis_user_pbs.append(((pico.Position_x[k] - cluster[i].ue.Position_x)**2 + (pico.Position_y[k] - cluster[i].ue.Position_y[a])**2)**(1/2))
    #         cluster[a].pbs.Position_inCiuster_x.append(pico.Position_x[dis_user_pbs.index(min(dis_user_pbs))])
    #         cluster[a].pbs.Position_inCiuster_y.append(pico.Position_y[dis_user_pbs.index(min(dis_user_pbs))])


    return para, cluster, macro, user

def Cluster_GenerateChannelZFbeamParameter(para, cluster, macro):
    for i in range(len(cluster)):  # Initialize list (size = # of users)
        cluster[i].channel.Dis_User_MBSs = []
        cluster[i].channel.Dis_User_PBSs = []
        cluster[i].channel.User_BSs = []
        cluster[i].channel.User_MBSs = []
        cluster[i].channel.User_PBSs = []
        cluster[i].channel.User_BSs_normalize = []
        cluster[i].channel.G = []
        cluster[i].channel.G_Hermitian = []
        cluster[i].channel.G_psuedo = []
        cluster[i].channel.identity_matrix = []
        cluster[i].channel.ZF_beam_weight = []
        cluster[i].channel.ZF_beam_weight_normalize = []
        for k in range(len(cluster[i].ue.Position_x)):
            cluster[i].channel.Dis_User_MBSs.append([])
            cluster[i].channel.Dis_User_PBSs.append([])
            cluster[i].channel.User_BSs.append([])
            cluster[i].channel.User_MBSs.append([])
            cluster[i].channel.User_PBSs.append([])
            cluster[i].channel.User_BSs_normalize.append([])
            cluster[i].channel.G.append([])
            cluster[i].channel.G_Hermitian.append([])
            cluster[i].channel.G_psuedo.append([])
            cluster[i].channel.identity_matrix.append([])
            cluster[i].channel.ZF_beam_weight.append([])
            cluster[i].channel.ZF_beam_weight_normalize.append([])
        for k in range(len(cluster[i].ue.Position_x)):
            for j in range(len(cluster)):  # Initially each user has 19 channelos
                cluster[i].channel.Dis_User_MBSs[k].append([])
                cluster[i].channel.Dis_User_PBSs[k].append([])
                cluster[i].channel.User_BSs[k].append([])
                cluster[i].channel.User_MBSs[k].append([])
                cluster[i].channel.User_PBSs[k].append([])
                cluster[i].channel.User_BSs_normalize[k].append([])
            for j in range(len(cluster)):  # Initialize user channel as an array
                cluster[i].channel.Dis_User_MBSs[k][j] = np.zeros(macro.numAntenna * len(cluster[j].mbs.Position_inCiuster_x))
                # cluster[i].channel.Dis_User_PBSs[k][j] = np.zeros(pico.numAntenna * cluster[j].pbs.Num)
                cluster[i].channel.User_MBSs[k][j] = np.zeros(macro.numAntenna * len(cluster[j].mbs.Position_inCiuster_x), 'complex')
                # cluster[i].channel.User_PBSs[k][j] = np.zeros(pico.numAntenna * cluster[j].pbs.Num, 'complex')
            cluster[i].channel.G[k] = np.zeros(((macro.numAntenna * len(cluster[i].mbs.Position_inCiuster_x)), 1), 'complex')
            cluster[i].channel.identity_matrix[k] = np.identity((macro.numAntenna * len(cluster[i].mbs.Position_inCiuster_x)))
            # cluster[i].channel.G[k] = np.zeros(((macro.numAntenna * cluster[i].mbs.Num + pico.numAntenna * cluster[i].pbs.Num), 1), 'complex')
            # cluster[i].channel.identity_matrix[k] = np.identity((macro.numAntenna * cluster[i].mbs.Num + pico.numAntenna * cluster[i].pbs.Num))

    return para, cluster

def Cluster_GenerateChannel(para, cluster, macro):
    # MBS (PBS) channel
    norm = 0
    for i in range(len(cluster)):
        for k in range(len(cluster[i].ue.Position_x)):
            for j in range(len(cluster)):
                for n in range(len(cluster[j].mbs.Position_inCiuster_x)):
                    cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna] = math.sqrt(math.pow((cluster[i].ue.Position_x[k]) - (cluster[j].mbs.Position_inCiuster_x[n]),2) + \
                                                                                           math.pow((cluster[i].ue.Position_y[k]) - (cluster[j].mbs.Position_inCiuster_y[n]),2))
                    cluster[i].channel.User_MBSs[k][j][n * macro.numAntenna] = math.sqrt(math.pow((cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna]),-para.PathLoss_exponent)) * \
                                                                             (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))
                    for q in range(1, macro.numAntenna):
                        cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna + q] = cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna]
                        cluster[i].channel.User_MBSs[k][j][n * macro.numAntenna + q] = math.sqrt(math.pow((cluster[i].channel.Dis_User_MBSs[k][j][n * macro.numAntenna + q]),-para.PathLoss_exponent)) * \
                                                                                       (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))
                # for n in range(cluster[j].pbs.Num):
                #     cluster[i].channel.Dis_User_PBSs[k][j][n * pico.numAntenna] = math.sqrt(math.pow((cluster[i].ue.ActualPosition_x[k]) - (cluster[j].pbs.Position_inCiuster_x[n]), 2) + \
                #                                                                             math.pow((cluster[i].ue.ActualPosition_y[k]) - (cluster[j].pbs.Position_inCiuster_y[n]), 2))
                #     cluster[i].channel.User_PBSs[k][j][n * pico.numAntenna] = math.sqrt(math.pow(1 + (cluster[i].channel.Dis_User_PBSs[k][j][n * pico.numAntenna] / para.PassLoss_reference_d),-para.PathLoss_exponent)) * \
                #                                                               (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))
                #     for q in range(1, pico.numAntenna):
                #         cluster[i].channel.Dis_User_PBSs[k][j][n * pico.numAntenna + q] = cluster[i].channel.Dis_User_PBSs[k][j][n * pico.numAntenna]
                #         cluster[i].channel.User_PBSs[k][j][n * pico.numAntenna + q] = math.sqrt(math.pow(1 + (cluster[i].channel.Dis_User_PBSs[k][j][n * pico.numAntenna + q] / para.PassLoss_reference_d), -para.PathLoss_exponent)) * \
                #                                                                         (np.random.normal(0,1) + 1j * np.random.normal(0, 1)) * (1 / math.sqrt(2))
                cluster[i].channel.User_BSs[k][j] = np.append(cluster[i].channel.User_MBSs[k][j],cluster[i].channel.User_PBSs[k][j])
                norm = np.linalg.norm(cluster[i].channel.User_BSs[k][j])
                cluster[i].channel.User_BSs_normalize[k][j] = cluster[i].channel.User_BSs[k][j] / norm

    return para, cluster

def Cluster_GenerateZFbeam(para, cluster, macro):
    for j in range(len(cluster)):  # G
        for i in range(len(cluster[j].ue.Position_x)):
            for k in range(len(cluster[j].ue.Position_x)):
                if i != k:
                    cluster[j].channel.G[i] = np.c_[cluster[j].channel.G[i], cluster[j].channel.User_BSs_normalize[k][j]]
            cluster[j].channel.G[i] = np.delete(cluster[j].channel.G[i], 0, 1)
            cluster[j].channel.G_Hermitian[i] = np.conj(cluster[j].channel.G[i]).transpose()
            cluster[j].channel.G_psuedo[i] = np.dot(inv(np.dot(cluster[j].channel.G_Hermitian[i], cluster[j].channel.G[i])), cluster[j].channel.G_Hermitian[i])
            cluster[j].channel.ZF_beam_weight[i] = np.dot((cluster[j].channel.identity_matrix[i] - np.dot(cluster[j].channel.G[i], cluster[j].channel.G_psuedo[i])), \
                                                          cluster[j].channel.User_BSs_normalize[i][j])
            norm = np.linalg.norm(cluster[j].channel.ZF_beam_weight[i])
            cluster[j].channel.ZF_beam_weight[i] = cluster[j].channel.ZF_beam_weight[i] / norm

    return para, cluster

def Cluster_ReceivedPower(para, cluster, macro):
    for i in range(len(cluster[0].ue.Position_x)):
        cluster[0].ue.received_signal_power.append([])
        cluster[0].ue.received_interference_power.append([])
        cluster[0].ue.received_SINR_power.append([])
        cluster[0].ue.received_capacity_power.append([])

    for i in range(len(cluster[0].ue.Position_x)):
        cluster[0].ue.received_signal_power[i] = 0
        cluster[0].ue.received_interference_power[i] = 0

    para.TotalSignal = 0
    for i in range(len(cluster[0].ue.Position_x)):  # signal_W
        # print(len(cluster[0].mbs.Position_inCiuster_x))
        cluster[0].ue.received_signal_power[i] = ((len(cluster[0].mbs.Position_inCiuster_x) * macro.power_W)/ len(cluster[0].ue.Position_x)) * \
                                                   math.pow(abs(np.dot(np.conj(cluster[0].channel.User_BSs[i][0]).transpose(),cluster[0].channel.ZF_beam_weight[i])), 2)


    para.TotalInterference = 0
    for k in range(len(cluster[0].ue.Position_x)):  # intereference_W
        for j in range(1, len(cluster)):
            # print(len(cluster[j].ue.Position_x))
            for i in range(len(cluster[j].ue.Position_x)):
                cluster[0].ue.received_interference_power[k] = ((len(cluster[j].mbs.Position_inCiuster_x) * macro.power_W ) / len(cluster[j].ue.Position_x)) * \
                                                                 math.pow(abs(np.dot(np.conj(cluster[0].channel.User_BSs[k][j]).transpose(),cluster[j].channel.ZF_beam_weight[i])), 2) \
                                                                 + cluster[0].ue.received_interference_power[k]

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

    if len(cluster[0].mbs.Position_inCiuster_x) == 0:
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
SimulationTime = 1

# Macro
macro = Macro()
macro.numAntenna = 5
#macro.Radius = 10                                                                      # Cell radius
macro.density = pow(10,-3)                                   # BS density
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
user.density = pow(10,-2)
print(user.density)

# Para
para = Para()
para.LoadingFactor = 0.4
para.SINRthreshold = [10]
para.selectedMBSnum = 1 # cluster size

# SimulationRegion
para.SimulationRegion = 100
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
    #user = USER()
    user = GenerateUserPosition(para, user)
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