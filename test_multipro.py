import multiprocessing
import re
from joblib import Parallel, delayed

import matplotlib.pyplot as plt
import numpy as np
import math



# dif = [1,2,3,4,5,6,7,8,9,10]
# cest = [0,4,8,12,19,26,34,39,40,56,78]
# taille = [2,3,4,5,6]

# for d in dif:
#     for c in cest:
#         for t in taille:
#             y = 1 - abs(math.log10(d))
#             y = (y * c) * t
#             print('dif=',d, 'CEST=',c, 'len=',t,'=',y)
# exit()

'''
population = ['(([CHALQG]H{4})|((N+)(D+)))|((A{5}|W{4})|(T{3}(A|G)))',
              '((H{3}N{2})((L|C)[^SREQPALHKC]))(((F|K)[NLKA])|((M|V)|(AA)))',
              '(([^CNLEGFMWKQDAI][GRDIWN])((W|G)|(E|A)))|(([^KFQLDSRAETHYVWNMIP]|(C|F))([SVWHQIPMEGT]V{2}))',
              'N|E',
              '(((F|D)|(K+))|([^PVMICGKFRDSAHEQLY]|[^YICHRTVSLKQEANWDGP]))|(((RV)(A+))|(L{3}(FT)))',
              '[^QNFSEWLCAIHYR]L',
              '(((MG)(L|C))((M|Q)|(A+)))|(([CMNVFKHAGTS](DT))|(C{4}|(T|I)))',
              '(ID)|K',
              '(([Q]|T{2})|((LH)|G{2}))(([^KDVAHLNGMYP]|(H+))|([PIHDNLTYAV][PRCATSDGL]))']

training_set = ['LWSDIKMKLKKT','KMGKLIGIPVLK','MKVAAAMAPKQV','AKCKVQSANVCK','KMWDWEQKKKWI',
                'WDWEQKKKWI--','RLPKRVQGNVEK','VTRMTIQVKGSK','TSKSKKRMTAKK','ICLKSQPICGID','LELKLGKRPMGW',
                'ARNRKKIMMRWI','PVNRLGKMSKNR','IGVLRSVKQTVR','NKQRRMLSRERS','NVVVQRRNHHTS','VGSVKSGNLRMR',
                'EPSNLPKGMNEK','QCAGWVQKRQIQ','KGKLDKDRNL--','KYTKTRKQSSKA','RRRRRRRRRRRR',
                'VCNRIEPLKPIL','ERQEEKIKKW--','CQLAWRPCAKAS','WFGLQRHLKKKD','RTTTRTTTRTTT','RTRTRTRTRTRT',
                'PVARKVVQICHP','CCWHNPKWRRTR','PRSWEVKEKETM','NFLRAQRQCQKQ','RGRGRGRGRGRG','RGKMPLRWMTRK',
                'HWSTCTRTRTLS','KSKSKSKSKSKS','SSDQDRDKWL--','KPCKWAGRACAK','RTTRTTRTTRTT','MWQMKWTRKTRE',
                'KPWHGCASRTKR','HIKWRLTKGTRT','VSLQCWELGPNK','HGRKWKRTKFDD','KRIIEDDQLE--','RRMVNRTITRMW',
                'GNKKNWRWYKNR','AQCCQHRKGYMN','KWVVRPRIRRLL','VSVVATGCVWET','TKTKTKTKTKTK','NSSNHSNNMPCQ','IRRWNDRIRITS',
                'TTTKTTTKTTTK','TTTTTKTTTTTK','PIKQIAWPIIEH','RKHHGWRWEQWK','KVLWRMPAQIIQ','WDRTSTRPSSVL','KEEVWLKWLI--',
                'MWVKGMKHKKMK','KSSKSSKSSKSS','KSSSKSSSKSSS','IIRSPICCVSRV','RSRSRSRSRSRS','GRKRGAIWKDTK','KHKHKHKHKHKH',
                'LSQQPRKRATWR','TTKTTKTTKTTK','DKVCKIQKRKWH','KKKKKKKKKKKK','RLWNSGEGRGEN','KGGGKGGGKGGG','QCRAGAMPAMYV',
                'KKRLHWIRWHCG','KKRKKHKKGKKP','GQRWLYKMKDSM','KGGKGGKGGKGG','LDHTWGKWGHQS','KKGKKKGKKHKK',
                'PGGVRSNDLLEV','DKRKIKQKMWWG','KGKGKGKGKGKG','CHLKDLRKMGLR','IKGMNIKMPTDQ','KKAKKKGKKHKK',
                'RRCQAQEFWLGA','KKGKKKGKKPKK','DSSSSSDSSSSS','KVIRYVVAPMKL','NAPWKHWRIINE','GNCPMKVCSPMG','NDISMCNKNNNW',
                'KDRTSKPKRPWC','SDGSKIKDRD--','KVRCLVEARPSW','VINKVISCPCVN','VINKVISNPCVN','IRTYLRKRNSTQ','TVSEPVMMVSVS',
                'GIFKTTKCKHNS','QRHDSHRHGLWL','HDDKNKESDD--','DTTTTTDTTTTT','LSNRRGREQYAG','DSSSDSSSDSSS','MAMADAAAPMNA',
                'NQYSNWNKNYK-','EMRQWKWMWENA','WWWKPKREDFMK','TTTTTTTTTTTT','RPPMLNVVRVVG','GGRVWEWNVAA-','SNHKMSECRGLR','ETTTTTETTTTT',
                'GPMPMNAKMKLC','ETNVRVKVVSES','RHRHRHRHRHRH','NNKCQVVAAFVM','HLVVSPRVSWGC','FNSNKITPTSNM','VPNIQVKGSK--',
                'VNLPMVMPNLRM','ETTETTETTETT','PVVYKTVIQCCD','VAWVMKAHVCTM','VNSDPSNGQMRD','DTTDTTDTTDTT',
                'DTTTDTTTDTTT','ETTTETTTETTT','NRVTESVRNVKM','QTATENSQMNSG','GLGNQHVVVLGV','NWRDCLSLIVPN','LLRLLGLVER--',
                'GLIEARAMQQCC','DSDSDSDSDSDS','QERRDDILWD--','DTDTDTDTDTDT','ETETETETETET','QTEHYENSARNS',]

def evaluate_fitness(individual):
    # Evaluation de la fitness de l'individu
    regex = re.compile(individual)
    fitness = 0
    for sequence in training_set:
        for match in re.finditer(regex, sequence): 
            motif = match.group()
            if len(motif)>= 3:
                fitness+=1

    return fitness

def parallel_ga(population):
    pool = multiprocessing.Pool()
    fitness_values = pool.map(evaluate_fitness, population)
    print(fitness_values)
    pool.close()
    pool.join()
    return fitness_values


list_fitness = parallel_ga(population)

dico = {}
for indi, fitness in zip(population, list_fitness):
    dico[indi] = fitness



# ou 
def parallel_ga(population):
    fitness_values = Parallel(n_jobs=-1)(delayed(evaluate_fitness)(individual) for individual in population)
    return fitness_values


# ou
def evaluate_fitness(individual, argument1, argument2):
    # Evaluation de la fitness de l'individu
    return fitness

def parallel_ga(population, argument1, argument2):
    pool = multiprocessing.Pool()
    fitness_function = functools.partial(evaluate_fitness, argument1=argument1, argument2=argument2)
    fitness_values = pool.map(fitness_function, population)
    pool.close()
    pool.join()
    return fitness_values

print(dico)
'''
# ls = [-50,-21,-3,2,14,26,58] #score
# lx = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24] # nbr regle
# lx=[]
# for i in range(24):
#     lx.append(i)

# ly1 = []
# ly2 = []
# ly3 = []
# ly4 = []
# ly5 = []
# ly10 = []
# std_dev = 3
# mean = 12

# for x in lx:
#     y1 = (1/(1 * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((x - mean) / 0.5) ** 2)
#     y2 = (1/(2 * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((x - mean) / 2) ** 2)
#     y3 = (1/(3 * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((x - mean) / 3) ** 2)
#     y4 = (1/(4 * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((x - mean) / 4) ** 2)
#     y5 = (1/(5 * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((x - mean) / 5) ** 2)
#     y10 = (1/(10 * math.sqrt(2 * math.pi))) * np.exp(-0.5 * ((x - mean) / 10) ** 2)
#     ly1.append(y1)
#     ly2.append(y2)
#     ly3.append(y3)
#     ly4.append(y4)
#     ly5.append(y5)
#     ly10.append(y10)


# r = []
# for y,s in ls:
# for y in ly5:
#     s = -25
#     z = s*y
#     x = s/y
#     # print(s,'*',y,'=', z)
#     # print(s,'/',y,'=', x)
#     r.append(z)
# plt.plot(lx, r, label='4')
# plt.show()

# exit()

# plt.plot(lx, ly1, label='1')
# plt.plot(lx, ly2, label='2')
# plt.plot(lx, ly3, label='3')
# plt.plot(lx, ly4, label='4')
# plt.plot(lx, ly5, label='5')
# plt.plot(lx, ly10, label='10')
# plt.legend()
# plt.show()


# def logistic(x):
#     return 1/(1 / (np.exp(-x)))

# def logistic_function(x, k, x0, L):
#     return 1/(L / (np.exp(-k*(x-x0))))

# x1 = np.linspace(0, 100, 100)
# x2 = np.linspace(0, 100, 100)
# x3 = np.linspace(0, 100, 100)
# x4 = np.linspace(0, 100, 100)
# x5 = np.linspace(0, 100, 100)

# y1 = logistic_function(1,  0.05, 0, 1)
# y2 = logistic_function(10, 0.05, 0, 1)
# y3 = logistic_function(25, 0.05, 0, 1)
# y4 = logistic_function(37, 0.05, 0, 1)
# y5 = logistic_function(68, 0.05, 0, 1)
# y6 = logistic_function(99, 0.05, 0, 1)

# print(y1,y2,y3,y4,y5, y6)

# plt.scatter(1, y1, label='1')
# plt.scatter(1, y2, label='10')
# plt.scatter(1, y3, label='25')
# plt.scatter(1, y4, label='37')
# plt.scatter(1, y5, label='68')
# plt.scatter(1, y6, label='99')
# plt.xlabel("X")
# plt.ylabel("Logistic(X)")
# plt.title("Graph of Logistic Function")
# plt.legend()
# plt.show()

# x = np.linspace(0, 100, num=100)
# y = logistic(x)

# plt.plot(x, y)

# plt.show()


l = ['a','b','c','d']
n = [4,6,7,2]

for i,j in zip(l,n):
    print(i, j)