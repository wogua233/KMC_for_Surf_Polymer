##author: Wang Zicong  2024.12.7##
from .MC_run import mc_run
import math
for T in range(300,601,10):
    kAA=math.exp((-1.49E-19)/((1.38E-23)*T))
    kAB=math.exp((-1.74E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    kBB=math.exp((-1.90E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    kBA=math.exp((-1.74E-19)/((1.38E-23)*T))
    kA=kAA+kAB
    kB=kBB+kBA
    r1=kAA/kAB
    r2=kBB/kBA
    k_stop=math.exp((-1.17E-19)/((1.38E-23)*T))*0.00005 #设定终止速率常数
    