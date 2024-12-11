import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# 定义高斯函数
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def cha_dian(chain_len_lsts):
    max1=max(chain_len_lsts)
    min1=min(chain_len_lsts)
    cha_dian=[]
    for i in range(min1,max1+1):
        cha_dian.append(i)
    return cha_dian

def gaussian_fit(xdata,ydata):# 拟合数据
    popt, pcov = curve_fit(gaussian, xdata, ydata,maxfev=100000)


    
    # 绘制原始数据和拟合曲线
    # plt.plot(xdata, ydata, 'b-', label='data')
    # plt.plot(xdata, gaussian(xdata, *popt), 'r-', label='fit')
    # plt.legend()
    # plt.show()
    return popt

# temperature_lst=[400,450,500,550,600],n_maxmono=15000,put_in=6e-8,k_termin_c=0.000008,k_termin_p=-1.20E-19
# T=400 [2, 7] [5590, 189]
# T=450 [2, 7, 12, 17] [1178, 1494, 159, 5]
# T=500 [2, 7, 12, 17, 22, 27, 32] [129, 248, 336, 287, 135, 40, 1]
# T=550 [2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62] [20, 37, 51, 53, 61, 61, 69, 66, 60, 42, 13, 6, 1]
# T=600 [2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97, 102, 107, 112, 117] [11, 15, 16, 18, 20, 17, 29, 21, 18, 24, 25, 19, 18, 17, 9, 12, 6, 8, 7, 7, 4, 5, 1, 1]
temperature_lst=[450,500,550,600]
chain_len_lsts=[[2, 7, 12, 17],[2, 7, 12, 17, 22, 27, 32],[2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62],[2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 92, 97, 102, 107, 112, 117]]
count_lsts=[[1178, 1494, 159, 5],[129, 248, 336, 287, 135, 40, 1],[20, 37, 51, 53, 61, 61, 69, 66, 60, 42, 13, 6, 1],[11, 15, 16, 18, 20, 17, 29, 21, 18, 24, 25, 19, 18, 17, 9, 12, 6, 8, 7, 7, 4, 5, 1, 1]]
for i in range(len(chain_len_lsts)):
    print(temperature_lst[i])
    popt=gaussian_fit(chain_len_lsts[i],count_lsts[i])
    bars=plt.bar(chain_len_lsts[i],count_lsts[i],label='T='+str(temperature_lst[i])+'K')
    bar_color = bars[i].get_facecolor()
    plt.plot(cha_dian(chain_len_lsts[i]), gaussian(cha_dian(chain_len_lsts[i]), *popt),'-')
plt.xlabel('Chain length')
plt.ylabel('Count')
plt.legend()
plt.show()
