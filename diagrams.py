##author: Wang Zicong  2024.12.7##
import numpy as np
from MC_run import *
import matplotlib.pyplot as plt

def inst_mcrun_to_T(temperature_lst,put_in=2e-10,n_maxmono=1000): # 实例化mc_run类成一个mc_run实例列表
    mc_run_lst=[]
    for T in temperature_lst:
        T=T # K
        p_in=put_in
        kAA=math.exp((-1.49E-19)/((1.38E-23)*T))
        kAB=math.exp((-1.74E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
        kBB=math.exp((-1.90E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
        kBA=math.exp((-1.74E-19)/((1.38E-23)*T))
        KA=kAA+kAB
        KB=kBB+kBA
        rA=kAA/KA
        rB=kBB/KB
        k_termin=math.exp((-1.20E-19)/((1.38E-23)*T))*0.0000005 #设定终止速率常数
        # 实例化mc_run类
        mc_run_lst.append(mc_run(n_maxmono=n_maxmono,p_in=p_in,K_A=KA,K_B=KB,rA=rA,rB=rB,K_termin=k_termin,seed=1))
    return mc_run_lst

def mean_chain_len_to_T(temperature_lst,n_maxmono=1000): # 温度对链长的影响
    mc_run_lst=inst_mcrun_to_T(temperature_lst,n_maxmono=n_maxmono)
    mean_chain_len_lst=[]
    for i in mc_run_lst:
        i.run(run_statis=True,run_config=False)
        # 输出结果
        mean_chain_len_lst.append((i.do_event.n_maxmono-i.do_event.n_unic_mono)/i.do_event.n_chains)
    plt.plot(temperature_lst, mean_chain_len_lst, marker='o')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Mean Chain Length')
    plt.title('Mean Chain Length vs Temperature')
    plt.grid(True)
    plt.show()

def chain_len_to_p_in(temperature_lst,p_in_lst,n_maxmono=1000): # 温度对链长与输入速率的关系
    num_plots = len(p_in_lst)
    num_rows = int(np.ceil(num_plots / 2))
    num_cols = 2 if num_plots > 1 else 1
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(10, 8))
    axes = axes.flatten() if num_plots > 1 else [axes]
    for i in range(len(p_in_lst)):
        mc_run_lst=inst_mcrun_to_T(temperature_lst,put_in=p_in_lst[i],n_maxmono=n_maxmono)
        mean_chain_len_lst=[]
        for j in mc_run_lst:
            j.run(run_statis=True,run_config=False)
            # 输出结果
            mean_chain_len_lst.append((j.do_event.n_maxmono-j.do_event.n_unic_mono)/j.do_event.n_chains)
        axes[i].plot(temperature_lst, mean_chain_len_lst, marker='o', label='p_in='+str(p_in_lst[i]))
        axes[i].legend()
        axes[i].set_xlabel('Temperature (K)')
        axes[i].set_ylabel('Mean Chain Length')
    for i in range(num_plots, num_rows * num_cols):
        fig.delaxes(axes[i])
    fig.tight_layout()
    plt.show()

def bonds_to_T(temperature_lst,n_maxmono=1000): # 温度对成键的影响
    mc_run_lst=inst_mcrun_to_T(temperature_lst,n_maxmono=n_maxmono)
    AA_ratio=[]
    AB_ratio=[]
    BB_ratio=[]
    for i in mc_run_lst:
        i.run(run_statis=True,run_config=False)
        # 输出结果
        all_bonds=i.do_event.A_A+i.do_event.A_B+i.do_event.B_B+i.do_event.B_A
        AA_ratio.append((i.do_event.A_A/all_bonds)*100)
        AB_ratio.append(((i.do_event.A_B+i.do_event.B_A)/all_bonds)*100)
        BB_ratio.append((i.do_event.B_B/all_bonds)*100)
    plt.scatter(temperature_lst, AA_ratio, label='AtA', color='blue')
    plt.scatter(temperature_lst, AB_ratio, label='AtB', color='red')
    plt.scatter(temperature_lst, BB_ratio, label='BtB', color='green')
    plt.legend()
    plt.title('bonds ratio vs Temperature')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Count (%)')
    plt.grid(True)
    plt.show()
    
def bonds_to_p_in(temperature_lst,p_in_lst,n_maxmono=1000): # 温度对成键的影响与输入速率的关系
    num_plots = len(p_in_lst)
    num_rows = int(np.ceil(num_plots / 2))
    num_cols = 2 if num_plots > 1 else 1
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(10, 8))
    axes = axes.flatten() if num_plots > 1 else [axes]
    for i in range(len(p_in_lst)):
        mc_run_lst=inst_mcrun_to_T(temperature_lst,put_in=p_in_lst[i],n_maxmono=n_maxmono)
        AA_ratio=[]
        AB_ratio=[]
        BB_ratio=[]
        for j in mc_run_lst:
            j.run(run_statis=True,run_config=False)
            all_bonds=j.do_event.A_A+j.do_event.A_B+j.do_event.B_B+j.do_event.B_A
            AA_ratio.append((j.do_event.A_A/all_bonds)*100)
            AB_ratio.append(((j.do_event.A_B+j.do_event.B_A)/all_bonds)*100)
            BB_ratio.append((j.do_event.B_B/all_bonds)*100)
        axes[i].plot(temperature_lst, AA_ratio, label='AtA', color='blue', linestyle='-')
        axes[i].plot(temperature_lst, AB_ratio, label='AtB', color='red', linestyle='--')
        axes[i].plot(temperature_lst, BB_ratio, label='BtB', color='green', linestyle=':')
        axes[i].legend()
        axes[i].set_xlabel('Temperature (K)')
        axes[i].set_ylabel('Count (%)')
        axes[i].set_title('p_in='+str(p_in_lst[i]))
    for i in range(num_plots, num_rows * num_cols):
        fig.delaxes(axes[i])
    fig.tight_layout()
    plt.show()

def chain_len_dist_to_T(temperature_lst,n_maxmono=1000): # 温度对链长分布的影响
    mc_run_lst=inst_mcrun_to_T(temperature_lst,n_maxmono=n_maxmono)
    chain_len_dict={}
    for i in range(len(mc_run_lst)):
        mc_run_lst[i].run(run_statis=True,run_config=True)
        chain_len_dict=mc_run_lst[i].config.count_distribution() # 得到链长分布字典 e.g. {2: 2, 8: 1}
        # 输出结果
        plt.bar(chain_len_dict.keys(), chain_len_dict.values(), label='T='+str(temperature_lst[i]))
    plt.legend()
    plt.title('Chain length distribution vs Temperature')
    plt.xlabel('Chain length')
    plt.ylabel('Count')
    plt.show()

def chain_len_to_maxmono(p_in,temperature=[500],n_maxmono_lst=range(100,10001,100)): # 采样单体数对链长的影响
    mean_chain_len_lst=[]
    for i in n_maxmono_lst:
        mc_run_lst=inst_mcrun_to_T(temperature,put_in=p_in,n_maxmono=i)
        for j in mc_run_lst:
            j.run(run_statis=True,run_config=False)
            # 输出结果
            mean_chain_len_lst.append((j.do_event.n_maxmono-j.do_event.n_unic_mono)/j.do_event.n_chains)
    plt.plot(n_maxmono_lst, mean_chain_len_lst, marker='o')
    plt.xlabel('Number of monomer')
    plt.ylabel('Mean Chain Length')
    plt.title('Mean Chain Length vs Number of monomers')
    plt.grid(True)
    plt.show()

def chain_len_to_T_to_maxmono(p_in,n_maxmono_lst,temperature=range(500,601,20)): # 采样单体数对链长对温度的影响
    for i in n_maxmono_lst:
        mc_run_lst=inst_mcrun_to_T(temperature,put_in=p_in,n_maxmono=i)
        mean_chain_len_lst=[]
        for j in mc_run_lst:
            j.run(run_statis=True,run_config=False)
            # 输出结果
            mean_chain_len_lst.append((j.do_event.n_maxmono-j.do_event.n_unic_mono)/j.do_event.n_chains)
        plt.plot(temperature, mean_chain_len_lst, marker='o', label='n_maxmono='+str(i))
    plt.legend()    
    plt.xlabel('Temperature (K)')
    plt.ylabel('Mean Chain Length')
    plt.title('Mean Chain Length vs Temperature')
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    temperature_lst=range(300,601,50)
    # mean_chain_len_to_T(temperature_lst,n_maxmono=1000)
    # bonds_to_T(temperature_lst,n_maxmono=1000)
# def generate_logarithmic_sequence(start, end, num_points):
#     # 计算对数范围
#     log_start = np.log10(start)
#     log_end = np.log10(end)
    
#     # 生成对数序列
#     log_sequence = np.linspace(log_start, log_end, num_points)
    
#     # 将对数序列转换为原始数值序列
#     sequence = 10 ** log_sequence
    
#     return sequence
    # p_lst=generate_logarithmic_sequence(2e-8,2e-13,5)
    # chain_len_to_p_in(temperature_lst,[2e-8,2e-9,2e-10,2e-11,2e-12,2e-13],n_maxmono=1000)

    bonds_to_p_in(temperature_lst,[2e-8,2e-9,2e-10,2e-11,2e-12,2e-13],n_maxmono=1000)

    # temperature_lst=range(500,601,20)
    # chain_len_dist_to_T(temperature_lst,n_maxmono=1000)

    # chain_len_to_maxmono(p_in=2e-9,n_maxmono_lst=range(100,30001,200))

    # chain_len_to_T_to_maxmono(p_in=2e-11,n_maxmono_lst=range(1000,200001,10000))