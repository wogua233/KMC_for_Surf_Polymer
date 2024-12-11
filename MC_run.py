##author: Wang Zicong  2024.12.7##
import os
import random

from matplotlib import pyplot as plt
from MC_config import run_config
from MC_do_event import do_event
import math

def cleanxyz():
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".xyz"):
                os.remove(os.path.join(root, file))
                print(f"Deleted file: {os.path.join(root, file)}")

class mc_run():
    def __init__(self,n_maxmono,p_in,K_A,K_B,rA,rB,K_termin,seed=1):
        # n_maxmono = 最大的单体数量
        # p_in = 总的净输入速率
        # K_A = A端发生增长的速率常数
        # K_B = B端发生增长的速率常数
        # rA = kAA/kA
        # rB = kBB/kB
        # K_termin = 终止速率常数
        # seed = 随机数种子

        self.n_maxmono = n_maxmono
        self.p_in = p_in
        self.K_A = K_A
        self.K_B = K_B
        self.rA=rA
        self.rB=rB
        self.K_termin = K_termin
        self.seed = seed
        # 初始化
        random.seed(self.seed)
        # 实例化一个do_event类
        self.do_event=do_event(self.n_maxmono,self.p_in,self.K_A,self.K_B,self.rA,self.rB,self.K_termin)
        # 实例化一个run_config类
        self.config=run_config(seed=self.seed)

    def build_P(self,P=[]): # 创建一个事件集合
        sum_P=sum(P)
        for i in range(len(P)):
            P[i]=P[i]/sum_P
        new_P=[sum(P[:i+1]) for i in range(len(P))]
        return new_P

    def random_select(self,P): # 归一化事件概率
        w=random.random()
        for i in range(len(P)):
            if w<=P[i]:
                # print(w)
                return i


    def run(self,run_statis=True,run_config=False):
        print('MC模拟开始！')
        # mc步数：
        # n_mcstep=0
        while self.do_event.mc_go_on() :
            # 构建事件概率组：
            P_event=self.do_event.events_gener() # 创建一个事件集合
            P_event=self.build_P(P_event) # 归一化事件概率
            # 随机选择一个事件得到一个此事件的idx：
            event_idx=self.random_select(P_event)
            # print('event_idx:',event_idx)
            # 进行这个事件do_event：
            self.do_event.do_event(event_idx)
            if run_statis:
                self.do_event.do_statis(event_idx)
            if run_config:
                self.config.run_config(event_idx)
            # 测试时输出：
            # print('event_idx:',event_idx)
            # self.do_event.testing_output()
            # if n_mcstep%1000==0: # 每step_interval步输出一次快照
            #     print('n_mcstep:',n_mcstep)
            #     print('cover:',self.do_event.n_in/self.n_maxmono)
            #     print('孤立单体数：',self.do_event.n_unic_mono)
            #     print('链数：',self.do_event.n_alive)
            #     print('平均链长:',self.do_event.n_in/(self.do_event.n_chains+self.do_event.n_unic_mono))
            #     print('A_A占比:',self.do_event.A_A/(self.do_event.A_A+self.do_event.A_B+self.do_event.B_B+self.do_event.B_A) if (self.do_event.A_A+self.do_event.A_B+self.do_event.B_B+self.do_event.B_A)!=0 else 'none')
            #     print('\n')
            # n_mcstep+=1
        # 输出结果
        if run_statis:
            self.do_event.mc_output()
        if run_config:
            self.config.input_unic_mono(self.do_event.n_unic_mono)
            self.config.output()
        print('MC模拟结束！')

    def run_snapshoot(self,step_interval=100): # 输出快照
        print('MC模拟开始！')
        # mc步数：
        n_mcstep=0
        while self.do_event.mc_go_on() :
            # 构建事件概率组：
            P_event=self.do_event.events_gener() # 创建一个事件集合
            P_event=self.build_P(P_event) # 归一化事件概率
            # 随机选择一个事件得到一个此事件的idx：
            event_idx=self.random_select(P_event)
            # 进行这个事件do_event：
            self.do_event.do_event(event_idx)
            self.do_event.do_statis(event_idx)
            self.config.run_config(event_idx)
            n_mcstep+=1
            # 快照输出：
            if n_mcstep%step_interval==0: # 每step_interval步输出一次快照
                self.config.input_unic_mono(self.do_event.n_unic_mono)
                self.config.output_snap_config(n_mcstep=n_mcstep)
        # 输出结果
        self.do_event.mc_output()
        self.config.input_unic_mono(self.do_event.n_unic_mono)
        self.config.output_snap_config(n_mcstep)
        print('MC模拟结束！')
    
    def run_to_cover(self,run_statis=True,run_config=False,interval=1000):
        cover=[]
        n_unic_mono=[]
        p_AA=[]
        chain_alive_site=[]
        mean_len=[]
        print('MC模拟开始！')
        # mc步数：
        n_mcstep=0
        while self.do_event.mc_go_on() :
            # 构建事件概率组：
            P_event=self.do_event.events_gener() # 创建一个事件集合
            P_event=self.build_P(P_event) # 归一化事件概率
            # 随机选择一个事件得到一个此事件的idx：
            event_idx=self.random_select(P_event)
            # print('event_idx:',event_idx)
            # 进行这个事件do_event：
            self.do_event.do_event(event_idx)
            if run_statis:
                self.do_event.do_statis(event_idx)
            if run_config:
                self.config.run_config(event_idx)
            if n_mcstep%interval==0: # 每step_interval步输出一次快照
                # print('n_mcstep:',n_mcstep)
                # print('cover:',self.do_event.n_in/self.n_maxmono)
                cover.append(self.do_event.n_in/self.n_maxmono)
                # print('孤立单体数：',self.do_event.n_unic_mono)
                n_unic_mono.append(self.do_event.n_unic_mono)
                chain_alive_site.append(self.do_event.n_chain_act_A+self.do_event.n_chain_act_B)
                # print('链数：',self.do_event.n_alive)
                # print('平均链长:',self.do_event.n_in/(self.do_event.n_chains+self.do_event.n_unic_mono))
                mean_len.append(self.do_event.n_in/(self.do_event.n_chains+self.do_event.n_unic_mono))
                # print('A_A占比:',self.do_event.A_A/(self.do_event.A_A+self.do_event.A_B+self.do_event.B_B+self.do_event.B_A) if (self.do_event.A_A+self.do_event.A_B+self.do_event.B_B+self.do_event.B_A)!=0 else 'none')
                p_AA.append(self.do_event.A_A/(self.do_event.A_A+self.do_event.A_B+self.do_event.B_B+self.do_event.B_A) if (self.do_event.A_A+self.do_event.A_B+self.do_event.B_B+self.do_event.B_A)!=0 else 0)
                # print('\n')
            n_mcstep+=1
        # 输出结果
        max_n_unic_mono=max(n_unic_mono)
        n_unic_mono=[unic_mono/max_n_unic_mono for unic_mono in n_unic_mono]
        max_chain_alive_site=max(chain_alive_site)
        chain_alive_site=[chain/max_chain_alive_site for chain in chain_alive_site]
        max_mean_len=max(mean_len)
        mean_len=[mean/max_mean_len for mean in mean_len]

        print('MC模拟结束！')
        fig,(ax1,ax2)=plt.subplots(1,2,sharex=True)
        ax1.plot(cover, n_unic_mono, label='n unic mono')
        ax1.legend()
        ax1.set_ylabel('n unic mono')
        ax2.plot(cover, chain_alive_site, label='chain\'s alive sites')
        ax2.plot(cover, mean_len, label='mean chain length')
        ax2.plot(cover, p_AA, label='AtA ratio')
        ax2.set_xlabel('coverage')
        ax2.set_ylabel('normalized value') 
        ax2.legend()
        plt.show()











if __name__ == '__main__':
    # for T in range(400,401,50):
    #     T=T # K
    #     p_in=6e-8
    #     kAA=math.exp((-1.49E-19)/((1.38E-23)*T))
    #     kAB=math.exp((-1.74E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    #     kBB=math.exp((-1.90E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    #     kBA=math.exp((-1.74E-19)/((1.38E-23)*T))
    #     KA=kAA+kAB
    #     KB=kBB+kBA
    #     rA=kAA/KA
    #     rB=kBB/KB
    #     k_termin=math.exp((-1.20E-19)/((1.38E-23)*T))*0.000008 #设定终止速率常数
    #     # 实例化mc_run类
    #     mc_run1=mc_run(n_maxmono=50000,p_in=p_in,K_A=KA,K_B=KB,rA=rA,rB=rB,K_termin=k_termin,seed=1)
    #     # 运行mc_run类
    #     mc_run1.run(run_statis=True,run_config=False)
    
    
    # # 快照输出：
    # T=400
    # p_in=2e-10
    # kAA=math.exp((-1.49E-19)/((1.38E-23)*T))
    # kAB=math.exp((-1.74E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    # kBB=math.exp((-1.90E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    # kBA=math.exp((-1.74E-19)/((1.38E-23)*T))
    # KA=kAA+kAB
    # KB=kBB+kBA
    # rA=kAA/KA
    # rB=kBB/KB
    # k_termin=math.exp((-1.20E-19)/((1.38E-23)*T))*0.0000005 #设定终止速率常数
    # # 实例化mc_run类
    # mc_run1=mc_run(n_maxmono=10000,p_in=p_in,K_A=KA,K_B=KB,rA=rA,rB=rB,K_termin=k_termin,seed=1)
    # # 运行mc_run类
    # mc_run1.run_snapshoot(step_interval=1000)
    # # cleanxyz()
    
    # 对cover统计
    T=600
    p_in=6e-8
    kAA=math.exp((-1.49E-19)/((1.38E-23)*T))
    kAB=math.exp((-1.74E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    kBB=math.exp((-1.90E-19)/((1.38E-23)*T))*math.exp((-0.3E-19)/((1.38E-23)*T))
    kBA=math.exp((-1.74E-19)/((1.38E-23)*T))
    KA=kAA+kAB
    KB=kBB+kBA
    rA=kAA/KA
    rB=kBB/KB
    k_termin=math.exp((-1.20E-19)/((1.38E-23)*T))*0.000008 #设定终止速率常数
    # 实例化mc_run类
    mc_run1=mc_run(n_maxmono=10000,p_in=p_in,K_A=KA,K_B=KB,rA=rA,rB=rB,K_termin=k_termin,seed=1)
    # 运行mc_run类
    mc_run1.run_to_cover(run_statis=True,run_config=False,interval=100)

    
    
    
    
    
    
    
    
    
    
    