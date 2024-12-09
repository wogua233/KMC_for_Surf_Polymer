##author: Wang Zicong  2024.12.7##
import os
import random
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
        while self.do_event.mc_go_on() :
            # 构建事件概率组：
            P_event=self.do_event.events_gener() # 创建一个事件集合
            P_event=self.build_P(P_event) # 归一化事件概率
            # 随机选择一个事件得到一个此事件的idx：
            event_idx=self.random_select(P_event)
            # 进行这个事件do_event：
            self.do_event.do_event(event_idx)
            if run_statis:
                self.do_event.do_statis(event_idx)
            if run_config:
                self.config.run_config(event_idx)
            # 测试时输出：
            # print('event_idx:',event_idx)
            # self.do_event.testing_output()
        # 输出结果
        if run_statis:
            self.do_event.mc_output()
        if run_config:
            self.config.input_unic_mono(self.do_event.n_unic_mono)
            self.config.output()
        print('MC模拟结束！')

    def run_snapshoot(self,step_interval=10000): # 输出快照
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
                self.config.output_snap_config(n_mcstep=n_mcstep)
        # 输出结果
        self.do_event.mc_output()
        self.config.input_unic_mono(self.do_event.n_unic_mono)
        self.config.output_snap_config(n_mcstep)
        print('MC模拟结束！')
    











if __name__ == '__main__':
    for T in range(600,601,50):
        T=T # K
        p_in=2e-10
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
        mc_run1=mc_run(n_maxmono=1000,p_in=p_in,K_A=KA,K_B=KB,rA=rA,rB=rB,K_termin=k_termin,seed=1)
        # 运行mc_run类
        mc_run1.run(run_statis=True,run_config=True)
    # 快照输出：
    T=400
    p_in=2e-10
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
    mc_run1=mc_run(n_maxmono=1000,p_in=p_in,K_A=KA,K_B=KB,rA=rA,rB=rB,K_termin=k_termin,seed=1)
    # 运行mc_run类
    mc_run1.run_snapshoot(step_interval=100)
    # cleanxyz()
    

    
    
    
    
    
    
    
    
    
    
    