##author: Wang Zicong  2024.12.7##
class do_event():
    def divide(self,a, b): # 保证分母不为0
        if b == 0:
            return 0
        else:
            return a / b
    def subt(self,a, b):
        if b == 0:
            return 0
        else:
            return a - b
        
    def __init__(self,n_maxmono,p_in,K_A,K_B,rA,rB,K_termin):
        # 读入常数
        self.p_in=p_in # p_in = 总的净输入速率
        self.n_maxmono=n_maxmono  # n_maxmono = 最大的单体数量
        self.K_termin=K_termin # K_termin = 终止速率常数
        self.K_A=K_A # K_A = A端发生增长的速率常数
        self.K_B=K_B # K_B = B端发生增长的速率常数
        self.AAtA=rA # rA = kAA/kA
        self.ABtA=1-rA 
        self.BBtB=rB # rB = kBB/kB
        self.BAtB=1-rB 
        # 初始化循环变量
        self.n_in=1 #总的输入单体数
        self.n_unic_mono=1 #总的独立单体数(没有加到链上的单体数)
        self.n_chain_act_A=0 #链上的总活性A数
        self.n_chain_act_B=0 #链上的总活性B数
        self.n_alive=self.n_chain_act_A+self.n_chain_act_B+2*self.n_unic_mono # 总的活性端数
        # 初始化统计变量
        self.n_chains=0 #链数
        self.A_A=0
        self.A_B=0
        self.B_B=0
        self.B_A=0
    def events_gener(self):
        # 加成
        self.P_chain_AaddA=self.n_chain_act_A*self.n_unic_mono*self.K_A*self.AAtA # 在链的A端加单体A
        self.P_chain_AaddB=self.n_chain_act_A*self.n_unic_mono*self.K_A*self.ABtA # 在链的A端加单体B
        self.P_chain_BaddB=self.n_chain_act_B*self.n_unic_mono*self.K_B*self.BBtB # 在链的B端加单体B
        self.P_chain_BaddA=self.n_chain_act_B*self.n_unic_mono*self.K_B*self.BAtB # 在链的B端加单体A
        self.P_monomer_AaddA=self.n_unic_mono*(self.n_unic_mono-1)*self.K_A*self.AAtA # 在单体A上加单体A
        self.P_monomer_AaddB=self.n_unic_mono*(self.n_unic_mono-1)*self.K_A*self.ABtA # 在单体A上加单体B
        self.P_monomer_BaddB=self.n_unic_mono*(self.n_unic_mono-1)*self.K_B*self.BBtB # 在单体B上加单体B
        self.P_monomer_BaddA=self.n_unic_mono*(self.n_unic_mono-1)*self.K_B*self.BAtB # 在单体B上加单体A
        # put_in
        # self.P_in=self.p_in*(1-self.n_in/self.n_maxmono) # 输入单体 (速率随剩余格点变化)
        self.P_in=self.p_in if self.n_in<self.n_maxmono else 0 # 输入单体 (速率不随剩余格点变化)
        # self.P_in=self.p_in *((self.n_maxmono-self.n_in)/self.n_maxmono)  # 输入单体 (速率随剩余格点变化且归一化)
        # 加H终止
        self.P_termin_chain_AH=self.n_chain_act_A*self.K_termin # 链A端被H终止
        self.P_termin_chain_BH=self.n_chain_act_B*self.K_termin # 链B端被H终止


        return [self.P_chain_AaddA,self.P_chain_AaddB,self.P_chain_BaddB,self.P_chain_BaddA,             # 0~3
                self.P_monomer_AaddA,self.P_monomer_AaddB,self.P_monomer_BaddB,self.P_monomer_BaddA,     # 4~7
                self.P_in,                                                                               # 8
                self.P_termin_chain_AH,self.P_termin_chain_BH,                                           # 9~10
                ]

    def do_event(self,event_idx): # 更新循环变量
        if event_idx==0: # 在链的A端加A
            self.n_unic_mono -= 1
            self.n_chain_act_A -= 1
            self.n_chain_act_B += 1
        elif event_idx==1: # 在链的A端加B
            self.n_unic_mono -= 1
        elif event_idx==2: # 在链的B端加B
            self.n_unic_mono -= 1
            self.n_chain_act_B -= 1
            self.n_chain_act_A += 1
        elif event_idx==3: # 在链的B端加A
            self.n_unic_mono -= 1
        elif event_idx==4: # 在单体A上加A
            self.n_unic_mono -= 2
            self.n_chain_act_B += 2
        elif event_idx==5: # 在单体A上加B
            self.n_unic_mono -= 2
            self.n_chain_act_A += 1
            self.n_chain_act_B += 1
        elif event_idx==6: # 在单体B上加B
            self.n_unic_mono -= 2
            self.n_chain_act_A += 2
        elif event_idx==7: # 在单体B上加A
            self.n_unic_mono -= 2
            self.n_chain_act_A += 1
            self.n_chain_act_B += 1
        elif event_idx==8: # 输入单体
            self.n_unic_mono += 1
            self.n_in += 1
        elif event_idx==9: # 链A端被H终止
            self.n_chain_act_A -= 1
        elif event_idx==10: # 链B端被H终止
            self.n_chain_act_B -= 1
        # 更新剩余总的活性端数
        self.n_alive=self.n_chain_act_A+self.n_chain_act_B+2*self.n_unic_mono 

    def do_statis(self,event_idx): # 更新统计变量
        if event_idx==0: # 在链的A端加A
            self.A_A+=1
        elif event_idx==1: # 在链的A端加B
            self.A_B+=1
        elif event_idx==2: # 在链的B端加B
            self.B_B+=1
        elif event_idx==3: # 在链的B端加A
            self.B_A+=1
        elif event_idx==4: # 在单体A上加A
            self.n_chains += 1
            self.A_A+=1
        elif event_idx==5: # 在单体A上加B
            self.n_chains += 1
            self.A_B+=1
        elif event_idx==6: # 在单体B上加B
            self.n_chains += 1
            self.B_B+=1
        elif event_idx==7: # 在单体B上加A
            self.n_chains += 1
            self.B_A+=1
        elif event_idx==8: # 输入单体
            pass
        elif event_idx==9: # 链A端被H终止
            pass
        elif event_idx==10: # 链B端被H终止
            pass
    
    def mc_go_on(self):
        # print('n_alive:',self.n_alive,'n_in::',self.n_in)
        if self.n_in==self.n_maxmono and (self.n_chain_act_A+self.n_chain_act_B)==0 and self.n_unic_mono <=1: # 当单体加入达到最大且没有链活性端剩余且单体最多剩1时才停止
            return False
        else:
            return True
        
    def mc_output(self):
        print('占位数：',self.n_in,'\n',
              '链数：',self.n_chains,'\n',
              '孤立单体数：',self.n_unic_mono,'\n',
              '总活性端数：',self.n_alive,'\n',
              '活性端A数：',self.n_chain_act_A,'\n',
              '活性端B数A：',self.n_chain_act_B,'\n',
              'AA数：',self.A_A,'\n',
              'AB数：',self.A_B+self.B_A,'\n',
              'BB数：',self.B_B,'\n')
    
    def testing_output(self):
        print('占位数：',self.n_in,'\n',
              '链数：',self.n_chains,'\n',
              '孤立单体数：',self.n_unic_mono,'\n',
              '总活性端数：',self.n_alive,'\n',
              '活性端A数：',self.n_chain_act_A,'\n',
              '活性端B数A：',self.n_chain_act_B,'\n',
              'AA数：',self.A_A,'\n',
              'AB数：',self.A_B+self.B_A,'\n',
              'BB数：',self.B_B,'\n')