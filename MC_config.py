##author: Wang Zicong  2024.12.7##
import random
from config_to_xyz import to_xyz

class run_config():
    def __init__(self,seed):
        random.seed(seed)
        self.config=[] # 记录链的具体构型
        self.n_unic_mono=0 # 记录独立的单体数
    def random_index_2d(self,lst): # 二维随机索引
        # print(lst)
        lst_one_dim=[]
        for head_tail in range(len(lst)):
            for chain_idx in lst[head_tail]:
                lst_one_dim.append((head_tail,chain_idx))
        # print(lst_one_dim)
        return random.choice(lst_one_dim)
    
    def find_ends(self,chain_str):  # 更新：chain_str
        find_string_idx = [[],[]] # 找到的字符串的位置[[头],[尾]]
        for i in range(len(self.config)): # 生成找到字符串的位置的列表（头尾会重复录入）
            if self.config[i].startswith(chain_str):
                find_string_idx[0].append(i)
                # print(self.config[i],'start with',chain_str)
            if self.config[i].endswith(chain_str):
                find_string_idx[1].append(i)
                # print(self.config[i],'end with',chain_str)
        # print(find_string_idx)
        pick_idx=self.random_index_2d(find_string_idx) # 随机选择一个符合的链
        return pick_idx
    
    def add_mono(self,pick_idx,add_str): # 更新：chain_str加add_str
        if add_str=='A':
            end_str='B'
        elif add_str=='B':
            end_str='A'
        if pick_idx[0]==0: # 头
            self.config[pick_idx[1]]=end_str+add_str+'-'+self.config[pick_idx[1]]
        elif pick_idx[0]==1: # 尾
            self.config[pick_idx[1]]=self.config[pick_idx[1]]+'-'+add_str+end_str
        else:
            print('pick_idx错误！')
            return None
        
    def add_H(self,pick_idx):
        if pick_idx[0]==0: # 头
            self.config[pick_idx[1]]='H'+'-'+self.config[pick_idx[1]]
        elif pick_idx[0]==1: # 尾
            self.config[pick_idx[1]]=self.config[pick_idx[1]]+'-'+'H'
        else:
            print('pick_idx错误！')
            return None
        
    def run_config(self,event_idx):
        if event_idx==0: # 在链的A端加A
            pick_idx=self.find_ends('A')
            self.add_mono(pick_idx,'A')
        elif event_idx==1: # 在链的A端加B
            pick_idx=self.find_ends('A') 
            self.add_mono(pick_idx,'B')
        elif event_idx==2: # 在链的B端加B
            pick_idx=self.find_ends('B')
            self.add_mono(pick_idx,'B')
        elif event_idx==3: # 在链的B端加A
            pick_idx=self.find_ends('B') 
            self.add_mono(pick_idx,'A')
        elif event_idx==4: # 在单体A上加A
            self.config.append('BA-AB')
        elif event_idx==5: # 在单体A上加B
            self.config.append('AB-AB')
        elif event_idx==6: # 在单体B上加B
            self.config.append('AB-BA')
        elif event_idx==7: # 在单体B上加A
            self.config.append('BA-BA')
        elif event_idx==8: # 输入单体
            pass 
        elif event_idx==9: # 链A端被H终止
            pick_idx=self.find_ends('A') 
            # print('在',self.config[pick_idx[1]],'的',pick_idx[0],'A端被H终止')
            self.add_H(pick_idx)
        elif event_idx==10: # 链B端被H终止
            pick_idx=self.find_ends('B') 
            # print('在',self.config[pick_idx[1]],'的',pick_idx[0],'B端被H终止')
            self.add_H(pick_idx) 

    def input_unic_mono(self,n_unic_mono):
        self.n_unic_mono=n_unic_mono
    def config_to_xyz(self):
        TO1=to_xyz(self.config,n_unic_mono=self.n_unic_mono)
        TO1.config_split()
        TO1.bonds_to_atoms()
        TO1.to_xyz()
        TO1.output_xyz()

    def output(self): 
        self.config_to_xyz()

    def output_snap_config(self,n_mcstep):
        TO2=to_xyz(self.config,n_unic_mono=self.n_unic_mono)
        TO2.config_split()
        TO2.bonds_to_atoms()
        TO2.to_xyz()
        TO2.output_snapxyz(n_mcstep)


    def count_distribution(self,bin=5): # 统计链长分布
        TO2=to_xyz(self.config)
        TO2.config_split() #  -> TO2.bonds
        distri_lst=[]
        for chain in TO2.bonds:
            distri_lst.append(len(chain)-1) # 不算-H AND -H
        count_dict = {}
        for chain_lenth in distri_lst:
            if chain_lenth in count_dict:
                count_dict[chain_lenth] += 1
            else:
                count_dict[chain_lenth] = 1
        max_length = max(count_dict.keys())
        normalized_dict={}
        if max_length <= bin :
            normalized_dict[max_length]=sum(count_dict.values())
        else:
            for i in range(0,max_length+1,bin):
                normalized_dict[i+bin//2]=0
                for j in range(i,i+bin):
                    if j in count_dict:
                        normalized_dict[i+bin//2]+=count_dict[j]
        return normalized_dict


if __name__ == '__main__':
    run_config1=run_config(seed=2)
    run_config1.config=['H-AB-BA-H','H-AB-AB-H','H-BA-AB-AB-BA-AB-BA-BA-BA-H','H-AB-BA-H','H-AB-AB-H','H-BA-AB-AB-BA-AB-BA-BA-BA-H','H-AB-BA-H','H-AB-AB-H','H-BA-AB-AB-BA-AB-BA-BA-BA-H']
    # find_ends=run_config1.find_ends('H')
    # print(run_config1.config[find_ends[1]])
    print(run_config1.count_distribution())