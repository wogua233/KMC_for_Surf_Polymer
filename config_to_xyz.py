##author: Wang Zicong  2024.12.7##
import pandas as pd


class to_xyz():
    def __init__(self,config,n_unic_mono=0,x_shift=0.6,y_shift=0.3):
        self.config=config
        self.bonds=[]
        self.atoms=[]
        self.elements=[]
        self.x=[]
        self.y=[]
        self.z=[]
        self.x_shift=x_shift
        self.y_shift=y_shift
        self.n_unic_mono=n_unic_mono
    def config_split(self): # 将config中的链拆分成键
        for i in range(len(self.config)):
            self.bonds.append([])
            for j in range(len(self.config[i])):
                if self.config[i][j]=='-':
                    self.bonds[i].append(self.config[i][j-1:j+2])
    def bonds_to_atoms(self):
        for i in range(len(self.bonds)): # 每一条链
            self.atoms.append([]) 
            for j in self.bonds[i]: # 每一条链的键
                if j[0]=='H' or j[-1]=='H': # 如果是H键
                    self.atoms[i].append('H')
                elif j[0]=='A' and j[-1]=='A': # 如果是A-A键
                    self.atoms[i].append('C')
                elif j[0]=='B' and j[-1]=='B': # 如果是B-B键
                    self.atoms[i].append('N')
                else: # 如果是A-B键
                    self.atoms[i].append('O')
    def to_xyz(self):
        for i in range(len(self.atoms)): # 每一条链
            for j in range(len(self.atoms[i])): # 每一条链的键
                self.x.append(i*self.x_shift)
                self.y.append(j*self.y_shift)
                self.z.append(0.)
                if self.atoms[i][j]=='H': # 如果是H键
                    self.elements.append('H')
                elif self.atoms[i][j]=='C': # 如果是C键
                    self.elements.append('C')
                elif self.atoms[i][j]=='N': # 如果是N键
                    self.elements.append('N')
                else: # 如果是O键
                    self.elements.append('O')
        if self.n_unic_mono>0: # 如果有单体
            for i in range(self.n_unic_mono): # 每一个单体
                self.x.append(len(self.atoms)*self.x_shift+i*self.x_shift)
                self.y.append(0.)
                self.z.append(0.)
                self.elements.append('Na')
    def output_xyz(self,filename='config'):
        data={'element':self.elements,'x':self.x,'y':self.y,'z':self.z}
        df=pd.DataFrame(data)
        length=len(self.elements)
        with open(f'{filename}.xyz', 'w') as f:
            description = f"{length}\n\n"
            f.write(description)
            df.to_csv(f,sep='\t',header=False,index=False,lineterminator='\n')
            f.close()

    def output_snapxyz(self,filename='snap_shoots',n_mcstep=0):
        data={'element':self.elements,'x':self.x,'y':self.y,'z':self.z}
        df=pd.DataFrame(data)
        length=len(self.elements)
        with open(f'{filename}_{n_mcstep}.xyz', 'w') as f:
            description = f"{length}\n\n"
            f.write(description)
            df.to_csv(f,sep='\t',header=False,index=False,lineterminator='\n')
            f.close()

    
        
        

if __name__=='__main__':
    config=['H-AB-BA-H','H-AB-H','H-BA-AB-AB-BA-AB-BA-BA-BA-H']
    TO1=to_xyz(config)
    TO1.config_split()
    print(TO1.bonds)
    TO1.bonds_to_atoms()
    print(TO1.atoms)
    TO1.to_xyz()
    TO1.output_xyz()