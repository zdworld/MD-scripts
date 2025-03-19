#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#用于读取非周期性的聚合物拓扑，输出周期性分子拓扑
import pandas as pd
import re
from io import StringIO
import os


# In[ ]:



class GromacsTop:
    #一个用于gmx拓扑文件读取和分析的类，仅支持itp文件
    def __init__(self,**kwargs):
        if not 'itp_name' in kwargs.keys():
            if len(kwargs)!=0:
                print('Error:no itp file input')
                return
            else:
                return
        else:
            top_file=kwargs['itp_name']
            
        if not top_file.endswith('.itp'):
            print('Error:topology file is not itp file')
            exit()
        else:
            with open(top_file,'r') as f:
                lines=f.readlines()
            self.__top_lines=[]
            for line in lines:
                # print(line)
                if not line.startswith(';'):
                    self.__top_lines.append(line.strip())
        self.__top_lines=[line for line in self.__top_lines if not line=='']
        
        #对读取的文件进行项目拆分
        self.terms=[]
        self.__terms_line_counter=[]
        for i,line in enumerate(self.__top_lines):
            if re.match(r'\[.*\]',line):
                term=line.replace('[ ','').replace(' ]','')
                self.terms.append(term)
                self.__terms_line_counter.append(i)
        self.__terms_line_counter.append(len(self.__top_lines))
    
        temp=[]
        for line in self.__top_lines:
            if ';' in line:
                temp.append(line.rsplit(';',1)[0])
            else:
                temp.append(line)
        self.__top_lines=temp
        # print(self.__top_lines)
        
    def top2df(self)->dict:
        #将拓扑文件的每个部分都转换为dataframe,首先规定各列的columns
        topology=dict()
        #按序将各部分转换
        for i in range(len(self.terms)):
            sep='\n'
            top_seris=StringIO(sep.join(self.__top_lines[self.__terms_line_counter[i]+1:self.__terms_line_counter[i+1]]))
            df=pd.read_csv(top_seris,sep=r'\s+',header=None)
            topology.update({self.terms[i]:df})
        self.topology=topology

        for key in self.topology.keys():
            if key=='moleculetype':
                columuns=['Name','nrexcl']
            elif key=='atoms':
                columuns=['atom_id','atom_type','resid','resname','atom_name','cgnr','charge','mass','typeB','chargeB','massB']
            elif key=='bonds':
                columuns=['atom1','atom2','funct','c0','c1','c2','c3']
            elif key=='pairs':
                columuns=['atom1','atom2','funct','c0','c1','c2','c3']
            elif key=='angles':
                columuns=['atom1','atom2','atom3','funct','c0','c1','c2','c3']
            elif key=='dihedrals':
                columuns=['atom1','atom2','atom3','atom4','funct','c0','c1','c2','c3','c4','c5']
            elif key=='impropers':
                columuns=['atom1','atom2','atom3','atom4','funct','c0','c1','c2','c3','c4','c5']
            elif key=='constraints':
                columuns=['atom1','atom2','funct','c0','c1','c2','c3']
            elif key=='exclusions':
                columuns=['atom','excl1','excl2','excl3','excl4','excl5','excl6','excl7','excl8','excl9','excl10','excl11']
            elif 'virtual' in key:
                columuns=['vs','funct','def1','def2','def3','def4','def5','def6','def7','def8','def9','def10','def11']
                
            self.topology[key].columns=columuns[0:self.topology[key].shape[1]]
                
    def write_itp(self,file_name:str='new.itp')->None:
        if os.path.exists(file_name):
            os.remove(file_name)
        for key in self.topology.keys():
            with open(file_name,'a') as file:
                file.write(f'[ {key} ]\n;')
            self.topology[key].to_csv(path_or_buf=file_name,sep='\t',header=True,index=False,index_label=';',mode='a',float_format='%.6f')
            with open(file_name,'a') as file:
                file.write('\n')


# In[ ]:

file_name=r'CEC236.itp'     #非周期性拓扑输入
output_path = r'CEC236-inf.itp' #周期性拓扑输出


unit_plus1=GromacsTop(itp_name=file_name)
unit_plus1.top2df()
unit_plus1.topology.keys()

unit_length=79              #周期性链每个残基的长度
pbc_length=8                #输入的残基数
pbc_length=pbc_length+1


# In[ ]:


#截取片段，修改得到周期性拓扑
pbc_dict=dict()
#moleculetype
df=unit_plus1.topology['moleculetype'].copy(deep=True)
# df['Name']='inf-CEC4'
pbc_dict.update({'moleculetype':df})

#atoms
df=unit_plus1.topology['atoms'][(unit_plus1.topology['atoms'].resid>=3) & (unit_plus1.topology['atoms'].resid<=pbc_length+1)].copy(deep=True)   #选中到
df.index=range(0,df.shape[0])

for i in range(df.shape[0]):
    df.at[i,'atom_id']=df.at[i,'atom_id']-unit_length
    df.at[i,'resid']=df.at[i,'resid']-2
    df.at[i,'cgnr']=df.at[i,'cgnr']-unit_length
pbc_dict.update({'atoms':df})

#bonds
df=unit_plus1.topology['bonds'].copy(deep=True)
df.insert(df.shape[1],'conmment',';')
df.insert(df.shape[1],'connect_res','1-1')
for i in range(df.shape[0]):
    atom_1=df.at[i,'atom1']
    atom_2=df.at[i,'atom2']
    res1=unit_plus1.topology['atoms'].query(f'atom_id=={atom_1}')['resid'].iloc[0]
    res2=unit_plus1.topology['atoms'].query(f'atom_id=={atom_2}')['resid'].iloc[0]
    df.at[i,'connect_res']=f'{res1}-{res2}'

select_col=[]
for i in range(df.shape[0]):
    # 将字符串格式的 `connect_res` 转换为数值列表
    connect_res = list(map(int, df.at[i, 'connect_res'].split('-')))
    # 判定条件：connect_res 包含 2 到 pbc_length+1 的范围，且不包含 1
    if any([j in connect_res for j in range(3, pbc_length + 2)]) & (2 not in connect_res):
        select_col.append(i)

# 根据选中行索引创建新 DataFrame，并重新索引
df_sel = df.iloc[select_col, :].copy()
df_sel.index = range(df_sel.shape[0])

for i in range(df_sel.shape[0]):
    # 将 `connect_res` 转换为数值列表，便于判定
    connect = list(map(int, df_sel.at[i, 'connect_res'].split('-')))
    
    if pbc_length + 2 in connect:
        if pbc_length + 2 == connect[0]:  # 如果残基6在 atom1 中
            atom_replace_id = df_sel.at[i, 'atom1']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom1'] = atom_in_res2
        if pbc_length + 2 == connect[1]:  # 如果残基6在 atom2 中
            atom_replace_id = df_sel.at[i, 'atom2']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom2'] = atom_in_res2    

for i in range(1,3):
    df_sel.loc[:,f'atom{i}']=df_sel.loc[:,f'atom{i}']-unit_length
    
pbc_dict.update({'bonds':df_sel})

#pairs
df=unit_plus1.topology['pairs'].copy(deep=True)
df.insert(df.shape[1],'conmment',';')
df.insert(df.shape[1],'connect_res','1-1')
for i in range(df.shape[0]):
    atom_1=df.at[i,'atom1']
    atom_2=df.at[i,'atom2']
    res1=unit_plus1.topology['atoms'].query(f'atom_id=={atom_1}')['resid'].iloc[0]
    res2=unit_plus1.topology['atoms'].query(f'atom_id=={atom_2}')['resid'].iloc[0]
    df.at[i,'connect_res']=f'{res1}-{res2}'

select_col=[]
for i in range(df.shape[0]):
    # 将字符串格式的 `connect_res` 转换为数值列表
    connect_res = list(map(int, df.at[i, 'connect_res'].split('-')))
    
    # 判定条件：connect_res 包含 2 到 pbc_length 的范围，且不包含 1
    if any(j in connect_res for j in range(3, pbc_length + 2)) & (2 not in connect_res):
        select_col.append(i)

# 根据选中行索引创建新 DataFrame，并重新索引
df_sel = df.iloc[select_col, :].copy()
df_sel.index = range(df_sel.shape[0])

# 替换残基6中的原子为残基2中对应原子
for i in range(df_sel.shape[0]):
    # 将 `connect_res` 转换为数值列表，便于判定
    connect = list(map(int, df_sel.at[i, 'connect_res'].split('-')))
    
    if pbc_length + 2 in connect:
        if pbc_length + 2 == connect[0]:  # 如果残基6在 atom1 中
            atom_replace_id = df_sel.at[i, 'atom1']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom1'] = atom_in_res2
        if pbc_length + 2 == connect[1]:  # 如果残基6在 atom2 中
            atom_replace_id = df_sel.at[i, 'atom2']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom2'] = atom_in_res2   
for i in range(1,3):
    df_sel.loc[:,f'atom{i}']=df_sel.loc[:,f'atom{i}']-unit_length
pbc_dict.update({'pairs':df_sel})


# In[ ]:


#angles
df=unit_plus1.topology['angles'].copy(deep=True)
df.insert(df.shape[1],'conmment',';')
df.insert(df.shape[1],'connect_res','1-1-1')
for i in range(df.shape[0]):
    atom_1=df.at[i,'atom1']
    atom_2=df.at[i,'atom2']
    atom_3=df.at[i,'atom3']
    res1=unit_plus1.topology['atoms'].query(f'atom_id=={atom_1}')['resid'].iloc[0]
    res2=unit_plus1.topology['atoms'].query(f'atom_id=={atom_2}')['resid'].iloc[0]
    res3=unit_plus1.topology['atoms'].query(f'atom_id=={atom_3}')['resid'].iloc[0]
    df.at[i,'connect_res']=f'{res1}-{res2}-{res3}'
select_col=[]
for i in range(df.shape[0]):
    # 将字符串格式的 `connect_res` 转换为数值列表
    connect_res = list(map(int, df.at[i, 'connect_res'].split('-')))
    
    # 判定条件：connect_res 包含 2 到 pbc_length 的范围，且不包含 1
    if any(j in connect_res for j in range(3, pbc_length + 2)) & (2 not in connect_res):
        select_col.append(i)

# 根据选中行索引创建新 DataFrame，并重新索引
df_sel = df.iloc[select_col, :].copy()
df_sel.index = range(df_sel.shape[0])

# 替换残基6中的原子为残基2中对应原子
for i in range(df_sel.shape[0]):
    # 将 `connect_res` 转换为数值列表，便于判定
    connect = list(map(int, df_sel.at[i, 'connect_res'].split('-')))
    
    if pbc_length + 2 in connect:
        if pbc_length + 2 == connect[0]:  # 如果残基6在 atom1 中
            atom_replace_id = df_sel.at[i, 'atom1']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom1'] = atom_in_res2
        if pbc_length + 2 == connect[1]:  # 如果残基6在 atom2 中
            atom_replace_id = df_sel.at[i, 'atom2']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom2'] = atom_in_res2 
        if pbc_length + 2 == connect[2]:  # 如果残基6在 atom3 中
            atom_replace_id = df_sel.at[i, 'atom3']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom3'] = atom_in_res2  
for i in range(1,4):
    df_sel.loc[:,f'atom{i}']=df_sel.loc[:,f'atom{i}']-unit_length
pbc_dict.update({'angles':df_sel})


# In[ ]:


#dihedrals
df=unit_plus1.topology['dihedrals'].copy(deep=True)
df.insert(df.shape[1],'conmment',';')
df.insert(df.shape[1],'connect_res','1-1-1-1')
for i in range(df.shape[0]):
    atom_1=df.at[i,'atom1']
    atom_2=df.at[i,'atom2']
    atom_3=df.at[i,'atom3']
    atom_4=df.at[i,'atom4']
    res1=unit_plus1.topology['atoms'].query(f'atom_id=={atom_1}')['resid'].iloc[0]
    res2=unit_plus1.topology['atoms'].query(f'atom_id=={atom_2}')['resid'].iloc[0]
    res3=unit_plus1.topology['atoms'].query(f'atom_id=={atom_3}')['resid'].iloc[0]
    res4=unit_plus1.topology['atoms'].query(f'atom_id=={atom_4}')['resid'].iloc[0]
    df.at[i,'connect_res']=f'{res1}-{res2}-{res3}-{res4}'
select_col=[]
for i in range(df.shape[0]):
    # 将字符串格式的 `connect_res` 转换为数值列表
    connect_res = list(map(int, df.at[i, 'connect_res'].split('-')))
    
    # 判定条件：connect_res 包含 2 到 pbc_length 的范围，且不包含 1
    if any(j in connect_res for j in range(3, pbc_length + 2)) & (2 not in connect_res):
        select_col.append(i)

# 根据选中行索引创建新 DataFrame，并重新索引
df_sel = df.iloc[select_col, :].copy()
df_sel.index = range(df_sel.shape[0])

# 替换残基6中的原子为残基2中对应原子
for i in range(df_sel.shape[0]):
    # 将 `connect_res` 转换为数值列表，便于判定
    connect = list(map(int, df_sel.at[i, 'connect_res'].split('-')))
    
    if pbc_length + 2 in connect:
        if pbc_length + 2 == connect[0]:  # 如果残基6在 atom1 中
            atom_replace_id = df_sel.at[i, 'atom1']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom1'] = atom_in_res2
        if pbc_length + 2 == connect[1]:  # 如果残基6在 atom2 中
            atom_replace_id = df_sel.at[i, 'atom2']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom2'] = atom_in_res2
        if pbc_length + 2 == connect[2]:  # 如果残基6在 atom3 中
            atom_replace_id = df_sel.at[i, 'atom3']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom3'] = atom_in_res2  
        if pbc_length + 2 == connect[3]:  # 如果残基6在 atom4 中
            atom_replace_id = df_sel.at[i, 'atom4']
            atom_name = unit_plus1.topology['atoms'].query(f"atom_id=={atom_replace_id}")['atom_name'].iloc[0]
            atom_in_res2 = unit_plus1.topology['atoms'].query(f"atom_name=='{atom_name}' & resid==3")['atom_id'].iloc[0]
            df_sel.at[i, 'atom4'] = atom_in_res2
for i in range(1,5):
    df_sel.loc[:,f'atom{i}']=df_sel.loc[:,f'atom{i}']-unit_length
pbc_dict.update({'dihedrals':df_sel})


# In[ ]:


pbc_top=GromacsTop()
pbc_top.topology=pbc_dict
pbc_top.write_itp(output_path)

