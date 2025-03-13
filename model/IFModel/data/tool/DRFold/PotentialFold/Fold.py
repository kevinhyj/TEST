#! /nfs/amino-home/liyangum/miniconda3/bin/python
import torch
torch.manual_seed(6) # my lucky number
import torch.autograd as autograd
import numpy as np 
np.random.seed(9) # my unlucky number
import random
random.seed(9)
import Cubic,Potential
import operations
import os,json,sys

import a2b,rigid
import torch.optim as opt
from torch.nn.parameter import Parameter
import torch.nn as nn
import math
from scipy.optimize import fmin_l_bfgs_b,fmin_cg,fmin_bfgs
from scipy.optimize import minimize
import lbfgs_rosetta

Scale_factor=1.0
USEGEO=False

def readconfig(configfile=''):
    config=[]
    expdir=os.path.dirname(os.path.abspath(__file__))
    if configfile=='':
        configfile=os.path.join(expdir,'lib','ddf.json')
    config=json.load(open(configfile,'r'))
    return config 

    
class Structure:
    def __init__(self,fastafile,geofile,saveprefix,coornpys):
        self.config=readconfig()
        self.seqfile=fastafile
        # try:
        #     self.geos = np.load(geofile,allow_pickle=True).item()
        # excpet:
        #    pass
        self.geo = np.load(geofile,allow_pickle=True).item()
        self.saveprefix=saveprefix
        self.seq=open(fastafile).readlines()[1].strip()
        self.L=len(self.seq)
        basenpy = np.load( os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib','base.npy'  )  )
        self.basex = operations.Get_base(self.seq,basenpy)
        othernpy = np.load( os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib','other2.npy'  )  )
        self.otherx = operations.Get_base(self.seq,othernpy)
        sidenpy = np.load( os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib','side.npy'  )  )
        self.sidex = operations.Get_base(self.seq,sidenpy)        
        self.txs=[]
        for coorfile in coornpys:
            self.txs.append( torch.from_numpy(np.load(coorfile)).double()     )
        
        self.init_mask()
        self.init_paras()
        self._init_fape()
    def _init_fape(self):
        self.tx2ds=[]
        for tx in self.txs:
            true_rot,true_trans   = operations.Kabsch_rigid(self.basex,tx[:,0],tx[:,1],tx[:,2])
            true_x2 = tx[:,None,:,:] - true_trans[None,:,None,:]
            true_x2 = torch.einsum('ijnd,jde->ijne',true_x2,true_rot.transpose(-1,-2))
            self.tx2ds.append(true_x2)
    def init_mask(self):
        halfmask=np.zeros([self.L,self.L])
        fullmask=np.zeros([self.L,self.L])
        for i in range(self.L):
            for j in range(i+1,self.L):
                halfmask[i,j]=1
                fullmask[i,j]=1
                fullmask[j,i]=1
        self.halfmask=torch.DoubleTensor(halfmask) > 0.5
        self.fullmask=torch.DoubleTensor(fullmask) > 0.5
        self.clash_mask = torch.zeros([self.L,self.L,22,22])
        for i in range(self.L):
            for j in range(i+1,self.L):
                self.clash_mask[i,j]=1
        # for i in range(self.L-1):
        #     self.clash_mask[i,i+1,5,0]=0
        #     self.clash_mask[i,i+1,0,5]=0
        for i in range(self.L):
             self.clash_mask[i,i,:6,7:]=1

        for i in range(self.L-1):
            self.clash_mask[i,i+1,:,0]=0
            self.clash_mask[i,i+1,0,:]=0
            self.clash_mask[i,i+1,:,5]=0
            self.clash_mask[i,i+1,5,:]=0

        self.side_mask = rigid.side_mask(self.seq)
        self.side_mask = self.side_mask[:,None,:,None] * self.side_mask[None,:,None,:]
        self.clash_mask = (self.clash_mask > 0.5) * (self.side_mask > 0.5)
        self.confimask_cc = torch.DoubleTensor(self.geo['cc'][:,:,-1]) < 0.5
        self.confimask_pp = torch.DoubleTensor(self.geo['pp'][:,:,-1]) < 0.5
        self.confimask_nn = torch.DoubleTensor(self.geo['nn'][:,:,-1]) < 0.5

        self.confimask_pccp =  torch.DoubleTensor(self.geo['pccp'][:,:,-1]) < 0.5
        self.dynamic_pccp   =  self.confimask_pccp*self.halfmask
        self.pccpi,self.pccpj = torch.where(self.dynamic_pccp > 0.5)
        self.dynamic_pccp_np=self.dynamic_pccp.numpy()

        self.confimask_cnnc =  torch.DoubleTensor(self.geo['cnnc'][:,:,-1]) < 0.5
        self.dynamic_cnnc   =  self.confimask_cnnc*self.halfmask
        self.cnnci,self.cnncj = torch.where(self.dynamic_cnnc > 0.5)
        self.dynamic_cnnc_np=self.dynamic_cnnc.numpy()

        self.confimask_pnnp=  torch.DoubleTensor(self.geo['pnnp'][:,:,-1]) < 0.5
        self.dynamic_pnnp   =  self.confimask_pnnp*self.halfmask
        self.pnnpi,self.pnnpj = torch.where(self.dynamic_pnnp > 0.5)
        self.dynamic_pnnp_np=self.dynamic_pnnp.numpy()


    def init_paras(self):
        self.cc_cs,self.cc_decs=Cubic.dis_cubic(self.geo['cc'],2,24,44)
        self.pp_cs,self.pp_decs=Cubic.dis_cubic(self.geo['pp'],2,30,56)
        self.nn_cs,self.nn_decs=Cubic.dis_cubic(self.geo['nn'],2,18,32)

        self.pccp_cs,self.pccp_decs=Cubic.torsion_cubic(self.geo['pccp'],-math.pi,math.pi,36)
        self.pccp_cs,self.pccp_decs=self.pccp_cs[self.dynamic_pccp_np],self.pccp_decs[self.dynamic_pccp_np]
        self.pccp_coe=torch.DoubleTensor(np.array([acs.c for acs in self.pccp_cs]))
        self.pccp_x = torch.DoubleTensor(np.array([acs.x for acs in self.pccp_cs]))

        self.cnnc_cs,self.cnnc_decs=Cubic.torsion_cubic(self.geo['cnnc'],-math.pi,math.pi,36)
        self.cnnc_cs,self.cnnc_decs=self.cnnc_cs[self.dynamic_cnnc_np],self.cnnc_decs[self.dynamic_cnnc_np]
        self.cnnc_coe=torch.DoubleTensor(np.array([acs.c for acs in self.cnnc_cs]))
        self.cnnc_x = torch.DoubleTensor(np.array([acs.x for acs in self.cnnc_cs]))


        self.pnnp_cs,self.pnnp_decs=Cubic.torsion_cubic(self.geo['pnnp'],-math.pi,math.pi,36)
        self.pnnp_cs,self.pnnp_decs=self.pnnp_cs[self.dynamic_pnnp_np],self.pnnp_decs[self.dynamic_pnnp_np]
        self.pnnp_coe=torch.DoubleTensor(np.array([acs.c for acs in self.pnnp_cs]))
        self.pnnp_x = torch.DoubleTensor(np.array([acs.x for acs in self.pnnp_cs]))
     

    def init_quat(self,ii):
        x = torch.rand([self.L,21])
        x[:,18:] = self.txs[ii].mean(dim=1)
        init_coor = self.txs[ii] 
        biasq=torch.mean(init_coor,dim=1,keepdim=True)
        q=init_coor-biasq
        m = torch.einsum('bnz,bny->bzy',self.basex,q).reshape([self.L,-1])

        x[:,:9] = x[:,9:18] = m
        x.requires_grad=True
        return x
        
    def init_quat_safe(self,ii):
        x = torch.rand([self.L,21])
        x[:,18:] = self.txs[ii].mean(dim=1)
        init_coor = self.txs[ii] 
        biasq=torch.mean(init_coor,dim=1,keepdim=True)
        q=init_coor-biasq + torch.rand([self.L,3,3])
        m = (torch.einsum('bnz,bny->bzy',self.basex,q) + torch.eye(3)[None,:,:])  .reshape([self.L,-1])

        x[:,:9] = x[:,9:18] = m
        x.requires_grad=True
        return x

    def compute_bb_clash(self,coor,other_coor):
        com_coor = torch.cat([coor,other_coor],dim=1)
        com_dis  = (com_coor[:,None,:,None,:] - com_coor[None,:,None,:,:]).norm(dim=-1)
        dynamicmask2_vdw= (com_dis <= 3.15) * (self.clash_mask)
        #vdw_dynamic=torch.nn.functional.softplus(3.15-com_dis[dynamicmask2_vdw])
        vdw_dynamic = Potential.LJpotential(com_dis[dynamicmask2_vdw],3.15)
        return vdw_dynamic.sum()*self.config['weight_vdw']

    def compute_full_clash(self,coor,other_coor,side_coor):
        com_coor = torch.cat([coor[:,:2],other_coor,side_coor],dim=1)
        com_dis  = (com_coor[:,None,:,None,:] - com_coor[None,:,None,:,:]).norm(dim=-1)
        dynamicmask2_vdw= (com_dis <= 2.5) * (self.clash_mask)
        #vdw_dynamic=torch.nn.functional.softplus(3.15-com_dis[dynamicmask2_vdw])
        vdw_dynamic = Potential.LJpotential(com_dis[dynamicmask2_vdw],2.5)
        return vdw_dynamic.sum()*self.config['weight_vdw']




    def compute_cc_energy(self,coor):
        min_dis,max_dis,bin_num = 2,24,44
        c_atoms=coor[:,1]
        upper_th=max_dis - ((max_dis-min_dis)/bin_num)*0.5
        lower_th=3.10
        cc_map=operations.pair_distance(c_atoms,c_atoms)
        dynamicmask_cc= (cc_map <= upper_th) * (self.confimask_cc) * (self.fullmask) * (cc_map >= 2.5)
        dynamicmask_cc_np=dynamicmask_cc.numpy()
        if dynamicmask_cc_np.sum()>1:
            self.cc_coe=torch.DoubleTensor(np.array([acs.c for acs in self.cc_cs[dynamicmask_cc_np]]))
            self.cc_x = torch.DoubleTensor(np.array([acs.x for acs in self.cc_cs[dynamicmask_cc_np]]))
            E_cb = Potential.cubic_distance(cc_map[dynamicmask_cc],self.cc_coe,self.cc_x,min_dis,max_dis,bin_num).sum()*self.config['weight_cc']*0.5
        else:
            E_cb =0 
        E_cb = E_cb + ((     (cc_map <= 2.5)* (self.fullmask) * (self.confimask_cc)      ).sum() * 5)*self.config['weight_cc']
        return E_cb
    def compute_pp_energy(self,coor):
        min_dis,max_dis,bin_num = 2,30,56
        p_atoms=coor[:,0]
        upper_th=max_dis - ((max_dis-min_dis)/bin_num)*0.5
        lower_th=3.10
        pp_map=operations.pair_distance(p_atoms,p_atoms)
        dynamicmask_pp= (pp_map <= upper_th) * (self.confimask_pp) * (self.fullmask) * (pp_map >= 2.5)
        dynamicmask_pp_np=dynamicmask_pp.numpy()
        if dynamicmask_pp_np.sum()>1:
            self.pp_coe=torch.DoubleTensor(np.array([acs.c for acs in self.pp_cs[dynamicmask_pp_np]]))
            self.pp_x = torch.DoubleTensor(np.array([acs.x for acs in self.pp_cs[dynamicmask_pp_np]]))
            E_cb = Potential.cubic_distance(pp_map[dynamicmask_pp],self.pp_coe,self.pp_x,min_dis,max_dis,bin_num).sum()*self.config['weight_pp']*0.5
        else:
            E_cb =0 
        E_cb = E_cb + ((     (pp_map <= 2.5)* (self.fullmask) * (self.confimask_pp)      ).sum() * 5)*self.config['weight_pp']
        return E_cb
    def compute_nn_energy(self,coor):
        min_dis,max_dis,bin_num = 2,18,32
        n_atoms=coor[:,-1]
        upper_th=max_dis - ((max_dis-min_dis)/bin_num)*0.5
        lower_th=3.10
        nn_map=operations.pair_distance(n_atoms,n_atoms)
        dynamicmask_nn= (nn_map <= upper_th) * (self.confimask_nn) * (self.fullmask) * (nn_map >= 2.5)

        dynamicmask_nn_np=dynamicmask_nn.numpy()
        if dynamicmask_nn_np.sum()>1:
            self.nn_coe=torch.DoubleTensor(np.array([acs.c for acs in self.nn_cs[dynamicmask_nn_np]]))
            self.nn_x = torch.DoubleTensor(np.array([acs.x for acs in self.nn_cs[dynamicmask_nn_np]]))
            E_cb = Potential.cubic_distance(nn_map[dynamicmask_nn],self.nn_coe,self.nn_x,min_dis,max_dis,bin_num).sum()*self.config['weight_nn']*0.5
        else:
            E_cb =0 
        E_cb = E_cb + ((     (nn_map <= 2.5)* (self.fullmask) * (self.confimask_nn)      ).sum() * 5)*self.config['weight_nn']
        return E_cb
    def compute_pccp_energy(self,coor):
        p_atoms=coor[:,0]
        c_atoms=coor[:,1]
        pccpmap=operations.dihedral( p_atoms[self.pccpi], c_atoms[self.pccpi], c_atoms[self.pccpj] ,p_atoms[self.pccpj]                  )
        neg_log = Potential.cubic_torsion(pccpmap,self.pccp_coe,self.pccp_x,36)
        return neg_log.sum()*self.config['weight_pccp']

    def compute_cnnc_energy(self,coor):
        n_atoms=coor[:,-1]
        c_atoms=coor[:,1]
        pccpmap=operations.dihedral( c_atoms[self.cnnci], n_atoms[self.cnnci], n_atoms[self.cnncj] ,c_atoms[self.cnncj]                  )
        neg_log = Potential.cubic_torsion(pccpmap,self.cnnc_coe,self.cnnc_x,36)
        return neg_log.sum()*self.config['weight_cnnc']

    def compute_pnnp_energy(self,coor):
        n_atoms=coor[:,-1]
        p_atoms=coor[:,0]
        pccpmap=operations.dihedral( p_atoms[self.pnnpi], n_atoms[self.pnnpi], n_atoms[self.pnnpj] ,p_atoms[self.pnnpj]                  )
        neg_log = Potential.cubic_torsion(pccpmap,self.pnnp_coe,self.pnnp_x,36)
        return neg_log.sum()*self.config['weight_pnnp']

    def compute_pcc_energy(self,coor):
        p_atoms=coor[:,1]
        c_atoms=coor[:,2]
        pccmap=operations.angle( p_atoms[self.pcci], c_atoms[self.pcci], c_atoms[self.pccj]                   )
        neg_log = Potential.cubic_angle(pccmap,self.pcc_coe,self.pcc_x,12)
        return neg_log.sum()*self.config['weight_pcc']

    def compute_fape_energy(self,coor,ep=1e-3,epmax=20):
        energy= 0
        for tx in self.tx2ds:
            px_mean = coor[:,[1]]
            p_rot   = operations.rigidFrom3Points(coor)
            p_tran  = px_mean[:,0]
            pred_x2 = coor[:,None,:,:] - p_tran[None,:,None,:] # Lx Lrot N , 3
            pred_x2 = torch.einsum('ijnd,jde->ijne',pred_x2,p_rot.transpose(-1,-2)) # transpose should be equal to inverse
            errmap=torch.sqrt( ((pred_x2 - tx)**2).sum(dim=-1) + ep )
            energy = energy + torch.sum(  torch.clamp(errmap,max=epmax)        )
        return energy * self.config['weight_fape']

    def compute_bond_energy(self,coor,other_coor):
        # 3.87
        o3 = other_coor[:-1,-2]
        p  = coor[1:,0]
        dis = (o3-p).norm(dim=-1)
        energy = ((dis-1.607)**2).sum()
        return energy * self.config['weight_bond']


    def compute_fape_energy_fromquat(self,x,coor,ep=1e-3,epmax=100):
        energy= 0
        p_rot,px_mean = a2b.Non2rot(x[:,:9],x.shape[0]),x[:,9:]
        pred_x2 = coor[:,None,:,:] - px_mean[None,:,None,:] # Lx Lrot N , 3
        pred_x2 = torch.einsum('ijnd,jde->ijne',pred_x2,p_rot.transpose(-1,-2)) # transpose should be equal to inverse
        #coor  = a2b.quat2b(x)
        for tx in self.tx2ds:
            # px_mean = coor[:,[1]]
            # p_rot   = operations.rigidFrom3Points(coor)
            # p_tran  = px_mean[:,0]

            errmap=torch.sqrt( ((pred_x2 - tx)**2).sum(dim=-1) + ep )
            energy = energy + torch.sum(  torch.clamp(errmap,max=epmax)        )
            #energy = energy + torch.sum(  errmap       )

        return energy * self.config['weight_fape']
    def energy(self,rama):
        #rama=torch.cat([rama,self.betas],dim=-1)
        
        coor=a2b.quat2b(self.basex,rama[:,9:])
        other_coor = a2b.quat2b(self.otherx,rama[:,9:])
        side_coor = a2b.quat2b(self.sidex,torch.cat([rama[:,:9],coor[:,-1]],dim=-1))

        #print(coor.shape,other_coor.shape,side_coor.shape)
        
        if self.config['weight_cc']>0:
            E_cc= self.compute_cc_energy(coor)
        else:
            E_cc=0
        if self.config['weight_pp']>0:
            E_pp= self.compute_pp_energy(coor)
        else:
            E_pp=0
        if self.config['weight_nn']>0:
            E_nn= self.compute_nn_energy(coor)
        else:
            E_nn=0

        if self.config['weight_pccp']>0:
            E_pccp= self.compute_pccp_energy(coor)
        else:
            E_pccp=0

        if self.config['weight_cnnc']>0:
            E_cnnc= self.compute_cnnc_energy(coor)
        else:
            E_cnnc=0

        if self.config['weight_pnnp']>0:
            E_pnnp= self.compute_pnnp_energy(coor)
        else:
            E_pnnp=0

        if self.config['weight_vdw']>0:
            E_vdw= self.compute_full_clash(coor,other_coor,side_coor)
        else:
            E_vdw=0

        if self.config['weight_fape']>0:
            E_fape= self.compute_fape_energy_fromquat(rama[:,9:],coor)
        else:
            E_fape=0
        if self.config['weight_bond']>0:
            E_bond= self.compute_bond_energy(coor,other_coor)
        else:
            E_bond=0
        return  E_vdw + E_fape + E_bond + E_pp + E_cc + E_nn + E_pccp + E_cnnc + E_pnnp


    def obj_func_grad_np(self,rama_):
        rama=torch.DoubleTensor(rama_)
        rama.requires_grad=True
        if rama.grad:
            rama.grad.zero_()
        f=self.energy(rama.view(self.L,21))*Scale_factor
        grad_value=autograd.grad(f,rama)[0]
        return grad_value.data.numpy().astype(np.float64)
    def obj_func_np(self,rama_):
        rama=torch.DoubleTensor(rama_)
        rama=rama.view(self.L,21)
        with torch.no_grad():
            f=self.energy(rama)*Scale_factor
            #print('score',f)
            return f.item()


    def foldning(self):
        minenergy=1e16
        count=0
        for tx in self.txs:
            #self.outpdb_coor(tx,self.saveprefix+f'.{str(count)}.pdb',energystr=str(0))
            count+=1
        minirama=None
        for ilter in range( len(self.txs)):
            try:
                rama=self.init_quat(ilter).data.numpy()
                self.outpdb(torch.DoubleTensor(rama) ,self.saveprefix+f'_{str(ilter)}'+'.pdb',energystr='')
                self.config=readconfig(os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib','vdw.json'))
                rama = fmin_l_bfgs_b(func=self.obj_func_np, x0=rama,  fprime=self.obj_func_grad_np,iprint=10)[0]
                rama = rama.flatten()
            except:
                rama=self.init_quat_safe(ilter).data.numpy()
                #self.outpdb(torch.DoubleTensor(rama) ,self.saveprefix+f'_init'+str(ilter)+'.pdb',energystr='')
                self.config=readconfig(os.path.join(os.path.dirname(os.path.abspath(__file__)),'lib','vdw.json'))
                rama = fmin_l_bfgs_b(func=self.obj_func_np, x0=rama,  fprime=self.obj_func_grad_np,iprint=10)[0]
                rama = rama.flatten()
                
            self.config=readconfig()
            self.config['weight_pp'] =1000 * self.config['weight_pp']
            self.config['weight_cc'] =1000 * self.config['weight_cc']
            self.config['weight_nn'] =1000 * self.config['weight_nn']
            self.config['weight_pccp'] =1000 * self.config['weight_pccp']
            self.config['weight_cnnc'] =1000 * self.config['weight_cnnc']
            self.config['weight_pnnp'] =1000 * self.config['weight_pnnp']
            for i in range(3):
                line_min = lbfgs_rosetta.ArmijoLineMinimization(self.obj_func_np,self.obj_func_grad_np,True,len(rama),120)
                lbfgs_opt = lbfgs_rosetta.lbfgs(self.obj_func_np,self.obj_func_grad_np)
                rama=lbfgs_opt.run(rama,256,lbfgs_rosetta.absolute_converge_test,line_min,8000,self.obj_func_np,self.obj_func_grad_np,1e-9)
            newrama=rama+0.0
            newrama=torch.DoubleTensor(newrama) 
            current_energy =self.obj_func_np(rama)
            #self.outpdb(newrama,self.saveprefix+f'_{str(ilter)}'+'.pdb',energystr=str(current_energy))
            if current_energy < minenergy:
                print(current_energy,minenergy)
                minenergy=current_energy
                self.outpdb(newrama,self.saveprefix+'.pdb',energystr=str(current_energy))


    def outpdb(self,rama,savefile,start=0,end=10000,energystr=''):
        #rama=torch.cat([rama.view(self.L,2),self.betas],dim=-1)
        coor_np=a2b.quat2b(self.basex,rama.view(self.L,21)[:,9:]).data.numpy()
        other_np=a2b.quat2b(self.otherx,rama.view(self.L,21)[:,9:]).data.numpy()
        shaped_rama=rama.view(self.L,21)
        coor = torch.FloatTensor(coor_np)
        side_coor_NP = a2b.quat2b(self.sidex,torch.cat([shaped_rama[:,:9],coor[:,-1]],dim=-1)).data.numpy()
        
        Atom_name=[' P  '," C4'",' N1 ']
        Other_Atom_name = [" O5'"," C5'"," C3'"," O3'"," C1'"]
        other_last_name = ['O',"C","C","O","C"]

        side_atoms=         [' N1 ',' C2 ',' O2 ',' N2 ',' N3 ',' N4 ',' C4 ',' O4 ',' C5 ',' C6 ',' O6 ',' N6 ',' N7 ',' N8 ',' N9 ']
        side_last_name =    ['N',      "C",   "O",   "N",   "N",   'N',   'C',   'O',   'C',   'C',   'O',   'N',    'N', 'N','N']

        base_dict = rigid.base_table()
        last_name=['P','C','N']
        wstr=[f'REMARK {str(energystr)}']
        templet='%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s'
        count=1
        for i in range(self.L):
            if self.seq[i] in ['a','g','A','G']:
                Atom_name = [' P  '," C4'",' N9 ']
                #atoms = ['P','C4']

            elif self.seq[i] in ['c','u','C','U']:
                Atom_name = [' P  '," C4'",' N1 ']
            for j in range(coor_np.shape[1]):
                outs=('ATOM  ',count,Atom_name[j],self.seq[i],'A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],0,0,last_name[j],'')
                #outs=('ATOM  ',count,Atom_name[j],'ALA','A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],1.0,90,last_name[j],'')
                #print(outs)
                if i>=start-1 and i < end:
                    wstr.append(templet % outs)

            for j in range(other_np.shape[1]):
                outs=('ATOM  ',count,Other_Atom_name[j],self.seq[i],'A',i+1,other_np[i][j][0],other_np[i][j][1],other_np[i][j][2],0,0,other_last_name[j],'')
                #outs=('ATOM  ',count,Atom_name[j],'ALA','A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],1.0,90,last_name[j],'')
                #print(outs)
                if i>=start-1 and i < end:
                    wstr.append(templet % outs)



                count+=1
            
        wstr='\n'.join(wstr)
        wfile=open(savefile,'w')
        wfile.write(wstr)
        wfile.close()
    def outpdb_coor(self,coor_np,savefile,start=0,end=1000,energystr=''):
        #rama=torch.cat([rama.view(self.L,2),self.betas],dim=-1)
        Atom_name=[' P  '," C4'",' N1 ']
        last_name=['P','C','N']
        wstr=[f'REMARK {str(energystr)}']
        templet='%6s%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s'
        count=1
        for i in range(self.L):
            if self.seq[i] in ['a','g','A','G']:
                Atom_name = [' P  '," C4'",' N9 ']
                #atoms = ['P','C4']

            elif self.seq[i] in ['c','u','C','U']:
                Atom_name = [' P  '," C4'",' N1 ']
            for j in range(coor_np.shape[1]):
                outs=('ATOM  ',count,Atom_name[j],self.seq[i],'A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],0,0,last_name[j],'')
                #outs=('ATOM  ',count,Atom_name[j],'ALA','A',i+1,coor_np[i][j][0],coor_np[i][j][1],coor_np[i][j][2],1.0,90,last_name[j],'')
                #print(outs)
                if i>=start-1 and i < end:
                    wstr.append(templet % outs)
                count+=1
            
        wstr='\n'.join(wstr)
        wfile=open(savefile,'w')
        wfile.write(wstr)
        wfile.close()

if __name__ == '__main__': 
    fastafile=sys.argv[1]
    geofile  =sys.argv[2]
    saveprefix=sys.argv[3]
    coornpys =sys.argv[4:]

    stru=Structure(fastafile,geofile,saveprefix,coornpys)

    
    stru.foldning()
