import torch
import torch.autograd as autograd
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import torchvision
import torchvision.transforms as transforms
import torchvision.datasets as datasets

import matplotlib.pyplot as plt
from matplotlib import cm

from skimage import io

import numpy as np 


class Statistical_Dataset(Dataset):
  def __init__(self,X,y):
    if X.shape[0] != y.shape[0]:
      print('First dimesion of X and y must be the same')
      return
      
    self.X = X
    self.y = y
  
  def __len__(self):
    return self.X.shape[0]
  
  def __getitem__(self,idx):
    Xi = self.X[idx]
    yi = self.y[idx]
    return Xi,yi

class MLP(nn.Module):
    def __init__(self, neural_layers, act_fun):
        super().__init__()
        
        assert not len(neural_layers) < 2, "len(neural_layes) must be larger than 2"
        assert not len(act_fun) < 1, "len(act_fun) must be larger than 1"
        assert not len(neural_layers) - len(act_fun) != 1, "len(neural_layes) - len(act_fun)  = 1"
   
        layers_list  = []
        for i in range(len(neural_layers)-1):
            layers_list.append(
                nn.Linear(neural_layers[i],neural_layers[i+1]).cuda()
            )
        self.layers = nn.ModuleList(layers_list)
        self.activation_fun = nn.ModuleList(act_fun)
        
    def forward(self, x):
        for index in range(len(self.layers)):
            x = self.layers[index](x)
            x = self.activation_fun[index](x)
        return x


class MINE:
    def __init__(self, net):
        self.mine_net = net
        self.optimizer = torch.optim.Adam(self.mine_net.parameters(), lr=0.0001)
        
    def mutual_information(self, joint, marginal, is_first_metric=True):
        t = self.mine_net(joint)
        if is_first_metric:
            et = torch.exp(self.mine_net(marginal))
            eps=1e-7
            mean_et = torch.log(torch.mean(et)+eps)
        else:
            et = torch.exp(self.mine_net(marginal)-1)
            mean_et = torch.mean(et)
        mi_lb = torch.mean(t) - mean_et
        return mi_lb, t, et
        
    def train_step(self, joint, marginal, ma_et, is_first_metric=True, ma_rate=0.01):
        # Distributions
        joint = torch.autograd.Variable(torch.FloatTensor(joint)).cuda()
        marginal = torch.autograd.Variable(torch.FloatTensor(marginal)).cuda()
        # Calculate MI
        mi_lb, t, et = self.mutual_information(joint, marginal, is_first_metric)
        # Unbiased with moving avg
        ma_et = (1-ma_rate)*ma_et + ma_rate*torch.mean(et)
        # Loss
        loss = -(torch.mean(t) - (1/ma_et.mean()).detach()*torch.mean(et))
        # Train
        self.mine_net.zero_grad()
        autograd.backward(loss)
        self.optimizer.step()
        return mi_lb.detach().cpu().numpy(), ma_et


def calcularTC(samplesjoint1, samplesjoint2, n_epochs, batch_size):
  dataset = Statistical_Dataset(samplesjoint1, samplesjoint2)
  dataloader = torch.utils.data.DataLoader(dataset, batch_size = batch_size)
  n_variables = samplesjoint1.shape[1]
  ma_et = 1
  ma_et_lista = [[1] for i in range(1, n_variables)]

  listaderedes = [MINE(MLP([1+i, 12, 1],[nn.ReLU(), nn.Identity()])) for i in range(1, n_variables)]

  for e in range(n_epochs) :
    tc = []

    infos_mutuas_batch = [[] for i in range(1, n_variables)]

    for i, (joint1,joint2) in enumerate(dataloader): #for del batch

      tc_sum_batch = 0

      for j in range(1, n_variables): #for de cada sumando de la info mutua. N variables son N-1 sumandos

        marginal1 = joint1[:,j].float()
        marginal1 = torch.unsqueeze(marginal1, 1).float()
        marginal2 = joint2[:,:j].float()

        marginal1resto = joint1[:,:j].float()
        marginal1cat = torch.cat([marginal1, marginal1resto], dim=1).float()

        mi, ma_et = listaderedes[j-1].train_step(marginal1cat,torch.cat([marginal1, marginal2], dim=1).float(), ma_et_lista[j-1][-1])
        ma_et_lista[j-1].append(ma_et)
      
        tc_sum_batch += mi

      tc.append(tc_sum_batch)

    if (e+1)%(n_epochs/10) == 0:
      print('Epoch {:3d} || TC : {:.3f}'.format(e+1,np.mean(tc)))

  return np.mean(tc)



def calcularDTC(samplesjoint1, samplesjoint2, n_epochs, batch_size):
  dataset = Statistical_Dataset(samplesjoint1, samplesjoint2)
  dataloader = torch.utils.data.DataLoader(dataset, batch_size = batch_size)
  n_variables = samplesjoint1.shape[1]
  ma_et = 1
  ma_et_lista2 = [[1] for i in range(1, n_variables-1)]
  ma_et_lista3 = [[1] for i in range(n_variables-1, 0, -1)]

  net = MLP([n_variables, 12, 1],[nn.ReLU(), nn.Identity()])

  mine1 = MINE(net)

  listaderedes1 = [MINE(MLP([n_variables, 12, 1],[nn.ReLU(), nn.Identity()])) for i in range(1, n_variables-1)]
  listaderedes2 = [MINE(MLP([i, 12, 1],[nn.ReLU(), nn.Identity()])) for i in range(n_variables-1, 0, -1)]

  for e in range(n_epochs):
    dtc = []

    for i, (joint1,joint2) in enumerate(dataloader):#(joint1, joint2)
      
      dtc_sum = 0

      marginal1max = joint1[:, n_variables-1]
      marginal1max = torch.unsqueeze(marginal1max, 1)
      marginal2max = joint2[:, :n_variables-1]

      mimax, ma_et = mine1.train_step(joint1[:,:].float(),torch.cat([marginal1max, marginal2max], dim=1).float(), ma_et)

      dtc_sum += mimax

      for j in range(1, n_variables-1):

        # I(X;Y,Z) - I(X;Z)
        ## JOINT 1
        # X
        marginal1 = joint1[:,j]
        marginal1 = torch.unsqueeze(marginal1, 1)
        # Y, Z
        marginaltriple = torch.cat([joint1[:,:j], joint1[:,j+1:]], dim=1)
        # Z
        marginaldoble = joint1[:,j+1:]

        ## JOINT 2
        # Y, Z
        marginal2triple = torch.cat([joint2[:,:j], joint2[:,j+1:]], dim=1)
        # Z
        marginal2doble = joint2[:,j+1:]
        

        mi1, ma_et2 = listaderedes1[j-1].train_step(torch.cat([marginal1, marginaltriple], dim=1).float(),torch.cat([marginal1, marginal2triple], dim=1).float(), ma_et_lista2[j-1][-1])
        mi2, ma_et3 = listaderedes2[j-1].train_step(torch.cat([marginal1, marginaldoble], dim=1).float(),torch.cat([marginal1, marginal2doble], dim=1).float(), ma_et_lista3[j-1][-1])

        ma_et_lista2[j-1].append(ma_et2)
        ma_et_lista3[j-1].append(ma_et3)

        infos_mutuas_batch_pos[j-1].append(mi1)
        infos_mutuas_batch_neg[j-1].append(mi2)
        
        dtc_sum += (mi1-mi2)
      dtc.append(dtc_sum)

    if (e+1)%(n_epochs/10) == 0:
      print('Epoch {:3d} || DTC : {:.3f}'.format(e+1,np.mean(dtc)))

  return np.mean(dtc)