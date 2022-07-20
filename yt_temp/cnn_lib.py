#!/usr/bin/env python
# coding: utf-8

import torch
import torch.nn as  nn
import torch.nn.functional as F
from torchvision.transforms import transforms
import pandas as pd
import time
from tqdm import tqdm
import numpy as np
import subprocess
import torch.optim as optim
import torchvision
import os
from PIL import Image
from visium_datasets import get_symbol_values, unparse_tile_name
from torch.utils.data import Dataset
from pathlib import Path
import csv

def calc_cor(outputs, labels):
    num_gene = outputs.shape[1]
    corR = []
    for i in range(num_gene):
        corR.append(np.corrcoef(outputs[:,i].to('cpu').detach().numpy(), labels[:,i].to('cpu').detach().numpy())[0,1])
    corR = np.array(corR)
    corR[np.isnan(corR)] = 0.0      
    print("corR: "+str(corR))

    return np.mean(corR)


def loss_function(outputs, labels):
    criterion = nn.SmoothL1Loss()
    num_gene = outputs.shape[1]
    loss = 0
    for i in range(num_gene):
        loss += criterion(outputs[:,i], labels[:,i]) / num_gene

    return loss


def build_resnet50(pretrained=True, num_features=1):
    # load resnet50 model
    net = torchvision.models.resnet50(pretrained=pretrained)
    # need to change the fc-layer for regression, since origin torch resnet model is for classification
    fc_inputs = net.fc.in_features
    # 添加全连接层
    net.fc = nn.Linear(fc_inputs, num_features)
    # train mode
    net.train()

    return net


def run_train(net, dataloader_dict, optimizer, num_epochs=100, device='cpu', outDir='./output', name=''):
    ### set device
    if str(device) != 'cpu':
        net.to(device)
        
    res_df = pd.DataFrame(columns=['train_loss','valid_loss','train_cor','valid_cor'])

    valid_loss_best = 1e+100;
    valid_cor_best = -1e+100;
    # epoch roop
    for epoch in range(num_epochs+1):
        print('Epoch {}/{}'.format(epoch, num_epochs))
        print('-------------')

        train_loss = 0
        valid_loss = 0
        train_cor = 0
        valid_cor = 0
        
        # train and valid roop per epoch
        for phase in ['train', 'valid']:
            net.train() if phase == 'train' else net.eval() # train or eval mode

            epoch_loss = 0.0  # sum of loss
            epoch_corrects = 0  # sum of corrects or correlation

            # skip trainging if epoch == 0
            if (epoch == 0) and (phase == 'train'): continue
                
            # extract minibatch from dataloader
            for inputs, labels in tqdm(dataloader_dict[phase]):
                # if GPU is avalable
                if str(device) != 'cpu':
                    inputs = inputs.to(device)
                    labels = labels.to(device)

                # initialize optimizer
                optimizer.zero_grad()

                # forward calculation
                with torch.set_grad_enabled(phase == 'train'):
                    outputs = net(inputs)
                    
                    # calculate loss
                    loss = loss_function(outputs, labels)
                    
                    # backpropagation if training phase
                    if phase == 'train': 
                        loss.backward()
                        optimizer.step()

                    # update sum of loss
                    epoch_loss += loss.item() * inputs.size(0)
                    
                    # update sum of cor
                    epoch_corrects += calc_cor(outputs, labels)


            # print loss
            epoch_loss = epoch_loss / len(dataloader_dict[phase].dataset)
            print('{} Loss: {:.4f}'.format(phase, epoch_loss))

            # train loss or valid loss
            if phase == 'train':
                train_cor = float(epoch_corrects)
                train_loss = epoch_loss
            else:
                valid_cor = float(epoch_corrects)
                valid_loss = epoch_loss
                
        ### append loss to DataFrame
        print('train_loss: {} valid_loss: {} train_cor: {} valid_cor: {}'.format(train_loss,valid_loss,train_cor,valid_cor))
        res_df = res_df.append([pd.Series([train_loss,valid_loss,train_cor,valid_cor],index=res_df.columns)], ignore_index=True)
        
        ### save training_loss.txt
        res_df.to_csv(outDir+"/training_loss_"+name+".txt", sep='\t', float_format='%.6f')

        ### save best model
        save_best = False
        
        if valid_cor_best < valid_cor:
            valid_cor_best = valid_cor
            save_best = True

        if save_best:
                subprocess.call(['rm','-r',outDir+'/model_'+name+'/'])
                subprocess.call(['mkdir',outDir+'/model_'+name+'/'])
                torch.save(net.state_dict(), outDir+"/model_"+name+"/model_"+str(epoch)+".pth")


class VisiumImagePredDataset(Dataset):
    def __init__(
        self,
        data_path: Path,
        count_mat_dir: str,
        raw_tile_dir: str,
        tile_size: tuple,
        symbol: str,
        
        #transform=None,
        # at least ToTensor is needed
        transform=transforms.Compose([transforms.ToTensor(),
                                     transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225])]),
        target_transform=None
    ):
        count_mat_dpath = data_path / count_mat_dir
        self.tile_labels = get_symbol_values(count_mat_dpath, symbol)
        self.tile_dpath = data_path / raw_tile_dir
        self.tile_size = tile_size
        self.symbol = symbol
        self.transform = transform
        self.target_transform = target_transform

    def __len__(self):
        return len(self.tile_labels)

    def __getitem__(self, idx):
        sample_id, tissue_id, barcode = self.tile_labels.index[idx]
        tile_name = unparse_tile_name(sample_id, tissue_id, barcode, self.tile_size)
        tile_fpath = self.tile_dpath / tile_name
        # image = read_image(str(tile_fpath))
        image = Image.open(str(tile_fpath))  # RGB, W, H, C
        #label = self.tile_labels.iloc[idx]
        label = self.tile_labels.iloc[idx].tolist()
        if self.transform:
            image = self.transform(image)
        if self.target_transform:
            label = self.target_transform(label)
        return image, label, tile_name


def get_prediction(net, data_dir, device='cpu', out_dir='./output', symbol='GFAP', name='', with_label=True):
    ### set device
    # if str(device) != 'cpu':
    #     net.to(device)
    
    net.eval()

    data_dir = Path(data_dir) 
    count_mat_dir = "assay_data" 
    raw_tile_dir = "Raw_Tiles" 
    tile_size = (256, 256) 
    symbol = symbol

    pred_data = VisiumImagePredDataset(data_dir, count_mat_dir, raw_tile_dir, tile_size, symbol)

    out_fdir = str(Path(out_dir) / name) +'_'+ symbol + '_predictions.csv'
    print('File Path: ',out_fdir)
    csv_file = open(str(out_fdir),'w')
    writer = csv.writer(csv_file)
    writer.writerow(['tile_name', 'label', 'prediction'])
    for i, data in enumerate(pred_data):
        print(i,'/', len(pred_data))
        image = data[0]; label = data[1]; tile_name = data[2]
        pred = net(image.unsqueeze(0))
        pred = np.round(pred.item(),3)
        writer.writerow([tile_name, label, pred])
    csv_file.close()

def parse_tile_name(tile_name: str):
    """Parse tile name into its components"""
    components = tile_name.split("-")
    sample_id = components[0]
    tissue_id = components[1]
    barcode = components[2] + '-' + components[3]
    size = components[4].removesuffix(".jpeg")
    size = size.split("x")
    size = tuple(int(dim) for dim in size)
    return {
        "sample_id": sample_id,
        "tissue_id": tissue_id,
        "barcode": barcode,
        "size": size,
    }