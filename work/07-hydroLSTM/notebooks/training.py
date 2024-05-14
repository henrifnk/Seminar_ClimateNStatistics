import lightning.pytorch as pl
from lightning.pytorch.callbacks import EarlyStopping, LearningRateMonitor
from lightning.pytorch.tuner import Tuner
# from pytorch_forecasting import TimeSeriesDataSet, TemporalFusionTransformer, Baseline
from pytorch_forecasting.metrics import MAE, SMAPE, PoissonLoss, QuantileLoss
import os
import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

project_dir = os.path.join('work','07-hydroLSTM')

meteo = pd.read_pickle(os.path.join(project_dir,'data','processed','meteo_processed.pkl'),compression = 'zip')

stream = pd.read_pickle(os.path.join(project_dir,'data','processed','stream_processed.pkl'),compression= 'zip')

stream = stream.iloc[:3652]
data = stream.join(meteo).drop(columns = ['Datum','date'])

class TimeseriesDataset(Dataset):   
    def __init__(self, X: np.ndarray, seq_len : int):
        super().__init__()
        self.X = torch.tensor(X).float()
        self.seq_len = seq_len

    def __len__(self):
        return self.X.__len__() - (self.seq_len)

    def __getitem__(self, index):
        return self.X[index:index+self.seq_len], self.X[index+self.seq_len]
    

class StreamFlowDataModule(pl.LightningDataModule):
    def __init__(self, X : np.array, num_workers = 0,seq_len = 364, batch_size=16):
        super().__init__()
        self.X = X
        self.seq_len = seq_len
        self.batch_size = batch_size
        self.num_workers = num_workers

    def prepare_data(self):
        self.X_val = self.X[0:int(len(self.X)*0.15)]
        self.X_train = self.X[int(len(self.X)*0.15)+1:int(len(self.X)*0.85)]
        self.X_test = self.X[int(len(self.X)*0.85)+1:]
        
    def train_dataloader(self):
        train_dataset = TimeseriesDataset(self.X_train, 
                                          seq_len=self.seq_len)
        
        train_loader = DataLoader(train_dataset,
                                  batch_size = self.batch_size, 
                                  shuffle = True, 
                                  num_workers = self.num_workers)
        return train_loader
    
    def val_dataloader(self):
        val_dataset = TimeseriesDataset(self.X_val, 
                                        seq_len=self.seq_len)
        
        val_loader = DataLoader(val_dataset, 
                                batch_size = self.batch_size, 
                                shuffle = False, 
                                num_workers = self.num_workers)
        return val_loader
    
dm = StreamFlowDataModule(data.to_numpy())

dm.prepare_data()
print(next(iter(dm.train_dataloader())))