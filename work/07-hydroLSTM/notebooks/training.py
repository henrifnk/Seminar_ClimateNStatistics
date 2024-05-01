import lightning.pytorch as pl
from lightning.pytorch.callbacks import EarlyStopping, LearningRateMonitor
from lightning.pytorch.tuner import Tuner
from pytorch_forecasting import TimeSeriesDataSet, TemporalFusionTransformer
import pandas as pd
import os

project_dir = os.path.join('work','07-hydroLSTM')

meteo = pd.read_pickle(os.path.join(project_dir,'data','processed','meteo_processed.pkl'),compression = 'zip')

stream = pd.read_pickle(os.path.join(project_dir,'data','processed','stream_processed.pkl'),compression= 'zip')

print(meteo.head())