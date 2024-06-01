import lightning.pytorch as pl
from lightning.pytorch.callbacks import EarlyStopping, LearningRateMonitor
from lightning.pytorch.tuner import Tuner
from pytorch_forecasting import TimeSeriesDataSet, TemporalFusionTransformer, AutoRegressiveBaseModel,AutoRegressiveBaseModelWithCovariates
from pytorch_forecasting.metrics import MAE,RMSE
import os
import pandas as pd # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np # type: ignore
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from pytorch_forecasting.models.nn import LSTM
from typing import Dict
from sklearn.preprocessing import StandardScaler
import random



class LSTMModel(AutoRegressiveBaseModel):
    def __init__(
        self,
        input_size:int,
        target: str,
        target_lags: Dict[str, Dict[str, int]],
        n_layers: int,
        hidden_size: int,
        dropout: float = 0.1,
        **kwargs,
    ):
        # arguments target and target_lags are required for autoregressive models
        # even though target_lags cannot be used without covariates
        # saves arguments in signature to `.hparams` attribute, mandatory call - do not skip this
        self.save_hyperparameters()
        # pass additional arguments to BaseModel.__init__, mandatory call - do not skip this
        super().__init__(**kwargs)

        # use version of LSTM that can handle zero-length sequences
        self.lstm = LSTM(
            hidden_size=self.hparams.hidden_size,
            input_size=self.hparams.input_size,
            num_layers=self.hparams.n_layers,
            dropout=self.hparams.dropout,
            batch_first=True,
        )
        self.output_layer = nn.Linear(self.hparams.hidden_size, 1)
        


    def encode(self, x: Dict[str, torch.Tensor]):
        # we need at least one encoding step as because the target needs to be lagged by one time step
        # because we use the custom LSTM, we do not have to require encoder lengths of > 1
        # but can handle lengths of >= 1
        assert x["encoder_lengths"].min() >= 1
        input_vector = x["encoder_cont"].clone()
        # lag target by one
        input_vector[..., self.target_positions] = torch.roll(
            input_vector[..., self.target_positions], shifts=1, dims=1
        )
        input_vector = input_vector[:, 1:]  # first time step cannot be used because of lagging

        # determine effective encoder_length length
        effective_encoder_lengths = x["encoder_lengths"] - 1
        # run through LSTM network
        _, hidden_state = self.lstm(
            input_vector, lengths=effective_encoder_lengths, enforce_sorted=False  # passing the lengths directly
        )  # second ouput is not needed (hidden state)
        return hidden_state

    def decode(self, x: Dict[str, torch.Tensor], hidden_state):
        # again lag target by one
        input_vector = x["decoder_cont"].clone()
        input_vector[..., self.target_positions] = torch.roll(
            input_vector[..., self.target_positions], shifts=1, dims=1
        )
        # but this time fill in missing target from encoder_cont at the first time step instead of throwing it away
        last_encoder_target = x["encoder_cont"][
            torch.arange(x["encoder_cont"].size(0), device=x["encoder_cont"].device),
            x["encoder_lengths"] - 1,
            self.target_positions.unsqueeze(-1),
        ].T
        input_vector[:, 0, self.target_positions] = last_encoder_target

        if self.training:  # training mode
            lstm_output, _ = self.lstm(input_vector, hidden_state, lengths=x["decoder_lengths"], enforce_sorted=False)

            # transform into right shape
            prediction = self.output_layer(lstm_output)
            prediction = self.transform_output(prediction, target_scale=x["target_scale"])

            # predictions are not yet rescaled
            return prediction

        else:  # prediction mode
            target_pos = self.target_positions
            test = input_vector[:, 0, target_pos]

            def decode_one(idx, lagged_targets, hidden_state):
                x = input_vector[:, [idx]]
                # overwrite at target positions
                x[:, 0, target_pos] = lagged_targets[-1]  # take most recent target (i.e. lag=1)
                lstm_output, hidden_state = self.lstm(x, hidden_state)
                # transform into right shape
                prediction = self.output_layer(lstm_output)[:, 0]  # take first timestep
                return prediction, hidden_state

            # make predictions which are fed into next step
            output = self.decode_autoregressive(
                decode_one,
                first_target=input_vector[:, 0, target_pos],
                first_hidden_state=hidden_state,
                target_scale=x["target_scale"],
                n_decoder_steps=input_vector.size(1),
            )

            # predictions are already rescaled
            return output

    def forward(self, x: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        hidden_state = self.encode(x)  # encode to hidden state
        output = self.decode(x, hidden_state)  # decode leveraging hidden state

        return self.to_network_output(prediction=output)
    
    @property
    def target_positions(self) -> torch.LongTensor:
        """
        Positions of target variable(s) in covariates.

        Returns:
            torch.LongTensor: tensor of positions.
        """
        # todo: expand for categorical targets
        return torch.tensor(
            [-1],
            device=self.device,
            dtype=torch.long,
        )


if __name__ == '__main__':

    stream = pd.read_pickle('/Users/lemarx/Documents/01_projects/Seminar_ClimateNStatistics/work/07-hydroLSTM/data/processed/stream_processed.pkl',compression= 'zip')
    meteo = pd.read_pickle('/Users/lemarx/Documents/01_projects/Seminar_ClimateNStatistics/work/07-hydroLSTM/data/processed/meteo_processed.pkl',compression = 'zip')

    stream['15207507_smooth'] = stream['15207507'].rolling(window=3).mean()


    stream = stream.iloc[:3652]
    meteo = meteo.shift(1).drop(columns = 'date')
    data = stream.join(meteo).drop(columns = 'Datum').dropna()
    data['groups'] = 0
    data['date'] = data.index


    max_pred_len = 2
    enc_len = 3
    training_cutoff = int(data["date"].max()*0.7) - max_pred_len
    val_cutoff = int(data["date"].max()*0.85) - max_pred_len

    training = TimeSeriesDataSet(
        data= data.iloc[: training_cutoff],
        target='15207507_smooth',#'15207507_smooth'
        group_ids=['groups'],#groups
        time_idx='date',#date
        max_encoder_length=enc_len,
        min_encoder_length=enc_len,
        min_prediction_length=max_pred_len,
        max_prediction_length=max_pred_len,
        time_varying_unknown_reals=['15207507_smooth'],
        time_varying_known_reals = data.columns[22:24].values.tolist(),
        target_normalizer = None,
        scalers = {col : None for col in data.columns[22:24]},
        # lags = {col : [i for i in range(1,max_pred_len)] for col in data.columns[20:-2]}, #lags = {'d2m_8200' : [i for i in range(1,max_pred_len)]}
        #data.columns[20:-1].values.tolist()
    )

    validation = TimeSeriesDataSet.from_dataset(training,data.iloc[training_cutoff: val_cutoff],min_prediction_idx=training_cutoff+365,stop_randomization=True)

    val_dataloader = validation.to_dataloader(train=False, batch_size=1, num_workers=1)

    train_dataloader = training.to_dataloader(train=True, batch_size=1, num_workers=1)




    trainer = pl.Trainer(
    accelerator="cpu",
    max_epochs=1,
    # clipping gradients is a hyperparameter and important to prevent divergance
    # of the gradient for recurrent neural networks
    gradient_clip_val=0.1,
    )

    model = LSTMModel.from_dataset(training,
                               input_size = 3,
                               learning_rate = 0.001,
                               n_layers=1,
                               hidden_size = 2,
                               loss = MAE(),
                               dropout = 0.2,
                               optimizer = 'adam',
                               reduce_on_plateau_patience = 10,
                               weight_decay = 0.001)

    trainer.fit(
        model,
        train_dataloaders=train_dataloader,
        val_dataloaders=val_dataloader,
    )
