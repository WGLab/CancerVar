import torch
import torch.nn as nn
import torch.nn.functional as F


class ensemble_CNN1d_5layer(nn.Module):
    def __init__(self, CNN_hparams):
        super(ensemble_CNN1d_5layer, self).__init__()
        
        self.in_channels  = CNN_hparams['in_channels']    # [5, 20, 40]
        self.out_channels = CNN_hparams['out_channels']   # [20, 40, 200]
        self.cnnkernel    = CNN_hparams['cnnkernel']      # 5
        self.cnnstride    = CNN_hparams['cnnstride']      # 1
        self.dropout      = CNN_hparams['dropout']        # 0.1
        self.num_classes  = CNN_hparams['num_classes']

        self.cnn1 = nn.Sequential(
            nn.Conv1d(self.in_channels[0], self.out_channels[0], self.cnnkernel, self.cnnstride),
            nn.BatchNorm1d(self.out_channels[0]),
            nn.Tanh(),
            nn.Dropout(p=self.dropout),
        )
        
        self.cnn2 = nn.Sequential(
            nn.Conv1d(self.in_channels[1], self.out_channels[1], self.cnnkernel, self.cnnstride),
            nn.BatchNorm1d(self.out_channels[1]),
            nn.Tanh(),
            nn.Dropout(p=self.dropout),
        )
        
        self.cnn3 = nn.Sequential(
            nn.Conv1d(self.in_channels[2], self.out_channels[2], self.cnnkernel, self.cnnstride),
            nn.BatchNorm1d(self.out_channels[2]),
#             nn.Tanh(), # 20210223-20:57
#             nn.Dropout(p=self.dropout), # 20210223-20:57
        )
        
        self.linear1 = nn.Sequential(nn.Linear(13000, 2, bias=True))


    def forward(self, x):
        x  = self.cnn1(x)
        x  = self.cnn2(x)
        x  = self.cnn3(x)
        x  = x.view(x.size(0), -1)
        x1 = self.linear1(x)
        return x, x1

def ensembleNet():
    
    CNN_hparams                   = {}
    CNN_hparams['in_channels']    = [1, 10, 40]
    CNN_hparams['out_channels']   = [10, 40, 200]
    CNN_hparams['cnnkernel']      = 3
    CNN_hparams['cnnstride']      = 1
    CNN_hparams['dropout']        = 0.3
    CNN_hparams['num_classes']    = 2
    
    return ensemble_CNN1d_5layer(CNN_hparams)



class evs_CNN1d_5layer(nn.Module):
    def __init__(self, CNN_hparams):
        super(evs_CNN1d_5layer, self).__init__()
        
        self.in_channels  = CNN_hparams['in_channels']    # [5, 20, 40]
        self.out_channels = CNN_hparams['out_channels']   # [20, 40, 200]
        self.cnnkernel    = CNN_hparams['cnnkernel']      # 5
        self.cnnstride    = CNN_hparams['cnnstride']      # 1
        self.dropout      = CNN_hparams['dropout']        # 0.1
        self.num_classes  = CNN_hparams['num_classes']

        self.cnn1 = nn.Sequential(
            nn.Conv1d(self.in_channels[0], self.out_channels[0], self.cnnkernel, self.cnnstride),
            nn.BatchNorm1d(self.out_channels[0]),
            nn.Tanh(),
            nn.Dropout(p=self.dropout),
        )
        
        self.cnn2 = nn.Sequential(
            nn.Conv1d(self.in_channels[1], self.out_channels[1], self.cnnkernel, self.cnnstride),
            nn.BatchNorm1d(self.out_channels[1]),
            nn.Tanh(),
            nn.Dropout(p=self.dropout),
        )
        
        self.cnn3 = nn.Sequential(
            nn.Conv1d(self.in_channels[2], self.out_channels[2], self.cnnkernel, self.cnnstride),
            nn.BatchNorm1d(self.out_channels[2]),
#             nn.Tanh(), # 20210223-20:57
#             nn.Dropout(p=self.dropout), # 20210223-20:57
        )
        
        self.linear1 = nn.Sequential(nn.Linear(8400, 2, bias=True))


    def forward(self, x):
        x  = self.cnn1(x)
        x  = self.cnn2(x)
        x  = self.cnn3(x)
        x  = x.view(x.size(0), -1)
        x1 = self.linear1(x)
        return x, x1

def evsNet():
    
    CNN_hparams                   = {}
    CNN_hparams['in_channels']    = [1, 10, 40]
    CNN_hparams['out_channels']   = [10, 40, 200]
    CNN_hparams['cnnkernel']      = 3
    CNN_hparams['cnnstride']      = 1
    CNN_hparams['dropout']        = 0.3
    CNN_hparams['num_classes']    = 2
    
    return evs_CNN1d_5layer(CNN_hparams)
