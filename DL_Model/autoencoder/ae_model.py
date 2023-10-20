
## Implementing an encoder transformer model
##
## a full transformer model is probably overkill
## going to be hard to get enough data to do the encoder-decoder training
## (plus GPUs are super expensive)
##
## First try at doing this is an autoencoder architecture that can
## compress the sequence of domains into a lower-dimensionality
## representation. We might be able to cluster these or use other
## models to analyze them w.r.t. disease or other outcomes.

## NN architecture
import torch
from torch import nn, Tensor
from torch.nn import Sequential, ModuleList, Linear, ReLU, Sigmoid
from collections import OrderedDict

## Model
class AutoEncoderModel(nn.Module):
  ## layersizes should include the input and output layers
  ## e.g. for (input) 3 -> 2 -> 2 -> 1 (output), pass [3,2,2,1]
  def __init__(self, layersizes_encoder: list[int], layersizes_decoder: list[int]):
    super().__init__()
    nl_enc = len(layersizes_encoder)
    nl_dec = len(layersizes_encoder)

    self.encoder = OrderedDict()
    self.decoder = OrderedDict()
    for i in range(nl_enc-1):
      self.encoder['Lin'+str(i)] = Linear(layersizes_encoder[i], layersizes_encoder[i+1])
      self.encoder['Relu'+str(i)] = ReLU()
    for i in range(nl_dec-1):
      self.decoder['Lin'+str(i)] = Linear(layersizes_decoder[i], layersizes_decoder[i+1])
      self.decoder['Relu'+str(i)] = ReLU() if i != nl_dec-2 else Sigmoid()

    self.encoder = Sequential(self.encoder)
    self.decoder = Sequential(self.decoder)

  # def init_weights(self, m) -> None:
  #   initrange = 0.1
  #   biasfill = 0.01
  #   if(isinstance(m, Linear)):
  #     nn.init.uniform_(m.weight, -initrange, initrange)
  #     m.bias.data.fill_(biasfill)

  def forward(self, x):
    encoded = self.encoder(x)
    decoded = self.decoder(encoded)
    return decoded

  def encode(self, x):
    return self.encoder(x)