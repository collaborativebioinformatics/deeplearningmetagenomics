# Quick Deep Learning Models

The goal of this was getting something working quickly, so these are very much off-the-shelf algorithms with no tuning.
It would probably be worthwhile to tune these models so they actually perform well.

# LLM

The LLM is located in the `transformer` folder. That folder has a `README.md` file written by the people at pytorch that
is very good. That document specifies all the commandline arguments.

To train a model, make a folder inside `transformer/data`. Within that folder, add three files: `train.txt, test.txt, valid.txt`.
These correspond to train, test, and validation data (respectively). Data is expected to be space separated words. similar to
like a paragraph or an essay. Words should be masked by replacing them with `<unk>` so that the model can learn to fill in gaps.

## Building a training set

Given a very long list of domain annotations, I built a train/test/validation set as follows:

```{r}
# R script to parse a character vector of domain accession numbers into train/test/validation data

my_domain_annots <- 'path/to/file.RData' # Assuming this has an object 'txt_seq'
outdir <- '.'
train_pct <- 0.8
test_pct <- 0.1
valid_pct <- 1 - train_pct - test_pct
pct_to_mask <- 0.05

load(my_domain_annots)

# txt_seq is a large character vector of length 1
# it contains space-separated accession numbers
txt_seq <- strsplit(txt_seq, ' ')[[1]]
orig_seq <- txt_seq

tr_end <- floor(length(txt_seq) * train_pct)
tst_start <- tr_end + 1L
tst_end <- tr_start + floor(length(txt_seq) * test_pct)
vd_start <- tst_end + 1L
vd_end <- length(txt_seq)

# mask the specified number of positions
mask_pos <- sample(seq_along(txt_seq), floor(length(txt_seq)*pct_to_mask))
txt_seq[mask_pos] <- '<unk>'

# build the testing sets
train <- txt_seq[seq(1,tr_end)]
test <- txt_seq[seq(tst_start, tst_end)]
validate <- txt_seq[seq(vd_start, vd_end)]

# save files before masking, just in case
cat(train, file=file.path(outdir, 'train.txt'))
cat(test, file=file.path(outdir, 'test.txt'))
cat(validate, file=file.path(outdir, 'valid.txt'))


# save off relevant info
save(orig_seq, mask_pos, tr_end, tst_start, tst_end, vd_start, vd_end, file='dataset_details.RData')
```

Then we can train a model with the following commandline call:
```
python3 main.py --data /path/to/newdata/folder --model Transformer
```

For a general LSTM, omit the `--model Transformer` argument.

# Autoencoder

The autoencoder is not quite as clean because I wrote more of it...

An example of running it with some random data is in `autoencoder/testing_randomdata.py`.

Building a model in python looks like this:

```{python}
import torch
from ae_model import AutoEncoderModel

# specify the input sizes of each layer and the output size of the final layer
# here we have 1000 => 100 => 80 => 50 => 10 => 50 => 80 => 100 => 1000
# note the smallest layer is duplicated (last entry of encoder, first of decoder)
encoder_layersizes = [1000, 100, 80, 50, 10]
decoder_layersizes = [10, 50, 80, 100, 1000]

# build the model
model = AutoEncoderModel(encoder_layersizes, decoder_layersizes)

# specify the loss function, example here is mean-squared error
loss_function = torch.nn.MSELoss()

# specify the optimizer, here we'll use Adam
# second argument is learning rate
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-8)

# print the model to make sure we configured it right
print(model)

## Training for a single sample
## imagine I have argument `values`, length 1000 vector

reconstructed = model(values)
loss = loss_function(reconstructed, values)

# zero out the optimizer and then update the model
optimizer.zero_grad()
loss.backward()
optimizer.step()


## Encoding a single sample
encoded = model.encode(values)
```

For testing on larger samples, look into using the `DataLoader` object. Autoencoders just need a single dataset, since
the goal is to learn to map inputs to inputs through a reduced dimensionality space. Thus, we don't need labeled data,
just a big set of data. An example of using a `DataLoader` object is in the example file.
