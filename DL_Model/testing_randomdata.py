# testing model.py
from ae_model import AutoEncoderModel
import torch
from torchvision import datasets
from torchvision import transforms
import matplotlib.pyplot as plt


## Example taken from online g4g tutorial -- trying to replicate results
# link: https://www.geeksforgeeks.org/implementing-an-autoencoder-in-pytorch/

# Transforms images to a PyTorch Tensor
tensor_transform = transforms.ToTensor()

# Download the MNIST Dataset
dataset = datasets.MNIST(root = "./data",
                         train = True,
                         download = True,
                         transform = tensor_transform)

# DataLoader is used to load the dataset
# for training
loader = torch.utils.data.DataLoader(dataset = dataset,
                                     batch_size = 32,
                                     shuffle = True)


## Encoder: 28*28 = 784 ==> 128 ==> 64 ==> 36 ==> 18 ==> 9
## Decoder: 9 ==> 18 ==> 36 ==> 64 ==> 128 ==> 784 ==> 28*28 = 784

encoder_layersizes = [784, 128, 64, 36, 18, 9]
decoder_layersizes = [9, 18, 36, 64, 128, 784]
model = AutoEncoderModel(encoder_layersizes, decoder_layersizes)
loss_function = torch.nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(),
                              lr=1e-3,
                              weight_decay=1e-8)

print(model)
#print(list(model.parameters())[0].clone())
epochs = 10
outputs = []
losses = []
for epoch in range(epochs):
  print("Epoch #", epoch+1)
  for (image, _) in loader:
    # Reshaping the image to (-1, 784)
    image = image.reshape(-1, 28*28)

    # Output of Autoencoder
    reconstructed = model(image)

    # Calculating the loss function
    loss = loss_function(reconstructed, image)

    # The gradients are set to zero,
    # the gradient is computed and stored.
    # .step() performs parameter update
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    # Storing the losses in a list for plotting
    losses.append(loss.detach().numpy())
  #print(list(model.parameters())[0].clone())
  outputs.append((epochs, image, reconstructed))

# Defining the Plot Style
plt.style.use('fivethirtyeight')
plt.xlabel('Iterations')
plt.ylabel('Loss')

# Plotting the last 100 values
#plt.plot([loss.detach().numpy() for loss in losses[-100:]])
plt.plot(losses[-100:])
plt.show()



for i, item in enumerate(image[-3:]):
  # Reshape the array for plotting
  item = item.reshape(-1, 28, 28)
  plt.imshow(item[0].detach().numpy())
  plt.show()
for i, item in enumerate(reconstructed[-3:]):
  item = item.reshape(-1, 28, 28)
  plt.imshow(item[0].detach().numpy())
  plt.show()