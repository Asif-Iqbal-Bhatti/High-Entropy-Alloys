#!/usr/bin/env python3
	
'''
#========================================================
# AUTHOR:: AsifIqbal -> @AIB_EM
# USAGE :: TRAINING ON A MATERIAL PROJECT DATABASE
#       :: USING PYTORCH
# https://learn.microsoft.com/en-us/windows/ai/windows-ml/tutorials/pytorch-analysis-data
# https://towardsdatascience.com/pytorch-tabular-binary-classification-a0368da5bb89
#========================================================
'''

from torch import nn, optim
import pandas as pd
import torch.nn.functional as F 
import matplotlib.pyplot as plt
import torch, numpy as np
from sklearn.datasets import make_circles
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler  
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import random_split, TensorDataset
from sklearn.metrics import confusion_matrix, classification_report

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu") 
print("The model will be running on", device, "device") 

#=============== READ A EXISTING DATASET
df = pd.read_csv('test.csv')
#X =  torch.tensor(df['TE_at'].values, dtype=torch.float32)
#y =  torch.tensor(df['Total_SRO'].values, dtype=torch.float32)
X = df.iloc[:, 0].to_numpy()
y = df.iloc[:, -1].to_numpy()
X= X.reshape(-1, 1)
y= y.reshape(-1, 1)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.30, random_state=1234)
print(len(X_train), len(X_test), len(X_train)+len(X_test))
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

#=============== VISUALIZE THE DATA
def testDataplot(X_train, y_train):
	fig, (train_ax, test_ax) = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(15, 5))
	train_ax.scatter(x=X_train, y=y_train)
	train_ax.set_title("Training Data")
	train_ax.set_xlabel("Feature #0")
	train_ax.set_ylabel("Feature #1")
	test_ax.scatter(x=X_test, y=y_test)
	test_ax.set_xlabel("Feature #0")
	test_ax.set_title("Testing data")
	plt.savefig('TE_vs_SRO.png', dpi=300,bbox_inches='tight')
testDataplot(X_train, y_train)

#=============== CONVERT DATA TO TORCH TENSORS
'''
# -----> MANUAL CONSTRUCTION
data = TensorDataset(X, y)
number_rows = len(X)  
test_split = int(number_rows*0.3)  
validate_split = int(number_rows*0.2) 
train_split = number_rows - test_split - validate_split     
train_set, validate_set, test_set = random_split(data, [train_split, validate_split, test_split])
'''

class Data(Dataset):
	def __init__(self, X_data, y_data):
			self.X_data = torch.tensor(X_data, dtype=torch.float32)
			self.y_data = torch.tensor(y_data, dtype=torch.float32)
	def __getitem__(self, index):
			return self.X_data[index], self.y_data[index]
	def __len__ (self):
			return len(self.X_data)
				
train_data = Data(X_train, y_train)			
test_data = Data(X_test, y_test)
#=============== INSTANTIATE TRAINING AND TEST DATA
train_batch_size = 100
train_dataloader = DataLoader(dataset=train_data, batch_size=train_batch_size, shuffle=False)
test_dataloader = DataLoader(dataset=test_data, batch_size=1)

input_dim = 1
output_dim = 1

class NeuralNetwork(nn.Module):
	def __init__(self, input_dim, output_dim):
			super(NeuralNetwork, self).__init__()
			self.layer_1 = nn.Linear(input_dim, 64)
			self.layer_2 = nn.Linear(64, 32)
			self.layer_3 = nn.Linear(32, output_dim)
			
	def forward(self, x):
			x1 = torch.relu(self.layer_1(x))
			x2 = torch.relu(self.layer_2(x1))
			x4 = self.layer_3(x2)
			return x4

model = NeuralNetwork(input_dim, output_dim)
model.to(device)
print((model))

#loss_fn = nn.CrossEntropyLoss()
loss_fn = nn.MSELoss()
optimizer = torch.optim.SGD(model.parameters(), lr=0.01)
#optimizer = torch.optim.Adam(model.parameters(), lr=0.3)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=4, gamma=0.9)
#scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, start_factor=0.5, total_iters=4)

def binary_acc(y_pred, y_test):
	correct_results_sum = (y_pred).sum().float()
	acc = correct_results_sum/y_test.shape[0]
	acc = torch.round(acc * 100)
	return acc

#========== START A LOOP
loss_values = []
for epoch in range(num_epochs:=501):
	running_train_loss = 0.0 
	running_accuracy = 0.0
	running_vall_loss = 0.0
	
	for X, y in train_dataloader:
		X, y = X.to(device), y.to(device)
		optimizer.zero_grad()
		pred = model(X)

		train_loss = loss_fn(pred, y)
		acc = binary_acc(pred, y)
		
		loss_values.append(train_loss.item())
		
		train_loss.backward()
		optimizer.step()
		
		running_train_loss += train_loss.item()
		running_accuracy += acc.item()
		train_loss_value = running_train_loss/len(train_dataloader)
	scheduler.step()	
	if epoch%40 == 0:
		print(f'Epoch: {epoch+0:03}, train_loss: {train_loss_value:.4f}, Acc: {running_accuracy/len(train_dataloader):.4f}, lr: {scheduler.get_last_lr()[0]:.4f}')
print("TRAINING COMPLETE ... ")

step = np.linspace(0, 100, len(loss_values))
fig, ax = plt.subplots(figsize=(8,8))
plt.plot(step, np.array(loss_values))
plt.title("Step-wise Loss")
plt.xlabel("Epochs")
plt.ylabel("Loss")
plt.savefig('output.png', dpi=300,bbox_inches='tight')

#========== TEST DATA ON A TRAINED MODEL
y_pred_list = []
model.eval()
with torch.no_grad():
	for X, y in test_dataloader:
		X = X.to(device)
		y_test_pred = model(X);
		y_test_pred = torch.sigmoid(y_test_pred)
		y_pred_list.append(y_test_pred)

y_pred_list = [a.squeeze().tolist() for a in y_pred_list]
print(len(y_pred_list)+len(X_train))

