import numpy as np
import pickle
import torch
import torch.nn as nn


## Part 1 - Model to learn the Q-function
# This part defines a simple model that learns a mapping between the canvas and the Q-values associated to each possible movements. 
# So it's a model with a size corresponding to the canvas size, and four outputs. 

class NLinearModels(nn.Module):
    def __init__(self, x_example, number_of_regressors=4):
        super(NLinearModels, self).__init__()
        in_features = x_example.reshape(-1).shape[0]
        self.linear = nn.Linear(in_features, number_of_regressors)
    
    def forward(self, x):
        x = x.reshape(x.shape[0], -1)
        return self.linear(x)

    def load(self):
        self.load_state_dict(torch.load('save_rl/weights.pt'))

    def save(self):
        torch.save(self.state_dict(), 'save_rl/weights.pt')


def train_on_batch(model, inputs, targets, criterion, optimizer):
    ### Simple helper function to train the model given a batch of inputs and targets, optimizes the model and returns the loss
    targets = torch.FloatTensor(targets)
    # Zero the parameter gradients
    optimizer.zero_grad()
    # Forward + backward + optimize
    outputs = model(inputs)
    loss = criterion(outputs, targets)
    loss.backward()
    optimizer.step()
    return loss.item()


## Part 2 - Experience Replay
## This part has to be read and understood in order to code the main.py file. 

class ExperienceReplay(object):
    """
    During gameplay all experiences < s, a, r, s’ > are stored in a replay memory. 
    During training, batches of randomly drawn experiences are used to generate the input and target for training.
    """
    def __init__(self, max_memory=100, discount=.9):
        """
        Setup
        max_memory: the maximum number of experiences we want to store
        memory: a list of experiences
        discount: the discount factor for future experience
        
        In the memory the information whether the game ended at the experience is stored seperately in a nested array
        [...,
        [experience, game_over],
        [experience, game_over],
        ...]
        """
        self.max_memory = max_memory
        self.memory = list()
        self.discount = discount

    def remember(self, experience, game_over):
        # Save an experience to memory
        self.memory.append([experience, game_over])
        # We don't want to store infinite memories, so if we have too many, we just delete the oldest one
        if len(self.memory) > self.max_memory:
            del self.memory[0]

    def get_batch(self, model, batch_size=10):    
        # How many experiences do we have?
        len_memory = len(self.memory)
        
        # Number of actions that can possibly be taken in the game (up, down, left, right)
        num_actions = 4
        
        # Dimensions of the game field
        env_dim = list(self.memory[0][0][0].shape)
        env_dim[0] = min(len_memory, batch_size)
        
        
        # We want to return an input and target vector with inputs from an observed state...
        inputs = torch.zeros(env_dim)
        #...and the target r + gamma * max Q(s’,a’)
        # Note that our target is a matrix, with possible fields not only for the action taken but also
        # for the other possible actions. The actions do not take the same values as the prediction to not affect them.
        Q = torch.zeros((inputs.shape[0], num_actions))
        
        # We randomly draw experiences to learn from
        for i, idx in enumerate(np.random.randint(0, len_memory,
                                                  size=inputs.shape[0])):
            """
            Here we load one transition <s, a, r, s’> from memory
            state_t: initial state s
            action_t: action taken a
            reward_t: reward earned r
            state_tp1: the state that followed s’
            """
#            idx = -1
            state, action_t, reward_t, state_tp1 = self.memory[idx][0]
            # We also need to know whether the game ended at this state
            game_over = self.memory[idx][1]

            # Add the state s to the input
            inputs[i:i+1] = state
            # First, we fill the target values with the predictions of the model.
            # They will not be affected by training (since the training loss for them is 0).
            model.eval()
            with torch.no_grad():
                Q[i] = model(state.unsqueeze(dim=0))[0]

                # If the game ended, the expected reward Q(s,a) should be the final reward r.
                # Otherwise the target value is r + gamma * max Q(s’,a’)
                
                # If the game ended, the reward is the final reward
                if game_over:  # if game_over is True
                    Q[i, action_t] = reward_t
                else:
                    # r + gamma * max Q(s’,a’)
                    next_round = model(state_tp1.unsqueeze(dim=0))[0]
                    Q[i, action_t] = reward_t + self.discount * torch.max(next_round)
        return inputs, Q

    def load(self):
        self.memory = pickle.load(open("save_rl/memory.pkl","rb"))
    
    def save(self):
        pickle.dump(self.memory,open("save_rl/memory.pkl","wb"))
