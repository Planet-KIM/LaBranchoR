import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics

from tensorflow.keras import Sequential
from tensorflow.keras.layers import Dense, Conv2D, LSTM, Dropout, Reshape
from tensorflow.keras.layers import Bidirectional, TimeDistributed
from tensorflow.keras import regularizers, optimizers
from tensorflow.keras.models import load_model

%matplotlib inline

X_train, X_valid, X_test = np.load('./X_train.npy'), np.load('./X_valid.npy'), np.load('./X_test.npy')
y_train, y_valid, y_test = np.load('./Y_train.npy'), np.load('./Y_valid.npy'), np.load('./Y_test.npy')
print(X_train.shape, X_valid.shape, X_test.shape)
print(y_train.shape, y_valid.shape, y_test.shape)

def matching(preds, true):
    total, match = 0, 0
    for p, t in zip(preds, true):
        if np.argmax(p) in np.nonzero(t)[0]:
            match += 1
        total += 1
    return match, total

K = 3
counts = np.zeros((2*K+1, 4))
for target, seq in zip(y_train, X_train):
    for bp in np.nonzero(target)[0]:
        if 0 > bp-K or bp+K+1 > seq.shape[0]: continue
        counts = counts + seq[bp-K: bp+K+1]
print(counts.T)

class ModelTrainer:
    def __init__(self, model):
        self.model = model
        self.train_auc = []
        self.train_match = []
        self.valid_auc = []
        self.valid_match = []
        
    def train(self, X_train, X_valid,
              y_train, y_valid, PATIENCE = 15, EPOCHS = 30):
        print(self.model.summary())
        for i in range(EPOCHS):
            self.model.fit(X_train, y_train, epochs = 1, verbose = 0, batch_size = 8)
            self._evaluate(X_train, X_valid, y_train, y_valid)
            if (i > PATIENCE
                and max(self.valid_match[-PATIENCE:])
                < max(self.valid_match)):
                break
            print(i, self.valid_match[-1], self.train_match[-1])
        self._plot_scores()
        print(max(self.valid_match), max(self.valid_auc))
    
    def predict(X):
        return self.model.predict(X)
                
    def _evaluate(self, X_train, X_valid, y_train, y_valid):
        valid_preds = self.model.predict(X_valid, verbose=0)
        train_preds = self.model.predict(X_train, verbose=0)
        self.valid_match += [matching(valid_preds, y_valid)[0]
                             / float(y_valid.shape[0])]
        self.train_match += [matching(train_preds, y_train)[0]
                             / float(y_train.shape[0])]
        self.valid_auc += [metrics.roc_auc_score(y_valid.flatten(),
                                                 valid_preds.flatten())]
        self.train_auc += [metrics.roc_auc_score(y_train.flatten(),
                                                 train_preds.flatten())]
    
    def _plot_scores(self):
        plt.plot(self.valid_match, label = 'Validation')
        plt.plot(self.train_match, label = 'Training')
        plt.ylabel('Match')
        plt.xlabel('epoch')
        plt.legend()
        plt.show()
        
        plt.plot(self.valid_auc, label = 'Validation')
        plt.plot(self.train_auc, label = 'Training')
        plt.ylabel('auROC')
        plt.xlabel('epoch')
        plt.legend()
        plt.show()
        
new_model = load_model('./2layer.h5')
ModelTrainer(model=new_model).train(X_train[:, :, :4], X_valid[:, :, :4],
                          y_train.reshape(-1, 70, 1), y_valid.reshape(-1, 70, 1))