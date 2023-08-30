import tensorflow as tf
from collections import defaultdict
import matplotlib.pyplot as plt

#################################################################################################
# LossHistory #
#################################################################################################
class LossHistory(tf.keras.callbacks.Callback):
    """ Records the history of the training per epoch and per batch """
    def on_train_begin(self, logs={}):
        self.epochs  = defaultdict(list) 
        self.batches = defaultdict(list) 
        self.pre_batch = 0

    def on_batch_end(self, batch, logs={}):
        self.batches['batch'].append(batch+self.pre_batch)
        for key,val in logs.items():
            self.batches[key].append(val)
        if 'lr' not in logs.keys():
            self.batches['lr'].append(tf.keras.backend.eval(self.model.optimizer.lr))

    def on_epoch_end(self, epoch, logs={}):
        self.epochs['epoch'].append(epoch)
        for key,val in logs.items():
            self.epochs[key].append(val)
        if 'lr' not in logs.keys():
            self.epochs['lr'].append(tf.keras.backend.eval(self.model.optimizer.lr))
        self.pre_batch = self.batches['batch'][-1] 

        

#################################################################################################
# PlotHistory #
#################################################################################################
def PlotHistory(history,params,outputName):
    """ Takes history from Keras training and makes loss plots (batch and epoch) and learning rate plots """
    #----- Figure -----#
    variables = sorted([key for key in history.epochs.keys() if 'val' not in key and 'val_'+key in history.epochs.keys()])
    variables += ["lr"]
    N = len(variables)
    fig, ax = plt.subplots(N,2,figsize=(12,N*2),sharex='col')
    plt.subplots_adjust(left    = 0.1,
                        right   = 0.6,
                        top     = 0.9,
                        bottom  = 0.1,
                        hspace  = 0.5,
                        wspace  = 0.4)

    #----- Batch Plots -----#
    for i,var in enumerate(variables):
        ax[i,0].plot(history.batches['batch'],history.batches[var],'k')
        ax[i,0].set_title(var)
        ax[i,0].set_xlabel('Batch')
        
    #----- Epoch Plots -----#
    for i,var in enumerate(variables):
        ax[i,1].plot(history.epochs['epoch'],history.epochs[var],'b')
        if 'val_'+var in history.epochs.keys():
            ax_twin = ax[i,1].twinx()
            ax_twin.plot(history.epochs['epoch'],history.epochs['val_'+var],'g')
            ax[i,1].set_ylabel("Training",color='b')
            ax[i,1].tick_params(axis='y', labelcolor='b')
            ax_twin.set_ylabel("Validation",color='g')
            ax_twin.tick_params(axis='y', labelcolor='g')
        ax[i,1].set_title(var)
        ax[i,1].set_xlabel('Epoch')

    #----- Print parameters -----#
    paramStr = "Parameters\n"
    for paramName in sorted(list(params.keys())):
        line = "- {} : ".format(paramName)
        if isinstance(params[paramName],(int,float,str)):
            value = str(params[paramName])
        else:
            value = params[paramName].__name__
        line += "{}\n".format(value)
        if len(line)>25:
            line = "{}:\n    {}".format(*line.split(':'))
        paramStr += line

    plt.gcf().text(0.7, 0.5, paramStr, fontsize=14)

    # Save #
    fig.savefig(outputName)
    print('Curves saved as %s'%outputName)
