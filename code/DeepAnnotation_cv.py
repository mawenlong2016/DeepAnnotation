# -*- coding: UTF-8 -*-
# DeepAnnotation for cross-validation
# python ./code/DeepAnnotation_cv.py --training_steps 100 --early_stop 10 --kfold 2 --learning_rate 0.01 --genotype ./species/Duroc/Genotype_train_example.h5 --phenotype ./species/Duroc/phenotype_train_example.txt --function ./species/Duroc/GeneFunc_example.h5 --metagene ./species/Duroc/metagene_example.h5 --codingSNP ./species/Duroc/codingSNP_flag_example.txt --work_path ./species/Duroc/Genotype_example --seed 1 --verbose True --save_prcocess True --save_earlymodel True --phe_flag 5 --batch_ratio 0.1 --max_gpu 0.9 --decay_steps 100 --decay_rate 0.9 --seed 1  --architecture fnn,fnn,fnn,fnn,fnn,fnn,fnn --unit 90 --regularizer l2 --regularizer_rate 0.1 --momentum 0.95 --dropout 0.1 --bias_architecture none,none,none --unit_ratio 1,2,4,5,6,7 --calculate -1

# import packages
import tensorflow as tf
import h5py
tf.compat.v1.disable_eager_execution()
import os
import sys
sys.setrecursionlimit(10000)
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
import datetime

# define main function
# Run main() function when called directly
# define the DeepAnnotation class
# Class representation of DeepAnnotation network model
class DeepAnnotation(object):
    # Initialize model
    def __init__(self, function, metagene, X_train, y_train, X_valid, y_valid, X_test, y_test,
                 codingSNP, batch_size, decay_steps, decay_rate,training_steps, early_stop,
                 learning_rate, architecture, unit, unit_ratio, regularizer, regularizer_rate, 
                 momentum, dropout, bias_architecture, verbose, save_prcocess, save_earlymodel):
        self.function = tf.convert_to_tensor(function)
        self.metagene = tf.convert_to_tensor(metagene)
        self.X_train = X_train
        self.X_valid = X_valid
        self.y_train = y_train
        self.y_valid = y_valid
        self.X_test = X_test
        self.y_test = y_test
        self.codingSNP = codingSNP
        self.batch_size = batch_size
        self.features = X_train.shape[2]
        self.term = function.shape[1]
        self.link = metagene.shape[1]
        self.early_stop = early_stop
        self.lr_decay_step = decay_steps
        self.lr_decay_rate = decay_rate
        self.training_steps = training_steps
        self.learning_rate = learning_rate
        self.architecture = architecture
        self.unit = unit
        self.unit_ratio = unit_ratio
        self.regularizer = regularizer
        self.regularizer_rate = regularizer_rate
        self.momentum = momentum
        self.dropout = dropout
        self.bias_architecture = bias_architecture
        self.verbose = verbose
        self.save_prcocess = save_prcocess
        self.save_earlymodel = save_earlymodel
        # Initialize training, validation, testing datasets
        self.initialize_dataset()
        # Define tensor for updating global step
        self.global_step = tf.compat.v1.train.get_or_create_global_step()
        # Build graph for metagene model
        self.build_model()       

    # Initialize training, validation, testing datasets
    def initialize_dataset(self):
        # Define iterator for training dataset
        tep_position = np.arange(0,self.X_train.shape[0],5).tolist()
        tep_position.append(self.X_train.shape[0])
        i=0
        self.TRdataset = tf.compat.v1.data.Dataset.from_tensor_slices((self.X_train[tep_position[i]:tep_position[i+1],0:self.X_train.shape[1],0:self.X_train.shape[2]], self.y_train[tep_position[i]:tep_position[i+1]]))
        for i in range(len(tep_position)-1)[1:len(range(len(tep_position)-1))]:
            trdataset = tf.compat.v1.data.Dataset.from_tensor_slices((self.X_train[tep_position[i]:tep_position[i+1],0:self.X_train.shape[1],0:self.X_train.shape[2]], self.y_train[tep_position[i]:tep_position[i+1]]))
            self.TRdataset = self.TRdataset.concatenate(trdataset)
        
        self.TRdataset = self.TRdataset.shuffle(self.batch_size * 5)
        self.TRdataset = self.TRdataset.batch(self.batch_size)
        self.TRdataset = self.TRdataset.prefetch(self.batch_size * 5)
        self.TRdataset = self.TRdataset.repeat()
        self.TRdataset = tf.compat.v1.data.make_initializable_iterator(self.TRdataset)
        
        # Define iterator for validation dataset
        tep_position = np.arange(0,self.X_valid.shape[0],5).tolist()
        tep_position.append(self.X_valid.shape[0])
        i=0
        self.VAdataset = tf.compat.v1.data.Dataset.from_tensor_slices((self.X_valid[tep_position[i]:tep_position[i+1],0:self.X_valid.shape[1],0:self.X_valid.shape[2]], self.y_valid[tep_position[i]:tep_position[i+1]]))
        for i in range(len(tep_position)-1)[1:len(range(len(tep_position)-1))]:
            vadataset = tf.compat.v1.data.Dataset.from_tensor_slices((self.X_valid[tep_position[i]:tep_position[i+1],0:self.X_valid.shape[1],0:self.X_valid.shape[2]], self.y_valid[tep_position[i]:tep_position[i+1]]))
            self.VAdataset = self.VAdataset.concatenate(vadataset)
        
        self.VAdataset = self.VAdataset.shuffle(self.batch_size * 5)
        self.VAdataset = self.VAdataset.batch(self.batch_size)
        self.VAdataset = self.VAdataset.prefetch(self.batch_size * 5)
        self.VAdataset = self.VAdataset.repeat()
        self.VAdataset = tf.compat.v1.data.make_initializable_iterator(self.VAdataset)
                

    # Specify session for model evaluations
    def set_session(self, sess):
        self.sess = sess        

    # Specify saver for model restore
    def set_saver(self, saver):
        self.saver = saver

    # Framework of deep neural metagene (DNN)
    def DNN(self, x, function, metagene, training=True, reuse=None):     
        if self.regularizer == "l1" or self.regularizer == "L1":
            regular = tf.keras.regularizers.l1(self.regularizer_rate)
        elif self.regularizer == "l2" or self.regularizer == "L2":
            regular = tf.keras.regularizers.l2(self.regularizer_rate)
        elif self.regularizer == "l1_l2" or self.regularizer == "L1_L2":
            regular = tf.keras.regularizers.l1_l2(l1=self.regularizer_rate,l2=self.regularizer_rate) 
        elif  self.regularizer == "none" or self.regularizer == "None" or self.regularizer == "NONE":
            regular = None

        if self.bias_architecture[0] == "none" or self.bias_architecture[0] == "None" or self.bias_architecture[0] == "NONE":
            nn_bias = True
        else:
            nn_bias = False

        # First neural metagene (nn1) layer 1 for SNPs
        # Visible to receive the genome annotations
        if self.architecture[0] == "fnn" or self.architecture[0] == "Fnn" or self.architecture[0] == "FNN":
            nn1_coding = tf.compat.v1.layers.dense(tf.compat.v1.layers.flatten(x[:,0:self.codingSNP,0:1]), 
                                            units=self.unit * self.unit_ratio[0], 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn1_coding")
            nn1_noncoding = tf.compat.v1.layers.dense(tf.compat.v1.layers.flatten(x[:,self.codingSNP:,0:1]), 
                                            units=self.unit * self.unit_ratio[0], 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn1_noncoding")
        else:
            nn1_coding = tf.compat.v1.layers.conv2d(tf.reshape(x[:,0:self.codingSNP,0:1],[-1, 1, 1, self.codingSNP]), 
                                            filters=self.unit * self.unit_ratio[0], 
                                            kernel_size=1, strides=1, 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn1_coding")  
            
            nn1_noncoding = tf.compat.v1.layers.conv2d(tf.reshape(x[:,self.codingSNP:,0:1],[-1, 1, 1, x.shape[1]-self.codingSNP]), 
                                            filters=self.unit * self.unit_ratio[0], 
                                            kernel_size=1, strides=1, 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn1_noncoding") 
            nn1_coding = tf.reshape(nn1_coding,  [-1,nn1_coding.shape[3]])
            nn1_noncoding = tf.reshape(nn1_noncoding,  [-1,nn1_noncoding.shape[3]])
            
        # First functional annotation bias neural metagene (gnn1) layer 2 for functional SNPs
        # Visible to receive the epigenomic and secondary structure annotations
        if self.bias_architecture[0] == "fnn" or self.bias_architecture[0] == "Fnn" or self.bias_architecture[0] == "FNN":
            gnn1_coding = tf.compat.v1.layers.dense(tf.compat.v1.layers.flatten(x[:,0:self.codingSNP,1:]), 
                                            units=self.unit * self.unit_ratio[0], 
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn1_coding")
            gnn1_noncoding = tf.compat.v1.layers.dense(tf.compat.v1.layers.flatten(x[:,self.codingSNP:,1:]), 
                                            units=self.unit * self.unit_ratio[0], 
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn1_noncoding")
            snn1_coding = tf.keras.layers.concatenate([nn1_coding, gnn1_coding], name = "DNN_snn1_coding", dtype="float64")
            snn1_noncoding = tf.keras.layers.concatenate([nn1_noncoding, gnn1_noncoding], name = "DNN_snn1_noncoding", dtype="float64")
        elif self.bias_architecture[0] == "cnn" or  self.bias_architecture[0] == "Cnn" or self.bias_architecture[0] == "CNN":
            gnn1_coding = tf.compat.v1.layers.conv2d(tf.reshape(x[:,0:self.codingSNP,1:],[-1,x.shape[2]-1,1,self.codingSNP]), 
                                            filters=self.unit * self.unit_ratio[0],
                                            kernel_size=(x.shape[2]-1,1), strides=(x.shape[2]-1,1),  
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn1_coding")
            gnn1_noncoding = tf.compat.v1.layers.conv2d(tf.reshape(x[:,self.codingSNP:,1:],[-1,x.shape[2]-1,1,x.shape[1]-self.codingSNP]), 
                                            filters=self.unit * self.unit_ratio[0],
                                            kernel_size=(x.shape[2]-1,1), strides=(x.shape[2]-1,1),  
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn1_noncoding")
            gnn1_coding = tf.reshape(gnn1_coding,  [-1,gnn1_coding.shape[3]]) 
            gnn1_noncoding = tf.reshape(gnn1_noncoding,  [-1,gnn1_noncoding.shape[3]]) 
            snn1_coding = tf.keras.layers.concatenate([nn1_coding, gnn1_coding], name = "DNN_snn1_coding", dtype="float64")
            snn1_noncoding = tf.keras.layers.concatenate([nn1_noncoding, gnn1_noncoding], name = "DNN_snn1_noncoding", dtype="float64")
        else:
            snn1_coding = tf.identity(nn1_coding, name = "DNN_snn1_coding")
            snn1_noncoding = tf.identity(nn1_noncoding, name = "DNN_snn1_noncoding")

        nn1_coding_bn = tf.compat.v1.layers.batch_normalization(snn1_coding, momentum=self.momentum, epsilon=1e-5,
                                                         training=training, reuse=reuse, name="DNN_nn1_coding_bn")
        nn1_coding_ac = tf.compat.v1.nn.relu(nn1_coding_bn, name="DNN_nn1_coding_ac")
        nn1_coding_dp = tf.compat.v1.layers.dropout(nn1_coding_ac, rate=self.dropout, training=training, name="DNN_nn1_coding_dp")
        
        nn1_noncoding_bn = tf.compat.v1.layers.batch_normalization(snn1_noncoding, momentum=self.momentum, epsilon=1e-5,
                                                         training=training, reuse=reuse, name="DNN_nn1_noncoding_bn")
        nn1_noncoding_ac = tf.compat.v1.nn.relu(nn1_noncoding_bn, name="DNN_nn1_noncoding_ac")
        nn1_noncoding_dp = tf.compat.v1.layers.dropout(nn1_noncoding_ac, rate=self.dropout, training=training, name="DNN_nn1_noncoding_dp")
        
        if self.bias_architecture[1] == "none" or self.bias_architecture[1] == "None" or self.bias_architecture[1] == "NONE":
            nn_bias = True
        else:
            nn_bias = False

        # Second neural metagene (nn2) layer 3 for output of adjusted regulatory variants and adjusted functional variants
        if self.architecture[1] == "fnn" or self.architecture[1] == "Fnn" or self.architecture[1] == "FNN":            
            nn2_coding = tf.compat.v1.layers.dense(nn1_coding_dp, units=self.unit * self.unit_ratio[1], 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn2_coding") 
            nn2_noncoding = tf.compat.v1.layers.dense(nn1_noncoding_dp, units=self.unit * self.unit_ratio[1], 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn2_noncoding")  
        else:
            nn2_coding = tf.compat.v1.layers.conv2d(tf.reshape(nn1_coding_dp,[-1, 1, 1, nn1_coding_dp.shape[1]]), 
                                             filters=self.unit * self.unit_ratio[1], 
                                             kernel_size=1, strides=1, 
                                             use_bias=nn_bias, kernel_regularizer= regular,
                                             reuse=reuse, name="DNN_nn2_coding") 
            nn2_coding = tf.reshape(nn2_coding,  [-1,nn2_coding.shape[3]])
            nn2_noncoding = tf.compat.v1.layers.conv2d(tf.reshape(nn1_noncoding_dp,[-1, 1, 1, nn1_noncoding_dp.shape[1]]), 
                                             filters=self.unit * self.unit_ratio[1], 
                                             kernel_size=1, strides=1, 
                                             use_bias=nn_bias, kernel_regularizer= regular,
                                             reuse=reuse, name="DNN_nn2_noncoding")   
            nn2_noncoding = tf.reshape(nn2_noncoding,  [-1,nn2_noncoding.shape[3]])
            
        # Second functional annotation bias neural metagene (gnn2) layer 3
        # Visible to receive the transcriptomic and gene function annotations
        if self.bias_architecture[1] == "fnn" or self.bias_architecture[1] == "Fnn" or self.bias_architecture[1] == "FNN":
            gnn2 = tf.compat.v1.layers.dense(tf.reshape(function,[1,-1]), 
                                            units=self.unit * self.unit_ratio[1], 
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn2")   
            snn2 = tf.keras.layers.concatenate([nn2_coding, nn2_noncoding, tf.tile(gnn2,[tf.shape(x)[0],1])], name = "DNN_snn2", dtype="float64")
        elif self.bias_architecture[1] == "cnn" or self.bias_architecture[1] == "Cnn" or self.bias_architecture[1] == "CNN":
            gnn2 = tf.compat.v1.layers.conv2d(tf.reshape(function,[-1,function.shape[0],1,function.shape[1]]), 
                                            filters=self.unit * self.unit_ratio[1],
                                            kernel_size=(function.shape[0],1), strides=(function.shape[0],1), 
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn2")  
            gnn2 = tf.reshape(gnn2,  [-1,gnn2.shape[3]])
            snn2 = tf.keras.layers.concatenate([nn2_coding, nn2_noncoding, tf.tile(gnn2,[tf.shape(x)[0],1])], name = "DNN_snn2", dtype="float64")
        else:
            snn2 = tf.identity(tf.keras.layers.concatenate([nn2_coding, nn2_noncoding],dtype="float64"), name = "DNN_snn2")        
        
        nn2_bn = tf.compat.v1.layers.batch_normalization(snn2, momentum=self.momentum, epsilon=1e-5,
                                                         training=training, reuse=reuse, name="DNN_nn2_bn")
        nn2_ac = tf.compat.v1.nn.relu(nn2_bn, name="DNN_nn2_ac")
        nn2_dp = tf.compat.v1.layers.dropout(nn2_ac, rate=self.dropout, training=training, name="DNN_nn2_dp")
        
        if self.bias_architecture[2] == "none" or self.bias_architecture[2] == "None" or self.bias_architecture[2] == "NONE":
            nn_bias = True
        else:
            nn_bias = False

        if self.architecture[2] == "fnn" or self.architecture[2] == "Fnn" or self.architecture[2] == "FNN":
            nn3 = tf.compat.v1.layers.dense(nn2_dp, units=self.unit * self.unit_ratio[2], 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn3")    
        else:
            nn3 = tf.compat.v1.layers.conv2d(tf.reshape(nn2_dp,[-1, 1, 1, nn2_dp.shape[1]]), 
                                            filters=self.unit * self.unit_ratio[2], 
                                            kernel_size=1, strides=1, 
                                            use_bias=nn_bias, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn3")   
            nn3 = tf.reshape(nn3,  [-1,nn3.shape[3]])
            
        # Third functional annotation bias neural metagene (gnn3) layer 4
        # Visible to receive the regulatory modules annotations
        if self.bias_architecture[2] == "fnn" or self.bias_architecture[2] == "Fnn" or self.bias_architecture[2] == "FNN":            
            gnn3 = tf.compat.v1.layers.dense(tf.reshape(metagene,[1,-1]), units=self.unit * self.unit_ratio[2], 
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn3")  
            snn3 = tf.keras.layers.concatenate([nn3, tf.tile(gnn3,[tf.shape(x)[0],1])],axis=1, name = "DNN_snn3", dtype="float64")
        elif self.bias_architecture[2] == "cnn" or self.bias_architecture[2] == "Cnn" or self.bias_architecture[2] == "CNN":
            gnn3 = tf.compat.v1.layers.conv2d(tf.reshape(metagene,[-1,metagene.shape[0],1,metagene.shape[1]]), 
                                            filters=self.unit * self.unit_ratio[2],
                                            kernel_size=(metagene.shape[0],1), strides=(metagene.shape[0],1),  
                                            use_bias=False, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_gnn3")  
            gnn3 = tf.reshape(gnn3,  [-1,gnn3.shape[3]])
            snn3 = tf.keras.layers.concatenate([nn3, tf.tile(gnn3,[tf.shape(x)[0],1])],axis=1, name = "DNN_snn3", dtype="float64")
        else:
            snn3 = tf.identity(nn3, name = "DNN_snn3")        

        nn3_bn = tf.compat.v1.layers.batch_normalization(snn3, momentum=self.momentum, epsilon=1e-5,
                                                             training=training, reuse=reuse, name="DNN_nn3_bn")
        nn3_ac = tf.compat.v1.nn.relu(nn3_bn, name="DNN_nn3_ac")
        nn3_dp = tf.compat.v1.layers.dropout(nn3_ac, rate=self.dropout, training=training, name="DNN_nn3_dp")
        
        # Fourth neural metagene (nn4) layer 4 for outputs of high-order annotations
        if self.architecture[3] == "fnn" or self.architecture[3] == "Fnn" or self.architecture[3] == "FNN":
            # Third neural metagene (nn3) layer for higher information of expected gene regulatory metagene
            nn4 = tf.compat.v1.layers.dense(nn3_dp, units=self.unit * self.unit_ratio[3], 
                                            use_bias=True, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn4")    
        else:
            nn4 = tf.compat.v1.layers.conv2d(tf.reshape(nn3_dp,[-1, 1, 1, nn3_dp.shape[1]]), 
                                            filters=self.unit * self.unit_ratio[3], 
                                            kernel_size=1, strides=1, 
                                            use_bias=True, kernel_regularizer= regular,
                                            reuse=reuse, name="DNN_nn4")    
            nn4 = tf.reshape(nn4,  [-1,nn4.shape[3]])        

        nn4_bn = tf.compat.v1.layers.batch_normalization(nn4, momentum=self.momentum, epsilon=1e-5,
                                                         training=training, reuse=reuse, name="DNN_nn4_bn")
        nn4_ac = tf.compat.v1.nn.relu(nn4_bn, name="DNN_nn4_ac")
        nn4_dp = tf.compat.v1.layers.dropout(nn4_ac, rate=self.dropout, training=training, name="DNN_nn4_dp")

        # Fifth neural metagene (nn5) layer 5 for outputs of high-order annotations
        # Expected to be visible to receive the functional annotations of large proteins
        # Invisible due to unaccessible proteome data
        nn5 = tf.compat.v1.layers.dense(nn4_dp, units=self.unit * self.unit_ratio[4], 
                                        use_bias=True, kernel_regularizer= regular,
                                        reuse=reuse, name="DNN_nn5")
        nn5_bn = tf.compat.v1.layers.batch_normalization(nn5, momentum=self.momentum, epsilon=1e-5,
                                                         training=training, reuse=reuse, name="DNN_nn5_bn")
        nn5_ac = tf.compat.v1.nn.relu(nn5_bn, name="DNN_nn5_ac")
        nn5_dp = tf.compat.v1.layers.dropout(nn5_ac, rate=self.dropout, training=training, name="DNN_nn5_dp")

        # Sixth neural metagene (nn6) layer 6 for outputs of high-order annotations
        # Expected to be visible to receive the functional annotations of lsmall metabolites
        # Invisible due to unaccessible metabolome data
        nn6 = tf.compat.v1.layers.dense(nn5_dp, units=self.unit * self.unit_ratio[5], 
                                        use_bias=True, kernel_regularizer= regular,
                                        reuse=reuse, name="DNN_nn6")
        nn6_bn = tf.compat.v1.layers.batch_normalization(nn6, momentum=self.momentum, epsilon=1e-5,
                                                         training=training, reuse=reuse, name="DNN_nn6_bn")
        nn6_ac = tf.compat.v1.nn.relu(nn6_bn, name="DNN_nn6_ac")
        nn6_dp = tf.compat.v1.layers.dropout(nn6_ac, rate=self.dropout, training=training, name="DNN_nn6_dp")
        
        # Final lineary neural metagene (fnn) layer 7 for summarizing higher information to final phenotype
        fnn = tf.compat.v1.layers.dense(nn6_dp, units=1,
                                        use_bias=True, kernel_regularizer= regular,
                                        reuse=reuse, name="DNN_fnn") 
        return fnn

    def LOSS(self, y_predict, y_real):        
        loss1 = tf.reduce_mean(tf.square(y_predict - y_real))
        if self.regularizer == "none" or self.regularizer == "None" or self.regularizer == "NONE":
            loss2 = tf.convert_to_tensor(0, dtype=tf.float64)
        else:
            loss2 = tf.compat.v1.losses.get_regularization_loss()
        loss3 = loss1+0.00001*loss2
        return loss3, loss1, loss2

    # Evaluate model on specified batch of data
    def evaluate_model(self, X_data, function, metagene, y_data, training=True, reuse=None):
        # Predicting phenotype from DNN with bottom features extraction and top convergence features
        pred_y = self.DNN(X_data, function, metagene, training=training, reuse=reuse)
        # Compute loss between predicted and real phenotype adjusted by regularization
        loss, real_loss, l2_loss = self.LOSS(pred_y, y_data)
        return loss, real_loss, l2_loss

    # Define graph for model
    def build_model(self):
        # Define placeholder for dataset handle (to select training or validation)
        self.dataset_handle = tf.compat.v1.placeholder(tf.string, shape=[], name='dataset_handle')
        self.iterator = tf.compat.v1.data.Iterator.from_string_handle(self.dataset_handle, self.TRdataset.output_types,
                                                                      self.TRdataset.output_shapes)
        self.X_data, self.y_data = self.iterator.get_next()
        self.y_data = tf.reshape(self.y_data, [-1, 1]) 

        # Define learning rate with exponential decay
        #################################
        self.learning_rt = tf.compat.v1.train.exponential_decay(self.learning_rate, self.global_step,
                                                                self.lr_decay_step, self.lr_decay_rate)
        ##################################################
        # Define placeholder for training status
        self.training = tf.compat.v1.placeholder(tf.bool, name='training')

        # # Compute predictions and loss for training/validation datasets
        self.loss, self.real_loss, self.l2_loss = self.evaluate_model(self.X_data, self.function, self.metagene, self.y_data, training=self.training)

        # Create separate lists for discriminator and generator variables
        DNN_vars = [var for var in tf.compat.v1.trainable_variables() if 'DNN_' in var.name]
        with tf.control_dependencies(tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.UPDATE_OPS)):
            self.DNN_optim = tf.compat.v1.train.AdamOptimizer(self.learning_rt) \
                .minimize(self.loss, var_list=DNN_vars, global_step=self.global_step)

        # Final predicted phenotype
        self.pred_y = self.DNN(self.X_data, self.function, self.metagene, training=False, reuse=True)

    # Train model
    def train(self, work_path, cv):
        # Define directory for saving models
        if not os.path.exists(os.path.join(work_path, 'cv_' + str(cv) + '/predict/')):
            os.makedirs(os.path.join(work_path, 'cv_' + str(cv) + '/predict/'))
        self.sess.run(tf.compat.v1.global_variables_initializer())
        self.sess.run(self.TRdataset.initializer)
        self.sess.run(self.VAdataset.initializer)        
        # Get handles for training and validation datasets
        self.training_handle, self.validation_handle = self.sess.run([self.TRdataset.string_handle(),
                                                                      self.VAdataset.string_handle()])        
        if self.save_prcocess == "true" or self.save_prcocess == "True" or self.save_prcocess == "TRUE" or self.save_prcocess == "t" or self.save_prcocess == "T":
            f=open(os.path.join(work_path, 'information_training_process.txt'),"a+") 
        count_flag = 0
        count_step = -1        
        vd_loss_check = 99999
        # Iterate through training steps
        while count_flag < self.training_steps:
            count_flag = count_flag + 1
            # Update global step
            step = tf.compat.v1.train.global_step(self.sess, self.global_step)            
            # Feed values for parameters
            fd = {self.dataset_handle: self.training_handle, self.training: True}
            fd_valid = {self.dataset_handle: self.validation_handle, self.training: False}
            loss, real_loss, l2_loss,_ = self.sess.run([self.loss, self.real_loss, self.l2_loss, self.DNN_optim], feed_dict=fd)                         
            if step < (self.training_steps - 1):
                if self.verbose == "true" or self.verbose == "True" or self.verbose == "TRUE" or self.verbose == "t" or self.verbose == "T":
                    print(
                    "Predict training Step %d:  %.5f [d_loss] %.5f [real_loss] %.5f [l_loss] "
                    % (step+1, loss, real_loss, l2_loss))
                if self.save_prcocess == "true" or self.save_prcocess == "True" or self.save_prcocess == "TRUE" or self.save_prcocess == "t" or self.save_prcocess == "T":
                    print(
                    "Predict training Step %d:  %.5f [d_loss] %.5f [real_loss] %.5f [l_loss] "
                    % (step+1, loss, real_loss, l2_loss), file=f) 
                vd_loss, vd_real_loss, vd_l2_loss= self.sess.run([self.loss, self.real_loss, self.l2_loss], feed_dict=fd_valid)
                if self.verbose == "true" or self.verbose == "True" or self.verbose == "TRUE" or self.verbose == "t" or self.verbose == "T":
                    print(
                        "Predict validation Step %d:  %.5f [d_loss] %.5f [real_loss] %.5f [l_loss] "
                        % (step+1, vd_loss, vd_real_loss, vd_l2_loss))
                if self.save_prcocess == "true" or self.save_prcocess == "True" or self.save_prcocess == "TRUE" or self.save_prcocess == "t" or self.save_prcocess == "T":
                    print(
                        "Predict validation Step %d:  %.5f [d_loss] %.5f [real_loss] %.5f [l_loss] "
                        % (step+1, vd_loss, vd_real_loss, vd_l2_loss), file=f)                
                if (vd_loss_check >= vd_loss) & (vd_loss > 0) :
                    vd_loss_check = vd_loss
                    count_step = 0   
                    if self.save_earlymodel == "true" or self.save_earlymodel == "True" or self.save_earlymodel == "TRUE" or self.save_earlymodel == "t" or self.save_earlymodel == "T":
                        # save early stopped model         
                        self.saver.save(sess=self.sess, save_path=os.path.join(work_path, 'cv_' + str(cv) + '/predict/best_model'),write_meta_graph=False)       
                else:
                    count_step = count_step + 1

                if (count_step == self.early_stop):
                    if self.verbose == "true" or self.verbose == "True" or self.verbose == "TRUE" or self.verbose == "t" or self.verbose == "T":
                        print("No improvement found within %d iterators, stopping optimization at %d steps."
                        % (self.early_stop, step+1))                    
                    if self.save_prcocess == "true" or self.save_prcocess == "True" or self.save_prcocess == "TRUE" or self.save_prcocess == "t" or self.save_prcocess == "T":
                        print("No improvement found within %d iterators, stopping optimization at %d steps."
                        % (self.early_stop, step+1), file=f)  
                    break
            else:
                if self.verbose == "true" or self.verbose == "True" or self.verbose == "TRUE" or self.verbose == "t" or self.verbose == "T":
                    print("The improvement was not increased within [%d] steps in total [%d] traning steps, a larger training step may need."
                        % (self.early_stop, self.training_steps))
                if self.save_prcocess == "true" or self.save_prcocess == "True" or self.save_prcocess == "TRUE" or self.save_prcocess == "t" or self.save_prcocess == "T":
                    print("The improvement was not increased within [%d] steps in total [%d] traning steps, a larger training step may need."
                        % (self.early_stop, self.training_steps), file=f)
        # save final model
        self.saver.save(sess=self.sess, save_path=os.path.join(work_path, 'cv_' + str(cv) + '/predict/final_model'),write_meta_graph=False)        
                 
    # Predict phenotype
    def predict_real(self, X_test):
        fd = {self.X_data: X_test}        
        real_predict = self.sess.run([self.pred_y], feed_dict=fd)        
        real_predict = np.array(real_predict)        
        real_predict = real_predict.reshape([real_predict.shape[1],1])
        return real_predict
   
def loadData_DeepAnnotation():
    '''
    the following parameters should be input：
    genotype: the SNP matrix
    phenotype: the interested trait values    
    ...other parameters continue
    '''
    tf.compat.v1.flags.DEFINE_string("genotype", "./data/genotype/Genotype_train.h5", 
                                     "genotype file")
    tf.compat.v1.flags.DEFINE_string("function", "./data/function/Function.h5", 
                                     "function annotation file")
    tf.compat.v1.flags.DEFINE_string("metagene", "./data/metagene/Network.h5", 
                                     "metagene annotation file")
    tf.compat.v1.flags.DEFINE_string("phenotype", "./data/phenotype/phenotype_train.txt",
                                     "phenotype file")
    tf.compat.v1.flags.DEFINE_string("codingSNP", "./data/variant/codingSNP_flag.txt", "coding SNP flag")
    tf.compat.v1.flags.DEFINE_integer("kfold", 5, "k-fold cross validation")
    tf.compat.v1.flags.DEFINE_integer("seed", 1, "set the random seed")
    tf.compat.v1.flags.DEFINE_float("batch_ratio", 0.1, "batch size ratio of the training set")
    tf.compat.v1.flags.DEFINE_float("max_gpu", 0.9, "set maximum gpu memory")
    tf.compat.v1.flags.DEFINE_string("work_path", "./DeepAnnotation", "work directory")    
    tf.compat.v1.flags.DEFINE_integer("decay_steps", 100, "set decay step")
    tf.compat.v1.flags.DEFINE_float("decay_rate", 0.9, "set decay rate")
    tf.compat.v1.flags.DEFINE_integer("training_steps", 100, "set total training step")
    tf.compat.v1.flags.DEFINE_integer("early_stop", 10, "set early stop steps")
    tf.compat.v1.flags.DEFINE_float("learning_rate", 0.01, 
                                     "learning rate: [0.1, 0.01, 0.001, 0.0001, 0.00001]")
    tf.compat.v1.flags.DEFINE_string("architecture", "cnn,fnn,fnn,fnn,fnn,fnn,fnn", 
                                     "base framework: [cnn, fnn]")
    tf.compat.v1.flags.DEFINE_integer("unit", 160, 
                                     "base feature number: [10, 30, 50, 100, 150, 250, 400, 500] followed by 2:3:6:3:2:1")
    tf.compat.v1.flags.DEFINE_string("unit_ratio", "2,3,6,6,3,2", 
                                     "unit ratio of each layer: [2:3:6:6:3:2, 1:2:3:3:2:1, 1:2:3:4:5:6]")
    tf.compat.v1.flags.DEFINE_string("regularizer", "l2", 
                                     "regularizer algorithm: [l1, l2, l1_l2]")
    tf.compat.v1.flags.DEFINE_float("regularizer_rate", 0.0001, 
                                     "regularizer rate: [0.1, 0.01, 0.001, 0.0001, 0.00001]")
    tf.compat.v1.flags.DEFINE_float("momentum", 0.9, 
                                     "batch normalization momentum: [1, 0.9, 0.8, 0.7, 0.6, 0.5]")
    tf.compat.v1.flags.DEFINE_float("dropout", 0.2, 
                                     "dropout rate: [0.3, 0.2, 0.1]")
    tf.compat.v1.flags.DEFINE_string("bias_architecture", "cnn,fnn,fnn", 
                                     "functional annotation bias framework: [cnn, fnn, none]")
    tf.compat.v1.flags.DEFINE_string("calculate", "0", 
                                     "calculate by cpu or gpu: [-1, 0, 1]")
    tf.compat.v1.flags.DEFINE_string("verbose", "True", 
                                     "training process displayed or not: [True, False]")
    tf.compat.v1.flags.DEFINE_integer("phe_flag", 5, "phenotype flag")
    tf.compat.v1.flags.DEFINE_string("save_prcocess", "True", 
                                     "save training process or not: [True, False]")
    tf.compat.v1.flags.DEFINE_string("save_earlymodel", "True", 
                                     "save early training model or not: [True, False]")
    tf.compat.v1.flags.DEFINE_string("output", "pred_phenotype.txt", 
                                     "the output name of phenotype prediction")
    FLAGS = tf.compat.v1.flags.FLAGS
    bias_architecture = str.split(FLAGS.bias_architecture,",")     
    f_genotye = h5py.File(FLAGS.genotype,'r')
    genotype = f_genotye['Genotype']
    if bias_architecture[1] != "none":
        f_function = h5py.File(FLAGS.function,'r')
        function = f_function['Function']
    else:
        function = np.array([[],[]])
    
    if bias_architecture[2] != "none":
        f_metagene = h5py.File(FLAGS.metagene,'r')
        metagene = f_metagene['Network']
    else:
        metagene = np.array([[],[]])
    phenotype = np.loadtxt(FLAGS.phenotype)
    codingSNP = np.int_(np.loadtxt(FLAGS.codingSNP))
    kfold = FLAGS.kfold
    seed = FLAGS.seed
    batch_ratio = FLAGS.batch_ratio
    max_gpu = FLAGS.max_gpu
    work_path = FLAGS.work_path    
    decay_steps = FLAGS.decay_steps
    decay_rate = FLAGS.decay_rate
    training_steps = FLAGS.training_steps
    early_stop = FLAGS.early_stop
    learning_rate = FLAGS.learning_rate
    architecture = str.split(FLAGS.architecture,",")
    unit = FLAGS.unit
    unit_ratio = np.float_(str.split(FLAGS.unit_ratio,","))
    regularizer = FLAGS.regularizer
    regularizer_rate = FLAGS.regularizer_rate
    momentum = FLAGS.momentum
    dropout = FLAGS.dropout
    calculate = FLAGS.calculate
    verbose = FLAGS.verbose
    phe_flag = FLAGS.phe_flag
    save_prcocess = FLAGS.save_prcocess
    save_earlymodel = FLAGS.save_earlymodel
    output = FLAGS.output
    hyper_parameters = str(learning_rate)+"_"+str(architecture)+"_"+str(unit)+"_"+str(unit_ratio) \
                       +"_"+regularizer+"_"+ str(regularizer_rate)+"_"+str(momentum)+"_" \
                       +str(dropout)+"_"+str(bias_architecture)
    if not os.path.exists(work_path):
        os.makedirs(work_path)
    return genotype, function, metagene, phenotype, codingSNP, kfold, seed, batch_ratio, \
           max_gpu, work_path, decay_steps, decay_rate, training_steps, early_stop, \
           learning_rate, architecture, unit, unit_ratio, regularizer, regularizer_rate, \
           momentum, dropout, bias_architecture, hyper_parameters, calculate,  \
           verbose, phe_flag, save_prcocess, save_earlymodel, output

# define the main() function
# Initialize and train DeepAnnotation deepannotation
def main():
    start1=datetime.datetime.now()
    # Define DeepAnnotation deepannotation parameters and options in dictionary of flags
    genotype, function, metagene, phenotype, codingSNP, kfold, seed, batch_ratio, \
    max_gpu, work_path, decay_steps, decay_rate, training_steps, early_stop, \
    learning_rate, architecture, unit, unit_ratio, regularizer, regularizer_rate, \
    momentum, dropout, bias_architecture, hyper_parameters, calculate, \
    verbose, phe_flag, save_prcocess, save_earlymodel, output \
    = loadData_DeepAnnotation()
    if verbose == "true" or verbose == "True" or verbose == "TRUE" or verbose == "t" or verbose == "T":
        print("Evaluation begins!")
        print("#################Architecture of deep learning framework######################")
        print("totally, there were: ",len(architecture)," base framework, and: ", len(bias_architecture), " functional annotation bias framework")
        print("learning rate: ",learning_rate)
        print("base architecture: ",architecture)
        print("base feature number: ", unit)
        print("feature number ratio: ", unit_ratio)
        print("regularizer algorithm: ", regularizer)
        print("regularizer rate: ", regularizer_rate)
        print("batch normalization momentum: ", momentum)
        print("dropout rate: ", dropout)
        print("functional annotation bias framework: ", bias_architecture)
        print("hyper_parameters: ", hyper_parameters)
        if calculate == "-1":
            print("calculated by cpu")
        else:
            print("calculated by gpu: ", calculate)        
        print("#################Architecture of deep learning framework######################")
    os.environ['CUDA_VISIBLE_DEVICES'] = calculate
    if save_prcocess == "true" or save_prcocess == "True" or save_prcocess == "TRUE" or save_prcocess == "t" or save_prcocess == "T":
        f=open(os.path.join(work_path, 'information_training_process.txt'),"a+")
        print("Evaluation begins!", file=f)
        print("#################Architecture of deep learning framework######################", file=f)
        print("totally, there were: ",len(architecture)," base framework, and: ", len(bias_architecture), " functional annotation bias framework", file=f)
        print("learning rate: ",learning_rate, file=f)
        print("base architecture: ",architecture, file=f)
        print("base feature number: ", unit, file=f)
        print("feature number ratio: ", unit_ratio, file=f)
        print("regularizer algorithm: ", regularizer, file=f)
        print("regularizer rate: ", regularizer_rate, file=f)
        print("batch normalization momentum: ", momentum, file=f)
        print("dropout rate: ", dropout, file=f)
        print("functional annotation bias framework: ", bias_architecture, file=f)
        print("hyper_parameters: ", hyper_parameters, file=f)
        if calculate == "-1":
            print("calculated by cpu", file=f)
        else:
            print("calculated by gpu: ", calculate, file=f)        
        print("#################Architecture of deep learning framework######################", file=f)
        print("#################Prediction performance of each cross validation######################", file=f)
        f.close()
    kfoldIndex = KFold(n_splits=kfold, shuffle=True, random_state=seed)
    TrainFlag = []
    TestFlag = []
    ValidFlag = []
    randomseed = seed
    cv = 1    
    pred_phenotype = np.zeros((phenotype.shape[0],2))
    pred_phenotype[:,0] = phenotype[:,phe_flag]    
    for trainIdx, test in kfoldIndex.split(genotype, phenotype):
        start=datetime.datetime.now()
        X_trainAll, y_trainAll, X_test, y_test \
            = genotype[trainIdx,:,:], phenotype[trainIdx, phe_flag], genotype[test,:,:], phenotype[test, phe_flag]
        X_train, X_valid, y_train, y_valid, train, valid \
            = train_test_split(X_trainAll, y_trainAll, trainIdx, test_size=1 / kfold, random_state=randomseed)
        randomseed = randomseed + 1
        TrainFlag.append(train)
        TestFlag.append(test)
        ValidFlag.append(valid)
        batch_size = int(len(X_train) * batch_ratio) + 1
        # Initialize model
        deepannotation = DeepAnnotation(function, metagene, X_train, y_train, X_valid, y_valid, X_test, y_test,
                        codingSNP, batch_size, decay_steps, decay_rate,training_steps, early_stop,
                        learning_rate, architecture, unit, unit_ratio, regularizer, regularizer_rate,
                         momentum, dropout, bias_architecture, verbose, save_prcocess, save_earlymodel)
        gpu_options = tf.compat.v1.GPUOptions(per_process_gpu_memory_fraction=max_gpu, allow_growth=True)
        config = tf.compat.v1.ConfigProto(gpu_options=gpu_options, allow_soft_placement=True)
        saver=tf.compat.v1.train.Saver()
        deepannotation.set_saver(saver)
        with tf.compat.v1.Session(config=config) as sess:    
            # Set model session
            deepannotation.set_session(sess)            
            # Train model
            deepannotation.train(work_path, cv)
            real_temp = deepannotation.predict_real(X_test)
            pred_phenotype[test,1] = real_temp[:,0]
            if verbose == "true" or verbose == "True" or verbose == "TRUE" or verbose == "t" or verbose == "T":
                print("predict correlation of [CV %d] is: %0.3f" %(cv, np.corrcoef(pred_phenotype[test,0],pred_phenotype[test,1])[0,1]))
            if save_prcocess == "true" or save_prcocess == "True" or save_prcocess == "TRUE" or save_prcocess == "t" or save_prcocess == "T":
                f=open(os.path.join(work_path, 'information_training_process.txt'),"a+")
                print("predict correlation of [CV %d] is: %0.3f" %(cv, np.corrcoef(pred_phenotype[test,0],pred_phenotype[test,1])[0,1]), file=f)
                f.close()
        tf.compat.v1.reset_default_graph()
        end=datetime.datetime.now()     
        if verbose == "true" or verbose == "True" or verbose == "TRUE" or verbose == "t" or verbose == "T":
            print("\n[ CV %d: trainin complete, time elapsed: %s ]\n" % (cv, end-start))         
        if save_prcocess == "true" or save_prcocess == "True" or save_prcocess == "TRUE" or save_prcocess == "t" or save_prcocess == "T":
            f=open(os.path.join(work_path, 'information_training_process.txt'),"a+")
            print("[CV %d]: trainin complete, time elapsed: %s" % (cv, end-start), file = f)
            f.close()
        cv = cv + 1    
    if verbose == "true" or verbose == "True" or verbose == "TRUE" or verbose == "t" or verbose == "T":
        print("final predict correlation of total is: %0.3f" %(np.corrcoef(pred_phenotype[:,0],pred_phenotype[:,1])[0,1]))
    if save_prcocess == "true" or save_prcocess == "True" or save_prcocess == "TRUE" or save_prcocess == "t" or save_prcocess == "T":
        f=open(os.path.join(work_path, 'information_training_process.txt'),"a+")
        print("final predict correlation of total is: %0.3f" %(np.corrcoef(pred_phenotype[:,0],pred_phenotype[:,1])[0,1]), file=f)
    np.savetxt(os.path.join(work_path,output), pred_phenotype, fmt='%0.3f')
    end1=datetime.datetime.now()
    if verbose == "true" or verbose == "True" or verbose == "TRUE" or verbose == "t" or verbose == "T":
        print("[ Training was complete, total time elapsed: %s ]\n" % (end1-start1))
    if save_prcocess == "true" or save_prcocess == "True" or save_prcocess == "TRUE" or save_prcocess == "t" or save_prcocess == "T":
        print("Training was complete, total time elapsed: %s" % (end1-start1), file = f)
        f.close()
    if verbose == "true" or verbose == "True" or verbose == "TRUE" or verbose == "t" or verbose == "T":
        print("#################Model training ends######################")
    if not os.path.exists(os.path.join(work_path,"TrainFlag.txt")):
        f = open(os.path.join(work_path,"TrainFlag-"+np.str(kfold)+".txt"),"w")
        for i in range(len(TrainFlag)):
            temp = TrainFlag[i]    
            for j in temp:
                print(j, end=',',file=f)
            print("",file=f)
        f.close()
        f = open(os.path.join(work_path,"TestFlag-"+np.str(kfold)+".txt"),"w")
        for i in range(len(TestFlag)):
            temp = TestFlag[i]    
            for j in temp:
                print(j, end=',',file=f)
            print("",file=f)
        f.close()
        f = open(os.path.join(work_path,"ValidFlag-"+np.str(kfold)+".txt"),"w")
        for i in range(len(ValidFlag)):
            temp = ValidFlag[i]    
            for j in temp:
                print(j, end=',',file=f)
            print("",file=f)
        f.close()

if __name__ == '__main__':
    main()


