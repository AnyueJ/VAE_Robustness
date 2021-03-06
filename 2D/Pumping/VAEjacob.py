import tensorflow as tf
import scipy.io as sio

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import time
from matplotlib import pyplot as mp

from tensorflow.contrib.slim import fully_connected as fc
from tensorflow.python.ops.parallel_for.gradients import jacobian
input_dim = 10000
class VAE8(object):

    def __init__(self, learning_rate=1e-4, batch_size=64, n_z=35,l2c = 0):
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.n_z = n_z
        self.l2c = l2c
        tf.reset_default_graph()
        self.build()

        self.sess = tf.InteractiveSession()
        self.sess.run(tf.global_variables_initializer())

    # Build the netowrk and the loss functions
    def build(self):
        self.x = tf.placeholder(
            name='x', dtype=tf.float32, shape=[None, input_dim])
        self.learning_rate = tf.placeholder(name = 'learning_rate',dtype = tf.float64,shape = None)
        self.is_training = tf.placeholder(name = 'is_training',dtype = tf.bool,shape = None)
        inp = tf.reshape(self.x,[-1,100,100,1])
        
        conv1 = tf.contrib.layers.conv2d(inp,8,[3,3],stride=2,
                                 padding='SAME',activation_fn=tf.nn.relu, scope = 'conv1')
        conv2 = tf.contrib.layers.conv2d(conv1,16,[3,3],stride=2,
                                 padding='SAME',activation_fn=tf.nn.relu, scope = 'conv2')
        conv3 = tf.contrib.layers.conv2d(conv2,32,[3,3],stride=2,
                                 padding='VALID',activation_fn=tf.nn.relu, scope = 'conv3')
        conv4 = tf.contrib.layers.conv2d(conv3,64,[3,3],stride=2,
                                 padding='SAME',activation_fn=tf.nn.relu, scope = 'conv4')
        conv5 = tf.contrib.layers.conv2d(conv4,128,[3,3],stride=2,
                                 padding='SAME',activation_fn=tf.nn.relu, scope = 'conv5')
        f0 = tf.reshape(conv5,[-1,3*3*128])
        #f1 =fc(f0, 100, scope='fc1',activation_fn=tf.nn.relu)
        
        self.z_mu = fc(f0, self.n_z, scope='enc_mu', 
                       activation_fn=None)
        self.z_log_sigma_sq = fc(f0, self.n_z, scope='enc_sigma', 
                                 activation_fn=None)
        eps = tf.random_normal(
            shape=tf.shape(self.z_log_sigma_sq),
            mean=0, stddev=1, dtype=tf.float32)
        self.z = self.z_mu + tf.sqrt(tf.exp(self.z_log_sigma_sq)) * eps

        # Decode
        # z -> x_hat
        dfc1 = fc(self.z, 3*3*128, scope='dfc1',activation_fn=tf.nn.relu)
        #dfc2 = fc(dfc1, 6*6*10, scope='dfc2',activation_fn=tf.nn.relu)
        inp = tf.reshape(dfc1,[-1,3,3,128])
        dconv0 = tf.layers.conv2d_transpose(inp,64,3,2,"SAME",activation=tf.nn.relu, name= 'dconv0')
        dconv1 = tf.layers.conv2d_transpose(dconv0,32,3,2,"SAME",activation=tf.nn.relu, name= 'dconv1')
        print(dconv1.shape)
        
        dconv2 = tf.layers.conv2d_transpose(dconv1,16,3,2,"VALID", activation=tf.nn.relu, name = 'dconv2')
        print(dconv2.shape)
        
        
        dconv3 = tf.layers.conv2d_transpose(dconv2,8,3,2,"SAME",  activation=tf.nn.relu, name = 'dconv3')
        print(dconv3.shape)
        
        dconv4 = tf.layers.conv2d_transpose(dconv3,1,3,2,"SAME",activation=tf.nn.sigmoid, name = 'dconv4')
        print(dconv3.shape)
        
        self.x_hat = tf.reshape(dconv4,[-1,dconv4.shape[2]*dconv4.shape[1]])
        epsilon = 1e-5
        recon_loss = -tf.reduce_sum(
            self.x * tf.log(epsilon+self.x_hat) + 
            (1-self.x) * tf.log(epsilon+1-self.x_hat), 
            axis=1
        )
        
        self.recon_loss = tf.reduce_mean(recon_loss)
        self.l2w = sum(tf.nn.l2_loss(var) for var in tf.trainable_variables() if not 'biases' in var.name)
        self.l2b = sum(tf.nn.l2_loss(var) for var in tf.trainable_variables() if not 'weights' in var.name)
        # Latent loss
        # KL divergence: measure the difference between two distributions
        # Here we measure the divergence between 
        # the latent distribution and N(0, 1)
        latent_loss = -0.5 * tf.reduce_sum(
            1 + self.z_log_sigma_sq - tf.square(self.z_mu) - 
            tf.exp(self.z_log_sigma_sq), axis=1)
        self.latent_loss = self.n_z*tf.reduce_mean(latent_loss)

        self.total_loss = self.recon_loss + 3*self.latent_loss + self.l2c*(self.l2w+self.l2b)
        self.train_op = tf.train.AdamOptimizer(
            learning_rate=self.learning_rate).minimize(self.total_loss)
        
        self.losses = {
            'recon_loss': self.recon_loss,
            'latent_loss': self.latent_loss,
            'total_loss': self.total_loss,
            'l2': self.l2w
        }        
        return

    # Execute the forward and the backward pass
    def run_single_step(self, x,learning_rate):
        is_training = True
        _, losses = self.sess.run(
            [self.train_op, self.losses],
            feed_dict={self.x: x, self.learning_rate:learning_rate,self.is_training : is_training}
        )
        return losses

    # x -> x_hat
    def reconstructor(self, x):
        is_training = False
        x_hat = self.sess.run(self.x_hat, feed_dict={self.x: x,self.is_training : is_training})
        return x_hat

    # z -> x
    def generator(self, z):
        is_training = False
        x_hat = self.sess.run(self.x_hat, feed_dict={self.z: z,self.is_training : is_training})
        return x_hat
    
    # x -> z
    def transformer(self, x):
        is_training = False
        z = self.sess.run(self.z, feed_dict={self.x: x,self.is_training : is_training})
        return z
    def caljacob(self,z):
        is_training = False
        jacob = self.sess.run(jacobian(self.x_hat,self.z), feed_dict={self.z:z,self.is_training : is_training})
        return jacob
model = VAE8(n_z = 16)
saver = tf.train.Saver()
saver.restore(model.sess, "tf_models/slimmodels/VAE16.ckpt")


x = sio.loadmat('x.mat')
x = np.array(x["x"])

VAEjacob = np.squeeze(model.caljacob(x))
VAEjacob = VAEjacob*(0.6+2.395)
model.sess.close()
sio.savemat('VAEjacob',{'VAEjacob':VAEjacob})
