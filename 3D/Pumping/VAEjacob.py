import tensorflow as tf
import scipy.io as sio

import numpy as np

import matplotlib.pyplot as plt
import time

from tensorflow.contrib.slim import fully_connected as fc
from tensorflow.python.ops.parallel_for.gradients import jacobian
input_dim = 192000
class VAE(object):

    def __init__(self, n_z=35, conv1n=4, wrecon=1, wlatent=1, l2c = 0):
        self.conv1n = conv1n
        self.n_z = n_z
        self.l2c = l2c
        self.wrecon = wrecon
        self.wlatent = wlatent
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
        inp = tf.reshape(self.x,[-1,40,120,40,1])
        
        conv1 = tf.contrib.layers.conv3d(inp,self.conv1n,3,stride=2,
                                 padding='SAME',activation_fn=tf.nn.relu, scope = 'conv1')
        print(conv1.shape)
        conv2 = tf.contrib.layers.conv3d(conv1,self.conv1n*2,3,stride=2,
                                 padding='SAME',activation_fn=tf.nn.relu, scope = 'conv2')
        print(conv2.shape)
        conv3 = tf.contrib.layers.conv3d(conv2,self.conv1n*4,3,stride=2,
                                 padding='same',activation_fn=tf.nn.relu, scope = 'conv3')
        print(conv3.shape)
        conv4 = tf.contrib.layers.conv3d(conv3,self.conv1n*8,3,stride=2,
                                 padding='valid',activation_fn=tf.nn.relu, scope = 'conv4')
        print(conv4.shape)
        f0 = tf.reshape(conv4,[-1,2*7*2*8*self.conv1n])
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
        dfc1 = fc(self.z, 2*7*2*8*self.conv1n, scope='dfc1',activation_fn=tf.nn.relu)
        #dfc2 = fc(dfc1, 6*6*10, scope='dfc2',activation_fn=tf.nn.relu)
        inp = tf.reshape(dfc1,[-1,2,7,2,64])
        print(inp.shape)
        dconv1 = tf.layers.conv3d_transpose(inp,self.conv1n*4,3,2,"valid",activation=tf.nn.relu, name= 'dconv1')
        print(dconv1.shape)
        
        dconv2 = tf.layers.conv3d_transpose(dconv1,self.conv1n*2,3,2,"same", activation=tf.nn.relu, name = 'dconv2')
        print(dconv2.shape)
        
        
        dconv3 = tf.layers.conv3d_transpose(dconv2,self.conv1n,3,2,"SAME",  activation=tf.nn.relu, name = 'dconv3')
        print(dconv3.shape)
        
        dconv4 = tf.layers.conv3d_transpose(dconv3,1,3,2,"SAME",activation=tf.nn.sigmoid, name = 'dconv4')
        print(dconv4.shape)
        
        self.x_hat = tf.reshape(dconv4,[-1,192000])
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

        self.total_loss = self.wrecon*self.recon_loss + self.wlatent*self.latent_loss + self.l2c*(self.l2w+self.l2b)
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
        jacob = self.sess.run(jacobian(self.x_hat,self.z), feed_dict={self.z:z})
        return jacob
    def calgrad(self,z,wellloc):
        is_training = False
        gradient = self.sess.run([tf.gradients(self.x_hat[:,n-1],(self.z))[0] for n in wellloc], feed_dict={self.z:z,self.is_training : is_training})
        return gradient
model = VAE(n_z = 32, conv1n = 8, wrecon = 1, wlatent = 1, l2c = 0.05)
saver = tf.train.Saver()
saver.restore(model.sess, "tf_models/VAE32.ckpt")



x = sio.loadmat('x.mat')
x = np.array(x["x"])

VAEjacob = np.squeeze(model.caljacob(x))
VAEjacob = VAEjacob*(-9)
model.sess.close()
sio.savemat('VAEjacob',{'VAEjacob':VAEjacob})
