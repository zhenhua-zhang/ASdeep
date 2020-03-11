"""Deep learning predictor module."""
import sys
import math

import numpy as np

try:
    from ased.zutils import split_train_test
except ImportError as ime:
    print("ImportError: {}".format(ime))


try:
    import tensorflow as tf
    from tensorflow import keras
    from tensorflow.keras import layers
    from tensorflow.keras import models
    from tensorflow.keras import losses
    from tensorflow.keras import metrics
    from tensorflow.keras import optimizers
except ImportError as ime:
    print("ImportError: {}".format(ime))
    print("[W]: It looks you miss TesforFlow or related dependencies, but you",
          "can still use the quantification function of ASE effect in this",
          "package, however, it's not possible to use the prediction fucntion",
          "which is based on the TensorFlow.")


class CNNModel(models.Model):
    """A CNN class."""
    def __init__(self, input_shape=(98, 8)):
        super(CNNModel, self).__init__(name='CNNModel')

        self.flatten = layers.Flatten()

        length = width = int(math.sqrt(input_shape[0] * input_shape[1]))
        self.reshape = layers.Reshape((length, width, 1))

        self.conv1 = layers.Conv2D(64, 3, activation="relu")
        self.bn1 = layers.BatchNormalization()
        self.mp1 = layers.MaxPooling2D((2, 2))

        self.conv2 = layers.Conv2D(32, 3, activation="relu")
        self.bn2 = layers.BatchNormalization()
        self.mp2 = layers.MaxPooling2D((2, 2))

        self.conv3 = layers.Conv2D(16, 3, activation="relu")
        self.bn3 = layers.BatchNormalization()
        self.mp3 = layers.MaxPooling2D((2, 2))

        self.d1 = layers.Dense(8, activation='relu')
        self.d2 = layers.Dense(3, activation='relu')

    def call(self, inputs):
        x = self.flatten(inputs)
        x = self.reshape(x)

        x = self.conv1(x)
        x = self.bn1(x)
        x = self.mp1(x)

        x = self.conv2(x)
        x = self.bn2(x)
        x = self.mp2(x)

        x = self.conv3(x)
        x = self.bn3(x)
        x = self.mp3(x)

        x = self.d1(x)
        return self.d2(x)


class DLPFactory:
    """A factory to train deep learning model on ASE quantification results.
    """
    model = None
    iptmtx_shape = None
    loss_obj, optim, train_loss, train_acc, test_loss, test_acc = [None] * 6

    def __init__(self):
        super(DLPFactory, self).__init__()
        self.train_dtst = None
        self.test_dtst = None

    def init(self, model=None, loss_obj=None, optim=None, train_loss=None,
             train_acc=None, test_loss=None, test_acc=None):
        """Init the factory.
        """
        self.model = model

        if loss_obj:
            self.loss_obj = loss_obj
        else:
            self.loss_obj = losses.SparseCategoricalCrossentropy(from_logits=True)

        if optim:
            self.optim = optim
        else:
            self.optim = optimizers.Adam()

        if train_loss:
            self.train_loss = train_loss
        else:
            self.train_loss = metrics.Mean(name="train_loss")

        if train_acc:
            self.train_acc = train_acc
        else:
            self.train_acc = metrics.SparseCategoricalAccuracy(name="train_acc")

        if test_loss:
            self.test_loss = test_loss
        else:
            self.test_loss = metrics.Mean(name="test_loss")

        if test_acc:
            self.test_acc = test_acc
        else:
            self.test_acc = metrics.SparseCategoricalAccuracy(name="test_acc")

        return self

    @tf.function
    def _pr_train_func(self, seqmts, labels):
        with tf.GradientTape() as tape:
            predictions = self.model(seqmts, training=True)
            loss = self.loss_obj(labels, predictions)

        gradients = tape.gradient(loss, self.model.trainable_variables)
        self.optim.apply_gradients(zip(gradients, self.model.trainable_variables))
        self.train_loss(loss)
        self.train_acc(labels, predictions)

    @tf.function
    def _pr_test_func(self, seqmts, labels):
        predictions = self.model(seqmts, training=False)
        loss = self.loss_obj(labels, predictions)
        self.test_loss(loss)
        self.test_acc(labels, predictions)

    def load_dataset(self, ts_path, gene_id="", test_prop=0.3, batch=32,
                     shuffle=True, shuf_buf_size=128):
        """Load dataset from saved Numpy file.
        """
        _dfarr = np.array([[0, 0]])
        dataset = np.load(ts_path, allow_pickle=True).get(gene_id, _dfarr)

        ntmtrx = np.vstack(dataset[:, 0]).astype(np.float32)
        dataset_size, length, width = ntmtrx.shape
        ntmtrx = ntmtrx.reshape((dataset_size, 1, length, width))

        max_sampl_num = int(math.sqrt(length / 2)) ** 2 * 2
        if max_sampl_num != length:
            print("[W]: reshape the input data for the product of length *",
                  "width have no integer square root solution", file=sys.stderr)
            ntmtrx = ntmtrx[:, :, :max_sampl_num, :]

        self.iptmtx_shape = ntmtrx.shape

        asepro = np.vstack(dataset[:, 1])[:, 1].astype(np.int32) + 1

        self.dataset = tf.data.Dataset.from_tensor_slices((ntmtrx, asepro))

        test_size = int(dataset_size * test_prop)
        self.test_dtst = self.dataset.take(test_size)

        self.train_dtst = self.dataset.skip(test_size)

        if shuffle:
            self.train_dtst = self.train_dtst.shuffle(shuf_buf_size)
            self.test_dtst = self.test_dtst.shuffle(shuf_buf_size)

        if batch:
            self.train_dtst = self.train_dtst.batch(batch)
            self.test_dtst = self.test_dtst.batch(batch)

        return self

    def train(self, epochs=50):
        if self.model is None:
            self.model = CNNModel(self.iptmtx_shape[2:])

        for _epoch in range(epochs):
            self.train_loss.reset_states()
            self.train_acc.reset_states()
            self.test_loss.reset_states()
            self.test_acc.reset_states()

            for train_seqmts, train_labels in self.train_dtst:
                self._pr_train_func(train_seqmts, train_labels)

            for test_seqmts, test_labels in self.test_dtst:
                self._pr_test_func(test_seqmts, test_labels)
            
            if _epoch % 10 == 0:
                print("Ep: {}, loss: {}, acc: {}, test loss: {}, test acc: {}" \
                      .format(_epoch, self.train_loss.result(),
                              self.train_acc.result() * 100,
                              self.test_loss.result(),
                              self.test_acc.result() * 100))

        return self

    def visulization(self):
        # TODO
        print("Plot some fig to show the results.")
        return self

    def predict(self):
        # TODO
        print("Predict for new dataset.")
        return self
    
    def save_model(self):
        print("Save the trained model into disk.")
        return self

    def cleanup(self):
        # TODO: a place holder, perhaps not useful
        print("Clean up all the shit.")


if __name__ == "__main__":
    print("[W]: This module should not be executed directly.", file=sys.stderr)
