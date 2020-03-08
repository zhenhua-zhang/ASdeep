"""Deep learning predictor module."""
import warnings
import numpy as np

# Supress some warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    import tensorflow as tf
    from tensorflow.keras.layers import BatchNormalization
    from tensorflow.keras.layers import Dense, Conv2D, Flatten
    from tensorflow.keras.models import Model
except ImportError as ime:
    print("ImportError: {}".format(ime))

class CNNModel(Model):
    """A CNN class."""
    def __init__(self):
        super(CNNModel, self).__init__(name='')

        self.conv1 = Conv2D(32, 3)
        self.bn1 = BatchNormalization()

        self.conv2 = Conv2D(32, 3)
        self.bn2 = BatchNormalization()

        self.flatten = Flatten()

        self.d1 = Dense(128, activation='relu')
        self.d2 = Dense(10, activation='relu')

    def call(self, input_tensor):
        x = self.conv1(input_tensor)
        x = self.bn1(x)
        x = tf.nn.relu(x)

        x = self.conv2(x)
        x = self.bn2(x)
        x = tf.nn.relu(x)

        x = self.flatten(x)
        x = self.d1(x)
        return self.d2(x)


class CNNFactory:
    """A factory to train CNN model on ASE results.
    """
    def __init__(self, model, loss_obj, optimizer, train_loss, train_acc,
                 test_loss, test_acc):
        super(CNNFactory, self).__init__()
        self.model = model

        self.loss_obj = loss_obj
        self.optimizer = optimizer

        self.train_loss = train_loss
        self.train_acc = train_acc

        self.test_loss = test_loss
        self.test_acc = test_acc

        self.train_dtst = None
        self.test_dtst = None

    @tf.function
    def _p_train_func(self, seqmts, labels):
        with tf.GradientTape() as tape:
            predictions = self.model(seqmts, training=True)
            loss = self.loss_obj(labels, predictions)
        gradients = tape.gradient(loss, self.model.trainable_variables)
        self.optimizer.apply_gradients(
            zip(gradients, self.model.trainable_variables)
        )
        self.train_loss(loss)
        self.train_acc(labels, predictions)

    @tf.function
    def _p_test_func(self, seqmts, labels):
        predictions = self.model(seqmts, training=False)
        loss = self.loss_obj(labels, predictions)
        self.test_loss(loss)
        self.test_acc(labels, predictions)

    def load_dataset(self, path, target_gene="", test_prop=0.3, batch=8,
                     shuffle=True, shuf_buf_size=100):
        """Load dataset from saved Numpy file.
        """
        dataset = np.load(path, allow_pickle=True)
        ntmtrx = dataset.get(target_gene)[1]
        asepro = dataset.get(target_gene)[0]

        dataset_size, *_ = ntmtrx.shape
        test_size = int(dataset_size * test_prop)

        self.dataset = tf.data.Dataset.from_tensor_slices((ntmtrx, asepro))
        if shuffle:
            self.dataset = self.dataset.shuffle(shuf_buf_size)

        if batch:
            self.dataset = self.dataset.batch(batch)

        self.test_dtst = self.dataset.take(test_size)
        self.train_dtst = self.dataset.skip(test_size)

        return self

    def train(self, epochs=5000):
        for _epoch in range(epochs):
            self.train_loss.reset_states()
            self.train_acc.reset_states()
            self.test_loss.reset_states()
            self.test_acc.reset_states()

            for train_seqmts, train_labels in self.train_dtst:
                self._p_train_func(train_seqmts, train_labels)

            for test_seqmts, test_labels in self.test_dtst:
                self._p_test_func(test_seqmts, test_labels)

        return self

    def visulization(self):
        # TODO
        print("Plot some fig to show the results.")
        return self

    def predict(self):
        # TODO
        print("Predict for new dataset.")
        return self

    def cleanup(self):
        # TODO: a place holder, perhaps not useful
        print("Clean up all the shit.")
