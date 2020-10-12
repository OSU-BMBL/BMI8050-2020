## TODO: Define the first convolutional layer

A. model.add(Conv1D(filters=32, kernel_size=12, input_shape=(train_features.shape[1], 2)))
B. model.add(Conv1D(filters=32, pool_size=12, input_shape=(test_features.shape[1], 4)))
C. model.add(Conv1D(filters=32, kernel_size=12, input_shape=(train_features.shape[1], 4)))
D. model.add(Conv1D(filters=32, pool_size=12, input_shape=(test_features.shape[1], 2)))

## TODO: Define the first max pooling layer

A. model.add(MaxPooling1D(pool_size=1))
B. model.add(MaxPooling1D(pool_size=4))
C. model.add(MaxPooling1D(pool_size=1, activation='relu'))
D. model.add(MaxPooling1D(pool_size=4, activation='relu'))
   
## TODO: Define the first dense layer

A. model.add(Dense(16, activation='relu'))
B. model.add(Dense(16, activation='softmax'))
C. model.add(Dense(16, activation='sigmoid'))
   
## TODO: Define the last Dense layer to output the classification probabilities

A. model.add(Dense(1, activation='sigmoid'))
B. model.add(Dense(1, activation='softmax'))
C. model.add(Dense(2, activation='sigmoid'))
D. model.add(Dense(2, activation='softmax'))


## TODO: Define the compile operation with your optimizer

A. model.compile(loss='adam', optimizer='binary_crossentropy', metrics=['binary_accuracy'])
B. model.compile(loss='adam', optimizer='binary_crossentropy', activation='softmax')
C. model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])
D. model.compile(loss='binary_crossentropy', optimizer='adam', activation='softmax')
   
## TODO: Use model.fit to train the CNN model, with the number of epochs, verbose type

A. model.fit(train_features, test_labels, epochs=50, verbose=1, validation_split=0.25)
B. model.fit(test_features, train_labels, epochs=50, verbose=1, validation_split=0.25)
C. model.fit(train_features, train_labels, epochs=50, verbose=1, validation_split=0.25)
D. model.fit(test_features, test_labels, epochs=50, verbose=1, validation_split=0.25)



## TODO: Use model.predict to train the CNN model

A. model.predict(train_features)
B. model.predict(test_features)
C. model.predict(train_labels)
D. model.predict(test_labels)

## Other questions
98% accuracy is great, can we improve the performance?
Add activation function to convolutional layer
What if we add more convolutional layers?
What if we use sigmoid in the output layer?
The more epoch the better?

Notes:
Sigmoid function is used for the two-class logistic regression, whereas the softmax function is used for the multiclass logistic regression.
The output of the softmax describes the probability (or if you may, the confidence) of our neural network that a particular sample belongs to a certain class

