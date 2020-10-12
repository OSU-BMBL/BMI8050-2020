## TODO: Define the first convolutional layer

1. model.add(Conv1D(filters=32, kernel_size=12, input_shape=(train_features.shape[1], 2)))
2. model.add(Conv1D(filters=32, pool_size=12, input_shape=(test_features.shape[1], 4)))
3. model.add(Conv1D(filters=32, kernel_size=12, input_shape=(train_features.shape[1], 4)))
4. model.add(Conv1D(filters=32, pool_size=12, input_shape=(test_features.shape[1], 2)))

## TODO: Define the first max pooling layer

1. model.add(MaxPooling1D(pool_size=1))
1. model.add(MaxPooling1D(pool_size=4))
2. model.add(MaxPooling1D(pool_size=1, activation='relu'))
2. model.add(MaxPooling1D(pool_size=4, activation='relu'))
   
## TODO: Define the fist dense layer

1. model.add(Dense(16, activation='relu'))
2. model.add(Dense(16, activation='softmax'))
3. model.add(Dense(16, activation='sigmoid'))
   
## TODO: Define the last Dense layer to output the classification probabilities. 

1. model.add(Dense(1, activation='sigmoid'))
2. model.add(Dense(1, activation='softmax'))
3. model.add(Dense(2, activation='sigmoid'))
4. model.add(Dense(2, activation='softmax'))


## TODO: Define the compile operation with your optimizer and learning rate of choice'''

1. model.compile(loss='adam', optimizer='binary_crossentropy', metrics=['binary_accuracy'])
2. model.compile(loss='adam', optimizer='binary_crossentropy', activation='softmax')
3. model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['binary_accuracy'])
4. model.compile(loss='binary_crossentropy', optimizer='adam', activation='softmax')
   
## TODO: Use model.fit to train the CNN model, with the number of epochs, vervose type.

1. model.fit(train_features, test_labels, epochs=50, verbose=1, validation_split=0.25)
2. model.fit(test_features, train_labels, epochs=50, verbose=1, validation_split=0.25)
3. model.fit(train_features, train_labels, epochs=50, verbose=1, validation_split=0.25)
4. model.fit(test_features, test_labels, epochs=50, verbose=1, validation_split=0.25)

An epoch is a single pass through the full training set
The more epoch the better?

## TODO: Use model.predict to train the CNN model, with the number of epochs, vervose type.

1. model.predict(train_features)
2. model.predict(test_features)
1. model.predict(train_labels)
2. model.predict(test_labels)

98% accuracy is great, can we imporve the performance?
Add activation function to convolutional layer

What if we add more convolutional layers?