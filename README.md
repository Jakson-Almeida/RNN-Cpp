
# Recurrent Neural Network Library (C++)

Welcome to the Recurrent Neural Network Library, an open-source project written entirely in C++ for implementing and testing Recurrent Neural Networks (RNNs) on resource-constrained environments like microcontrollers (e.g., Arduino and ESP32). This library was designed to be lightweight and efficient, making it ideal for embedded systems with limited memory.

## Key Features

Pure C++ Implementation: All linear algebra and matrix calculations were implemented from scratch, without the use of third-party libraries.

Dynamic Memory Management: The library automatically deallocates unused memory, optimizing its use in devices with limited resources.

## Customizable Architecture:

Supports N inputs and M outputs.

Can be configured as a deep network with up to N layers.

For each layer, M temporal (recurrent) layers can be added, enabling memory and recurrent functionalities.


##Supported Activation Functions:

Sigmoid

ReLU


## External Training: 

The training process is external to the network. An example neural network was trained using the Particle Swarm Optimization (PSO) algorithm.


## Example Usage

Here’s an example of how to set up, configure, and use the library:

```c++
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Neural_Network.h" // Include your library header

float *out;
int deepLearning[]   = {2, 10, 10, 2};
int temporalLayers[] = {2, 1, 0, 0};
float vet[] = {1, -1};

int main()
{
    int TAM = sizeof(deepLearning) / sizeof(int);
    int tamVet = sizeof(vet) / sizeof(float);
    Neural_Network neuron(deepLearning, temporalLayers, TAM);
    Matrix weightsAndBias(neuron.getVetValuesLength(), 1);
    Matrix vetEntrance(tamVet, 1);
    
    for(int i = 0; i < weightsAndBias.getNumLinhas(); i++)
        weightsAndBias.set(i, 0, i * 0.3 - 3);
    
    neuron.updateValues(weightsAndBias);
    vetEntrance.updateValues(vet, tamVet);
    neuron.updateEntrance(&vetEntrance);
    neuron.feedforward();
    out = neuron.getOutputVet();

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "First Output:" << std::endl;
    std::cout << "v1 = " << out[0] << std::endl;
    std::cout << "v2 = " << out[1] << std::endl << std::endl;

    std::cout << "Second Output:" << std::endl;
    vet[0] = 0;
    vet[1] = 0;
    vetEntrance.updateValues(vet, tamVet);
    neuron.updateEntrance(&vetEntrance);
    neuron.feedforward();
    std::cout << "v1 = " << neuron.getOutput(0) << std::endl;
    std::cout << "v2 = " << neuron.getOutput(1) << std::endl;

    return 0;
}
```


Installation

Clone the repository into your project:

git clone https://github.com/your-username/your-repository.git

Include the necessary header files in your C++ project, and you’re ready to go.

Requirements

This library was built using only the following standard C++ headers:

#include <iostream>
#include <iomanip>
#include <cmath>

It does not depend on any external libraries, ensuring portability and simplicity.

Applications

Embedded Systems: The library is optimized for use on microcontrollers such as Arduino and ESP32.

Lightweight AI Prototyping: Test recurrent neural networks in environments with limited resources.

Custom AI Architectures: Easily define custom deep and recurrent layers for a wide variety of use cases.


Author

This library was developed by Jakson-Almeida, with all linear algebra and matrix operations implemented from scratch. The design ensures flexibility and efficiency for low-memory environments.

Contribution

Contributions are welcome! Feel free to open issues or submit pull requests to improve the library. Let’s make lightweight neural networks accessible to everyone.

License

This project is licensed under the MIT License.