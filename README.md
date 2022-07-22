# PAD_HE

# Private-Anomaly-Detection
# Setup

This project requires <a href="https://github.com/microsoft/SEAL/releases/tag/v3.6.2"> the SEAL Version 3.6.2 </a>.
The examples.cpp 
<a href="https://github.com/microsoft/SEAL/releases/tag/v3.6.2"> The SEAL Version 3.6.2 </a> should be built and run with the mentioned dependencies in the <a href="https://github.com/microsoft/SEAL"> the SEAL directory </a>, i.e., "/SEAL/native/examples/" directory should be built and run. 

# Building and Running

Once it is done, the committed programming files is downloaded and put into the "/SEAL/native/examples/". That is, the content of the "/SEAL/native/examples/" is filled with the programming files.). Then, following is done:

````
$ cmake --build build 
$ cd build/bin
$ ./sealexamples
````

The theoretically client and server are evaluated as two separate parties. However, the code contains both functionalities in the same main function, which could be separated programmatically if needed/required.   


