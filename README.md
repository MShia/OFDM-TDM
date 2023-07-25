OFDM/TDM with MMSE-FDE Readme

This code implements the paper titled "Performance of OFDM/TDM with MMSE-FDE Using Pilot-assisted Channel Estimation". It provides an implementation of various functions required for the transmitter and receiver operations in an OFDM/TDM system with MMSE-FDE.
Files

The code is organized into three files:

    main.c: This file contains the main function from where the other functions are called.
    ofdmCA.h: This header file includes the function declarations for all the functions defined in ofdm.c.
    ofdmCA.c: This file contains the implementations of the functions declared in ofdm.h.

Transmitter Functions

The following functions are declared in ofdmCA.h and implemented in ofdmCA.c. They are used for the transmitter operations.
QPSKModulation

Function Declaration:
void QPSKModulation(int* InputArray, complex* symbols);

This function takes binary data as input and maps it to QPSK symbols.
addPilotSequence

Function Declaration:
void addPilotSequence(complex* symbols, complex* symbolswithPilot);

The addPilotSequence function adds a pilot sequence to the modulated symbols at position K-1.
divideIntoTimeSlots

Function Declaration:
void divideIntoTimeSlots(complex* symbols);

The divideIntoTimeSlots function divides the modulated symbols into time slots and stores the output in tSlot_symbols.
IFFT

Function Declaration:
void IFFT(complex* input, int size);

The IFFT function performs an N-point Inverse Fast Fourier Transform (IFFT) on the input.
GenFrame

Function Declaration:
void GenFrame(complex tSlot_symbols[][Nm], complex* gThFrame);

The GenFrame function generates the gth frame as Sg in the paper's signal.
addCyclicPrefix

Function Declaration:
void addCyclicPrefix(complex* gThFrame, complex* gThFramewithGI);

The addCyclicPrefix function adds a guard interval (GI) to the gth frame and stores the output in gThFramewithGI.
softLimiting

Function Declaration:
void softLimiting(complex* gThFrameWithGI, complex* HPAout);

The softLimiting function applies soft limiting to the gth frame passed through the High Power Amplifier (HPA).

CHANNEL FUNCTIONS

Function Declaration:
complex complexMultiply(complex a, complex b);

The complexMultiply function performs complex multiplication of two complex numbers.
convolution

Function Declaration:
void convolution(complex* ChImpRes, int ChLength, complex* TxOut, int TxOutLen);

The convolution function performs convolution between complex channel coefficients and the transmitted signal.
removeGI



RECEIVER FUNCTIONS
The following functions are also declared in ofdmCA.h and implemented in ofdmCA.c. They are used for the receiver operations. The receiver side is not complete - specially the channel estimation part has to be added to achieve the full functiontionality.   

complexMultiply


Function Declaration:
void removeGI(complex* ChOut, complex* RxSymbWOgi, complex* gi_symbols, int input_len);

The removeGI function removes the guard interval (GI) from the received symbols.
FFT

Function Declaration:
void FFT(complex* input, int size);

The FFT function performs Fast Fourier Transform (FFT) on the input.


/* After getting the R(k) - MMSE-FDE is to be performed which could not be finalized yet
    As in Eq(4) - R(k) = S(k)H(k) + N(k) - Here the channel response and Noise have to be estimated 
    and then Eq(6) i.e. R^(k) = w(k)R(k) is used for MMSE-FDE. 
    Arbitrary channel response h(t) is generated to bypass the channel estimation step and complete
    validate the overall flow of the code  - The section is being commented out though   
    */ 
mmseFDE

Function Declaration:
void mmseFDE(complex R[], complex H[], int R_len, double noise_variance);

The mmseFDE function implements Minimum Mean Square Error (MMSE) Frequency Domain Equalization (FDE) using the received signal, channel coefficients, and noise variance.
demodulateQPSK

Function Declaration:
void demodulateQPSK(complex* Rxsymbols, int numSymbols, int* demodulatedData);

The demodulateQPSK function demodulates QPSK symbols to obtain the corresponding binary data.

Usage

To use this code, you need to include the ofdmCA.h header file in your main program and call the required functions based on your application requirements. Make sure to provide the necessary inputs and handle the outputs as per your specific use case.

Please refer to the function descriptions and comments within the code for more detailed information about each function's usage and inputs/outputs.
Note:

This code is an implementation of the paper "Performance of OFDM/TDM with MMSE-FDE Using Pilot-assisted Channel Estimation". It provides a framework for OFDM/TDM systems and includes functions for the transmitter, channel modeling and receiver operations. You may need to modify and adapt the code to suit your specific requirements and system configurations.
References

If you use this code or refer to the implementation, please consider citing the original paper:

Paper Title: "Performance of OFDM/TDM with MMSE-FDE Using Pilot-assisted Channel Estimation"
Authors: Haris Gacanin and Fumiyuki Adachi


Please consult the paper for more information about the underlying algorithms and techniques used in this code.