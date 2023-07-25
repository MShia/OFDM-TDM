#ifndef _ofdmCA_h_
#define _ofdmCA_h_

#define Nc 256      // Number of subcarriers 
#define K 16        // K time slots 
#define Nm Nc/K     // Number of subcarriers in Kth time slot 
#define PI	3.14159265358979323846264338327950288
#define GI Nm       // Nm Guard Interval 
#define beta 0.9     // HPA saturation level
#define P 1         // Transmit Power Coefficient
#define NFFT 16      // # FFT points - FFT window size 
#define No 0.0001    // Noise power 
#define Tc 256       // FFT sampling period
#define ChLen 8      // length of h - channel impulse resp
#define convLen (Nc+GI) + ChLen -1 // length of channel covolution 


//Defining complex structure with Re and Im 
typedef struct {
    double Re;
    double Im;
} complex;


// Variable declarations  
/*Transmitter*/
complex symbols[Nc];             // arry  : QPSK symbols  
complex symbolswithPilot[Nc];    // array : with pilot symbols added
complex tSlot_symbols[K][Nm];    // array : K time slots with Nm subcarriers 
complex output_IFFT[NFFT];       // array : IFFT out 
complex gThFrame[Nc];            // array : Sg signal - gth frame with out CyclicPrefix/GI  
complex gThFrameWithGI[Nc + GI]; // array : gth frame with GI 
complex HPAout[Nc+GI];           // array : HPA output soft limit (HPA) function
/*Channel*/
complex ChOut[convLen];          // array : channel output 
/*Receiver*/
complex RxSymbWOgi[convLen - GI];// array : Removal of GI 
complex gi_symbols[GI];          // array : CP/ GI symols 
complex RxFFT[Nc];               // array : FFT out 


/*
     _______    
    |__   __|    \    /
       | |        \  /
       | |         \/
       | |         /\
       | |        /  \
       | |       /    \
       |_|      /      \
    */
/*Function Declaration for transmitter*/

//Data Mod  - Takes in binary data and maps to QPSK symbols 
void QPSKModulation(int* InputArray, complex* symbols); 
// Add pilot sequence to the mod symbols at K-1
void addPilotSequence(complex *symbols, complex *symbolswithPilot); 
// TDM - Slots - stores output in tSlot_symbols 
void divideIntoTimeSlots(complex *symbols);
// Nm - Point IFFT 
void IFFT(complex* input, int size);
// Generation of gth Frame as Sg in papersignal
void GenFrame(complex tSlot_symbols[][Nm], complex *gThFrame);
// Addition of Guard interval as GI in paper - output stored in gThFramewithGI
void addCyclicPrefix (complex *gThFrame, complex *gThFramewithGI);
// Passing the gthFrame through HPA - soft Limiting 
void softLimiting(complex *gThFrameWithGI, complex *HPAout);

/*
      CCCCCC     H     H
     C      C    H     H
    C            HHHHHHH
   C             H     H
   C             H     H
    C      C     H     H
     CCCCCC      H     H
*/

//complex Multiply function used in implementing Convolution Operation
complex complexMultiply(complex a, complex b);
//Convolution function takes complex channel coefficients, its length , Received signal(HPA out) and its length
//The output is stored in ChOut  
void convolution(complex *ChImpRes, int ChLength, complex *TxOut, int TxOutLen); 
/*
      RRRRR    X       X
     R     R    X     X
     R      R    X   X
     RRRRRRRR     XXX
     R     R     X   X
     R      R   X     X
     R       R X       X
*/
// removeGI function takes the output from 
void removeGI(complex *ChOut, complex *RxSymbWOgi, complex *gi_symbols, int input_len);
// FFT implementation - takes complex input and FFT size 
void FFT(complex* input, int size);
//mmseFDE implementation - inputs R(k), H(k) and length and nosie variance 
void mmseFDE(complex R[], complex H[], int R_len, double noise_variance); 

// QPSK symbol demapping to bits 
void demodulateQPSK(complex* Rxsymbols, int numSymbols, int* demodulatedData); 

#endif