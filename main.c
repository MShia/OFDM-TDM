#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "ofdmCA.h"

int main(void) {
    
// Channel Related 
    // Seed the random number generator
    srand(time(NULL));
    
    //***************** Generate Input stream of 1's and 0's**************
    //********************************************************************
    int InputArray[2 * (Nc)];
    for (int i = 0; i < 2 * Nc; i++) {
        InputArray[i] = rand() % 2; // Generates 0 or 1
    }
    printf("The Input bit Stream \n");
    for (int x = 0; x < 2*Nc; x++){
        printf (" %d ", InputArray[x])  ;
    }
    int inputsize = sizeof(InputArray) / sizeof(InputArray[0]); 

    // Arbitrary Channel Impulse Response 
    complex ChImpRes[ChLen] = {
        {0.5, 0.2},
        {0.8, -0.4},
        {0.3, -0.6},
        {-0.1, 0.7},
        {-0.9, 0.3},
        {-0.6, -0.2},
        {-0.3, 0.1},
        {0.4, -0.5}
    };
    
   
    /*
     _______     _______
    |__   __|   |__   __|
       | |        \  /
       | |         \/
       | |         /\
       | |        /  \
       | |       /    \
       |_|      /______\
    */

    // Maps the binary InputArray to QPSK symbols stores in symbols
    QPSKModulation(InputArray, symbols); 
    int modsymsize = sizeof(symbols) / sizeof(symbols[0]);

    
    //Adds Pilot sequence at K-1 slot and stores in symbolswithPilot
   addPilotSequence(symbols,symbolswithPilot);


    //Creates the K TDM slots and stores in tSlot_symbols
    divideIntoTimeSlots(symbols); 
    
    
   //kth OFDM signal with Nm subcarriers
   // Computes Nm- point IFFT for each time slot Generates gth OFDM/TDM frame and 
   //stores output in gThFrame - without Guard Interval / Cyclic Prefix as Eq(2)    
    GenFrame(tSlot_symbols,gThFrame);
    
    //Adding Guard Interval or CyclicPrefix to the Gth Frame
    //and stores the frame with GI in gThFramewithGI
    addCyclicPrefix(gThFrame, gThFrameWithGI); 

    // gth signal frame with GI is then passed to softlimiting function which implements HPA characteristics
    // as in Eq(3) Completion of Fig2.a Transmitter
    softLimiting(gThFrameWithGI, HPAout);
   
    
/*
      CCCCCC     H     H
     C      C    H     H
    C            HHHHHHH
   C             H     H
   C             H     H
    C      C     H     H
     CCCCCC      H     H
*/
  
    // The channel effects are implemented in convolution function
    // Takes output of softlimit func and arbitrary Ch impulse response as input output is stored in ChOut
    convolution(ChImpRes, ChLen, HPAout, Nc+GI);
   

 
      
/*
      RRRRR    X       X
     R     R    X     X
     R      R    X   X
     RRRRRRRR     XXX
     R     R     X   X
     R      R   X     X
     R       R X       X
*/
    int sizeafterchannel = sizeof(ChOut) / sizeof(ChOut[0]); // computes size of channel output


    // At receiver side the frames received include the GI payload which is removed 
    // Takes ChOut as input, RxSymbWOgi stores the data without GI, and gi_symbols are removed GI symbols
    removeGI(ChOut, RxSymbWOgi, gi_symbols, sizeafterchannel);
    int RxSymbGILen = sizeof(RxSymbWOgi) / sizeof(RxSymbWOgi[0]); 
    



    // Get R(k) the frequncy domain of Received signal(without GI)
    FFT(RxSymbWOgi, NFFT);

    for (int x =0 ; x < Nc; x++){
        RxFFT[x] = RxSymbWOgi[x];    // copies output to RxFFT 
    }

   
    /* After getting the R(k) - MMSE-FDE is to be performed which could not be finalized yet
    As in Eq(4) - R(k) = S(k)H(k) + N(k) - Here the channel response and Noise have to be estimated 
    and then Eq(6) i.e. R^(k) = w(k)R(k) is used for MMSE-FDE. 
    Arbitrary channel response h(t) is generated to bypass the channel estimation step and complete
    validate the overall flow of the code  - The section is being commented out though   
    */ 
	/*
    int X = 16;
	int Y =16; 

	FILE *file;
	complex  h[X * Y];

	// Open the file
    file = fopen("channel_response.txt", "r");
    if (file == NULL) {
        printf("Error opening the file.\n");
        return -1;
    }

    
    // Read the channel response data from the file
    for (int i = 0; i < X * Y; i++) {
        double real, imag;
        if (fscanf(file, "%lf %lf", &real, &imag) != 2) {
            printf("Error reading data from the file.\n");
            return -1;
        }
        h[i].Re = real; 
        h[i].Im = imag;
    }

    // Close the file
    fclose(file);

    // taking FFT of channel response 
    FFT(h, NFFT);   // H(k)

    // Calls mmsEFDE function which takes h(H(k)) , RxFFT (R(k)) and computes the weights
    double n_var = (2*No)/ Tc ; 
    mmseFDE(RxFFT, h, Nc, n_var);   
    int lenghtMMSE = sizeof(RxFFT) / sizeof(RxFFT[0]); 
    */

       
    printf("\n"); 


    // After equalization (not completed) - R^(k) is converted back to time domain
    IFFT(RxFFT, NFFT);
    int size_RxFFT = sizeof(RxFFT) / sizeof(RxFFT[0]); 
  
  
    // And Demodulation is the final step where symobls are converted back to 1s n 0's
    int demodulatedData[2 * Nc]; 
    demodulateQPSK(RxFFT, size_RxFFT, demodulatedData); 
    printf("The Demodulated BitStream \n"); 
    for (int x = 0; x < 2* size_RxFFT; x++){
        printf (" %d ", demodulatedData[x]) ;
    }

   




return 0;
}


