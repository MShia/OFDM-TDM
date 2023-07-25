#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include  "ofdmCA.h"

/*
     _______     _______
    |__   __|   |__   __|
       | |        \  /
       | |         \/
       | |         /\
       | |        /  \
       | |       /    \
       |_|      /______\
    */ // Transmitter Function definitions 

// Data Mod - QPSK mapping
void QPSKModulation(int *InputArray, complex *symbols) {
    for (int i = 0; i < (Nc); i++) {
        int bit1 = InputArray[2 * i];
        int bit2 = InputArray[2 * i + 1];

        symbols[i].Re = (bit1 == 0) ? -1.0 : 1.0;
        symbols[i].Im = (bit2 == 0) ? -1.0 : 1.0;

       
    }
}

// Add Pilot Sequence Chu-pilot sequence p(i)= cos(PI i^2/Nc) + jsin(PI i^2/Nc)
void addPilotSequence(complex *symbols, complex *symbolswithPilot) {
int i = Nc - Nm;
while (i < Nc) {
    symbolswithPilot[i].Re = cos(PI * pow(i,2) / Nc);
    symbolswithPilot[i].Im = sin(PI * pow(i,2) / Nc);
    i++;
}
}

// Dividing Nc subcarriers to K time slots each 
void divideIntoTimeSlots(complex *symbols) {
    int symbolIndex = 0;
    for (int slot = 0; slot < K; slot++) {
        for (int subcarrier = 0; subcarrier < Nm; subcarrier++) {
            tSlot_symbols[slot][subcarrier] = symbols[symbolIndex];
            symbolIndex++;
        }
    }
}


void GenFrame(complex tSlot_symbols[][Nm], complex *gThFrame){
    for (int i = 0; i < K ; i++){
        IFFT(tSlot_symbols[i],  NFFT) ;  
    }
    complex Temp[Nc];
    for (int slot = 0; slot < K; slot++) {
        for (int subcarrier = 0; subcarrier < Nm; subcarrier++) {
            Temp[subcarrier + slot * Nm] = tSlot_symbols[slot][subcarrier];
            gThFrame[subcarrier + slot * Nm].Re = (sqrt(2) / sqrt(Nm)) * Temp[subcarrier + slot * Nm].Re;
            gThFrame[subcarrier + slot * Nm].Im = (sqrt(2) / sqrt(Nm)) * Temp[subcarrier + slot * Nm].Im;
        }
    }
}


void addCyclicPrefix (complex *gThFrame, complex *gThFramewithGI) {
    complex cyclicPrefix[GI];

    // Extract cyclic prefix samples from the end of S_g
    for (int i = Nc - GI, j = 0; i < Nc; i++, j++) {
        cyclicPrefix[j] = gThFrame[i];
    }

    // Place cyclic prefix samples at the start of S_g
    for (int i = Nc + GI - 1, j = Nc - 1; j >= 0; i--, j--) {
        gThFramewithGI[i] = gThFrame[j];
    }

    // Copy cyclic prefix samples to the start of S_g
    for (int i = 0; i < GI; i++) {
        gThFramewithGI[i] = cyclicPrefix[i];
    }
}

void softLimiting(complex *gThFrameWithGI, complex *HPAout) {
    //complex output;
    double mag[Nc + GI]; 
    for (int m= 0 ; m < Nc +GI; m++){
         mag[m] = sqrt(gThFrameWithGI[m].Re * gThFrameWithGI[m].Re + gThFrameWithGI[m].Im * gThFrameWithGI[m].Im);
        if (mag[m] < beta) {
        HPAout[m].Re = sqrt(P)*gThFrameWithGI[m].Re;
        HPAout[m].Im = sqrt(P)*gThFrameWithGI[m].Im;
    }   
        else {
        double scaleFactor = beta / mag[m];
        HPAout[m].Re = sqrt(P)*(gThFrameWithGI[m].Re * scaleFactor);
        HPAout[m].Im = sqrt(P)*(gThFrameWithGI[m].Im * scaleFactor);
    }
          
    }
    
        return ;
}

/*
      CCCCCC     H     H
     C      C    H     H
    C            HHHHHHH
   C             H     H
   C             H     H
    C      C     H     H
     CCCCCC      H     H
*/


complex complexMultiply(complex a, complex b) {
    complex result;
    result.Re = a.Re * b.Re - a.Im * b.Im;
    result.Im = a.Re * b.Im + a.Im * b.Re;
    return result;
}

double* generateAWGN(int length) {
    double* noise = (double*)malloc(length * sizeof(double));
    double noise_var = (2*No)/Tc;
    if (noise == NULL) {
        printf("Memory allocation failed!\n");
        return NULL;
    }

    srand(time(NULL));  // Seed the random number generator

    for (int i = 0; i < length; i += 2) {
        double u1 = (double)rand() / RAND_MAX;  // Random number between 0 and 1
        double u2 = (double)rand() / RAND_MAX;  // Random number between 0 and 1

        double z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
        double z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * PI * u2);

        noise[i] = sqrt(noise_var/2)*z1;
        if (i + 1 < length) {
            noise[i + 1] = sqrt(noise_var/2)*z2;
        }
    }

    return noise;
}


void convolution(complex *ChImpRes, int ChLength, complex *TxOut, int TxOutLen) {
    int ConvLen = TxOutLen + ChLength - 1;
    double *noise; 
    noise = generateAWGN(ConvLen); 
  
    for (int i = 0; i < ConvLen; i++) {
        ChOut[i].Re = 0.0;
        ChOut[i].Im = 0.0;

        for (int j = 0; j < TxOutLen; j++) {
            if (i - j >= 0 && i - j < ChLength) {
                complex mult = complexMultiply(ChImpRes[j], TxOut[i - j]);
                ChOut[i].Re += mult.Re + noise[i];
                ChOut[i].Im += mult.Im + noise[i];
            }
        }
    }
}


/*
      RRRRR    X       X
     R     R    X     X
     R      R    X   X
     RRRRRRRR     XXX
     R     R     X   X
     R      R   X     X
     R       R X       X
*/

void removeGI(complex *ChOut, complex *RxSymbWOgi, complex *gi_symbols, int input_len) {
    int output_len = input_len - GI;

    for (int i = 0; i < output_len; i++) {
        RxSymbWOgi[i].Re = HPAout[i + GI].Re;
        RxSymbWOgi[i].Im = HPAout[i + GI].Im;
    }

    for (int i = 0; i < GI; i++) {
        gi_symbols[i].Re = HPAout[i].Re;
        gi_symbols[i].Im = HPAout[i].Im;
    }
}

void mmseFDE(complex R[], complex H[], int R_len, double noise_variance) {
    // Perform MMSE-FDE equalization
    for (int k = 0; k < R_len; k++) {
        double H_real = H[k].Re;
        double H_imag = H[k].Im;
        double denominator = pow(H_real, 2) + pow(H_imag, 2) + 2 * noise_variance;
        complex W;
        W.Re = H_real / denominator; 
        W.Im = - H_imag / denominator;
        R[k].Re *= W.Re;
        R[k].Im *= W.Im;
    }
}


void demodulateQPSK(complex* Rxsymbols, int numSymbols, int* demodulatedData) {
    for (int i = 0; i < numSymbols; i++) {
        // Extract the real (I) and imaginary (Q) components of the symbol
        double I = Rxsymbols[i].Re;
        double Q = Rxsymbols[i].Im;

        // Determine the binary values based on the I and Q components
        demodulatedData[2 * i] = (I >= 0) ? 1 : 0;
        demodulatedData[2 * i + 1] = (Q >= 0) ? 1 : 0;
    }
}

void IFFT(complex* input, int size) {
    if (size <= 1) {
        // Base case: do nothing
        return;
    }

    // Divide the input into even and odd indices
    complex even[size / 2];
    complex odd[size / 2];
    for (int i = 0; i < size / 2; i++) {
        even[i] = input[2 * i];
        odd[i] = input[2 * i + 1];          
    }

    // Recursive calls to IFFT
    IFFT(even, size / 2);
    IFFT(odd, size / 2);

    // Combine the results of the recursive calls
    for (int k = 0; k < size / 2; k++) {
        double angle = 2 * PI * k / size;
        double cosAngle = cos(angle);
        double sinAngle = sin(angle);

        complex t;
        t.Re = cosAngle * odd[k].Re - sinAngle * odd[k].Im;
        t.Im = cosAngle * odd[k].Im + sinAngle * odd[k].Re;

        input[k].Re = even[k].Re + t.Re;
        input[k].Im = even[k].Im + t.Im;
        input[k + size / 2].Re = even[k].Re - t.Re;
        input[k + size / 2].Im = even[k].Im - t.Im;
    }
}


void FFT(complex* input, int size) {
    if (size <= 1) {
        // Base case: do nothing
        return;
    }

    // Divide the input into even and odd indices
    complex even[size / 2];
    complex odd[size / 2];
    for (int i = 0; i < size / 2; i++) {
        even[i] = input[2 * i];
        odd[i] = input[2 * i + 1];          
    }

    // Recursive calls to IFFT
    IFFT(even, size / 2);
    IFFT(odd, size / 2);

    // Combine the results of the recursive calls
    for (int k = 0; k < size / 2; k++) {
        double angle = -2 * PI * k / size;
        double cosAngle = cos(angle);
        double sinAngle = sin(angle);

        complex t;
        t.Re = cosAngle * odd[k].Re - sinAngle * odd[k].Im;
        t.Im = cosAngle * odd[k].Im + sinAngle * odd[k].Re;

        input[k].Re = even[k].Re + t.Re;
        input[k].Im = even[k].Im + t.Im;
        input[k + size / 2].Re = even[k].Re - t.Re;
        input[k + size / 2].Im = even[k].Im - t.Im;
    }
}
