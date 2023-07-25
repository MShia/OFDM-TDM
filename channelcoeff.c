#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define NUM_SUBCARRIERS 64
#define NUM_PATHS 4

int main() {
    FILE *file;
    complex double h[NUM_SUBCARRIERS * NUM_PATHS];

    // Open the file
    file = fopen("channel_response.txt", "r");
    if (file == NULL) {
        printf("Error opening the file.\n");
        return -1;
    }

    // Read the channel response data from the file
    for (int i = 0; i < NUM_SUBCARRIERS * NUM_PATHS; i++) {
        double real, imag;
        if (fscanf(file, "%lf %lf", &real, &imag) != 2) {
            printf("Error reading data from the file.\n");
            return -1;
        }
        h[i] = real + I * imag;
    }

    // Close the file
    fclose(file);

    // Use the channel response data in your C code

    return 0;
}
