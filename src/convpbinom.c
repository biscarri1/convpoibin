#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <R.h>

//used to make fft computations more readable
#define REAL 0
#define IMAG 1

    
//This function carries out FFT based convolution on a two sequences. One called the "kernel" and the other called the "signal"
//The work flow is 1. transform the two sequences to fourier domain 2. multiply the two resulting sequences 3. transform back 
//For more detail see: https://en.wikipedia.org/wiki/Convolution_theorem
double fft_convolution(double kernel[], int kernellen, double signal[], int signallen, double result[]){

    int i;
    int j;

    int convolution_length = kernellen + signallen - 1; //the length of the resulting convolution aka our final result

    //We have to add 0's at the end of the arrays (called 0 padding) so that the convolution will work
    double padded_signal[convolution_length];
    double padded_kernel[convolution_length];
    int k;
    
    
    for (k=0; k< convolution_length; k++){
        padded_signal[k] = 0;
        padded_kernel[k] = 0;
    }

    
    //The for loops below fill out the values of padded_signal and padded_kernel with the values of 
    //signal and kernel
    for (k=0; k< signallen; k++){
        padded_signal[k] = signal[k];
    }

    for (k=0; k< kernellen; k++){
        padded_kernel[k] = kernel[k];
    }

    
    //allocate memory for fft_multi, which will hold the result of the fourier domain multiplcation of the two sequences
    fftw_complex *fft_multi; 
    fft_multi = fftw_malloc(sizeof (fftw_complex)*convolution_length);


    //allocate memory for the result of doing the fft on the signal and the fft on the kernel
    fftw_complex *signal_fft, *kernel_fft; 
    signal_fft = fftw_malloc(sizeof(fftw_complex)* (convolution_length));
    kernel_fft = fftw_malloc(sizeof(fftw_complex)* (convolution_length));

    //create the fft plans, used by the fftw3 library
    fftw_plan plan1 = fftw_plan_dft_r2c_1d(convolution_length,
                                             padded_signal,
                                             signal_fft,
                                             FFTW_ESTIMATE);

    fftw_plan plan2 = fftw_plan_dft_r2c_1d(convolution_length,
                                         padded_kernel,
                                         kernel_fft,
                                         FFTW_ESTIMATE);

    //execute the plans, this is what actually does the fft and transforms the signal and kernel
    fftw_execute(plan1);
    fftw_execute(plan2);

    //Now carry out the multiplcation of the sequences in the fourier domain. Note that we normalize by the length.
    for(i=0 ; i < convolution_length; i++){

        fft_multi[i][REAL] = (signal_fft[i][REAL]*kernel_fft[i][REAL] - signal_fft[i][IMAG]*kernel_fft[i][IMAG])/convolution_length;
        fft_multi[i][IMAG] = (signal_fft[i][REAL]*kernel_fft[i][IMAG] + signal_fft[i][IMAG]*kernel_fft[i][REAL])/convolution_length;

    }

    //Transform fft_multi back from the fourier domain to get the final result
    fftw_plan plan3 = fftw_plan_dft_c2r_1d(convolution_length,fft_multi,result,FFTW_ESTIMATE);
    fftw_execute(plan3);


    //The below code destroys the plans and frees up memory
    fftw_destroy_plan(plan3);
    fftw_destroy_plan(plan2);
    fftw_destroy_plan(plan1);

    fftw_free(kernel_fft);
    fftw_free(signal_fft);
    fftw_free(fft_multi);


}


//This function calculates the Poisson binomial density by means of directly calculating the convolution definition of
//the distribution of a sum of random variables. See: https://en.wikipedia.org/wiki/Convolution_of_probability_distributions
double direct_convolution_local(double probs[], int probslen, double result[]){

  int i,j;
  int oldlen = 2; // length of old kernel
  double signal[2];
  double t,tmp;


  // initialize (old kernel)
  result[0] = 1-probs[0];
  result[1] = probs[0];

  // loop through all other probs
  for(i=1; i < probslen; i++){

    // set signal
    signal[0] = probs[i];
    signal[1] = 1-probs[i];

    // initialize result and calculate the two edge cases
    result[oldlen] = signal[0] * result[oldlen-1];

    t = result[0];
    result[0] = signal[1]*t;
      
    //calculate the interior cases
    for(j=1; j < oldlen; j++){
      tmp=result[j];
      result[j] = signal[0] * t + signal[1] * result[j];
      t=tmp;
    }
      
    oldlen++;
  }

}




//The same as direct_convolution_local except the variables are pointers so that it can be called directly from R.
void direct_convolution(double probs[], int *probslen, double result[]){

  int i,j;
  int oldlen = 2; // length of old kernel
  double signal[2];
  double t, tmp;

  // initialize (old kernel)
  result[0] = 1-probs[0];
  result[1] = probs[0];

  // loop through all other probs
  for(i=1; i < *probslen; i++){

    // set signal
    signal[0] = probs[i];
    signal[1] = 1-probs[i];
   
    // calculate last element of result
    result[oldlen] = signal[0] * result[oldlen-1];
           
    t = result[0];//signal[1] * result[0];
      
    result[0] = signal[1]*t;
      
    for(j=1; j < oldlen; j++){
      tmp=result[j];
      result[j] = signal[0] * t + signal[1] * result[j];
      t=tmp;
    }
      
      oldlen++;

  }

}



//This functions calculates the Poisson binomial density via a binary tree structured FFT convolution approach.
//First, the function partitions the inputted probabilities into different groups, making them as equal in size as possible.
//Then, it applies direct_convolution_local to each of the individual groups.
//Finally, it takes pairs of groups and convolves them together using fft_convolution, always operating on pairs of equal size
//It continues this process until only a single sequence is left. 
void tree_convolution(double p_vec[], int *p_vec_len, int *num_splits, double result_vec[]){

  int number_splits = *num_splits+0;

  //see how many times p_vec_len is divided by num_splits
  int leftover = *p_vec_len % number_splits; 

  int total_splits = log(number_splits)/log(2);

  int base = *p_vec_len/number_splits;//see how many points will be in each bin. Note that 
  //int division in C is like a floor function, so, for example if p_vec_len/num_splits is 
  //45.78 then base will equal 45. This means that we will be missing some points.
  //To remedy this problem, the next two for loops are used.  
    
  //create an array containing the number of points in each bin
  int num_in_each[number_splits];
  int i;
  for (i = 0; i < number_splits; i++){
    num_in_each[i] = base;
  }

  //This iterates through the num_splits array and keeps adding 1 to each element, in a circular
  //fashion, leftover amount of times. This is basically distributing the remainder points evenly
  //across the num_splits array.
  for (i = 0; i < leftover; i++){
    num_in_each[i % number_splits] += 1;
  }

    //Do the initial regular convolution on each bin
    int memory_index_tracker = 0;
    int left_memory_index_tracker = 0;
    int right_memory_index_tracker = 0;

    for (i = 0; i < number_splits; i++){
        direct_convolution_local(p_vec+memory_index_tracker, num_in_each[i], result_vec+memory_index_tracker+i);
        memory_index_tracker += num_in_each[i];
    }

    int temp_num_splits = number_splits/2;

      // add 1 to each num_in_each to account for the fact that the length of each bin increases by 1
    // after doing the regular convolution.
    for(i = 0; i < number_splits; i++){
        num_in_each[i] += 1;
    }


    int res_tracker = 0;
    
    //Take each pair and apply fft_convolution to them, and continue doing so until we only have one sequence, which is the result
    while (total_splits > 0){

        memory_index_tracker = 0;
        left_memory_index_tracker = 0;
        right_memory_index_tracker = 0;
        res_tracker = 0;
        for (i = 0; i < temp_num_splits; i++){

            left_memory_index_tracker += num_in_each[2*i];


            fft_convolution(result_vec+memory_index_tracker,num_in_each[2*i],
                            result_vec+left_memory_index_tracker,num_in_each[2*i+1],
                            result_vec+res_tracker);

            left_memory_index_tracker += num_in_each[2*i+1];


            memory_index_tracker += num_in_each[2*i]+num_in_each[2*i+1];

            res_tracker = memory_index_tracker - (i+1);

        }

        for (i = 0; i < temp_num_splits; i++){
            num_in_each[i] = num_in_each[2*i]+num_in_each[2*i+1]-1;
        }

        temp_num_splits = temp_num_splits/2;

        total_splits = total_splits - 1;

    }


  }