#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <R.h>

#define REAL 0
#define IMAG 1

double fft_convolution(double kernel[], int kernellen, double signal[], int signallen, double result[]){

    int i;
    int j;

    int convolution_length = kernellen + signallen - 1; //the length of the resulting convolution

    //We have to add 0's at the end of the arrays (called 0 padding) so that the convolution will work
    double padded_signal[convolution_length];
    double padded_kernel[convolution_length];
    int k;

    for (k=0; k< convolution_length; k++){
        padded_signal[k] = 0;
        padded_kernel[k] = 0;
    }


    for (k=0; k< signallen; k++){
        padded_signal[k] = signal[k];
    }

    for (k=0; k< kernellen; k++){
        padded_kernel[k] = kernel[k];
    }


    fftw_complex *fft_multi; 
    fft_multi = fftw_malloc(sizeof (fftw_complex)*convolution_length);



    fftw_complex *signal_fft, *kernel_fft; //the results of the fft for the signal and for the kernel
    signal_fft = fftw_malloc(sizeof(fftw_complex)* (convolution_length));
    kernel_fft = fftw_malloc(sizeof(fftw_complex)* (convolution_length));

    //create the fft plans
    fftw_plan plan1 = fftw_plan_dft_r2c_1d(convolution_length,
                                             padded_signal,
                                             signal_fft,
                                             FFTW_ESTIMATE);

    fftw_plan plan2 = fftw_plan_dft_r2c_1d(convolution_length,
                                         padded_kernel,
                                         kernel_fft,
                                         FFTW_ESTIMATE);

    //execute them, this is what actually does the fft
    fftw_execute(plan1);
    fftw_execute(plan2);


    for(i=0 ; i < convolution_length; i++){

        fft_multi[i][REAL] = (signal_fft[i][REAL]*kernel_fft[i][REAL] - signal_fft[i][IMAG]*kernel_fft[i][IMAG])/convolution_length;
        fft_multi[i][IMAG] = (signal_fft[i][REAL]*kernel_fft[i][IMAG] + signal_fft[i][IMAG]*kernel_fft[i][REAL])/convolution_length;

    }


    fftw_plan plan3 = fftw_plan_dft_c2r_1d(convolution_length,fft_multi,result,FFTW_ESTIMATE);
    fftw_execute(plan3);




    fftw_destroy_plan(plan3);
    fftw_destroy_plan(plan2);
    fftw_destroy_plan(plan1);

    fftw_free(kernel_fft);
    fftw_free(signal_fft);
    fftw_free(fft_multi);


}


double direct_convolution_local(double probs[], int probslen, double result[]){

  int i,j;
  int oldlen = 2; // length of old kernel
  int currlen = 0; // length of current result (new kernel)
  double signal[2];
  double t,tmp;
  int k;

  // initialize (old kernel)
  result[0] = 1-probs[0];
  result[1] = probs[0];

  // loop through all other probs
  for(i=1; i < probslen; i++){

    // set signal
    signal[0] = probs[i];
    signal[1] = 1-probs[i];
    
    // create current kernel
    currlen = oldlen + 1;

    // initialize result
    result[currlen-1] = signal[0] * result[oldlen-1];
    oldlen++;
    t = result[0];//signal[1] * result[0];
    result[0] = signal[1]*t;
    for(j=1; j < currlen-1; j++){
      tmp=result[j];
      result[j] = signal[0] * t + signal[1] * result[j];
      t=tmp;
    }

  }

}




//The same as direct_convolution_local except it is formatted so that it can be called directly from R
void direct_convolution(double probs[], int *probslen, double result[]){

  int i,j;
  int oldlen = 2; // length of old kernel
  double signal[2];
  double t;
  int currlen;
  double tmp;

  // initialize (old kernel)
  result[0] = 1-probs[0];
  result[1] = probs[0];

  // loop through all other probs
  for(i=1; i < *probslen; i++){

    // set signal
    signal[0] = probs[i];
    signal[1] = 1-probs[i];
   
    // create current kernel
    currlen = oldlen + 1;

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



void tree_convolution(double p_vec[], int *p_vec_len, int *num_splits, double result_vec[]){

  int number_splits = *num_splits+0;

  int leftover = *p_vec_len % number_splits; 

  //int final_size = *p_vec_len + number_splits; 

  int total_splits = log(number_splits)/log(2);

  //double result_vec[final_size]; 

  int base = *p_vec_len/number_splits;  
  int num_in_each[number_splits];
  int i;
  for (i = 0; i < number_splits; i++){
    num_in_each[i] = base;
  }


  for (i = 0; i < leftover; i++){
    num_in_each[i % number_splits] += 1;
  }


    int memory_index_tracker = 0;
    int left_memory_index_tracker = 0;
    int right_memory_index_tracker = 0;

    for (i = 0; i < number_splits; i++){
        direct_convolution_local(p_vec+memory_index_tracker, num_in_each[i], result_vec+memory_index_tracker+i);
        memory_index_tracker += num_in_each[i];
    }

    int temp_num_splits = number_splits/2;


    for(i = 0; i < number_splits; i++){
        num_in_each[i] = num_in_each[i] += 1;
    }


    int res_tracker = 0;
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