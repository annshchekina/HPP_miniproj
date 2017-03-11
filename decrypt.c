#include <stdio.h>
#include <stdlib.h>
#include "file_utils.h"
#include "utils.h"
#include "funcs.h"
#include "timings.h"

int main(int argc, char** argv) {
  if(argc != 4) {
    printf("Please give 3 arguments: input_file_name output_file_name n_threads\n");
    return -1;
  }

  const char* input_file_name = argv[1];
  const char* output_file_name = argv[2];
  int n_threads = atoi(argv[3]);
  printf("input_file_name = '%s'\n", input_file_name);
  printf("output_file_name = '%s'\n", output_file_name);
  printf("n_threads = %d\n", n_threads);

  // Get size of input file
  int fileSize = get_file_size(input_file_name);
  if(fileSize <= 0) {
    printf("Error: (fileSize <= 0)\n");
    return -1;
  }
  // Allocate buffer for input data
  unsigned char* buf0 = (unsigned char*)malloc(fileSize*sizeof(unsigned char));
  // Read input file
  if(read_file_into_buffer(input_file_name, buf0, fileSize) != 0) {
    printf("Error in read_file_into_buffer for input_file_name = '%s'.\n", input_file_name);
    return -1;
  }
  printf("Input file read OK, now starting computation...\n");

  double computationStartTime = get_wall_seconds();

  double Func1_time = 0, Func1_Func2_time = 0, Func2_time = 0, Func3_time = 0;
  double characteristicNumber = ComputeNumber(buf0, fileSize, 3, 
  	&Func1_time, &Func1_Func2_time, &Func2_time, &Func3_time);
  printf("Func1() time: %f\n Func2() time: %f\n Func3() time: %f\n",
    Func1_time + Func1_Func2_time, Func2_time - Func1_Func2_time, Func3_time);

  double start;
  
  unsigned char* buf1 = (unsigned char*)malloc(fileSize*sizeof(unsigned char));
  start = get_wall_seconds();
  for(int i = 0; i < fileSize; i++)
	  buf1[i] = ShuffleBitsA(buf0[i]);
  printf("ShuffleBitsA() time: %f\n", get_wall_seconds() - start);

  unsigned char* buf2 = (unsigned char*)malloc(fileSize*sizeof(unsigned char));
  start = get_wall_seconds();
  for(int i = 0; i < fileSize; i++)
	  buf2[i] = ShuffleBitsB(buf1[i]);
  printf("ShuffleBitsB() time: %f\n", get_wall_seconds() - start);
  
  unsigned char* buf3 = (unsigned char*)malloc(fileSize*sizeof(unsigned char));
  int fileSize_even = fileSize % 2 == 0 ? fileSize : fileSize - 1; 
  buf3[fileSize_even] = buf2[fileSize_even];
  for(int i = 0; i < fileSize_even - 1; i += 2)
  {
      buf3[i] = buf2[i+1];
      buf3[i+1] = buf2[i];
  }

  unsigned char* buf4 = (unsigned char*)malloc(fileSize*sizeof(unsigned char));
  start = get_wall_seconds();
  for(int i = 0; i < fileSize; i++)
	  buf4[i] = ShuffleBitsC(buf3[i]);
  printf("ShuffleBitsC() time: %f\n", get_wall_seconds() - start);

  unsigned char* buf5 = (unsigned char*)malloc(fileSize*sizeof(unsigned char));
  for(int i = 0; i < fileSize; i++)
      buf5[i] = buf4[fileSize - 1 - i];
  
  double computationTimeTaken = get_wall_seconds() - computationStartTime;
  printf("Computation took %f wall seconds.\n", computationTimeTaken);

  // Write output file
  if(write_buffer_to_file(output_file_name, buf5, fileSize) != 0) {
    printf("Error in write_buffer_to_file for output_file_name = '%s'.\n", output_file_name);
    return -1;
  }
  printf("Done. Output file '%s' written OK, %d bytes. characteristicNumber = %15.10f\n", 
	output_file_name, fileSize, characteristicNumber);

  free(buf0);
  free(buf1);
  free(buf2);
  free(buf3);
  free(buf4);
  free(buf5);
  
  return 0;
}
