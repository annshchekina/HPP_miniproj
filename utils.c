#include "utils.h"
#include <stdio.h>
#include <stdlib.h>

int bits[8], bitsB[8];

inline unsigned char ShuffleBitsA(unsigned char c) {
  for(int i = 0; i < 8; i++)
    bits[i] = (c >> i) & 0x01;
  // Swap bits pairwise, starting at 0
  for(int i = 0; i < 8; i+=2) {
    // Swap bits i and i+1
    int tmp = bits[i];
    bits[i] = bits[i+1];
    bits[i+1] = tmp;
  }
  unsigned char result = 0;
  for(int i = 0; i < 8; i++) {
    if(bits[i])
      result |= (0x01 << i);
  }
  return result;
}

inline unsigned char ShuffleBitsB(unsigned char c) {
  for(int i = 0; i < 8; i++)
    bits[i] = (c >> i) & 0x01;
  // Swap bits pairwise, starting at 1
  for(int i = 1; i < 7; i+=2) {
    // Swap bits i and i+1
    int tmp = bits[i];
    bits[i] = bits[i+1];
    bits[i+1] = tmp;
  }
  // Create bitsB list where bits from bitsA are placed in reverse order 
  for(int i = 0; i < 8; i++)
    bitsB[i] = bits[7-i];
  unsigned char result = 0;
  for(int i = 0; i < 8; i++) {
    if(bitsB[i])
      result |= (0x01 << i);
  }
  return result;
}

inline unsigned char ShuffleBitsC(unsigned char c) {
  for(int i = 0; i < 8; i++)
    bits[i] = (c >> i) & 0x01;
  // Swap bits 3 and 6
  int tmp1 = bits[3];
  bits[3] = bits[6];
  bits[6] = tmp1;
  // Swap bits 1 and 7
  int tmp2 = bits[1];
  bits[1] = bits[7];
  bits[7] = tmp2;
  unsigned char result = 0;
  for(int i = 0; i < 8; i++) {
    if(bits[i])
      result |= (0x01 << i);
  }
  return result;
}

void PrintBufferContents(const unsigned char* buf, int n) {
  for(int i = 0; i < n; i++)
    printf("%d ", (int)buf[i]);
  printf("\n");
}
