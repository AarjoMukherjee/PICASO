#include "Common.h"
float Offset[96], Slope[96], Quad[96], Cube[96], SqrtCoeff[96];
float TimeBinWidthPs = 200.0;
float TimeScale = 5.0;
float TimeOffsetBin = 2000.0;
float EnergyScale = 1.0;
int ECal = 0;
int rootWrite = 0;
int storeTrace = 0;
