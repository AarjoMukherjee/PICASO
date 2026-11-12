#include "Common.h"
#include "Calibration.h"
#include <cstdio>
#include <iostream>

void GainMatch(int Param, int ECal, char *gainfile)
{
    // defaults
    for(int i=0;i<Param;i++)
    {
        Offset[i]    = 0.0;
        Slope[i]     = 1.0;
        Quad[i]      = 0.0;
        Cube[i]      = 0.0;
        SqrtCoeff[i] = 0.0;
    }

    if(!ECal)
        return;

    FILE *fp = fopen(gainfile,"r");

    if(!fp)
    {
        printf("Cannot open gain calibration file: %s\n", gainfile);
        return;
    }

    int ch = 0;

    while(ch < Param)
    {
        double off, slope, quad, cube, sqrtc;

        int nread =
            fscanf(fp,
                   "%lf %lf %lf %lf %lf",
                   &off,
                   &slope,
                   &quad,
                   &cube,
                   &sqrtc);

        if(nread != 5)
            break;

        Offset[ch]    = off;
        Slope[ch]     = slope;
        Quad[ch]      = quad;
        Cube[ch]      = cube;
        SqrtCoeff[ch] = sqrtc;

        ch++;
    }

    fclose(fp);

    printf("Loaded calibration for %d channels\n", ch);
}
