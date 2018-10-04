
#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

// Global Variables
int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
int masky[3][3] = {{1, 2, 1}, {0, 0, 0}, {-1, -2, -1}};
double ival[256][256], maxival;

main(argc, argv) int argc;
char **argv;
{
    // Declare Variables
    int i, j, p, q, mr, sum1, sum2;
    double loThreshold, hiThreshold;
    FILE *fo1, *fo2, *fo3, *fp1, *fopen();
    char *foobar;

    //Open main file for reading
    argc--;
    argv++;
    foobar = *argv;
    fp1 = fopen(foobar, "rb");

    // Read in the high threshhold value
    argc--;
    argv++;
    foobar = *argv;
    hiThreshold = atoi(foobar);

    // Read in the low threshold value.
    argc--;
    argv++;
    foobar = *argv;
    loThreshold = atoi(foobar);

    // Create file for output for magnitude
    argc--;
    argv++;
    foobar = *argv;
    fo1 = fopen(foobar, "wb");

    // Create file for output for high threshold
    argc--;
    argv++;
    foobar = *argv;
    fo2 = fopen(foobar, "wb");

    // Create file for output for low threshold
    argc--;
    argv++;
    foobar = *argv;
    fo3 = fopen(foobar, "wb");

    // Output the headers for all three files
    fprintf(fo1, "P5\n256 256\n255\n"); 
    fprintf(fo2, "P5\n256 256\n255\n");  
    fprintf(fo3, "P5\n256 256\n255\n"); 

    // Get all the pixels from the file
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            pic[i][j] = getc(fp1);
        }
    }

    // Set mask radius due to matrix size. Traverse through the elements by using a scanning convulution.
    mr = 1;
    for (i = mr; i < 256 - mr; i++)
    {
        for (j = mr; j < 256 - mr; j++)
        {
            sum1 = 0;
            sum2 = 0;
            for (p = -mr; p <= mr; p++)
            {
                for (q = -mr; q <= mr; q++)
                {
                    sum1 += pic[i + p][j + q] * maskx[p + mr][q + mr];
                    sum2 += pic[i + p][j + q] * masky[p + mr][q + mr];
                }
            }
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
        }
    }

    // Find the max gradient and perform a sqrt of squares. Use this value to do the scaling.
    maxival = 0;
    for (i = mr; i < 256 - mr; i++)
    {
        for (j = mr; j < 256 - mr; j++)
        {
            ival[i][j] = sqrt((double)((outpicx[i][j] * outpicx[i][j]) +
                                       (outpicy[i][j] * outpicy[i][j])));
            if (ival[i][j] > maxival)
                maxival = ival[i][j];
        }
    }

    //First loop to calculate the magnitude of the picture
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            ival[i][j] = (ival[i][j] / maxival) * 255; // Scale it down to fit within the max limit of 255
            fprintf(fo1, "%c", (char)((int)(ival[i][j])));
        }
    }

    // Second loop to calculate the hithreshold picture
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {

            if (ival[i][j] > hiThreshold)
                fprintf(fo2, "%c", ((char)(255))); // If element is greater than threshold, set it to max value

            else
                fprintf(fo2, "%c", ((char)(0))); // Else, set the value to the lowest possible
        }
    }

    // Third loop to calculate the low threshold picture
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {

            if (ival[i][j] > loThreshold)
                fprintf(fo3, "%c", (char)(255)); // If element is greater than threshold, set it to max value

            else
                fprintf(fo3, "%c", (char)(0)); // Else, set the value to the lowest possible
        }
    }
}