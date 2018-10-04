#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//Global Variables
#define PICSIZE 256
#define MAXMASK 100

int pic[PICSIZE][PICSIZE];
int final[PICSIZE][PICSIZE];
int peaks[PICSIZE][PICSIZE];
double x_mask[MAXMASK][MAXMASK];
double y_mask[MAXMASK][MAXMASK];
double conv[PICSIZE][PICSIZE];
double tanDir[PICSIZE][PICSIZE];

// Prototype Functions
void peakFunction(double conv[PICSIZE][PICSIZE], double tanDir, int i, int j, FILE *input);
int histogramFunction(double percent, double conv[PICSIZE][PICSIZE], double maxival);
void hysteresisFunction(double low, double high);

// Begin Main Function
int main(argc, argv) int argc;
char **argv;
{
    // Declare Variables
    int i, j, p, q, mr, centx, centy;
    double y_maskval, x_maskval, sum1, sum2, sig, maxival, percent, hi, lo;
    char *foobar, shiftImage;

    FILE *fo1, *fo2, *fo3, *fp1, *fopen();

    // Argument command parameters
    argc--;
    argv++;
    foobar = *argv;
    fp1 = fopen(foobar, "rb"); // Open file for reading

    argc--;
    argv++;
    foobar = *argv;
    fo1 = fopen(foobar, "wb"); // Output pgm file for magnitude

    argc--;
    argv++;
    foobar = *argv;
    fo2 = fopen(foobar, "wb"); // Output pgm file for high peaks

    argc--;
    argv++;
    foobar = *argv;
    fo3 = fopen(foobar, "wb"); // Output pgm file for final edges

    argc--;
    argv++;
    foobar = *argv;
    sig = atof(foobar); // Read in the sig value and convert to integer

    argc--;
    argv++;
    foobar = *argv;
    percent = atof(foobar); // Read in the input value and convert to integer

    // Output the headers for all three files
    fprintf(fo1, "P5\n256 256\n255\n");
    fprintf(fo2, "P5\n256 256\n255\n");
    fprintf(fo3, "P5\n256 256\n255\n");

    mr = (int)(sig * 3);   // p.code , Mask Radius set
    centx = (MAXMASK / 2); // p.code, x Centroid found
    centy = (MAXMASK / 2); // p.code, y Centroid found

    // Read in the header file and store in a buffer to shift image back to original position.
    for (i = 0; i < 15; i++)
    {
        shiftImage = getc(fp1);
    }

    // Get all the pixels from the file
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            pic[i][j] = getc(fp1);
        }
    }
    // Perform the Guassian Equations for first order partial derivative.
    for (p = -mr; p <= mr; p++)
    {
        for (q = -mr; q <= mr; q++)
        {
            x_maskval = q * (exp(-1.0 * (((p * p) + (q * q)) / (2.0 * (sig * sig))))); // Find the X derivative
            x_mask[p + centy][q + centx] = x_maskval;                                  // X mask is the partial

            y_maskval = p * (exp(-1.0 * (((p * p) + (q * q)) / (2.0 * (sig * sig))))); // Find the y derivative
            y_mask[p + centy][q + centx] = y_maskval;                                  // Y mask is the partial
        }
    }

    maxival = 0; // Set max value to zero

    // Traverse through the elements by using a scanning convulution. Double up as we have two outputs
    for (i = mr; i <= 255 - mr; i++)
    {
        for (j = mr; j <= 255 - mr; j++)
        {
            sum1 = 0;
            sum2 = 0;
            for (p = -mr; p <= mr; p++)
            {
                for (q = -mr; q <= mr; q++)
                {
                    sum1 += pic[i + p][j + q] * x_mask[p + centy][q + centx];
                    sum2 += pic[i + p][j + q] * y_mask[p + centy][q + centx];
                }
            }

            conv[i][j] = sqrt(sum1 * sum1 + sum2 * sum2); // Perform sqrt of squares

            // Defensive case, since we can't divide by zero, we give it a very small value.
            if (sum1 == 0.0)
            {
                sum1 = 0.000001;
            }

            tanDir[i][j] = sum2 / sum1; // The slope, the tanDir value, used for high peaks

            // Perform the scaling
            if (conv[i][j] > maxival)
            {
                maxival = conv[i][j];
            }
        }
    }
    // Output the magnitude image and rescale
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            conv[i][j] = (conv[i][j] / maxival) * 255;     // Scale it down to fit within the max limit of 255
            fprintf(fo1, "%c", (char)((int)(conv[i][j]))); // Output magnitude image
        }
    }

    // Perform the High Peak Algorithm to search for all high peaks in the image
    for (i = 0; i < 256; i++)
    {
        for (j = 0; j < 256; j++)
        {
            peakFunction(conv, tanDir[i][j], i, j, fo2); // Check to see if it's a max peak
        }
    }

    // Find threshold values
    hi = histogramFunction(percent, conv, maxival);
    lo = 0.35 * hi;

    hysteresisFunction(lo, hi); // Begin checking the double threshold

    // Output the final picture
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            if (final[i][j] == 1)
            {
                fprintf(fo3, "%c", ((char)255));
            }
            else
                fprintf(fo3, "%c", ((char)0));
        }
    }
}

// Begin peak function to find max peaks
void peakFunction(double conv[PICSIZE][PICSIZE], double tanDir, int i, int j, FILE *input)
{
    double magLeft, magRight; // Temporary variables to use to check final pixel

    if ((tanDir <= .4142) && (tanDir > -.4142))
    {
        magLeft = conv[i][j - 1];
        magRight = conv[i][j + 1];
    }
    else if ((tanDir <= 2.4142) && (tanDir > .4142))
    {
        magLeft = conv[i - 1][j - 1];
        magRight = conv[i + 1][j + 1];
    }
    else if ((tanDir <= -.4142) && (tanDir > -2.4142))
    {
        magLeft = conv[i + 1][j - 1];
        magRight = conv[i - 1][j + 1];
    }
    else // Use as final else for the positions below and above, since it's hard to calculate, we make it last.
    {
        magLeft = conv[i - 1][j];
        magRight = conv[i + 1][j];
    }

    // Check results from if-else statement and flag pixel as peak or not.
    if ((conv[i][j] > magLeft) && (conv[i][j] > magRight))
    {
        fprintf(input, "%c", ((char)255)); // Mark pixel as high value
        peaks[i][j] = 1;                   // Flag pixel as possible true high
    }
    else
    {
        fprintf(input, "%c", ((char)0)); // Mark pixel as low value
        peaks[i][j] = 0;                 // Flag as false
    }
} // End peak function

// Begin function to find the double threshholds
void hysteresisFunction(double low, double high)
{
    int i, j, p, q; // Declare Counters
    int flag = 1;   // Status Flag

    // Check for peaks and high values
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            if (peaks[i][j] == 1 && conv[i][j] > high)
            {
                final[i][j] = 1; // if peak next to high value, make it high.
                peaks[i][j] = 0;
            }
            else if (conv[i][j] < low) // Else, make it zero
            {
                final[i][j] = 0;
                peaks[i][j] = 0;
            }
        }
    }

    // Now we need to check for peaks adjacent to high values
    while (flag == 1)
    {
        flag = 0;

        for (i = 0; i < PICSIZE; i++)
        {
            for (j = 0; j < PICSIZE; j++)
            {
                if (peaks[i][j] == 1) // If peak is true
                {
                    for (p = -1; p <= 1; p++) // Check magnitudes of adjacent peaks
                    {
                        for (q = -1; q <= 1; q++)
                        {
                            if (final[i + p][j + q] == 1) // If adjacent peak is high
                            {
                                peaks[i][j] = 0;
                                final[i][j] = 1; // Make peak a final peak
                                flag = 1;        // Set to 1 and start over
                            }
                        }
                    }
                }
            }
        }
    }

    return;
} // End function to find the double threshholds

//Begin histrogram function
int histogramFunction(double percent, double conv[PICSIZE][PICSIZE], double maxival)
{
    double cutoff, pixel, histogram[PICSIZE], areaOfTops;
    int i, j;

    areaOfTops = 0.0; // Initialize areaoftops

    cutoff = percent * PICSIZE * PICSIZE; // Set cutoff point to desired percent

    // Initialize histogram array
    for (i = 0; i < PICSIZE; i++)
    {
        histogram[i] = 0;
    }

    // Add the number of values at each magnitude to the histogram
    for (i = 0; i < PICSIZE; i++)
    {
        for (j = 0; j < PICSIZE; j++)
        {
            conv[i][j] = (conv[i][j] / maxival) * 255; // Scale it down to fit within the max limit of 255
            (histogram[(int)conv[i][j]])++;
        }
    }

    // Check for cutoff point
    for (i = PICSIZE - 1; i > 1; i--)
    {
        areaOfTops += histogram[i];

        if (areaOfTops > cutoff)
        {
            break;
        }
    }
    return i; // return hi value
} // End histogram function