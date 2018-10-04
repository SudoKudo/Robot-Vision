# How to run
# Make sure exe file is in the same directory

gcc sobel.c -o sobel -lm
 
./sobel face05.pgm 100 35 mag.pgm high.pgm low.pgm
