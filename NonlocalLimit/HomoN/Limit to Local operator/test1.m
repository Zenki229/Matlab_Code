clear;clc;close all
kerfun = @kernel;
Stiff_nonlocal_row_free(0.2,0.5,1,kerfun)