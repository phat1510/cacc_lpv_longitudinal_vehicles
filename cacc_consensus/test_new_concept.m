clc
clear
close all

A = [0 0 0 0;
    0 0 0 0 ;
    0 1 0 -1;
    0 0 0 0];
m = 1500;
h = 1;
B = 1/m*[1 0 0 0;
    0 1 0 0;
    0 0 -h 0;
    0 0 1 0];

Co = ctrb(A,B);

unco = length(A) - rank(Co)