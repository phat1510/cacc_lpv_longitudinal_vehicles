clc
clear



numerator = 1;
denominator = [0.5,1];
sys = tf(numerator,denominator)


step(sys)
