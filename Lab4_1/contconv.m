function [y, t] = contconv(x1, x2, t1, t2, dt)

    %This function computes an approximation to continuous-time 
    %convolution. It has 5 inputs: the two signals to convolute, two
    %starting times for the samples and the spacing of the samples.
    
    y=conv(x1,x2)*dt;                           %MatLab Built-in Function
    t = (0 : length(y)-1)*dt + t1 + t2;         %Time vector

end


