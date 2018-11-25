function [y] = eightpskmap(x1,x2,x3)

    %This function maps 3 bits from {0,1} to exp(i*j*2*pi/8) with i=0,1,...,7

    len = length(x1);                           %Length of input vector
    y = zeros(1,len);                           %Preallocating vector

    for i = 1:len                               %Loop size
        bits = [x1(i),x2(i),x3(i)];             %Bits vector
        if bits == [0,0,0]                      %Mapping 000 to i=0
            y(i) = exp(0*1j*2*pi/8);
        elseif bits == [0,0,1]                  %Mapping 001 to i=1
            y(i) = exp(1*1j*2*pi/8);
        elseif bits == [0,1,0]                  %Mapping 010 to i=2
            y(i) = exp(2*1j*2*pi/8);
        elseif bits == [0,1,1]                  %Mapping 011 to i=3
            y(i) = exp(3*1j*2*pi/8);
        elseif bits == [1,0,0]                  %Mapping 100 to i=4
            y(i) = exp(4*1j*2*pi/8);
        elseif bits == [1,0,1]                  %Mapping 101 to i=5
            y(i) = exp(5*1j*2*pi/8);
        elseif bits == [1,1,0]                  %Mapping 110 to i=6
            y(i) = exp(6*1j*2*pi/8);
        elseif bits == [1,1,1]                  %Mapping 111 to i=7
            y(i) = exp(7*1j*2*pi/8);
        end                                     %End of conditions
    end                                         %End of loop

end
