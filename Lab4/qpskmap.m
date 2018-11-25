function [y] = qpskmap(x1, x2)

    %This function maps 2 bits from {0,1} to {-1-j,-1+j,+1-j,+1+j}

    len = length(x1);                           %Length of input vector
    y = zeros(1,len);                           %Preallocating vector

    for i = 1:len                               %Loop size
        bits = [x1(i),x2(i)];                   %Bits vector
        if bits == [0,0]                        %Mapping 00 to -1-j
            y(i) = -1 -1j;
        elseif bits == [0,1]                    %Mapping 01 to -1+j
            y(i) = -1 +1j;
        elseif bits == [1,0]                    %Mapping 10 to +1-j
            y(i) = +1 -1j;
        elseif bits == [1,1]                    %Mapping 11 to +1+j
            y(i) = +1 +1j;
        end                                     %End of conditions
    end                                         %End of loop

end
