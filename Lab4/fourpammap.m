function [y] = fourpammap(x1, x2)

    %This function maps 2 bits from {0,1} to {-3,-1,+1,+3}

    len = length(x1);                           %Length of input vector
    y = zeros(1,len);                           %Preallocating vector

    for i = 1:len                               %Loop size
        bits = [x1(i),x2(i)];                   %Bits vector
        if bits == [0,0]                        %Mapping 00 to -3
            y(i) = -3;
        elseif bits == [0,1]                    %Mapping 01 to -1
            y(i) = -1;
        elseif bits == [1,0]                    %Mapping 10 to +1
            y(i) = +1;
        elseif bits == [1,1]                    %Mapping 11 to +3
            y(i) = +3;
        end                                     %End of conditions
    end                                         %End of loop

end
