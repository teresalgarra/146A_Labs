function [y] = sixteenqammap(x1, x2, x3, x4)

    %This function maps 4 bits from {0,1} to
    %{-3-3j,-3-1j,-3+1j,-3+3j,-1-3j,-1-1j,-1+1j,-1+3j,+1-3j,+1-1j,+1+1j,+1+3j,+3-3j,+3-1j,+3+1j,+3+3j}

    len = length(x1);                           %Length of input vector
    y = zeros(1,len);                           %Preallocating vector

    for i = 1:len                               %Loop size
        bits = [x1(i),x2(i),x3(i),x4(i)];       %Bits vector
        if bits == [0,0,0,0]                    %Mapping 0000 to -3-3j
            y(i) = -3 -3j;
        elseif bits == [0,0,0,1]                %Mapping 0001 to -3-1j
            y(i) = -3 -1j;
        elseif bits == [0,0,1,0]                %Mapping 0010 to -3+1j
            y(i) = -3 +1j;
        elseif bits == [0,0,1,1]                %Mapping 0011 to -3+3j
            y(i) = -3 +3j;
        elseif bits == [0,1,0,0]                %Mapping 0100 to -1-3j
            y(i) = -1 -3j;
        elseif bits == [0,1,0,1]                %Mapping 0101 to -1-1j
            y(i) = -1 -1j;
        elseif bits == [0,1,1,0]                %Mapping 0110 to -1+1j
            y(i) = -1 +1j;
        elseif bits == [0,1,1,1]                %Mapping 0111 to -1+3j
            y(i) = -1 +3j;
        elseif bits == [1,0,0,0]                %Mapping 1000 to +1-3j
            y(i) = +1 -3j;
        elseif bits == [1,0,0,1]                %Mapping 1001 to +1-1j
            y(i) = +1 -1j;
        elseif bits == [1,0,1,0]                %Mapping 1010 to +1+1j
            y(i) = +1 +1j;
        elseif bits == [1,0,1,1]                %Mapping 1011 to +1+3j
            y(i) = +1 +3j;
        elseif bits == [1,1,0,0]                %Mapping 1100 to +3-3j
            y(i) = +3 -3j;
        elseif bits == [1,1,0,1]                %Mapping 1101 to +3-1j
            y(i) = +3 -1j;
        elseif bits == [1,1,1,0]                %Mapping 1110 to +3+1j
            y(i) = +3 +1j;
        elseif bits == [1,1,1,1]                %Mapping 1111 to +3+3j
            y(i) = +3 +3j;
        end                                     %End of conditions
    end                                         %End of loop

end
