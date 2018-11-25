function [y1, y2] = fourpammapinverse(x)

    %This function maps 2 bits from {-3,-1,+1,+3} to {0,1}

    len = length(x);                            %Length of input vector
    y1 = zeros(1,len);                          %Preallocating vector
    y2 = zeros(1,len);                          %Preallocating vector

    for i = 1:len                               %Loop size
        if (x(i) <= -2)                         %Mapping -3 to 00
            y1(i) = 0;
            y2(i) = 0;
        elseif (x(i) > -2) && (x(i) <= 0)       %Mapping -1 to 01
            y1(i) = 0;
            y2(i) = 1;
        elseif (x(i) > 0) && (x(i) <= 2)        %Mapping +1 to 10
            y1(i) = 1;
            y2(i) = 0;
        elseif (x(i) > 2)                       %Mapping +3 to 11
            y1(i) = 1;
            y2(i) = 1;
        end                                     %End of conditions
    end                                         %End of loop

end
