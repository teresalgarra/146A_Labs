function [y] = bpskmapinverse(x)

    %This function maps 1 bit from {-1,+1} to {0,1}

    len = length(x);                            %Length of input vector
    y = zeros(1,len);                           %Preallocating vector

    for i = 1:len                               %Loop size
        if real(x(i)) < 0                       %Mapping -1 to 0
            y(i) = 0;
        elseif real(x(i)) > 0                   %Mapping 1 to 1
            y(i) = 1;
        end                                     %End of conditions
    end                                         %End of loop

end
