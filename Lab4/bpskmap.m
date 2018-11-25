function [y] = bpskmap(x)

    %This function maps 1 bit from {0,1} to {-1,+1}

    len = length(x);                            %Length of input vector
    y = zeros(1,len);                           %Preallocating vector

    for i = 1:len                               %Loop size
        if x(i) == 0                            %Mapping 0 to -1
            y(i) = -1;
        elseif x(i) == 1                        %Mapping 1 to 1
            y(i) = 1;
        end                                     %End of conditions
    end                                         %End of loop

end
