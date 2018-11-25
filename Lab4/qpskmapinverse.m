function [y1, y2] = qpskmapinverse(x)

    %This function maps 2 bits from {-1-j,-1+j,+1-j,+1+j} to {0,1}

    len = length(x);                            %Length of input vector
    y1 = zeros(1,len);                          %Preallocating vector
    y2 = zeros(1,len);                          %Preallocating vector

    for i = 1:len                                   %Loop size
        if (real(x(i))<=0) &&  (imag(x(i))<=0)      %Mapping -1-j to 00
            y1(i) = 0;
            y2(i) = 0;
        elseif (real(x(i))<=0) &&  (imag(x(i))>0)   %Mapping -1+j to 01
            y1(i) = 0;
            y2(i) = 1;
        elseif (real(x(i))>0) &&  (imag(x(i))<=0)   %Mapping +1-j to 10
            y1(i) = 1;
            y2(i) = 0;
        elseif (real(x(i))>0) &&  (imag(x(i))>0)    %Mapping +1+j to 11
            y1(i) = 1;
            y2(i) = 1;
        end                                         %End of conditions
    end                                             %End of loop

end
