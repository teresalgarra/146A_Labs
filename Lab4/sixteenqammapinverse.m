function [y1, y2, y3, y4] = sixteenqammapinverse(x)

    %This function maps 4 bits from {-3-3j,-3-1j,-3+1j,-3+3j,-1-3j,-1-1j,-1+1j,-1+3j,+1-3j,+1-1j,+1+1j,+1+3j,+3-3j,+3-1j,+3+1j,+3+3j}
    %to {0,1}

    len = length(x);                                                              %Length of input vector
    y1 = zeros(1,len);                                                            %Preallocating vector
    y2 = zeros(1,len);                                                            %Preallocating vector
    y3 = zeros(1,len);                                                            %Preallocating vector
    y4 = zeros(1,len);                                                            %Preallocating vector

    for i = 1:len                                                                 %Loop size
        if (real(x(i))<=-2)&&(imag(x(i))<=-2)                                     %Mapping -3-3j to 0000
            y1(i) = 0;
            y2(i) = 0;
            y3(i) = 0;
            y4(i) = 0;
        elseif (real(x(i))<=-2)&&(imag(x(i))>-2)&&(imag(x(i))<=0)                 %Mapping -3-1j to 0001
            y1(i) = 0;
            y2(i) = 0;
            y3(i) = 0;
            y4(i) = 1;
        elseif (real(x(i))<=-2)&&(imag(x(i))>0)&&(imag(x(i))<=2)                  %Mapping -3+1j to 0010
            y1(i) = 0;
            y2(i) = 0;
            y3(i) = 1;
            y4(i) = 0;
        elseif (real(x(i))<=-2)&&(imag(x(i))>2)                                   %Mapping -3+3j to 0011
            y1(i) = 0;
            y2(i) = 0;
            y3(i) = 1;
            y4(i) = 1;
        elseif (real(x(i))>-2)&&(real(x(i))<=0)&&(imag(x(i))<=-2)                 %Mapping -1-3j to 0100
            y1(i) = 0;
            y2(i) = 1;
            y3(i) = 0;
            y4(i) = 0;
        elseif (real(x(i))>-2)&&(real(x(i))<=0)&&(imag(x(i))>-2)&&(imag(x(i))<=0) %Mapping -1-1j to 0101
            y1(i) = 0;
            y2(i) = 1;
            y3(i) = 0;
            y4(i) = 1;
        elseif (real(x(i))>-2)&&(real(x(i))<=0)&&(imag(x(i))>0)&&(imag(x(i))<=2)  %Mapping -1+1j to 0110
            y1(i) = 0;
            y2(i) = 1;
            y3(i) = 1;
            y4(i) = 0;
        elseif (real(x(i))>-2)&&(real(x(i))<=0)&&(imag(x(i))>2)                   %Mapping -1+3j to 0111
            y1(i) = 0;
            y2(i) = 1;
            y3(i) = 1;
            y4(i) = 1;
        elseif (real(x(i))>0)&&(real(x(i))<=2)&&(imag(x(i))<=-2)                  %Mapping +1-3j to 1000
            y1(i) = 1;
            y2(i) = 0;
            y3(i) = 0;
            y4(i) = 0;
        elseif (real(x(i))>0)&&(real(x(i))<=2)&&(imag(x(i))>-2)&&(imag(x(i))<=0)  %Mapping +1-1j to 1001
            y1(i) = 1;
            y2(i) = 0;
            y3(i) = 0;
            y4(i) = 1;
        elseif (real(x(i))>0)&&(real(x(i))<=2)&&(imag(x(i))>0)&&(imag(x(i))<=2)   %Mapping +1+1j to 1010
            y1(i) = 1;
            y2(i) = 0;
            y3(i) = 1;
            y4(i) = 0;
        elseif (real(x(i))>0)&&(real(x(i))<=2)&&(imag(x(i))>2)                    %Mapping +1+3j to 1011
            y1(i) = 1;
            y2(i) = 0;
            y3(i) = 1;
            y4(i) = 1;
        elseif (real(x(i))>2)&&(imag(x(i))<=-2)                                   %Mapping +3-3j to 1100
            y1(i) = 1;
            y2(i) = 1;
            y3(i) = 0;
            y4(i) = 0;
        elseif (real(x(i))>2)&&(imag(x(i))>-2)&&(imag(x(i))<=0)                   %Mapping +3-1j to 1101
            y1(i) = 1;
            y2(i) = 1;
            y3(i) = 0;
            y4(i) = 1;
        elseif (real(x(i))>2)&&(imag(x(i))>0)&&(imag(x(i))<=2)                    %Mapping +3+1j to 1110
            y1(i) = 1;
            y2(i) = 1;
            y3(i) = 1;
            y4(i) = 0;
        elseif (real(x(i))>2)&&(imag(x(i))>2)                                     %Mapping +3+3j to 1111
            y1(i) = 1;
            y2(i) = 1;
            y3(i) = 1;
            y4(i) = 1;
        end                                                                       %End of conditions
    end                                                                           %End of loop

end
