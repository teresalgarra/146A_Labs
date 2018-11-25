function [y1, y2, y3] = eightpskmapinverse(x)

    %This function maps 3 bits from exp(i*j*2*pi/8) with i=0,1,...,7 to {0,1}

    len = length(x);                                            %Length of input vector
    y1 = zeros(1,len);                                          %Preallocating vector
    y2 = zeros(1,len);                                          %Preallocating vector
    y3 = zeros(1,len);                                          %Preallocating vector

    for i = 1:len                                               %Loop size
        if (angle(x(i))>=-pi/8)&&(angle(x(i))<=pi/8)            %Mapping i=0 to 000
            y1(i) = 0;
            y2(i) = 0;
            y3(i) = 0;
        elseif (angle(x(i))>pi/8)&&(angle(x(i))<=3*pi/8)        %Mapping i=1 to 001
            y1(i) = 0;
            y2(i) = 0;
            y3(i) = 1;
        elseif (angle(x(i))>3*pi/8)&&(angle(x(i))<=5*pi/8)      %Mapping i=2 to 010
            y1(i) = 0;
            y2(i) = 1;
            y3(i) = 0;
        elseif (angle(x(i))>5*pi/8)&&(angle(x(i))<=7*pi/8)      %Mapping i=3 to 011
            y1(i) = 0;
            y2(i) = 1;
            y3(i) = 1;
        elseif (angle(x(i))>7*pi/8)||(angle(x(i))<=-7*pi/8)     %Mapping i=4 to 100
            y1(i) = 1;
            y2(i) = 0;
            y3(i) = 0;
        elseif (angle(x(i))>-7*pi/8)&&(angle(x(i))<=-5*pi/8)    %Mapping i=5 to 101
            y1(i) = 1;
            y2(i) = 0;
            y3(i) = 1;
        elseif (angle(x(i))>-5*pi/8)&&(angle(x(i))<=-3*pi/8)    %Mapping i=6 to 110
            y1(i) = 1;
            y2(i) = 1;
            y3(i) = 0;
        elseif (angle(x(i))>-3*pi/8)&&(angle(x(i))<=-pi/8)      %Mapping i=7 to 111
            y1(i) = 1;
            y2(i) = 1;
            y3(i) = 1;
        end                                                     %End of conditions
    end                                                         %End of loop

end
