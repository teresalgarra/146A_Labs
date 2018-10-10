function [x] = signalx(t)

    %This funtion takes an array of values and gives an output vector x
    %with the values of our function evaluated at those time points.

    N = length(t);                              %Getting the length of the array
    x = zeros(N,1);                             %Creating the vector x
    j = 1;                                      %Iterator for the loop       

    for z = t                                   %Taking the loop through every value of t
        if z < -3                               %First condition
            value = 0;                          %First value
        elseif (z >= -3) && (z <= -1)           %Second condition
            value = 2.*exp(z+2);                %Second value
        elseif (z >= -1) && (z <= 4)            %Third condition
            value = 2.*exp(-z).*cos(2*pi*z);    %Third value
        elseif z > 4                            %Fourth condition
            value = 0;                          %Fourth value
        end                                     %End of the conditions

        x(j) = value;                           %Assigning the corresponding value to x
        j = j+1;                                %Moving forward with the iterator
    end                                         %End of the loop
end