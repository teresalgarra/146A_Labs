function [rc,time_axis] = srrc_pulse(a,m,l)

    %a is the excess bandwidth
    %m is the number of samples per T
    %length is the number of Ts taken
    length_os = floor(l*m);                     %Number of samples
    z = cumsum(ones(length_os,1))/m;            %Time vector
    N1 = 4*a*cos(pi*(1+a)*z);                   %Numerator Term 1
    N2 = pi*(1-a)*sinc((1-a)*z);                %Numerator Term 2
    D = (pi.*(1-16.*a.^2.*z.^2));               %Denominator
    zerotest = m/(4*a);                         %Location of zeros
    n = 1/4/a + 1e-7;                           %Value for the zeros
    if (zerotest == floor(zerotest))            %Zeros loop
        N1(zerotest) = 4*a*cos(pi*(1+a)*n);     %Numerator Term 1
        N2(zerotest) = pi*(1-a)*sinc((1-a)*n);  %Numerator Term 2
        D(zerotest) = (pi.*(1-16.*a.^2.*n.^2)); %Denominator
    end                                         %End of the loop
    G = (N1+N2)./D;                             %One side of peak
    rc = [flipud(G);1;G];                       %Peak reflection
    time_axis = [flipud(-z);0;z];               %Time vector
end