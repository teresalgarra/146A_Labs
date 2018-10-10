function [y, t] = contconv(x1, x2, t1, t2, dt)

    %This function computes an approximation to continuous-time 
    %convolution. It has 5 inputs: the two signals to convolute, two
    %starting times for the samples and the spacing of the samples.
    
    l1 = length(x1);                                %Length of signal 1
    l2 = length(x2);                                %Length of signal 2
    lt = l1 + l2 -1;                                %Both lengths
    
    %Making both vectors the same length
    
    h1 = zeros(1, lt);
    h1(1 : l1) = x1;

    h2 = zeros(1, lt);
    h2(1 : l2) = x2;
    y = zeros(lt, 1);                              %Declaring the conv signal
    x1 = fliplr(x1);
    
    for n = 1 : lt                                 %Convolution loop 1
        for k = 1 : n                             %Convolution loop 2
            v(k)=h1(k)*h2(n-k+1);
        end
        y(n)=sum(v)/100;
    end    
    
    t = (0 : length(y)-1)*dt + t1 + t2;                 %Vector length

end


