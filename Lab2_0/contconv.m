function [y, t] = contconv(x1, x2, t1, t2, dt)

    %This function computes an approximation to continuous-time 
    %convolution. It has 5 inputs: the two signals to convolute, two
    %starting times for the samples and the spacing of the samples.
%     
%     l1 = length(x1);                            %Length of signal 1
%     l2 = length(x2);                            %Length of signal 2
%     lt = l1 + l2 -1;                            %Both lengths
%     
    y=conv(x1,x2)*dt;
    
%     h1 = zeros(1, lt);                          %Making both vectors 
%     h1(1 : l1) = x1;                            %the same length
%     h2 = zeros(1, lt);                          %Making both vectors
%     h2(1 : l2) = x2;                            %the same length
%     
%     y = zeros(1, lt);                           %Declaring the conv signal
%     h1 = fliplr(h1);                            %Reversing x1
%     
%     for n = 1 : lt                              %Convolution loop 1
%         for k = 1 : n                           %Convolution loop 2
%             y(n)=y(n)+h1(k)*h2(n-k+1)*dt;       %Convolution formula
%         end
%     end    
    
    t = (0 : length(y)-1)*dt + t1 + t2;         %Time vector

end


