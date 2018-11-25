function [r] = randbit(len_in)
    
    %This function generates random bits taking values in {0,1}
    
    r = randi([0,1],len_in,1);                     %Random coefficients

end