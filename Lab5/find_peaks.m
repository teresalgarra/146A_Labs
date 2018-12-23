function [peaks,peak_indices] = find_peaks(row_vector)
    A = [0 row_vector 0];
    j = 1;
    for i=1:length(A)-2
        temp=A(i:i+2);
        if(max(temp)==temp(2))
            peaks(j) = row_vector(i);

            peak_indices(j) = i;
            j = j+1;
        end
    end
end