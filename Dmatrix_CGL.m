%% D Matrix for CGL nodes
% 
function D = Dmatrix_CGL(xk)

    N = length(xk) - 1;
    D = zeros(N+1, N+1);

    for j = 1:N+1
        for k = 1:N+1
        if (j == k) && (k == 1)
            D(k,j) = (2*(N^2)+1)/6;
        elseif (j == k) && (k == N+1)
            D(k,j) = -(2*(N^2)+1)/6;
        elseif (j == k) 
           D(k,j) = -xk(k)/(2*(1-xk(k)^2));
        elseif ne(j,k)
            if (j==1) || (j==N+1)
                if (k == 1) || (k==N+1)
                    D(k,j) = (((-1)^(j+k))/(xk(k)-xk(j)));
                else 
                    D(k,j) = 0.5*(((-1)^(j+k))/(xk(k)-xk(j)));
                end
            elseif (k == 1) || (k==N+1)
                D(k,j) = 2*(((-1)^(j+k))/(xk(k)-xk(j)));
            else 
                D(k,j) = (((-1)^(j+k))/(xk(k)-xk(j)));
            end 
        end 
        end
    end
    