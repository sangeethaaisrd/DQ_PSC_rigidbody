function w = cc_quad_weights(N)
%CC_QUAD_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
w = zeros(N+1,1);
w_middle = zeros(N-1,1);

if mod(N,2)==0
    w(1) = 1/(N^2-1);
else
    w(1) = 1/(N^2);
end
w(end) = w(1);

for s=1:floor(N/2)
    for j=0:floor(N/2)
        sum = (1/(1-4*(j^2)))*cos(2*pi*j*s/N);
        if j==0 || j==floor(N/2)
            sum = sum/2;
        end
        w_middle(s) = w_middle(s)+sum;
    end
    w_middle(s) = (4/N)*w_middle(s);
    w_middle(N-s) = w_middle(s);
end

w = [w(1);w_middle;w(end)];
end

