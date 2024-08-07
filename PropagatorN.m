function[ p ] = PropagatorN(N,lambda,area,z)
p = zeros(N,N);
for ii = 1:N
    for jj = 1:N
        alpha = lambda*( ii - N/2-1)/area;
        beta = lambda*( jj - N/2-1)/area;
        if((alpha^2 + beta^2)<=1)
        p( ii , jj ) = exp(-2*pi*1i*z*sqrt(1 - alpha^2 - beta^2)/lambda);
        end  % if
    end
end 