function k=KL_div(D1n,D2n)                
    eta = .00001 * ones(size(D1n));
    D1n = D1n+eta;
    D2n = D2n+eta;
    temp = D1n.*log(D1n./D2n);
    temp(isnan(temp)) = 0;
    k = sum(temp);
end