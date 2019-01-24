    function M = SKLD(D1 , D2 , lambda)
        if strcmp(class(D1),'double') == 0;D1 = double(D1);end
        if strcmp(class(D2),'double') == 0;D2 = double(D2);end
        D1n = D1./kron(sum(D1),ones(size(D1,1),1));
        D2n = D2./kron(sum(D2),ones(size(D2,1),1));
        k1 = KL_div(D1n,D2n);
        k2 = KL_div(D2n,D1n);
        k = (k1 + k2)/2;
        M = exp(-lambda * k);
    end