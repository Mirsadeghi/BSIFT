    function Pt = TransformPoint(P1 , H)
        P1 = [P1;ones(1,size(P1,2))];
        Ptmp = H * P1;
        Pt = [Ptmp(1,:)./Ptmp(3,:) ; Ptmp(2,:)./Ptmp(3,:)];
    end