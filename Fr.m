function Fr=Fr(A,Z)
    function R=R(x,A,Z)
        R=sqrt(1+A.^2-2.*A.*cos(x)+Z.^2);
    end
Fr1=integral(@(x) (R(x,A,Z)+cos(x).*log(R(x,A,Z)+A-cos(x))).*cos(x),0,pi);
Fr2=-pi/4.*(1+log(4))+integral(@(x) (R(x,A,Z)+cos(x).*log((R(x,A,Z)+1+A)./(R(x,A,Z)+1-A))).*cos(x),0,pi);
%瑕玷考虑，三种情况，两种积分；
if Z>0
    Fr=Fr1;
elseif Z==0
    if A>1
        Fr=Fr1;
    elseif A==1
        Fr=-1.042260020194811;
    else
        Fr=Fr2;
    end
else
    Fr=Fr1;
end
end