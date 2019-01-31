function Bm=Bm(a1,a2,b,N,I)
J=N.*I./(2*b.*(a2-a1));
rc = a1-1e-4;
zc = 0;
Bcr = Br(a1,a2,b,J,rc,zc);
Bcz = Bz(a1,a2,b,J,rc,zc);
Bc = sqrt(Bcr.^2+Bcz.^2);
re = (2*a1+a2)/3;
ze = b+1e-4;
Ber = Br(a1,a2,b,J,re,ze);
Bez = Bz(a1,a2,b,J,re,ze);
Be = sqrt(Ber.^2+Bez.^2);
    if Bc>=Be
    Bm = Bc;
    else
    Bm = Be;
    end
end

function Br=Br(a1,a2,b,J,r,z)
mu=4*pi*1e-7;
w=mu*J.*r./(2*pi);
A1=a1./r;
A2=a2./r;
Z1=(-b-z)./r;
Z2=(b-z)./r;
if r==0
    Br=0;
else
Br=w.*(Fr(A1,Z1)-Fr(A1,Z2)-Fr(A2,Z1)+Fr(A2,Z2));
end
end

function Bz=Bz(a1,a2,b,J,r,z)
mu=4*pi*1e-7;
w=mu*J.*r./(2*pi);
A1=a1./r;
A2=a2./r;
Z1=(-b-z)./r;
Z2=(b-z)./r;
    if r==0
    Bz=mu.*J./2.*((b+z).*log((a2+sqrt(a2.^2+(b+z).^2))./(a1+sqrt(a1.^2+(b+z).^2)))+...
        (b-z).*log((a2+sqrt(a2.^2+(b-z).^2))./(a1+sqrt(a1.^2+(b-z).^2))));
    else
    Bz =w.*(Fz(A1,Z1)-Fz(A1,Z2)-Fz(A2,Z1)+Fz(A2,Z2));
    end
end

function Fr=Fr(A,Z)
    function R=R(x,A,Z)
        R=sqrt(1+A.^2-2*A.*cos(x)+Z.^2);
    end
Fr1=integral(@(x) (R(x,A,Z)+cos(x).*log(R(x,A,Z)+A-cos(x))).*cos(x),0,pi);
Fr2=-pi/4*(1+log(4))+integral(@(x) (R(x,A,Z)+cos(x).*log((R(x,A,Z)+1+A)./(R(x,A,Z)+1-A))).*cos(x),0,pi);
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

function Fz=Fz(A,Z)
Z1=abs(Z);
int1 = integral(@(x) Z1.*log(sqrt(1+A.^2-2*A.*cos(x)+Z1.^2)+A-cos(x)),0,pi);
int2 = integral(@(x) -cos(x).*log(sqrt(1+A.^2-2*A.*cos(x)+Z1.^2)+Z1),0,pi);
int3 = integral(@(x) -sin(x).*atan((Z1.*(A-cos(x)))./(sqrt(1+A.^2-2*A.*cos(x)+Z1.^2).*sin(x))),0+1e-9,pi-1e-9);
if A>1
    wA=-pi/(2.*A);
else
    wA=-pi/2.*A;
end
Fz1=int1+wA+int2+int3;
if Z>0
    Fz=Fz1;
     elseif Z==0 
     Fz=0; 
else
    Fz=-Fz1;    
end
end
 
 