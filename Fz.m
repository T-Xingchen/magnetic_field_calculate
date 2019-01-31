function Fz=Fz(A,Z)
Z1=abs(Z);
int1 = integral(@(x) Z1.*log(sqrt(1+A.^2-2*A.*cos(x)+Z1.^2)+A-cos(x)),0,pi);
int2 = integral(@(x) -cos(x).*log(sqrt(1+A.^2-2.*A.*cos(x)+Z1.^2)+Z1),0,pi);
int3 = integral(@(x) -sin(x).*atan((Z1.*(A-cos(x)))./(sqrt(1+A.^2-2*A.*cos(x)+Z1.^2).*sin(x))),0+1e-9,pi-1e-9);
if A>1
    wA=-pi./(2.*A);
else
    wA=-pi./2.*A;
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