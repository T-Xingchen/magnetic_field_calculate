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