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