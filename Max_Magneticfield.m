
function []=Max_Magneticfield()
a1=input('请输入螺线管内半径（m）： ');
a2=input('请输入螺线管外半径(m)： ');
b=input('请输入螺线管半高(m)： ');
N=input('请输入螺线管匝数： ');
I=input('请输入螺线管电流(A)： ');
J=N.*I./(2*b.*(a2-a1));
    if b<0.5
        b_vec=0:0.0001:b;
        len=length(b_vec);
        B=zeros(len,1);
        for i=1:len
            B_r = Br(a1,a2,b,J,a1,b_vec(i));
            B_z = Bz(a1,a2,b,J,a1,b_vec(i));
            B(i) = sqrt(B_r.^2+B_z.^2);
        end
    [val,pos]=max(B);
    fprintf("最大磁场点为:（%f,%f）和（%f,-%f）\n最大磁场为：B_max=%fT\n",a1,(pos-1)/10000,a1,-(pos-1)/10000,val)
    else                    %如果螺线管较长则增大步长，提高速度
        b_vec=0:0.001:b;        
        len=length(b_vec);
        B=zeros(len,1);
        for i=1:len
            B_r = Br(a1,a2,b,J,a1,b_vec(i));
            B_z = Bz(a1,a2,b,J,a1,b_vec(i));
            B(i) = sqrt(B_r.^2+B_z.^2);
        end
    [val,pos]=max(B);
    fprintf("最大磁场点为:（%f,%f）和（%f,-%f）\n最大磁场为：B_max=%fT\n",a1,(pos-1)/1000,a1,-(pos-1)/1000,val)
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
 
 