clc
clear
format long
m=1;
k=zeros(1,length([-0.07:0.001:0.07]));
p=0.025:0.005:0.06;              %b
q=-0.07:0.001:0.07;              %zÖá×ø±ê

for j=1:length(p)
for i=1:length(q)
    a=Bz(0.025,0.02599,p(j),10^8,0.025,q(i));
    b=Br(0.025,0.02599,p(j),10^8,0.025,q(i));
    k(j,i)=sqrt(a.^2+b.^2);
    %k(m)=a;
end
end
[X,Y]=meshgrid(p,q);
plot3(X,Y,k)
