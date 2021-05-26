
close all
clc
tet6=(-20).*pi/180;
tet1=pi;
tet66=tet6+pi/2;
r66=3.4/10;
r2=2/10;
r1=5.3/10;
r4=2/10;
h4=0.6/10;
r5=4/10;
r2g=r2/2;
r5g=r5/2;
r4g=0.088;
r3g=4/10;
m2=1.2;
m3=4.8;
m4=1.5;
m5=2.4;
m6=0.64;
I2=m2.*(r2)^2/12;
I3=m3.*(2.*r3g)^2/12;
I4=(0.3.*((((h4)^2/12)+0.02^2/4)-(((h4)^2/12)+0.015^2/4))+0.3.*0.088^2)+(1.2.*0.18^2/1);
I5=m5.*(r5)^2/12;
I6=m6.*(((h4)^2/12)+0.02^2/4);
F6=1;
a2=0;
tet2=0:pi/200:2*pi;
w2=1;
xx=zeros(4,length(tet2)+1);
vv=zeros(4,length(tet2)+1);
aa=zeros(4,length(tet2)+1);
FF=zeros(14,length(tet2)+1);
TE=zeros(1,length(tet2));
FFF=zeros(1,length(tet2));
xxx=zeros(1,length(tet2));
vvv=zeros(1,length(tet2));
aaa=zeros(1,length(tet2));
xx(1,1)=r1-r2;
xx(2,1)=(-23).*pi/180;
xx(3,1)=9.2;
xx(4,1)=0;
for i=1:length(tet2)
%x1=r3
%x2=tet5
%x3=r6
%x4=tet3
eqns=@(x)[r1.*cos(tet1)+r2.*cos(tet2(i))+x(1).*cos(x(4));
    r1.*sin(tet1)+r2.*sin(tet2(i))+x(1).*sin(x(4));
    x(3).*cos(tet6)+r66.*cos(tet66)+r1.*cos(tet1)-r4.*cos(x(4)+pi/2)-r5.*cos(x(2));
    x(3).*sin(tet6)+r66.*sin(tet66)+r1.*sin(tet1)-r4.*sin(x(4)+pi/2)-r5.*sin(x(2))];
x0=[xx(1,i) xx(2,i) xx(3,i) xx(4,i)];
x=fsolve(eqns,x0);
xx(1,i+1)=x(1);
xx(2,i+1)=x(2);
xx(3,i+1)=x(3);
xx(4,i+1)=x(4);
%v1=w3
%v2=V3
%v3=V6
%v4=w5
A=[-x(1).*sin(x(4)) cos(x(4)) 0 0;
    x(1).*cos(x(4)) sin(x(4)) 0 0;
    r4.*sin(x(4)+pi/2) 0 cos(tet6) r5.*sin(x(2));
    -r4.*cos(x(4)+pi/2) 0 sin(tet6) -r5.*cos(x(2))];
Y=[r2.*w2.*sin(tet2(i));-r2.*w2.*cos(tet2(i));0;0];
v=A\Y;
vv(1,i+1)=v(1);
vv(2,i+1)=v(2);
vv(3,i+1)=v(3);
vv(4,i+1)=v(4);
%a1=a3
%a2=a5
%a3=A3
%a4=A6
AA=[-x(1).*sin(x(4)) 0 cos(x(4)) 0;
    x(1).*cos(x(4)) 0 sin(x(4)) 0;
    r4.*sin(x(4)+pi/2) r5.*sin(x(2)) 0 cos(tet6);
    -r4.*cos(x(4)+pi/2) -r5.*cos(x(2)) 0 sin(tet6)];
YY=[r2.*a2.*sin(tet2(i))+r2.*((w2)^2).*cos(tet2(i))+x(1).*((v(1))^2).*cos(x(4))+2.*v(2).*v(1).*sin(x(4));
    -r2.*a2.*cos(tet2(i))+r2.*((w2)^2).*sin(tet2(i))+x(1).*((v(1))^2).*sin(x(4))-2.*v(2).*v(1).*cos(x(4));
    -r4.*((v(1))^2).*cos(x(4)+pi/2)-r5.*((v(4))^2).*cos(x(2));
    -r4.*((v(1))^2).*sin(x(4)+pi/2)-r5.*((v(4))^2).*sin(x(2))];
a=(AA)\YY;
aa(1,i+1)=a(1);
aa(2,i+1)=a(2);
aa(3,i+1)=a(3);
aa(4,i+1)=a(4);
    Fi2y=m2.*r2g.*(((w2)^2).*sin(tet2(i))-(a2).*cos(tet2(i)))-10.*m2;
    Fi2x=m2.*r2g.*(((w2)^2).*cos(tet2(i))+(a2).*sin(tet2(i)));
Ti2=-I2.*a2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fi3y=m3.*(x(1)-r3g).*(((v(1))^2).*sin(x(4))-(a(1)).*cos(x(4)))-m3.*a(3).*sin(x(4))-2.*m3.*v(2).*v(1).*cos(x(4))-10.*m3;
    Fi3x=m3.*(x(1)-r3g).*(((v(1))^2).*cos(x(4))+(a(1)).*sin(x(4)))-m3.*a(3).*cos(x(4))+2.*m3.*v(2).*v(1).*sin(x(4));
Ti3=-I3.*a(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fi4y=m4.*(r4-r4g).*(-((v(1))^2).*sin(x(4)+pi/2)+(a(1)).*cos(x(4)+pi/2))-10.*m4;
    Fi4x=-m4.*(r4-r4g).*(((v(1))^2).*cos(x(4)+pi/2)+(a(1)).*sin(x(4)+pi/2));
Ti4=-I4.*a(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fi5y=m5.*r5g.*(-((v(4))^2).*sin(x(2))+(a(2)).*cos(x(2)))-10.*m5;
    Fi5x=-m5.*r5g.*(((v(4))^2).*cos(x(2))+(a(2)).*sin(x(2)));
Ti5=-I5.*a(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fi6x=-m6.*a(4).*cos(tet6);
Fi6y=-m6.*a(4).*sin(tet6)-10.*m6;
F6x=F6.*abs(cos(tet6+pi));
F6y=F6.*abs(sin(tet6+pi));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%F12x=F1
%F32x=F2
%F43=F3
%F54x=F4
%F65x=F5
%F16=F6
%F12y=F7
%F32y=F8
%T43=F9
%F54y=F10
%F65y=F11
%T2=F13
%F14x=F14
%F14y=F12
AAA=[1 1 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 1 0 0 0 0 0 0;
    0 -r2.*sin(tet2(i)) 0 0 0 0 0 r2.*cos(tet2(i)) 0 0 0 0 1 0;
    0 1 cos(x(4)+pi/2) 0 0 0 0 0 0 0 0 0 0 0;
    0 0 sin(x(4)+pi/2) 0 0 0 0 1 0 0 0 0 0 0;
    0 0 x(1).*sin(x(4)+pi/2).*cos(x(4))+x(1).*cos(x(4)+pi/2).*sin(x(4)) 0 0 0 0 0 1 0 0 0 0 0;
    0 0 cos(x(4)+pi/2) 1 0 0 0 0 0 0 0 0 0 1;
    0 0 sin(x(4)+pi/2) 0 0 0 0 0 0 1 0 1 0 0;
    0 0 0 -r4.*sin(x(4)+pi/2) 0 0 0 0 1 r4.*cos(x(4)+pi/2) 0 0 0 0;
    0 0 0 1 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 1 0 0 0;
    0 0 0 0 r5.*sin(x(2)) 0 0 0 0 0 r5.*cos(x(2)) 0 0 0;
    0 0 0 0 1 cos(tet6+pi/2) 0 0 0 0 0 0 0 0;
    0 0 0 0 0 sin(tet6+pi/2) 0 0 0 0 1 0 0 0];
YYY=[-Fi2x;
    -Fi2y;
    -(Fi2y.*r2g.*cos(tet2(i))-Fi2x.*r2g.*sin(tet2(i))+Ti2);
    -Fi3x;
    -Fi3y;
    -(r3g.*Fi3y.*cos(x(4))+r3g.*Fi3x.*sin(x(4))+Ti3);
    -Fi4x;
    -Fi4y;
    -(r4g.*Fi4y.*cos(x(4)+pi/2)-r4g.*Fi4x.*sin(x(4)+pi/2)+Ti4);
    -Fi5x;
    -Fi5y;
    -(r5g.*Fi5y.*cos(x(2))+r5g.*Fi5x.*sin(x(2))+Ti5);
    -Fi6x+F6x;
    -Fi6y-F6y];
F=AAA\YYY;
for k=1:14
    FF(k,i+1)=F(k);
end
% Fi2=sqrt(Fi2x^2+Fi2y^2);
% Fi3=sqrt(Fi3x^2+Fi3y^2);
% Fi4=sqrt(Fi4x^2+Fi4y^2);
% Fi5=sqrt(Fi5x^2+Fi5y^2);
% Fi6=sqrt(Fi6x^2+Fi6y^2);
E22=(r2-r2g).*w2.*(-sin(tet2(i)).*Fi2x+cos(tet2(i)).*Fi2y)+(Ti2).*w2;
E3=(x(1)-r3g).*v(1).*(-sin(x(4)).*Fi3x+cos(x(4)).*Fi3y)+(Ti3).*v(1)+v(2).*(cos(x(4)).*Fi3x+sin(x(4)).*Fi3y);
E4=(r4-r4g).*v(1).*(sin(x(4)+pi/2).*Fi4x-cos(x(4)+pi/2).*Fi4y)+(Ti4).*v(1);
E5=(r5-r5g).*v(4).*(sin(x(2)).*Fi5x-cos(x(2)).*Fi5y)+(Ti5).*v(4);
E6=v(3).*(cos(tet6).*(Fi6x)+sin(tet6).*(Fi6y))-v(3).*F6;
Etotal=E3+E4+E5+E6+E22;
TE(i)=-Etotal/(w2);
end
clc
for j=1:length(tet2)
    xxx(j)=xx(1,j+1);
end
plot(tet2,xxx)
xlabel('teta 2 (rad)')
ylabel('r3 (m)')
figure
for j=1:length(tet2)
    xxx(j)=xx(2,j+1);
end
plot(tet2,xxx)
xlabel('teta 2 (rad)')
ylabel('teta 5 (rad)')
figure
for j=1:length(tet2)
    xxx(j)=xx(3,j+1);
end
plot(tet2,xxx)
xlabel('teta 2 (rad)')
ylabel('r6 (m)')
figure
for j=1:length(tet2)
    xxx(j)=xx(4,j+1);
end
plot(tet2,xxx)
xlabel('teta 2 (rad)')
ylabel('teta 3 (rad)')
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(tet2)
    vvv(j)=vv(1,j+1);
end
plot(tet2,vvv)
xlabel('teta 2 (rad)')
ylabel('w3 (rad/s)')
figure
for j=1:length(tet2)
    vvv(j)=vv(2,j+1);
end
plot(tet2,vvv)
xlabel('teta 2 (rad)')
ylabel('V3 (m/s)')
figure
for j=1:length(tet2)
    vvv(j)=vv(3,j+1);
end
plot(tet2,vvv)
xlabel('teta 2 (rad)')
ylabel('V6 (m/S)')
figure
for j=1:length(tet2)
    vvv(j)=vv(4,j+1);
end
plot(tet2,vvv)
xlabel('teta 2 (rad)')
ylabel('w5 (rad/s)')
figure
%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(tet2)
    aaa(j)=aa(1,j+1);
end
plot(tet2,aaa)
xlabel('teta 2 (rad)')
ylabel('a3 (rad/s^2)')
figure
for j=1:length(tet2)
    aaa(j)=aa(2,j+1);
end
plot(tet2,aaa)
xlabel('teta 2 (rad)')
ylabel('a5 (rad/s^2)')
figure
for j=1:length(tet2)
    aaa(j)=aa(3,j+1);
end
plot(tet2,aaa)
xlabel('teta 2 (rad)')
ylabel('A3 (m/s^2)')
figure
for j=1:length(tet2)
    aaa(j)=aa(4,j+1);
end
plot(tet2,aaa)
xlabel('teta 2 (rad)')
ylabel('A6 (m/s^2)')
figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(tet2)
    FFF(j)=FF(1,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F12x (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(7,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F12y (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(2,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F32x (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(8,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F32y (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(3,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F43 (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(9,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('T43 (N.m)')
figure
for j=1:length(tet2)
    FFF(j)=FF(4,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F54x (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(10,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F54y (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(5,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F65x (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(11,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F65y (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(6,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F16 (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(14,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F14x (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(12,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('F14y (N)')
figure
for j=1:length(tet2)
    FFF(j)=FF(13,j+1);
end
plot(tet2,FFF)
xlabel('teta 2 (rad)')
ylabel('T2 (N.m)')
figure
plot(tet2,TE)
xlabel('teta 2 (rad)')
ylabel('TE2 (N.m)')