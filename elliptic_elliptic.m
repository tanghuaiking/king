%% PSF���� ������λ��
%% ����PSF
for nn=1:10
for mm=1:10
dx=nn*0.1;
dy=mm*0.1;
bx=10+dx;%��˹����������
by=10+dy;%��˹����������
sigmax=1;
sigmay=1;
syms x y
z=50000*exp(-(((x-bx)^2)/(2*sigmax^2)+((y-by)^2)/(2*sigmay^2)));
%z=(sin(sqrt((x-bx)^2+(y-by)^2))/sqrt((x-bx)^2+(y-by)^2))^2;
%% ������ ��ⱳ����Ϊ100���Ӳ���,���������밵����Ϊ5
noise=zeros(20);
for noix=1:20
    for noiy=1:20
noise(noix,noiy)=poissrnd(105);%%0;%;
    end
end         %;round(5*rand(20));%
%% pixel ���ĸ���Ϊ(0,0),(0,1),(1,0),(1,1)
v_all=zeros(20,20);
for ix=1:20
x1=ix;x2=x1+1;
for iy=1:20
y1=iy;y2=y1+1;
V=round(int(int(z,y,y1,y2),x,x1,x2));%(int(int(z,y,y1,y2),x,x1,x2));%round(int(int(z,y,y1,y2),x,x1,x2));
V_d=double(V);
v_all(ix,iy)=poissrnd(V_d);
end
end
v_all=v_all+noise;
%delete(gcp('nocreate'))
%surf(v_all);
%colorbar
%mesh(v_all);
%% ��Ȧ
v_psf=zeros(20,20);
for ix=1:20
for iy=1:20
   if v_all(ix,iy)>max(max(noise))
      v_psf(ix,iy)=v_all(ix,iy)-105; 
   end
end
end
lx=0;
for ix=1:20
lxi=(ix+0.5)*sum(v_psf(ix,:));
lx=lxi+lx;
end
ly=0;
for iy=1:20
lyi=(iy+0.5)*sum(v_psf(:,iy));
ly=lyi+ly;
end
cx=lx/sum(sum(v_psf));
% abs(bx-cx)
cy=ly/sum(sum(v_psf));
% abs(by-cy)
%% ��������
Q11=0;
Q22=0;
Q12=0;
for ix=1:20
for iy=1:20
   if v_psf(ix,iy)>0
      deltax=(floor(cx)-ix);
      deltay=(floor(cy)-iy);
      Q11=v_psf(ix,iy)*(deltax)^2+Q11; 
      Q22=v_psf(ix,iy)*(deltay)^2+Q22;
      Q12=v_psf(ix,iy)*(deltax)*(deltay)+Q12;    
   end
end
end
T=Q11+Q22;
e_plus(nn,mm)=(Q11-Q22)/T;
e_x(nn,mm)=2*Q12/T;

%% 
gx=1.5:1:20.5;
gy=1.5:1:20.5;
gridx=1.001:0.001:20;
gridy=1.001:0.001:20;
V_3d= interp2(gx',gy,v_psf,gridx',gridy,'cubic');
height_psf=max(max(v_psf));
[r,c]=find(abs(V_3d-height_psf/2)<0.2 );
% ���Բ׶���߷���
x=[c/1000,r/1000];
% F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,2).^2+p(3)*x(:,1)+p(4)*x(:,2)+p(5);
% p0=[1 1 1 1 1];
F=@(p,x)p(1)*x(:,1).^2+p(2)*x(:,2).^2+p(3)*x(:,1).*x(:,2)+p(4)*x(:,1)+p(5)*x(:,2)+p(6);
p0=[1 1 1 1 1 1];
% warning off
% ���ϵ������С���˷���
p=nlinfit(x,zeros(size(x,1),1),F,p0);
%  plot(x(:,1),x(:,2),'*');
% hold on;
% xmin=min(x(:,1));
% xmax=max(x(:,1));
% ymin=min(x(:,2));
% ymax=max(x(:,2));
% % ��ͼ
% ezplot(@(x,y)F(p,[x,y]),[-1+xmin,1+xmax,-1+ymin,1+ymax]);
% title('�������');
% legend('������','�������')
b=1;
a=sqrt(p(2)/p(1));
phi_semi=0;%��ֵ������������Զ��ȷ��phi�Ƕ�
e_plus2(nn,mm)=(a^2-b^2)/(a^2+b^2)*cos(2*phi_semi);
e_x2(nn,mm)=(a^2-b^2)/(a^2+b^2)*sin(2*phi_semi);

end
end






