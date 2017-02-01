clear all
%EULER method
x0=0;
x02=0;
b=0
load 'variables.mat'
% x0=10;
% x02=11;
% b=16;

nsh=500 
bact=1:20

dt=1
T=80000
a=(1:dt:T);
u=12;

F=zeros(length (bact)+u,1);
f=zeros(length (bact)+u,1);
D=zeros(length (bact)+u,1);
d=zeros(length (bact)+u,1);
E=zeros(length (bact)+u,1);
e=zeros(length (bact)+u,1);
DNA1=zeros(length (bact),1);

DNA2=zeros(length (bact),1);
DNA2(9)=8;
DNA2(10)=10;
DNA2(11)=8;


matF=zeros(round(length(a)/nsh),length (bact)+u); %integral number needed
matf=zeros(round(length(a)/nsh),length (bact)+u);
matD=zeros(round(length(a)/nsh),length (bact)+u);
matd=zeros(round(length(a)/nsh),length (bact)+u);
matE=zeros(round(length(a)/nsh),length (bact)+u);
mate=zeros(round(length(a)/nsh),length (bact)+u);
matDNA1=zeros(round(length(a)/nsh),length (bact)+u);
matDNA2=zeros(round(length(a)/nsh),length (bact)+u);
n=0;

gama=0.00008;
gama2=0.008;

m=1;
t=0;
step=16000;
maxbact=length (bact)+u+1;
while length(bact)<maxbact && t<T;
    t=t+1;
 
  if t>m*step;
     m=m+1;
     bact=1: (length(bact)+2 )
  else
  for x=bact %need to add border conditions
      if x==min(bact)
          
          dF(x)=DifF*2*(F(x+1)-F(x));
          df(x)=Diff*2*(f(x+1)-f(x));
          dD(x)=DifD*2*(D(x+1)-D(x));
          dd(x)=Difd*2*(d(x+1)-d(x));
          dE(x)=DifE*2*(E(x+1)-E(x));
          de(x)=Dife*2*(e(x+1)-e(x));
          
      elseif x==max(bact)
         dF(x)=DifF*2*(F(x-1)-F(x));
         df(x)=Diff*2*(f(x-1)-f(x)); 
         dD(x)=DifD*2*(D(x-1)-D(x));
         dd(x)=Difd*2*(d(x-1)-d(x));
         dE(x)=DifE*2*(E(x-1)-E(x));
         de(x)=Dife*2*(e(x-1)-e(x));
      else
        dF(x)=DifF*(F(x-1)+F(x+1)-2*F(x));
        df(x)=Diff*(f(x-1)+f(x+1)-2*f(x));
        dD(x)=DifD*(D(x-1)+D(x+1)-2*D(x));
        dd(x)=Difd*(d(x-1)+d(x+1)-2*d(x));
        dE(x)=DifE*(E(x-1)+E(x+1)-2*E(x));
        de(x)=Dife*(e(x-1)+e(x+1)-2*e(x));
      end
%         U=Gauss([1:20]',x0,b);
%         U2=Gauss([1:20]',x02,b);
% %       
        if length (bact)==24
            x0=17;
            x02=10;
            b=16;
            
            
        end
        if length (bact)>25
     
      U=Gauss([1:maxbact-1]',x0,b).*D*8;
      U2=Gauss([1:maxbact-1]',x02,b).*D*8;
      
      Force=diff(U);
      FT=sum(Force);
     
      
       
     
      Force2=diff(U2);
      FT2=sum(Force2);
      
      
    if x0<0
          x0=1;
      elseif x0>32
          x0=32;
      else
          
        x0=x0+FT2*gama;
      end
       
      if x02<0
          x02=1;
      elseif x02>32
          x02=32;
      else
          
        x02=x02+FT2*gama;
      end
       end
      noiseF= F(x)*(1+randi([-100 100],1,1)*0.0001);
       prodF=roF *f(x)*(((F(x)*F(x))+sigF)/(1+kF*F(x)*F(x)));
       degF=miuF*F(x)+miuDF*D(x)*F(x);
       F(x)=noiseF+(prodF-degF+dF(x))*dt;
      
       
       noisef=f(x)*(1+randi([-100 100],1,1)*0.0001);
       prodf=(sigf-roF*f(x)*(F(x).^2+sigF))/(1+kF*F(x).^2) ;
       degf=miuf*f(x);
       f(x)=noisef+(prodf-degf+df(x))*dt;
       
       noiseD=D(x)*(1+randi([-100 100],1,1)*0.0001);
       prodD=roD*d(x)*(D(x).^2 +sigD);
       degD=miuD*D(x)+miuDE*D(x)*E(x);
       D(x)=noiseD+(prodD-degD+dD(x))*dt;
       
       noised=d(x)*(1+randi([-100 100],1,1)*0.0001);
       prodd=sigd-(roD*d(x)*((D(x).^2)+sigD)) ;
       degd=miud*d(x);
       d(x)= noised+(prodd-degd+dd(x))* dt; 
      
       noiseE=E(x)*(1+randi([-100 100],1,1)*0.0001);
       prodE=roE*e(x)*((D(x)* ((E(x).^2)+sigE))/((1+kDE*D(x).^2)*(1+kE*E(x).^2)));
       degE=miuE*E(x);
       E(x)=noiseE+(prodE-degE+ dE(x))*dt;
       
       
       noisee=e(x)*(1+randi([-100 100],1,1)*0.0001);
       prode=sige-(roE*e(x)*((D(x)* ((E(x).^2) +sigE))/((1+kDE*(D(x).^2))*(1+kE*(E(x).^2)))));
       dege=miue*e(x);
       e(x)= noisee+(prode-dege+de(x))*dt;

      if mod(t,nsh)==0
      n=n+1;
      matF(n,:)=F;
      matf(n,:)=f;
      matD(n, :)=D;
      matd(n, :)=d;
      matE(n, :)=E;
      mate(n, :)=e;
      matDNA1(n,:)=Gauss([1:maxbact-1]',x0,b);
%      matDNA1(n, :)=DNA1;
      matDNA2(n, :)=Gauss([1:maxbact-1]',x02,b);
      end
      
          
  end
 
  end
end

% imshow(matF,[0,max(max(matF))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])%Change the size of the plot
% 
% 
% subplot(2,3,1)
% imshow(matF,[0,max(max(matF))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])
% 
% subplot(2,3,2)
% imshow(matD,[0,max(max(matD))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])
% 
% subplot(2,3,3)
% imshow(matE,[0,max(max(matE))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])
% 
% subplot(2,3,4)
% imshow(matf,[0,max(max(matf))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])
% 
% subplot(2,3,5)
% imshow(matd,[0,max(max(matd))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])
% 
% subplot(2,3,6)
% imshow(mate,[0,max(max(mate))])
% colormap(jet)
% set(gca,'DataAspectRatio',[1,10,1])


%Display F,D,E together
c_fig(:,:,1)=matF;
c_fig(:,:,2)=matD;
c_fig(:,:,3)=matE;
figure(1);
imshow (c_fig)
set(gca,'DataAspectRatio',[1,10,1])

figure(2)
subplot(1,3,1)
imshow(matDNA1,[0,max(max(matDNA1))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(1,3,2)
imshow(matDNA2,[0,max(max(matDNA2))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(1,3,3)
imshow(matD,[0,max(max(matD))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

