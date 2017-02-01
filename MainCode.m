clear all
%EULER method

load 'variables.mat'

nsh=500 
bact=1:20

dt=1
T=43000
a=(1:dt:T);

F=zeros(length (bact),1);
f=zeros(length (bact),1);
D=zeros(length (bact),1);
d=zeros(length (bact),1);
E=zeros(length (bact),1);
e=zeros(length (bact),1);



matF=zeros(round(length(a)/nsh),length (bact)); %integral number needed
matf=zeros(round(length(a)/nsh),length (bact));
matD=zeros(round(length(a)/nsh),length (bact));
matd=zeros(round(length(a)/nsh),length (bact));
matE=zeros(round(length(a)/nsh),length (bact));
mate=zeros(round(length(a)/nsh),length (bact));
n=0;
for t=1:length(a);
  %(2)
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
      
      plot(bact, matF(n,:) )
      matf(n,:)=f;
      matD(n, :)=D;
      matd(n, :)=d;
      matE(n, :)=E;
      mate(n, :)=e;
      end
 
  end
 
end
imshow(matF,[0,max(max(matF))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])%Change the size of the plot


subplot(2,3,1)
imshow(matF,[0,max(max(matF))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(2,3,2)
imshow(matD,[0,max(max(matD))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(2,3,3)
imshow(matE,[0,max(max(matE))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(2,3,4)
imshow(matf,[0,max(max(matf))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(2,3,5)
imshow(matd,[0,max(max(matd))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])

subplot(2,3,6)
imshow(mate,[0,max(max(mate))])
colormap(jet)
set(gca,'DataAspectRatio',[1,10,1])


%Display F,D,E together
c_fig(:,:,1)=matF;
c_fig(:,:,2)=matD;
c_fig(:,:,3)=matE;
figure;
imshow (c_fig)
set(gca,'DataAspectRatio',[1,10,1])
