close all;
clc;close all;clear;

#if using bezier surface put B=1 or put B=0;

#for Bezier surface
B=1;
if B==1
  [px,py,pz]=Bezier();
endif
#for surface from STL
if B==0
  [px,py,pz]=stl();
endif
#Plotting surface and surface normal
figure;
surf(px,py,pz)
hold on;
surfnorm(px,py,pz)
hold off;
[Nx,Ny,Nz]=surfnorm(px,py,pz);
view(3); box;  
title('\bf Surface and surface normal')

#pre-processing for creating dataset for K means
PX = reshape(px,size(px)(1)*size(px)(2),1);
PY = reshape(py,size(py)(1)*size(py)(2),1);
PZ = reshape(pz,size(pz)(1)*size(pz)(2),1);
NX = reshape(Nx,size(Nx)(1)*size(Nx)(2),1);
NY = reshape(Ny,size(Ny)(1)*size(Ny)(2),1);
NZ = reshape(Nz,size(Nz)(1)*size(Nz)(2),1);

#feature normalization
Normal=[NX NY NZ]; 
Normals=Normal./sqrt(sum(Normal.*Normal,2));
point=[PX,PY,PZ];
points=(point-min(point))./(max(point)-min(point));
Dataset=[points Normals];

#Kmeans clustering and deciding optimum number of clusters
K = 1;            # INitial K
max_iters = 50;   # Max iterations for K  means
naio=30;         #Threshold value of NAI
outlying=5000;
O=25;
while(outlying>O)
  NAI=[];
  #initiate centroids
  initial_centroids = kMeansInitCentroids(Dataset, K);
  #Run K-Means
  [centroids, idx] = runkMeans(Dataset, initial_centroids, max_iters);
  for k =1:size(PX)(1)
    h = idx(k);
    H=Dataset(k,4:6);
    H=H/sqrt(H*H');    #unit normal
    L=centroids(h,4:6);
    L=L/sqrt(L*L');    #unit centroid normal
    nai = acosd(H*L'); #Angle between normal at point and normal at centroid 
    NAI=[NAI;nai];
  endfor
  Checknai=NAI-naio;  
  outlying=sum(Checknai>0);
  printf('\nNo of points lying outside the threshold NAI of cluster\n')
  outlying
  printf('Total number of clusters formed\n')
  K
  if outlying>O
    K=K+1;
    printf('still some of the points have normals outside the threshold NAI limit, taking more clusters\n')
  endif
  
#plotting clusters
  figure;
  colors =['r','b','g','c','m','k'];
  colors=[colors,colors,colors,colors,colors];
  for k =1:K
    f=find(idx==k);
    scatter3(PX(f),PY(f),PZ(f),colors(k))
    hold on;
  endfor 
endwhile
figure;
colors =['r','b','g','c','m','k'];
colors=[colors,colors,colors,colors,colors];
for k =1:K
  f=find(idx==k);
  scatter3(PX(f),PY(f),PZ(f),colors(k))
  hold on;
endfor

#plotting points in 2D
figure;
colors =['r','b','g','c','m','k'];
colors=[colors,colors,colors,colors,colors];
for k =1:K
  f=find(idx==k);
  scatter(PX(f),PY(f),colors(k))
  hold on;
endfor

#Locating points NBG/COBG/CBG
NBG=[];
COBG=[];
CBG=[];
Four=[];
for k =1:size(PX)(1)
  d=Dataset(:,1:2)-Dataset(k,1:2);
  D=diag(d*d');
  closest=sort(D)(1:5);
  index=[];
  for i=1:5
    index=[index;find(D==closest(i))];
  endfor
  NC=size(unique(idx(index)))(1);
  if NC==1;
    NBG=[NBG;k];
  endif
  if NC==2;
    COBG=[COBG;k,unique(idx(index))'];
  endif
  if NC==3;
    CBG=[CBG;k,unique(idx(index))'];
  endif
  if NC==4;
    Four=[Four;k,unique(idx(index))'];
  endif
endfor
Boundary=[COBG(:,1);CBG(:,1)];
OBG=[];
x=[max(PX),min(PX)];  
y=[max(PY),min(PY)];  
for k =1:size(PX)(1)
  if PX(k)>=x(1)-0.01||PX(k)==x(2)||PY(k)>=y(1)-0.01||PY(k)==y(2);
    OBG=[OBG k];
  endif;  
endfor
 
figure;
hold on;
scatter(PX(COBG),PY(COBG),50,"filled")
hold on;
scatter(PX(CBG),PY(CBG),50,"filled")
hold on;
scatter(PX(OBG),PY(OBG),50,"filled")
hold off;
title('Locating points shared by NBG/COBG/three or four clusters')
legend('COBG','CBG','OBG')

#Identifying non machinable region
U=[];
R=0.08;
L0=1*R;
newidx=idx;
for r=1:K
  cluster=[];
  cluster=intersect(union(OBG,Boundary),[1:size(PX,1)](idx==r));
  U=[U;Dataset(cluster,1:3),r*ones(size(cluster,1),1)];
  I=setdiff((1:size(PX,1))(idx==r),cluster);
  P=intersect(Boundary,[1:size(PX,1)](idx==r));
  inside=[Dataset(I',1:2),I'];
  C=[];
  M=[];
  F=[];
  Inside=union(inside(:,3),cluster);
  Inside=[Dataset(Inside,1:2),Inside];
  for m =1:size(inside,1)
    t=Dataset(P,1:2)-inside(m,1:2);
    T=sqrt(diag(t*t'));
    p1=T==sort(T)(1);
    p2=T==sort(T)(2);
    e=Dataset(P(p1),1:2)-Dataset(P(p2),1:2);
    E=sqrt(e*e');
    if sum(T<R)>1&E>L0&E<2*R
        C=[C;inside(m,3)];
        M=[P(p1);P(p2)];
        x1= Dataset(M(1),1);
        y1= Dataset(M(1),2);
        x2= Dataset(M(2),1);
        y2= Dataset(M(2),2);
        [xc]=centroids(r,1);
        [yc]=centroids(r,2);
        g=yc-y1-(y2-y1)/(x2-x1)*(xc-x1);
        for i =1:size(Inside,1)
          x=Inside(i,1);
          y=Inside(i,2);
          v=y-y1-(y2-y1)/(x2-x1)*(x-x1);
          if v*g<0
            F=[F;Inside(i,3)];
            for b=1:size(F,1)              
              L=centroids(:,1:3)-Dataset(F(b),1:3);
              l=sqrt(diag(L*L'));
              I=l==sort(l)(1);
              newidx(F(b))=find(I==1);
            endfor  
          endif
        endfor
    endif
  endfor
  figure;
  scatter(PX(P),PY(P),'b')
  hold on;
  scatter(PX(F),PY(F),'r');
  axis([0,1,0,1]); 
  hold on;
endfor

#plotting re-clustered surface
figure;
colors =['r','b','g','c','m','k'];
colors=[colors,colors,colors,colors,colors];
for k =1:K
  f=find(newidx==k);
  scatter(PX(f),PY(f),colors(k))
  hold on;
endfor