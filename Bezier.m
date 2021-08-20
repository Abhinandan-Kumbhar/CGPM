function [px,py,pz]=Bezier()
  % Designs Bezier bicubic surface from 16 control points
  close all;clear;clc
  pkg load io windows
  p = xlsread('bzdata.xlsx','Sheet1');
  f = p(:,1);
  g = p(:,2);
  h = p(:,3);
  cpx = reshape(f,[4,4])';
  cpy = reshape(g,[4,4])';
  cpz = reshape(h,[4,4])';
  ux = 0:0.02:1;
  vx = 0:0.02:1;
  [m,n] = size(ux);
  for i = 1:n
      for j = 1:n
       u = ux(i);
       v = vx(j);
      #bu = [(1-u)^6 6*u*(1-u)^5 15*u^2*(1-u)^4 20*u^3*(1-u)^3 15*u^4*(1-u)^2 6*u^5*(1-u)^1 u^6];
      #bv = [(1-v)^6 6*v*(1-v)^5 15*v^2*(1-v)^4 20*v^3*(1-v)^3 15*v^4*(1-v)^2 6*v^5*(1-v)^1 v^6];
       bu=[(1-u)^3 3*u*(1-u)^2 3*u*u*(1-u) u^3];
       bv=[(1-v)^3 3*v*(1-v)^2 3*v*v*(1-v) v^3];
      px(i,j)= bu*cpx*bv';
      py(i,j)= bu*cpy*bv';
      pz(i,j)= bu*cpz*bv';
      end
  end 
  
  #Reshping the X,Y,Z matrices in to vectors for plotting

   xP=[]; yP=[]; zP=[];
   xP =horzcat(xP, reshape(cpx,1,[])); 
   yP =horzcat(yP, reshape(cpy,1,[]));
   zP =horzcat(zP, reshape(cpz,1,[])); 
endfunction
