
clear; clc;
load dates;
load names;
svf = load('svfmeansF20A1.txt');
svy = load('svymeansF20A1.txt');
%svf(622,:) = [];
%svy(622,:) = [];
load ferrorsF20A1;

%type fbeta.csv
%a = readmatrix('fbeta2.csv');
%a(1,:) = [];
%type ybeta.csv
%b = readmatrix('ybeta2.csv');
%b(1,:) = [];


% Compute objects from predictors
h   = 12;
fb  = sparse(fbetas);
thf = [svf(1,:).*(1-svf(2,:));svf(2,:);svf(3,:).^2];
xf  = svf(4:end-3,:);
gf  = svf(end-3+1:end,:);
[evf,phif] = compute_uf(xf,thf,fb,h);

% Compute uncertainty
T = 618;
N = 132;
ut    = zeros(T+1,N,h);
for i = 1:N
    tic;
    yb  = sparse(ybetas(i,:));
    thy = [svy(1,i).*(1-svy(2,i));svy(2,i);svy(3,i).^2];
    xy  = svy(4:end-3,i);
    ut(:,i,:) = compute_uy(xy,thy,yb,py,evf,phif);
    fprintf('Series %d, Elapsed Time = %0.4f \n',i,toc);
end
gy = svy(end-3+1:end,:);
save utF20A1 dates ut
save gewekeF20A1 dates gy gf names