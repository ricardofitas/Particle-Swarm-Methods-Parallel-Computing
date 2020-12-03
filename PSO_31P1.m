% application of Hybridization of Particle Swarm Optimisation (PSO) 
% with Genetic Algorithm - Parallel swarms - Parallel Off. Computing
% Ricardo Fitas, 2020

%% initial code
try
    myCluster = parcluster('local');
    parpool(myCluster.NumWorkers)
end
main_PSOst = evalin('base','main_PSOst');
if main_PSOst.selectfunc == 1
    main_PSOst = Part1(main_PSOst);
end
if main_PSOst.selectfunc == 2
    main_PSOst = Part2(main_PSOst);
end
if main_PSOst.selectfunc == 3
    main_PSOst = Part3(main_PSOst);
end
assignin('base','main_PSOst',main_PSOst)

%% functions
function main_PSOst = Part1(main_PSOst)
%% handle of variables
X_min = main_PSOst.X_min;
X_max = main_PSOst.X_max;
P = main_PSOst.P;
D = main_PSOst.D;

%% code
main_PSOst.X= round(rand(D,P).*(X_max-X_min)+X_min);
main_PSOst.xpbest = main_PSOst.X;
end

function main_PSOst = Part2(main_PSOst)
%% handle of variables
K = main_PSOst.K;
X = main_PSOst.X;
P = main_PSOst.P;
XG = main_PSOst.XG;
YG = main_PSOst.YG;
ZG = main_PSOst.ZG;
f1 = main_PSOst.f1;
iter = main_PSOst.iter;
detaux = str2func(string(main_PSOst.KDSstr{1}));
user_push1 = evalin('base','user_push1');
%% code
parfor ii = 1:P
    assignin('base','user_push1','user_push1');
    [Scherwinkel,XT,YT,ZT]=detaux(K, X(:,ii), XG, YG, ZG);
    f1(ii) = max(max(abs(Scherwinkel)));
    Scherwinkelmat(:,:,ii) = Scherwinkel;
    XTmat(:,:,ii) =  XT;
    YTmat(:,:,ii) =  YT;
    ZTmat(:,:,ii) =  ZT;
    
end
%% return of handle
if iter == 1
    [main_PSOst.gbest,k] = min(f1);
    main_PSOst.xgbest = X(:,k);
    main_PSOst.XT = XTmat(:,:,k);
    main_PSOst.YT = YTmat(:,:,k);
    main_PSOst.ZT = ZTmat(:,:,k);
    main_PSOst.Scherwinkel = Scherwinkelmat(:,:,k);
end
main_PSOst.f1 = f1;
main_PSOst.XTmat = XTmat;
main_PSOst.YTmat = YTmat;
main_PSOst.ZTmat = ZTmat;
main_PSOst.Scherwinkelmat = Scherwinkelmat;
end

function main_PSOst = Part3(main_PSOst)
%% handle of variables
K = main_PSOst.K;
X = main_PSOst.X;
P = main_PSOst.P;
XG = main_PSOst.XG;
YG = main_PSOst.YG;
ZG = main_PSOst.ZG;
XT = main_PSOst.XT;
YT = main_PSOst.YT;
ZT = main_PSOst.ZT;
Scherwinkel = main_PSOst.Scherwinkel;
f1 = main_PSOst.f1;
X_min = main_PSOst.X_min;
X_max = main_PSOst.X_max;
P = main_PSOst.P;
D = main_PSOst.D;
iter = main_PSOst.iter;
v = main_PSOst.v;
main_PSOst.omega = main_PSOst.o_max - iter*...
    (main_PSOst.o_max-main_PSOst.o_min)/(main_PSOst.maxiter);
c1 = main_PSOst.c1;
c2 = main_PSOst.c2;

%% code
K_K= 2/abs(2-c1-c2-sqrt((c1+c2)^2-4*(c1+c2)));
main_PSOst = Part2(main_PSOst);
f1 = main_PSOst.f1;
pbest = main_PSOst.pbest;
xpbest = main_PSOst.xpbest;
g = (abs(pbest)> abs(f1));
pbest(g)=f1(g);
xpbest(:,g) = X(:,g);
gbest = main_PSOst.gbest;
xgbest = main_PSOst.xgbest;

gbest2 = min(abs(pbest));
if gbest>gbest2
    gbest = gbest2;
    string1 = "Found a new gbest at iteration "+ string(iter)+": "+ string(gbest);
    message2 ="----- Notification from MATLAB: main_PSOsts PC ------" + newline + string1 + newline + "Kind regards,"+ newline + "Ricardo Fitas";
    message = "MATLAB: Notification of new gbest";
    run('sendemail.m')
    [~,k]=min(pbest);
    xgbest = X(:,k);
    Scherwinkel = main_PSOst.Scherwinkelmat(:,:,k);
    XT = main_PSOst.XTmat(:,:,k);
    YT = main_PSOst.YTmat(:,:,k);
    ZT = main_PSOst.ZTmat(:,:,k);
end
Xa = X(:, 1:round(size(X,2)/2));
Xb = X(:, (round(size(X,2)/2)+1):end);
Xb = gen_alg(round(Xb-X_min),f1((round(size(X,2)/2)+1):end),pc,pm,X_max-X_min,0);
r1 = rand(1,round(size(X,2)/2));
r2 = rand(1,round(size(X,2)/2));
v(:,1:round(size(X,2)/2)) = K*(v(:,1:round(size(X,2)/2))+c1.*r1.*...
    (xpbest(:,1:round(size(X,2)/2))...
    -Xa) + c2.*r2.*(xgbest-Xa));
Xa = Xa + v(1:round(size(X,2)/2));
X = [Xa, Xb+X_min];
X = round(min(X_max, max(X_min,X)));
rr = rand(1,P);
[~, coord] = sort(rr);
X= X(:,coord); v=v(:,coord); pbest = pbest(coord);xpbest=xpbest(:,coord);


%% return of handle
main_PSOst.X = X;
main_PSOst.v = v;
main_PSOst.gbest =gbest;
main_PSOst.pbest =pbest;
main_PSOst.xpbest =xpbest;
main_PSOst.xgbest = xgbest;
main_PSOst.XT = XT;
main_PSOst.YT = YT;
main_PSOst.ZT = ZT;
main_PSOst.Scherwinkel = Scherwinkel;
end



function [X] = gen_alg(X,f1,pc,pm,X_max,X_min)

tam_max = max([abs(X_max),abs(X_min)]); % detection of N = 2^i
Nn = ceil(log(tam_max)/log(2)); % calculation of i = log_2(N)
vect_X = zeros(size(X,2),Nn*size(X,1)); 
for ii = 1:size(X,2)
    for jj=1:size(X,1)
        evalX = X(jj,ii);
        for kk = Nn:-1:1
            vect_X(ii,kk+Nn*(jj-1))=mod(evalX,2);
            evalX = floor(evalX/2);
        end
    end
end
%% reproduction 
F = (1+f1).^(-1);
P = cumsum(F./sum(F));
R = rand(1,size(F,2));
new_vect = vect_X;
[rank_R, Count] = deal(zeros(1,size(F,2)));
for ii = 1:size(X,2)
    rank_R(ii) = size(X,2) - nnz(find(R(ii)-P<=0))+1;
    Count(rank_R(ii)) = Count(rank_R(ii)) + 1;
    vect_X(ii,:) = new_vect(rank_R(ii),:);
end
%% crossover
rP2 = rand(1,floor(size(X,2)/2));
rP2(rP2>pc)=1;rP2(rP2<=pc)=0;
iaux = 1;
for ii =1:2:size(vect_X,1)
    if rP2(iaux)==0
        rP = rand(1,size(vect_X,2));
        rP(rP>pc)=1;rP(rP<=pc)=0;
        for jj=1:size(vect_X,2)
            if rP(jj)==0
                aux= vect_X(ii,jj);
                vect_X(ii,jj)=vect_X(ii+1,jj);
                vect_X(ii+1,jj)=aux;
            end
        end
    end
end
%% mutation and final vector calculation
soma = zeros(size(X));
for ii = 1:size(vect_X,1)
    for jj= 1:size(vect_X,2)
        if rand(1,1) < pm
            if vect_X(ii,jj) == 0
                vect_X(ii,jj) = 1;
            else
                vect_X(ii,jj) = 0;
            end
        end
        kk = Nn - 1 - mod(jj-1,Nn);
        soma(1+floor((jj-1)/Nn),ii)=soma(1+floor((jj-1)/Nn),ii) + vect_X(ii,jj)*(2^(kk));
    end
end
X = soma;
end
