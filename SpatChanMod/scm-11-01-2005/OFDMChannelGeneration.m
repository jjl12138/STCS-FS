clc;
clear;

seed = 1;
rand('state',seed);
randn('state',seed);
% T: length of pilot
% i: number of time block
% N: number of Rx in BS
% S: number of RF chains
% ch: change of support 
%% generate channel matrix.
% num_link=1; num_path=6; num_subpath=20;
%H=scm(scmparset,linkparset(1),antparset);
N=256;
M=1;
scmpar=scmparset;
scmpar.NumBsElements=N;
scmpar.NumMsElements=M;
% scmpar.NumPaths=5; % default: 6
scmpar.NumPaths=6;
scmpar.NumSubPathsPerPath=20;
scmpar.NumTimeSamples=100;
scmpar.Scenario='urban_macro'; % default: 'none'
%scmpar.Scenario='urban_micro'; % default: 'none'
%scmpar.Scenario='suburban_macro'; % default: 'none'
%scmpar.BsUrbanMacroAS='fifteen'; % default: 'eight'
%scmpar.XpdIndependentPower='yes'; % default: 'no'
linkpar=linkparset;

% num_path=5;
% [H,delays,output]=scm(scmpar,linkpar,antparset);
% H1=permute(H,[2,1,3,4]);
%NUMRB=100;Nsc=12; % bandwidth=20M Hz
NUMRB=100;Nsc=5; % bandwidth=20M Hz
Ntot=NUMRB*Nsc;
%Ntot=NUMRB;
subcarrierSpacing=15e3;
Tb= 1/subcarrierSpacing;   % symbol time 
Nfft=2^(ceil(log2(Ntot)));
sampletime=1/subcarrierSpacing/Nfft;

%Velocity=0.1;
% num_path=path;
%scmpar.NumPaths=num_path; % default: 6
%linkpar.MsVelocity=Velocity; % default: 10 m/s
% linkpar.MsDirection=10;
% linkpar.MsBsDistance=50;
% linkpar.ThetaBs=0;
% linkpar.ThetaMs=0;
[H,delays,output]=scm(scmpar,linkpar,antparset);

H1=H;
%H1=permute(H,[2,1,3,4]);

sorted_delay=round(delays/sampletime);             
channel_matrix_size = size(H1);
channel_matrix_size(3) = max(sorted_delay)+1;
channel_out = zeros(channel_matrix_size);
for tap_i = 1:channel_matrix_size(3)
    tap_positions = find(sorted_delay == tap_i-1);
    if sum(tap_positions)>0
        channel_out(:,:,tap_i,:) =sum(H1(:,:,tap_positions,:),3);
    end
end





H_fft_large = fft(channel_out,Nfft,3);
% Eliminate guardband
H_fft= H_fft_large(:,:,[Nfft-Ntot/2+1:Nfft 2:Ntot/2+1],:);
%TMP = [Nfft-Ntot/2+1:Nfft 2:Ntot/2+1];


%% Full Frequency (maybe 512?) Channel Matrix
% figure(1);
% H = zeros(N, Nfft);
% for i = 1:Nfft
%     Hnew = H_fft_large(:,:,i);
%     Hnew = Hnew.';
%     H(:,i) = Hnew;
% end
% F = dftmtx(N)/sqrt(N);
% H_Ang_Fre = F'*H;
% %plot(abs(H_Ang_Fre))
% subplot(121);
% [x,y] = meshgrid(1:1:Nfft, 1:1:N);
% mesh(x,y, abs(H_Ang_Fre(:,1:1:Nfft)));
% title('(a)');
% ylabel('Number of antennas ( angular domain )');
% xlabel({'Number of subcarries';'( frequency domain )'});
% axis([1,Nfft,1,N]);
% colorbar;
% set(gca,'FontSize',20);
% 
% subplot(122);
% Fdel = dftmtx(Nfft)/sqrt(Nfft);
% H_Ang_Del = H_Ang_Fre*Fdel';
% %H_Ang_Del = ifft(H_Ang_Fre,Nfft,2);
% [x,y] = meshgrid(1:1:Nfft, 1:1:N);
% mesh(x,y, abs(H_Ang_Del(:,1:1:Nfft)));
% title('(b)');
% ylabel('Number of antennas ( angular domain )');
% xlabel({'Number of subcarries'; '( delay domain )'});
% axis([1,Nfft,1,N]);
% colorbar;
% set(gca,'FontSize',20);
% 
% set(gcf,'position',[0, 0, 1000, 800]);   

%% Pilot Frequency Channel Matrix
% figure(2);
% L = 64;
% Hp = zeros(N,L);
% for i = 1:L
%     Hnew = H_fft_large(:,:,(i-1)*8+1);
%     Hnew = Hnew.';
%     Hp(:,i) = Hnew;
% end
% F = dftmtx(N)/sqrt(N);
% Hp_Ang_Fre = F'*Hp;
% %var(Hp_Ang_Fre)
% 
% subplot(121);
% [x,y] = meshgrid(1:1:L, 1:1:N);
% surf(x,y, abs(Hp_Ang_Fre(:,1:1:L)));
% title('(a)');
% ylabel('Number of antennas ( angular domain )');
% xlabel({'Number of pilot carriers'});
% axis([1,L,1,N]);
% colorbar;
% shading interp;
% %set(gca,'FontSize',20);
% 
% subplot(122);
% Fdel = dftmtx(L)/sqrt(L);
% %Hp_Ang_Del = Hp_Ang_Fre*Fdel';
% Hp_Ang_Del = ifft(Hp_Ang_Fre,L,2);
% [x,y] = meshgrid(1:1:L, 1:1:N);
% %mesh(x,y, abs(Hp_Ang_Del(:,1:1:L)));
% surf(x,y, abs(Hp_Ang_Del(:,1:1:L)));
% title('(b)');
% ylabel('Number of antennas ( angular domain )');
% xlabel({'Number of delay taps'});
% axis([1,L,1,N]);
% colorbar;
% shading interp;
% % set(gca,'FontSize',20);
% % 
% % set(gcf,'position',[0, 0, 1000, 800]); 

%% Pilot Frequency Channel Matrix Generation
figure(3);
L = 64;
SIM = scmpar.NumTimeSamples;
Hp = zeros(N,L,SIM);
for n = 1:SIM
    for i = 1:L
        Hnew = H_fft_large(:,:,(i-1)*8+1,n);
        Hnew = Hnew.';
        Hp(:,i,n) = Hnew;
    end
end
F = dftmtx(N)/sqrt(N);
Hp_Ang_Fre = zeros(N,L,SIM);
for n = 1:SIM
    Hp_Ang_Fre(:,:,n) = F'*Hp(:,:,n);
end
%Hp_Ang_Fre = F'*Hp;
for n = 1:SIM
    tmp = Hp_Ang_Fre(:,:,n);
    Hp_Ang_Fre(:,:,n) = Hp_Ang_Fre(:,:,n)/sqrt(var(tmp(:)));
end
%var(Hp_Ang_Fre)

nSim = 1;
subplot(121);
[x,y] = meshgrid(1:1:L, 1:1:N);
surf(x,y, abs(Hp_Ang_Fre(:,1:1:L,nSim)));
title('(a)');
ylabel('Bin index in the angular domain');
xlabel({'Index of the pilot carriers'});
axis([1,L,1,N]);
colorbar;
shading interp;
set(gca,'FontSize',20);

subplot(122);
Fdel = dftmtx(L)/sqrt(L);
%Hp_Ang_Del = Hp_Ang_Fre*Fdel';

% have additional coefficients using ifft, might be inaccurate
Hp_Ang_Del = ifft(Hp_Ang_Fre,L,2);

% accurate one
Hp = zeros(N,L,SIM);
maxL = 22;
for n = 1:SIM
    for i = 1:maxL
        Hnew = channel_out(:,:,i,n);
        Hnew = Hnew.';
        Hp(:,i,n) = Hnew;
    end
end
Hp_Ang_Del = zeros(N,L,SIM);
for n = 1:SIM
    Hp_Ang_Del(:,:,n) = F'*Hp(:,:,n);
end

[x,y] = meshgrid(1:1:L, 1:1:N);
%mesh(x,y, abs(Hp_Ang_Del(:,1:1:L)));
surf(x,y, abs(Hp_Ang_Del(:,1:1:L,nSim)));
title('(b)');
ylabel('Bin index in the angular domain');
xlabel({'Index of the delay taps'});
axis([1,L,1,N]);
colorbar;
shading interp;
set(gca,'FontSize',20);
% 
set(gcf,'position',[0, 0, 1000, 800]);

% figure(5);
% HpTmp = zeros(256,100);
% for i = 1:100
%     HpTmp(:,i) = Hp_Ang_Fre(:,1,i);
% end
% plot(abs(HpTmp));

%% Sparse Cut
% Hp_Ang_Del1 = Hp_Ang_Del;
% index = find(abs(Hp_Ang_Del1)<0.2);
% Hp_Ang_Del1(index) = 0;
% Hp_Ang_Fre1 = fft(Hp_Ang_Del1,L,2);
% 
% figure(4);
% subplot(121);
% surf(x,y, abs(Hp_Ang_Fre1(:,1:1:L,nSim)));
% title('(a)');
% ylabel('Number of antennas ( angular domain )');
% xlabel({'Number of pilot carriers'});
% axis([1,L,1,N]);
% colorbar;
% shading interp;
% subplot(122);
% surf(x,y, abs(Hp_Ang_Del1(:,1:1:L,nSim)));
% title('(b)');
% ylabel('Number of antennas ( angular domain )');
% xlabel({'Number of delay taps'});
% axis([1,L,1,N]);
% colorbar;
% shading interp;
% 
save Hp_Ang_Fre Hp_Ang_Fre;
save Hp_Ang_Del Hp_Ang_Del;




