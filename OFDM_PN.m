function SER=OFDM_PN(f3dB,SNR,pc,u)
% pc: Whether apply compensation algorithm
% u: Phase Noise Compensation Order

%clear all; SNR=30; f3dB = 50;pc='True'; u=20;

% OFDM parameters
K=64;                   %No. of subcarriers
M=64;                   %QAM modulation level
CP=16;                  %No. of cyclic prefix samples
NoSym=1000;             %No. of symbols
NoBits=log2(M)*NoSym;   %No. of bits per branch
BW=20e6;                %Band Width
start=4; gap=8;         %Pilot


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generated random bits of 0's and 1's
B = randi([0 1],NoBits,K);

% 64-QAM modulator
X = qammod(B,M,'UnitAveragePower',true,'InputType','bit');
%scatterplot(X(:,10));

% Adding pilots to samples
for s=1:NoSym
    Xp(s,:)=X(s,:);
    samples=X(s,:);
    pilot=(1+1j)/abs(1+1j);
    for i=start:gap:K
        samples(i)=pilot;
        Xp(s,i)=pilot;
    end
    % IFFT
    samples_IFFT = ifft(samples,K);
    samples_CP = [samples_IFFT(K-CP+1:K) samples_IFFT];
    % Transmitted symbol
    OFDMsymbols(s,:) = samples_CP;
end

% Demodulation
XI_hat=qamdemod(Xp,M,'UnitAveragePower',true,'OutputType','integer');

% All symbols
OFDMSeries =  reshape(OFDMsymbols,1,[]);
OFDMSeries = OFDMSeries./sqrt(mean(abs(OFDMSeries).^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multipath
h = [1];
channel_output = conv(OFDMSeries,h);

% AWGN
SNR_dB = SNR;
SNR_lin = 10^(SNR_dB/10);
noisePower = sum(abs(h).^2)./SNR_lin;
wNoise = sqrt(0.5*noisePower).*randn(1,length(channel_output))+1j*sqrt(0.5*noisePower).*randn(1,length(channel_output));

% Phase Noise
for s = 1:NoSym
    phi(1)=0;
    for n=1:(K+CP)
        Ts=1/BW;
        w(n) = normrnd(0,sqrt(4*pi*f3dB*Ts));
        phi(n+1) = phi(n) + w(n);
        PN(s,n)=exp(1j*phi(n));
    end
end
pNoise =  reshape(PN,1,[]);

% Received signal
y = channel_output.*pNoise + wNoise;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove cyclic prefix
samples_parallel = reshape(y(1:NoSym*(K+CP)),NoSym,[]);
samples_noCP = samples_parallel(:,CP+1:end);
channel_output_parallel = reshape(channel_output(1:NoSym*(K+CP)),NoSym,[]);
channel_output_noCP = channel_output_parallel(:,CP+1:end); 

% FFT output
for s=1:NoSym
    Y(s,:)=fft(samples_noCP(s,:),K);
    A(s,:)=fft(channel_output_noCP(s,:),K);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase Noise compensation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(pc,'True'))
    
    % Coefficients    
    for m = 1:NoSym
        % Am
        for lk=(1+u):(K-u)
            Am(lk-u,1)= A(m,lk);
            for i=1:u
                Am(lk-u,2*i)= A(m,lk-i);
                Am(lk-u,2*i+1)= A(m,lk+i);
            end
        end
        
        % RImIm
        for k=1:K
            for l=1:K
                Temp(k,l)=exp(-0.5*abs(k-l)*4*pi*f3dB*Ts);
            end
        end
        RImIm=(1/K^2)*fliplr(fft2(Temp));
        

        % RJmJm
        for n=1:2*u+1
            for p=1:2*u+1
                RJmJm(n,p)=RImIm(K/2-u-1+n,K/2-u-1+p);
            end
        end
        %RJmJm=(1/K^2)*fft2(Temp,2*u+1,2*u+1);
        
        % REmEm
        %{
        Sum=0;
        for i=1:K
            if (i<=(K/2-u))||(i>(K/2+u+1))
                Sum=Sum+RImIm(i,i);
            end
        end
        %}
        REmEm = (4*pi*f3dB*Ts)*eye(K-2*u);
        
        % M
        MM=RJmJm*conj(Am')*(Am*RJmJm*conj(Am')+REmEm)';
        Jm=(MM*Y(m,(1+u):(K-u))');
        Um=fliplr(conj(Jm));
        
        % Circular Convolution
        Y(s,(1+u):(K-u))=(Am*Um)';
    end
    
    
    for s=1:NoSym
        %Y(s,:)=cconv(Y(s,:),fft(conj(PN(s,CP+1:end)),K),K);
    end
end


% Channel estimation
for s=1:NoSym    
    for i=start:gap:K
        channelEstimate(i)=pilot'*Y(s,i)./power(abs(pilot),2);
        channelEstimate((i-(gap/2)+1):(i+(gap/2)))=channelEstimate(i);
    end
    % Symbol estimation
    for i=1:K
        X_hat(s,i) = channelEstimate(i)'*Y(s,i)./power(abs(channelEstimate(i)),2);
    end

end

% Demodulation
RI_hat=qamdemod(X_hat,M,'UnitAveragePower',true,'OutputType','integer');

% Constellation diagram
%scatterplot(X_hat(:,64));
%axis([-1.2 1.2 -1.2 1.2]);

% SER measurements
SER = length(find(RI_hat~=XI_hat)) / (NoSym*K);
%disp(SER);
end
