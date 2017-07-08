for f3dB=[50 100 150]
    SER=[];
    SER_PC=[];
    SNR=[];
    
    for snr=0:30
        ser=OFDM_PN(f3dB,snr,'False',9);
        serr=OFDM_PN(f3dB,snr,'True',20);
        SER=[SER ser];
        SER_PC=[SER_PC serr];
        SNR=[SNR snr];
    end
    
    figure
    semilogy(SNR,SER,'^-',SNR,SER_PC,'^-');
    xlabel('SNR');
    ylabel('SER');
    title(['SER vs. SNR with f3dB = ',num2str(f3dB)]);
    legend('no Phase Noise Compensation','with Phase Noise Compensation')
    %legend('without Compensation','ideal compensation')
    grid on;
end
