clear all
fs=128; % sampling frequency
NumSubj=1;  %number of subject
ch=32;  % number of channels
for p=1:NumSubj
    % loading the data
    add=['H:\Post-Doc\DATA\Emotion\all subjects\p1\v' num2str(p) '.mat'];
    EEG1=importdata(add);
    
    % decomposition into brain frequency bands 
    % NOTE: needs installation of EEGLAB
    [delta,filtwts] = eegfilt(EEG1,fs,1,4);  %delta (1-4 Hz)
    [theta,filtwts] = eegfilt(EEG1,fs,4,8);  %theta (4-8 Hz)   
    [alpha,filtwts] = eegfilt(EEG1,fs,8,16);  %alpha (8-16 Hz)
    [beta,filtwts] = eegfilt(EEG1,fs,16,32);  %beta (16-32 Hz)
    [gamma,filtwts] = eegfilt(EEG1,fs,32,45);  %gamma (32-45 Hz)

    % ploting time-frequency spectrum
    [insamp,insfrq]=inst_frq_amp(beta,fs);  % extraction of frequncy and amplitude of a narrow band signal
    ptitle='time-frequency spectrum';
    PlotSpectra(insfrq,insamp,ptitle);
    
    % Phase-phase synchronization matrix
    %for ch=1:32
    %sig1=[delta(ch,:);theta(1,:);alpha(1,:);beta(1,:)];
    insphase=inst_phase(beta);          % extraction of instant phases
    sync=phasesync(insphase);           % connectivity matrix
    
    % global efficiency
    GE=efficiency_wei(sync)
    
    % local efficiency
    LE=efficiency_wei(sync,1)
    figure;
    mx=max(LE);mn=min(LE);
    topoplot(LE,'eeglab_chan32.locs','electrodes','off','colormap','hot','maplimits',[mn mx]);
    colorbar
    title('Local efficiency')
    
    %DFA : scaling exponent
    figure;    
    for i=1:ch %calculation for each channel
        x=beta(i,:);
        [y1,y2]=DFA_main(x);
        D(i)=y1;
        A(i)=y2;
    end
    mx=max(A);mn=min(A);
    figure;
    topoplot(A,'eeglab_chan32.locs','electrodes','off','colormap','hot','maplimits',[mn mx]);
    colorbar
    title('Scaling Exponent')

   %Entropy
    for i=1:ch %calculation for each channel
        x=[delta(i,:);theta(i,:);alpha(i,:);beta(i,:);gamma(i,:)];
        En(i)=FFTentropy(x);
    end
    mx=max(En);mn=min(En);
    figure;
    topoplot(En,'eeglab_chan32.locs','electrodes','off','colormap','hot','maplimits',[mn mx]);
    colorbar
    title('Entropy')
    
    %Fractal Dimention higuchin & petrosian algorithm
    kmax=13;
    for i=1:ch %calculation for each channel
        x=beta(i,:);
        [d,LM]=higuchin(x,kmax);
        higDim(i)=d;
        
        d=petrosian(x);
        petDim(i)=d;
    end

end

    
    
    
    
    
    
    
    
    
    
    
    
    
% add=['H:\Post-Doc\DATA\Emotion\part1\part1\p01\v' num2str(vd)];
% save(add,'EEG1')
