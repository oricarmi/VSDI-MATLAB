function analyzeTS(Z,fs)
   if nargin<2
       fs = [];
   end
   if isempty(fs)
       fs = 100; % default fs
   end
   t = linspace(0,(length(Z)-1)/fs,length(Z));
   [pxx,f] = periodogram(Z,[],[],fs); 
   figure; 
   subplot 311; plot(t,Z); xlabel('time [sec]'); ylabel('amplitude'); title('Z vs time');
   subplot 312;[r,lags] = xcorr(Z,Z); plot(lags,r); xlabel('lag [samples]'); ylabel('amplitude'); title('Z sample autocorrelation');
   subplot 313;plot(f,pxx); xlabel('freq [Hz]'); ylabel('amplitude'); title('Z spectrum');
end

