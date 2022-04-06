%{
/*-----------------------------------------------------------------------------------------------------+
| DAQ_module01.m         |Single DAQ Acquisitor. Includes frequency estimator through ipDTF, Harmonic  |
|                        |finder, and signal RMS, Harmonic RMS, THD, Average value and Power Spectrum  |
|                        |calculator                                                                   |
+------------------------------------------------------------------------------------------------------+
| Authors: Joao Barreiros C. Rodrigues nº99968, Francisco Simplício nº99940, Inês Castro nº99962       |
|          LEEC-IST                                                                                    |
| Date: 03 April 2022                                                                                  |
+-----------------------------------------------------------------------------------------------------*/
%}

%{
------------------
|DAQ INIT ENGAGED|
------------------
%}
clear all;
sampling_freq=25000;
N_aq=1000;
daq_storage= daq("ni");
daq_storage.Rate=sampling_freq; %transfer sampling frequency to daq_storage object parameter Rate
ch=addinput(daq_storage,"Dev1","ai0","Voltage");
ch.Range=[-2,2];
aqvec=read(daq_storage, N_aq, "OutputFormat", "Matrix"); % store data read from daq as a vector
szvec=length(aqvec);

%{
---------------------
|DAQ INIT TERMINATED|
---------------------
%}


%{
% example
N_aq=30000;
sampling_freq=40000;
amplitude=1;
n=0:1:N_aq-1;
t=n/sampling_freq;
%%sn = 2*sawtooth(2*pi*54*t);
aqvec=amplitude*square(2*pi*600*t);
%%sn=sawtooth(2*pi*frequencia*t);
%}

%{
-------------------------------       
|FREQUENCY ESTIMATOR SEQUENCE |      
-------------------------------      
%}


delta_freq=sampling_freq/N_aq;


spctrvec=fft(aqvec);
szspctr=length(spctrvec);

peak = 0;
for i=2:szspctr
	if(peak < spctrvec(i)) %Find fft peak without using MATLAB bloat functions%
		peak = spctrvec(i);
		peakidx = i;
    end
end
neighboridx=0;
if(spctrvec(peakidx-1) > spctrvec (peakidx+1)) %Find suitable neighbor to use in ipDTF formula%
	neighboridx=peakidx-1;
	L = neighboridx;
else
	L = peakidx;
   	neighboridx=peakidx+1;
end

Upk = real(spctrvec(L)); %Apply ipDTF formula%
Une = real(spctrvec(L+1));
Vpk = imag(spctrvec(L));
Vne = imag(spctrvec(L+1));

n=2*pi/N_aq;

Kopt=((Vne-Vpk)*sin(n*peakidx)+(Une-Upk)*cos(n*peakidx))/(Une-Upk);
Zpk=Vpk*(Kopt-cos(n*peakidx)/sin(n*peakidx)) + Upk;
Zne=Vne*(Kopt-cos(n*neighboridx)/sin(n*neighboridx)) + Une;
lambda=(1/n)*acos((Zne*cos(n*neighboridx)-Zpk*cos(n*peakidx))/(Zne-Zpk));
est_freq=lambda*delta_freq;


%{
--------------------------       
|AVG CALCULATOR SEQUENCE |      
--------------------------      
%}

Nsampoper= sampling_freq / est_freq;
Nperiods= N_aq / Nsampoper;
Ncperiods= floor(Nperiods);
N_rounded = round (Ncperiods * Nsampoper);

summation = 0;

for i = 1:N_rounded
    summation = summation + aqvec(i);
end

avg = summation/N_rounded;


%{
--------------------------       
|RMS CALCULATOR SEQUENCE |      
--------------------------      
%}

summation = 0;

for i=1:(N_rounded-1)
  tmp = aqvec(i) * aqvec(i);
  summation = summation+tmp;
end

rms = sqrt(summation/N_rounded);
%{
------------------------------
|HARMONIC CALCULATOR SEQUENCE|
------------------------------
%}

N_harmonics=11; %see max number of harmonics%
harmonicfreqvec(1)=peakidx;
harmonicvoltvec(1)=abs(peak);
lastharmonicfreq=peakidx;
for i=2:N_harmonics
        nextharmonicfreq=round(lastharmonicfreq+peakidx); %see cast%
        harmonicfreqvec(i)=nextharmonicfreq;
        if(nextharmonicfreq>length(spctrvec))
                break;
        end
        harmonicvoltvec(i)=abs(spctrvec(nextharmonicfreq)/N_aq);
        lastharmonicfreq=nextharmonicfreq;
end

harmonicvecsize=length(harmonicvoltvec);
%{
--------------------------
|HRMS CALCULATOR SEQUENCE |
--------------------------
%}

for i=1:harmonicvecsize
  harmonicrms(i) = abs(harmonicvoltvec(i))/sqrt(2)
end

%{
--------------------------
|THD CALCULATOR SEQUENCE |
--------------------------
%}

for i=2:harmonicvecsize
  tmp = harmonicvoltvec(i) * harmonicvoltvec(i);
  summation = summation+tmp;
end
thd = sqrt(summation)/harmonicvoltvec(1);
thddB=20*log10(thd);

%{
--------------------------       
|GRAPH PLOTTING SEQUENCE |      
--------------------------      
%}
subplot(2,1,1);

n=0:1:N_aq-1;
x1=n/sampling_freq;

plot(x1, aqvec);
xlim([0 5/est_freq]); 
xlabel('Sampler \Delta time / s');
ylabel('Sampled signal amplitude / V');
title(sprintf('Sampled signal \n Estimated original signal frequency: %f Hz\n Average value: %f  RMS: %f\n Number of acquisitions: %d Sampler Frequency:%f\n Range:[-2,2]', est_freq, avg,rms,N_aq, sampling_freq));

subplot(2,1,2);

spectral_resol = sampling_freq/N_aq;

power_vec = 2*abs(spctrvec)/sqrt(2);
x2 = 0:spectral_resol:spectral_resol*(N_aq/2-1);
plot(x2,20*log10(power_vec(1:N_aq/2)/N_aq), '.-'); %log scale

xlabel('Frequency / Hz');
ylabel('Logarithmic power / dBV^2');
title(sprintf('Single-sided Logarithmic signal power'));



