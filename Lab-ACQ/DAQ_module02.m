%{
/*-----------------------------------------------------------------------------------------------------+
| DAQ_module02.m         |                                                                             |
|                        |                                                                             |
|                        |                                                                             |
+------------------------------------------------------------------------------------------------------+
| Authors: Joao Barreiros C. Rodrigues nº99968, Francisco Simplício nº99940, Inês Castro nº99962       |
|          LEEC-IST                                                                                    |
| Date: 03 April 2022                                                                                  |
+-----------------------------------------------------------------------------------------------------*/
%}

N_it=200;
for i=1:N_it

	%{
	------------------
	|DAQ INIT ENGAGED|
	------------------
	%}
	sampling_freq=;
	N_aq=;
	delta_freq=sampling_freq/N_aq;
	daq_storage= daq("ni");
	daq_storage.Rate=sampling_freq; %transfer sampling frequency to daq_storage object parameter Rate
	addinput(daq_storage,"Dev1","ai0","Voltage"); 
	aqvec=read(daq_storage, N_aq, "OutputFormat", "Matrix"); % store data read from daq as a vector
	szvec=length(aqvec);

	%{
	---------------------
	|DAQ INIT TERMINATED|
	---------------------
	%}

	%{
	-----------------------------------
	|SIGNAL POWER CALCULATOR SEQUENCE |
	-----------------------------------

	%}
	spctrvec=fft(aqvec);
	tmp= 2*abs(spctrvec)/sqrt(2);
	spower_vec = tmp + spower_vec;
end

spower_vec=spower_vec/N_it


%{
------------------------------
|HARMONIC CALCULATOR SEQUENCE|
------------------------------
%}

N_harmonics=; %see max number of harmonics%
harmonicfreqvec(1)=peakidx;
harmonicvoltvec(1)=peak;
lastharmonicfreq=peakidx;
for i=2:N_harmonics
        nextharmonicfreq=lastharmonicfreq+est_freq; %see cast%
        harmonicfreqvec(i)=nextharmonicfreq;
        harmonicvoltvec(i)=spctrvec(nextharmonicfreq);
        lastharmonicfreq=nextharmonicfreq;
end

harmonicvecsize=lenght(harmonicvec);
%{
--------------------------
|HRMS CALCULATOR SEQUENCE |
--------------------------
%}

for i=2:harmonicvecsize
  harmonicrms = abs(harmonicvoltvec(i))/sqrt(2)
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


%{
----------------------------
|NOISE CALCULATOR SEQUENCE |
----------------------------
%}
power_vecsize = length(spower_vec);
tmp_power=20*log10(spower_vec);

for i=1:spower_vecsize/2
	summation=10^(tmp_power(i)/10)+summation;
		
end
noiserms=sqrt(summation);

%{
----------------------------
|SINAD CALCULATOR SEQUENCE |
----------------------------
%}
for i=2:spower_vecsize
        if(powerpeak < spower_vec(i)) %Find fft peak without using MATLAB bloat functions%
                powerpeakidx = i;
    end
end
for i=1:spower_vecsize
	if(i ~= powerpeak)
		summation = summation+ spower_vec(i);
	end
end

avgspower=summation/(spower_vecsize-1)
sinad = 20*log10(spower_vec(powerpeak))-20*log10(avgspower);

%{
---------------------------
|ENOB CALCULATOR SEQUENCE |
---------------------------
%}
enob= (sinad-1.76)/6.02;

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
xlabel('Sampler \delta time / s');
ylabel('Sampled signal amplitude / V');
title(sprintf('Last Sampled signal \n \n  RMS: %f Noise RMS: %f\n SINAD: %f ENOB: %f\n  Number of acquisitions: %d\n Sampler Frequency:%f\n Range:',rms,noiserms,sinad,enob, N_aq, sampling_freq);

subplot(2,1,2);

spectral_resol = sampling_freq/N_aq;

x2 = 0:spectral_resol:spectral_resol*(N_aq/2-1);
plot(x2,20*log10(spower_vec(1:N_aq/2)/N_aq), '.-'); %log scale

xlabel('Frequency / Hz');
ylabel('Logarithmic power / dBV^2');
title(sprintf('Single-sided Logarithmic signal power'));



