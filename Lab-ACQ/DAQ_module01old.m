%{
/*-----------------------------------------------------------------------------------------------------+
| DAQ_module01.m         |Acquired signal frequency estimator and Average value, RMS and signal calcu- |
|                        |lator                                                                        |
|                        |                                                                             |
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
-------------------------------       
|FREQUENCY ESTIMATOR SEQUENCE |      
-------------------------------      
%}
spctrvec=fft(aqvec);
szspctr=length(spctrvec);
peak=0;
for i=2:szspctr
	if(peak < spctrvec(i)) %Find fft peak without using MATLAB bloat functions%
		peak = spctrvec(i);
		peakidx = i;
	end

if(spctrvec(peakidx-1) > spctrvec (peakidx+1)) %Find suitaible neighbor to use in ipDTF formula%
	neighboridx=peakidx-1;
	neighbor=spctrvec(peakidx-1);
else
	neighboridx=peakidx+1;
	neighbor=spctrvec(peakidx+1);
end

Upk = real(peak); %Apply ipDTF formula%
Une = real(neighbor);
Vpk = imag(peak);
Vne = imag(neighbor);
n=2*pi/szvec;
Kopt=((Vne-Vpk)*sin(n*peakidx)+(Une-Upk)*cos(n*peakidx))/(Une-Upk);
Zpk=Vpk*(Kop-cos(n*peakidx)\sin(n*peakidx)) + Upk;
Zne=Vne*(Kop-cos(n*neighboridx)\sin(n*neighboridx)) + Une;
lambda=(1/n)*ascos((Zne*cos(n*neighboridx)-Zpk*cos(n*peakidx))\(Zne-Zpk));
est_freq=lambda*delta_freq;


%{
--------------------------       
|AVG CALCULATOR SEQUENCE |      
--------------------------      
%}

Nsampoper= sampling_freq / est_freq;
Nperiods= N_aq / Nsampoper;
Ncperiods= floor(Nperiods);
avg = round (Ncperiods * Nsampoper);


%{
--------------------------       
|RMS CALCULATOR SEQUENCE |      
--------------------------      
%}

for i=1:avg-1
  tmp = aqvec(i) * aqvec(i);
  summation = summation+tmp;
end
rms = sqrt(summation/avg);


%{
--------------------------       
|GRAPH PLOTTING SEQUENCE |      
--------------------------      
%}
subplot(2,1,N_aq);
x1 = linspace(0,(1/delta_freq));
xlabel('Sampler \delta time / s');
ylabel('Sampled signal amplitude / V');
title(sprintf('Sampled signal \n Estimated original signal frequency:', est_freq, 'Hz \n', 'Average value:', avg, ' RMS:', rms '\n Number of acquisitions:', N_aq, '\n Sampler Frequency:', sampling_freq, 'Range:');
plot(x1, aqvec);

for i=1:szspctr
	tmp= (spctrvec(i)*spctrvec(i))/2;
	powervec = 10 * log (tmp);
end

x2=linspace(0, szspctr)
xlabel('Frequency / Hz');
ylabel('Logarithmic power / dBV^2');
title(sprintf('Single-sided Logarithmic signal power');
plot(x2, powervec);



