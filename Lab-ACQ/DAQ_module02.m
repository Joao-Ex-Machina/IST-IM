%{
/*-----------------------------------------------------------------------------------------------------+
| DAQ_module02.m         |Acquired signal frequency estimator and Harmonics, RMS, THD, noise, SINAD and|    
|                        |ENOB calculator                                                              |
|                        |                                                                             |
+------------------------------------------------------------------------------------------------------+
| Authors: Joao Barreiros C. Rodrigues nº99968, Francisco Simplício nº99940, Inês Castro nº99962       |
|          LEEC-IST                                                                                    |
| Date: 05 April 2022                                                                                  |
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
szspctr=lenght(spctrvec);

for i=1:szspctr
	if(peak < spctrvec(i)) %Find fft peak without using MATLAB bloat functions%
		peak = spctrvec(i);
		peakidx = i;
	end

if(spctrvec(peakidx-1) > spctrvec (peakidx+1)) %Find sutaible neighbor to use in ipDTF formula%
	neighboridx=peakidx-1;
	neighbor=spctrvec(peakidx-1);
else
	neighboridx=peakidx+1;
	neighbor=spctrvec(peakfreq+1);
end

Upk = real(peak); %Apply ipDTF formula%
Une = real(neighbor);
Vpk = imag(peak);
Vne = imag(neighbor);
n=2*pi/szvec;
Kopt=((Vne-Vpk)*sin(n*peakidx)+(Une-Upk)*cos(n*peakidx))/(Une-Upk);
Zpk=Vpk*(Kop-cos(n*peakidx)\sin(n*peakidx)) + Upk;
Zne=Vne*(Kop-cos(n*neighboridx)\sin(n*neighboridx)) + Une;
lambda=(1/n)*arccos((Zne*cos(n*neighboridx)-Zpk*cos(n*peakidx))\(Zne-Zpk));
est_freq=lambda*delta_freq;

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
|RMS CALCULATOR SEQUENCE |      
--------------------------      
%}

for i=2:harmonicvecsize
  tmp = harmonicvoltvec(i) * harmonicvoltvec(i);
  tmp=tmp/sqrt(2);
  summation = summation+tmp;
end
tmp2=harmonicvoltvec(1)*harmonicvoltvec(1);
rms = sqrt(tmp2+summation);

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
--------------------------       
|GRAPH PLOTTING SEQUENCE |      
--------------------------      
%}
subplot(2,1,1);
x1 = linspace(0,(1/delta_freq));
xlabel('Sampler \delta time / s');
ylabel('Sampled signal amplitude / V');
title(sprintf('Sampled signal \n q, 'Hz \n', 'THD:', thd, ' RMS:', rms '\n Number of acquisitions:', N_aq, '\n Sampler Frequency:', sampling_freq, 'Range:');
plot(x, aqvec);

for i=1:szspctr
	tmp= (spctrvec(i)*spctrvec(i))/2;
	powervec = 10 * log (tmp);
	end

x2=linspace(0, szspctr)
xlabel('Frequency / Hz');
ylabel('Logarithmic power / dBV^2');
title(sprintf('Single-sided Logarithmic signal power');
plot(x, aqvec);



