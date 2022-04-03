%{
/*-----------------------------------------------------------------------------------------------------+
| DAQ_module01.m         |Acquired signal frequency estimator and RMS, apmain.c                        |
|                        |                                                                             |
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
est_freq=lambda*szspctr;

%{
--------------------------       
|RMS CALCULATOR SEQUENCE |      
--------------------------      
%}

for i=1:szvec
  tmp = aqvec(i) * aqvec (i);
  summation = summation+tmp;
end

rms = sqrt(summation / szvec);

%{
------------------------------   
|HARMONIC CALCULATOR SEQUENCE|  
------------------------------  
%}




