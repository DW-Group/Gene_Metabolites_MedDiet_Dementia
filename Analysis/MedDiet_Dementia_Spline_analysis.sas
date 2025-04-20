%include "/scripts/data_adrd_nhs_full_baseline.sas"; 

%let covbase = &hiedu_ &qnSES_ &actcc_ &smkcat_ htn depression antidep diabetes &bmic_ &husbedu_ fhdem &phmsstatus_ marry live_alone; 

%LGTPHCURV9(data=nhs1data, exposure=amedcon, case=adcase, model=cox, time=tad,
	strata=agecon period,
	adj= &covbase, 
	refval=0,
	NK=3,
	select=1,
    	outplot=JPEG,
    	pwhich=SPLINE,  
    	footer=NONE,
    	GRAPHTIT=NONE,
   	plot=4,
   	e=T,     
    	ci=2,    
   	displayx=F,  
   	klines=F,   
   	plotdec=F,
    plotdata=amed_dementia_spline_baseline.txt);
run;



 

























