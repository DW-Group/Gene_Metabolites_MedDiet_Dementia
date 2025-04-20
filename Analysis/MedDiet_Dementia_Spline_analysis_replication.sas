%include "/scripts/data_adrd_hpfs_full_baseline.sas"; 

%let covbase = &proff_ &qnSES_ &actcc_ &smkcat_ htn depression antidep diabetes &bmic_ fhdem marry live_alone; 

%LGTPHCURV9(data=hpfsdata, exposure=amedcon, case=adcase, model=cox, time=tad,
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
    plotdata=amed_dementia_spline_baseline_hpfs.txt);
run;



 

























