%include "/scripts/data_hpfs_adrd_full_baseline_r1.sas"; 

data hpfsdatas; 
    set hpfsdata end=_end_;  
    
    if apoe4grp ne .; 
    
    if apoe4grp in (1,2,3) then apoe_3cat=0;
    else if apoe4grp in (4,5) then apoe_3cat=1;
    else if apoe4grp in (6) then apoe_3cat=2;
    else apoe_3cat=.; 

    %indic3(vbl=apoe_3cat, prefix=apoe_3cat, min=1, max=2, reflev=0, missing=.,usemiss=0);
    run;

proc freq; tables apoe_3cat; run;

%let covapoe = &apoe_3cat_ &platforms_ PC1 PC2 PC3 PC4  ; 
%let covprs  = PGScon &platforms_ PC1 PC2 PC3 PC4  ; 
%let covbase = &proff_  &qnSES_ &actcc_ &smkcat_   htn depression antidep  diabetes &bmic_  fhdem  marry live_alone; 

/*HPFS*/

%macro m2(dataset, stra_str, stra, adj1, adj2);

%mphreg9(data=hpfsdatas, event=adcase, time=tad, 
         timevar=t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18 t20,
         id=id, tvar=period, where=&stra, modopt=%quote(maxiter=30),
         agevar=agecon,
         qret= irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 irt18 irt20,
         cutoff=cutoff,
         dtdx=dtdx_ad, dtdth=dtdth,labels=F,

         model1 = amedcon &adj1, 
         model2 = amedcon &adj2, 
         outdat=&dataset);

data &dataset;
set &dataset;
   
   RR=put(HazardRatio,4.2)|| ' (' ||put(LCL,4.2)|| ', ' ||put(UCL,4.2)|| ')';
   beta=put(Estimate,9.6);
   see =put(StdErr,9.6);

   strata="&stra_str.";
   cohorts="hpfs";
   keep variable cohorts strata modelno RR ProbChisq beta see;
   
   run; 

%mend;

%m2(dataset=apoe_3cat0hpfs, stra_str=apoe_3cat0, stra= %quote(apoe_3cat = 0), adj1= &covbase, adj2= &covbase &covprs);
%m2(dataset=apoe_3cat1hpfs, stra_str=apoe_3cat1, stra= %quote(apoe_3cat = 1), adj1= &covbase, adj2= &covbase &covprs);
%m2(dataset=apoe_3cat2hpfs, stra_str=apoe_3cat2, stra= %quote(apoe_3cat = 2), adj1= &covbase, adj2= &covbase &covprs);

data subgroup_hpfs;
set
apoe_3cat0hpfs apoe_3cat1hpfs apoe_3cat2hpfs
;
run;

proc sort data=subgroup_hpfs; by cohorts strata modelno; run;

PROC EXPORT DATA= subgroup_hpfs
            OUTFILE="/results/amed_APOE4_PRS_subgroup_hpfsfull_baseline.csv"
            DBMS=csv REPLACE;
RUN;
