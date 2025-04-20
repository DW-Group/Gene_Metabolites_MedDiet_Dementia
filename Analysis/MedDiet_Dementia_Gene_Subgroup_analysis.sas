%include "/scripts/data_nhs_adrd_full_baseline_r1.sas"; 

data nhs1datas; 
    set nhs1data end=_end_;  
    
    if apoe4grp ne .; 
    
    if apoe4grp in (1,2,3) then apoe_3cat=0;
    else if apoe4grp in (4,5) then apoe_3cat=1;
    else if apoe4grp in (6) then apoe_3cat=2;
    else apoe_3cat=.; 

    %indic3(vbl=apoe_3cat, prefix=apoe_3cat, min=1, max=2, reflev=0, missing=.,usemiss=0);
    run;

proc freq; tables apoe_3cat; run;

%let covapoe = &apoe_3cat_ &platforms_ PC1 PC2 PC3 PC4  ; 
%let covprs = PGScon &platforms_ PC1 PC2 PC3 PC4  ; 
%let covbase = &hiedu_  &qnSES_ &actcc_ &smkcat_   htn depression antidep  diabetes &bmic_  &husbedu_   fhdem   &phmsstatus_  marry live_alone; 

/*NHS1*/

%macro m2(dataset, stra_str, stra, adj1, adj2);

%mphreg9(data=nhs1datas, event=adcase, time=tad, 
         timevar=t80 t82 t84 t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18 t20,
         id=id, tvar=period, where=&stra, modopt=%quote(maxiter=30),
         agevar=agecon,
         qret= irt80 irt82 irt84 irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 irt18 irt20,
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
   cohorts="nhs1";
   keep variable cohorts strata modelno RR ProbChisq beta see;
   
   run; 

%mend;

%m2(dataset=apoe_3cat0nhs1, stra_str=apoe_3cat0, stra= %quote(apoe_3cat = 0), adj1= &covbase, adj2= &covbase &covprs);
%m2(dataset=apoe_3cat1nhs1, stra_str=apoe_3cat1, stra= %quote(apoe_3cat = 1), adj1= &covbase, adj2= &covbase &covprs);
%m2(dataset=apoe_3cat2nhs1, stra_str=apoe_3cat2, stra= %quote(apoe_3cat = 2), adj1= &covbase, adj2= &covbase &covprs);
%m2(dataset=qPGS002280_scaled0nhs1, stra_str=PGS0022800, stra= %quote(qPGS002280_scaled = 0), adj1= &covbase, adj2= &covbase &covapoe);
%m2(dataset=qPGS002280_scaled1nhs1, stra_str=PGS0022801, stra= %quote(qPGS002280_scaled = 1), adj1= &covbase, adj2= &covbase &covapoe);
%m2(dataset=qPGS002280_scaled2nhs1, stra_str=PGS0022802, stra= %quote(qPGS002280_scaled = 2), adj1= &covbase, adj2= &covbase &covapoe);
%m2(dataset=qPGS000334_scaled0nhs1, stra_str=PGS0003340, stra= %quote(qPGS000334_scaled = 0), adj1= &covbase, adj2= &covbase);
%m2(dataset=qPGS000334_scaled1nhs1, stra_str=PGS0003341, stra= %quote(qPGS000334_scaled = 1), adj1= &covbase, adj2= &covbase);
%m2(dataset=qPGS000334_scaled2nhs1, stra_str=PGS0003342, stra= %quote(qPGS000334_scaled = 2), adj1= &covbase, adj2= &covbase);

data subgroup_nhs1;
set
apoe_3cat0nhs1 apoe_3cat1nhs1 apoe_3cat2nhs1   
qPGS002280_scaled0nhs1 qPGS002280_scaled1nhs1 qPGS002280_scaled2nhs1
qPGS000334_scaled0nhs1 qPGS000334_scaled1nhs1 qPGS000334_scaled2nhs1
;
run;

proc sort data=subgroup_nhs1; by cohorts strata modelno; run;

PROC EXPORT DATA= subgroup_nhs1
            OUTFILE="/results/amed_APOE4_PRS_subgroup_nhsfull_baseline.csv"
            DBMS=csv REPLACE;
RUN;
