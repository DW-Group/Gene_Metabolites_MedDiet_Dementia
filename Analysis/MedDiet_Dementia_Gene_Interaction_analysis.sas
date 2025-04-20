%include "/scripts/data_nhs_adrd_full_baseline.sas"; 

data nhs1datas; 
    set nhs1data end=_end_; 
    
    if apoe4grp ne .; 
    
    if apoe4grp in (1,2,3) then apoe_3cat=0;
    else if apoe4grp in (4,5) then apoe_3cat=1;
    else if apoe4grp in (6) then apoe_3cat=2;
    else apoe_3cat=.; 

    %indic3(vbl=apoe_3cat, prefix=apoe_3cat, min=1, max=2, reflev=0, missing=., usemiss=0);

    int_apoe3cat1amed = apoe_3cat1*amedcon;
    int_apoe3cat2amed = apoe_3cat2*amedcon;

    run;

proc freq; tables apoe_3cat; run;

%let covgene = &platforms_ PC1 PC2 PC3 PC4  ; 
%let covapoe = &apoe_3cat_ &platforms_ PC1 PC2 PC3 PC4  ; 
%let covprs = PGScon &platforms_ PC1 PC2 PC3 PC4  ; 
%let covapoeprs = &apoe_3cat_ PGScon &platforms_ PC1 PC2 PC3 PC4  ; 
%let covbase = &hiedu_  &qnSES_ &actcc_ &smkcat_   htn   depression antidep  diabetes &bmic_  &husbedu_   fhdem   &phmsstatus_  marry live_alone; 

/*Model*/

%mphreg9(data=nhs1datas, event=adcase, time=tad,  
         timevar=t80 t82 t84 t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18 t20,
         id=id, tvar=period, modopt=%quote(maxiter=30),
         agevar=agecon, 
         qret= irt80 irt82 irt84 irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 irt18 irt20,   
         cutoff=cutoff,
         dtdx=dtdxdiab, dtdth=dtdth,labels=F,

         model1    = amedcon ,
         model2    = amedcon &covbase,
         model3    = amedcon &covbase &covapoeprs,

         model4    = amedcon &apoe_3cat_ int_apoe3cat1amed int_apoe3cat2amed &covgene,
         model5    = amedcon &apoe_3cat_ int_apoe3cat1amed int_apoe3cat2amed &covbase &covgene,
         model6    = amedcon &apoe_3cat_ int_apoe3cat1amed int_apoe3cat2amed &covbase &covprs,

         outdat=int_nhs1);
         run;

data int_nhs1; set int_nhs1; cohorts='nhs1'; factors='main_interaction';run;

PROC EXPORT DATA= int_nhs1
            OUTFILE="/results/amed_APOE4_interaction_nhsfull_baseline.csv"
            DBMS=csv REPLACE;
RUN;



















