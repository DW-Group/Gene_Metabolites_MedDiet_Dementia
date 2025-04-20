options mautosource sasautos='/usr/local/channing/sasautos';

PROC IMPORT DATAFILE="/data/merged_data_baseline_nhs_2023_for_sas_all_imputed.csv"
    OUT=all_imputed
    DBMS=CSV 
    REPLACE;
    GETNAMES=YES; 
    GUESSINGROWS=MAX; 
RUN;

PROC IMPORT DATAFILE="/data/merged_data_baseline_nhs_2023_for_sas_APOE4_noncarrier_imputed.csv"
    OUT=APOE4_noncarrier_imputed
    DBMS=CSV 
    REPLACE;
    GETNAMES=YES; 
    GUESSINGROWS=MAX; 
RUN;

PROC IMPORT DATAFILE="/data/merged_data_baseline_nhs_2023_for_sas_APOE4_carrier_imputed.csv"
    OUT=APOE4_carrier_imputed
    DBMS=CSV 
    REPLACE;
    GETNAMES=YES; 
    GUESSINGROWS=MAX; 
RUN;

%let med_list = HMDB0000462 HMDB0000658 HMDB0001348 HMDB0003282 HMDB0011103 HMDB0011214 HMDB0029377;
%let covar_nongene_list = bmi hiedu3 husbedu2_3_m phmsstatus3_4 fhdem actcon nSES smkk2 smkk3_4_5 dep_antidep htn sysbp2 sysbp3 sysbpm hchol alcocon daykcal;

title2 'All_imputed';
%mediate(data=all_imputed, 
exposure=AMED_avg, 
intermed=&med_list,
id=id, 
time=tad, 
event=adcase, 
intmiss=F, 
surv=T,
covars=agemo &covar_nongene_list);

title2 'APOE4_noncarrier_imputed';
%mediate(data=APOE4_noncarrier_imputed, 
exposure=AMED_avg, 
intermed=&med_list,
id=id, 
time=tad, 
event=adcase, 
intmiss=F, 
surv=T,
covars=agemo &covar_nongene_list);

title2 'APOE4_carrier_imputed';
%mediate(data=APOE4_carrier_imputed, 
exposure=AMED_avg, 
intermed=&med_list,
id=id, 
time=tad, 
event=adcase, 
intmiss=F, 
surv=T,
covars=agemo &covar_nongene_list);