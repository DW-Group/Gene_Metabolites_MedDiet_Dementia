libname hppdcase '/proj/hppars/hppar0g/case/';
filename nhstools '/proj/nhsass/nhsas00/nhstools/sasautos'; 
filename channing '/usr/local/channing/sasautos/'; 
libname readfmt '/proj/nhsass/nhsas00/formats/'; 
options mautosource sasautos=(channing nhstools); 
options fmtsearch=(readfmt);

************************************************************************************************************************************************
***********************************************************(1)   READ IN OUTCOME DATA   *******************************************************
************************************************************************************************************************************************;
 
/****** READ IN TICS variables and covariates already derived (create a dataset called "all")******/ 

******select the cognitive population (currently, includes all 23565 who were considered eligible,
      but many have never completed a cognitive interview)
      ------ to select just those who have completed a baseline cognitive assessment;

%include'/proj/nhalzs/nhalz00/CODE/nhscogfn_covariates01.sas'; 
 %include '/udd/hpypl/review/macro/avzgen6.sas';
 %include '/udd/hpypl/review/macro/fuavzgen6.sas';
 %include '/udd/hpypl/review/macro/avzgen4.sas';
 %include '/udd/hpypl/review/macro/fuavzgen4.sas';

/*** global, from 6 tests (using MEAN and SD from baseline scores ) ***/
 %avzgen6(all, score, animal, story1, story2, words2, digbac, ntotal); /*** First cognitive assessment ***/
 %fuavzgen6(all, score, animal, story1, story2, words2, digbac, scoreb, animalb, story1b, story2b, words2b, digbacb, ntotalb); /*** Second ***/ 
 %fuavzgen6(all, score, animal, story1, story2, words2, digbac, scorec, animalc, story1c, story2c, words2c, digbacc, ntotalc); /*** Third ***/ 
 %fuavzgen6(all, score, animal, story1, story2, words2, digbac, scored, animald, story1d, story2d, words2d, digbacd, ntotald); /*** Fourth***/ 


/*** verbal, from 4 tests (using MEAN and SD from baseline scores ) ***/
 %avzgen4(all, words, words2, story1, story2, verbn); /*** First cognitive assessment ***/ 
 %fuavzgen4(all, words, words2, story1, story2, wordsb, words2b, story1b, story2b, verbnb); /*** Second ***/  
 %fuavzgen4(all, words, words2, story1, story2, wordsc, words2c, story1c, story2c, verbnc); /*** Third ***/ 
 %fuavzgen4(all, words, words2, story1, story2, wordsd, words2d, story1d, story2d, verbnd); /*** Fourth ***/ 

/*** count those who do not have at least one complete battery ***/
data all;
set all;
if ntotal=. AND ntotalb=. AND ntotalc=. AND ntotald=. then cogm=1;
else cogm=0;
proc freq;
tables cogm;
where score ^= .;
proc sort; 
by id; 
proc freq;
tables intdt;
run;

/*proc freq data=all;
tables phy    anti    hbp   hchol   mi    db   alz   dprs  eft  asp  ibu   alc   alco  smk   nhor   bmi  vite  vitc   mvit;
run; 
proc means n nmiss min mean std max;run;*/


************************************************************************************************************************************************
*********************************************************** (2) READ IN EXPOSURE DATA   ********************************************************
************************************************************************************************************************************************;

/* Read in R */

/* tics baseline irt ranged from 1142-1219, 1995-2000, so prepare 4 expore level to match the tics baseline */

%include "scripts/food_nhs.sas"; 

data foodsnhs;
    set foodsnhs;

*cumulative average;
trmeat9094 = mean (of trmeat90d trmeat94d);
prmeat9094 = mean (of prmeat90d prmeat94d);
urmeat9094 = mean (of urmeat90d urmeat94d);
poultr9094 = mean (of poultr90d poultr94d);
fishal9094 = mean (of fishal90d fishal94d);
regegg9094 = mean (of regegg90d regegg94d);
hdairy9094 = mean (of hdairy90d hdairy94d);
ldairy9094 = mean (of ldairy90d ldairy94d);
tdairy9094 = mean (of tdairy90d tdairy94d);
soypro9094 = mean (of soypro90d soypro94d);
nutsal9094 = mean (of nutsal90d nutsal94d);
legume9094 = mean (of legume90d legume94d);
whgrns9094 = mean (of whgrns90d whgrns94d);
fruits9094 = mean (of fruits90d fruits94d);
vegeal9094 = mean (of vegeal90d vegeal94d);
rfgrns9094 = mean (of rfgrns90d rfgrns94d);
alcoho9094 = mean (of alco90n   alco94n);
calori9094 = mean (of calor90n  calor94n);
legsoy9094 = sum(0,soypro9094,legume9094);
nutleg9094 = mean (of nutleg90d nutleg94d);
decoff9094 = mean (of dcaf_s90  dcaf_s94);
coffee9094 = mean (of coff_s90  coff_s94);
sodium9094 = mean (of na90a     na94a);
 
run;

************************************************************************************************************************************************
*********************************************************** (3) READ IN OTHER COVARIATES *******************************************************
************************************************************************************************************************************************;

%ahei2010_8410(); proc contents; run; 

data ahei2010_8494;
set ahei2010_8410;
ahei_nomt_84=ahei2010_84-ahei2010_rmtI84;
ahei_nomt_86=ahei2010_86-ahei2010_rmtI86;
ahei_nomt_90=ahei2010_90-ahei2010_rmtI90;
ahei_nomt_94=ahei2010_94-ahei2010_rmtI94;
ahei_nomt_8494= mean (of ahei_nomt_84 ahei_nomt_86 ahei_nomt_90 ahei_nomt_94);

proc means n nmiss min mean std max nolabels;
var ahei2010_84     ahei2010_86     ahei2010_90     ahei2010_94
    ahei2010_rmtI84 ahei2010_rmtI86 ahei2010_rmtI90 ahei2010_rmtI94
    ahei_nomt_84    ahei_nomt_86    ahei_nomt_90    ahei_nomt_94   
    ahei_nomt_8494;
run;

    %include '/udd/hpypl/review/dietscore/ahei80_mc.sas'; 
    data ahei80; set diet80l; ahei_nomt_80=nAHEI80a-meataI80; keep id nAHEI80a ahei_nomt_80 ssb80; proc sort; by id; proc means n nmiss min mean std max nolabels; var ahei_nomt_80 ssb80; run;

    %n80_nts(keep=calor80n alco80n );   
    %n84_nts(keep=calor84n alco84n );   
    %n86_nts(keep=calor86n alco86n );  
    %n90_nts(keep=calor90n alco90n );   
    %n94_nts(keep=calor94n alco94n );    
 
    %n84_ant(keep=aofib84a);   
    %n86_ant(keep=aofib86a);  
    %n90_ant(keep=aofib90a);   
    %n94_ant(keep=aofib94a);  

%nses7614 ; run;

data nses7614;
   set nses7614;
   keep id division: /* Census division */
    mdinc: /* Median family income */
    mdvhs: /* Median family home value */
    pcolled: /* per over 25 with college or more */
    pfaminterest: /* per families receiving interest dividends or rent income */
    pohse: /* per occupied housing units */
    pwht: /* per population that is white */
    ppov: /* per population living in poverty */
    region: /* Census region */
    pfwchbf: /* per families with children headed by single female */
    phs: /* per over 25 with only high school */
    pnohs: /* per over 25 with no high school */
    pkid: /* per population under 5 */
    pold: /* per population over 65 */
    popd: /* population density num per sq km */
    ;
run;


%nur92 (keep= marry92 alone92 );  

%nur96 (keep= marry96 alone96 );
 
%nur00 (keep= marry00 alone00 ); 

libname hppdcase '/proj/hppars/hppar0g/case/';
/**** PD case in NHS between 1976-2012 ****/
data nhspdcase; set hppdcase.pd12nur; dtdx_pd=dtdx; drop dtdx status;
run; 
proc sort nodupkey; by id; run;

**READ IN APOE DATA FROM COGNITIVE STUDY OF BUCCAL CELL DNA**;
data apoe4dat;
 infile '/proj/nhblds/nhbld00/endpoints/cogfx/LABCODES/lab312/cog.98.312.results.APOE_3937' lrecl=11 recfm=d;
input
@1              id                  6.
@9              apoe4               2.
;

label
              id              = 'ID NUMBER'
              apoe4           = 'APOE4 0=wild,1=het,2=var,99=missing'
;
proc sort; by id;

data apoe2dat;
infile '/proj/nhblds/nhbld00/endpoints/cogfx/LABCODES/lab312/cog.98.312.results.APOE_4075' lrecl=11 recfm=d;
input
@1              id                  6.
@9              apoe2               2.
;

label
                id              = 'ID NUMBER'
                apoe2           = 'APOE2 0=wild,1=het,2=var,99=missing'
;

proc sort; by id;


data apoecell; merge apoe4dat apoe2dat; by id; run;

data apoecell; set apoecell;
if apoe4=0 and apoe2=0 then apoegenotype='e3e3';
else if apoe4=0 and apoe2=1 then apoegenotype='e2e3';
else if apoe4=0 and apoe2=2 then apoegenotype='e2e2';
else if apoe4=1 and apoe2=0 then apoegenotype='e3e4';
else if apoe4=1 and apoe2=1 then apoegenotype='e2e4';
else if apoe4=2 and apoe2=0 then apoegenotype='e4e4';

if id=258586 then apoegenotype=.;
 run;

**READ IN APOE DATA FROM RENAL STUDY**;
data gen;
infile '/proj/nhalzs/nhalz00/CODE/APOEDATA/apoe_from_renal_study.txt';
input
@1  id      6.
@8  apoety $4.
;
run;

proc sort; by id;

data gen; set gen; by id; if first.id and apoety not in (' ');run;


data gen; set gen;
if id=258586  then apoety=.;
run;


/* APOE imputed GWAS data genotype, file made by Hongyan Huang */
proc import datafile='/udd/sthoh/share/APOE/APOE_genotype_17Jul2017.txt' out=apoe_meas_imp2017 dbms=tab replace;
  getnames=yes;
run;


data apoe_nhs;
  set apoe_meas_imp2017;
  if study='HPFS' then delete;
  if study='NHS2' then delete;
  if study = 'PHS' then delete;
run; proc sort nodupkey; by id; run;


data apoe_nhs2017;
merge apoe_nhs(in=a) apoecell(in=b) gen(in=c); by id;
ingwas=a; incell=b; inrenal=c; run;

    %der7620 (keep=	can76 can78 can80 can82 can84 can86 can88 can90 can92 can94  can7694); 
              can7694=0; if  can76=1 or can78=1 or can80=1 or can82=1 or can84=1 or can86=1 or can88=1 or can90=1 or can92=1 or can94=1 then can7694=1;
              Run;

   
data alldat;
      merge  all (in=a)  foodsnhs
             ahei2010_8494 ahei80 n80_nts n84_nts n86_nts n90_nts n94_nts n84_ant n86_ant n90_ant n94_ant nses7614 nur92 nur96 nur00
             nhspdcase apoe_nhs2017 der7620
   ; 
   by id;
     exrec=1;
     if first.id and a then exrec=0;

aofib9094=mean (of aofib90a aofib94a);
calor9094=mean (of calor90n calor94n);
ahei_nomt_9094= mean (of ahei_nomt_90    ahei_nomt_94);
ssb9094       = mean (of ahei2010_ssb90  ahei2010_ssb94);

/**************************************************************
         Individual SES  -  marital status & living alone
**************************************************************/  
live_alone=alone92;
if alone92 ne 1 then live_alone=0;

marry=marry92;
if marry ne 1 then marry=0; 


* Neighborhood SES data ;
 
   zperpov=1-ppov10_94; *invert to make increases good;
   zsinglefemale=1-pfwchbf10_94; *invert to make increases good;
   zincome=mdinc10_94;
   zhomeval=mdvhs10_94;
   zpcollege=pcolled10_94;
   zphigh=phs10_94;
   zpnohigh=pnohs10_94;
   zpfamint=pfaminterest10_94;
   zpkids=pkid10_94;
   zpolds=pold10_94;
   zpocchse=pohse10_94;
   zpwhite=pwht10_94;

run;

proc standard data=alldat mean=0 std=1 out=alldat2;
     var zincome zhomeval zpcollege zpfamint zpocchse zpwhite zperpov zsinglefemale;
run; 

data n1data;
   set alldat2 end=_end_;
   nSES=zincome+zhomeval+zpcollege+zpfamint+zpocchse+zpwhite+zperpov;  
   

%beginex();

*******************************************************************************************
*                                      (4)  Baseline exclusion                            *
*******************************************************************************************;

*Get counts of starting N*;
	%exclude(id gt 0 or id le 0, nodelete=t); 
*Missing ids*;
	%exclude(id le 0);         
*Replicate ids or not in master file (cognition data)*;
	%exclude(exrec ne 0);     
*Missing first cognitive interview*;
	%exclude(score eq .);
*Incomplete battery - SELECT POPULATION WITH AT LEAST ONE COMPLETE BATTERY DURING FOLLOW-UP (including baseline)*;
	%exclude(cogm eq 1);
*Get counts of those left*;
	%exclude(id gt 0, nodelete=t); 
*Missing age at baeline *;
        %exclude(intage eq .); 
*Less than 70 (in master file)*;
	%exclude(intage lt 70);  

%output();

%endex();                    /* MUST do AFTER all exclusions */
run;


*******************************************************************************************
*                                   (5) baseline covariates                               *
*******************************************************************************************;

proc means n nmiss min mean std max median data=n1data;
var ahei_nomt_9094 aofib9094 calor9094 phy nSES trmeat9094 prmeat9094 urmeat9094;
run;
 


data nhsdata;
set n1data;
if ahei_nomt_8494=. then ahei_nomt_8494=48; *nmiss=281, median=48;
if aofib8494=. then aofib8494=18.55;        *nmiss=281, median=18.55;
if phy=. then phy=9.3;                      *nmiss=1148, median=9.3;
if nSES=. then nSES=-0.642;                 *nmiss=204, median=-0.642;
run;

     
%pctl9(data=nhsdata,
    varlist=  trmeat9094 ahei_nomt_9094 calor9094 ssb9094 fruits9094 vegeal9094  whgrns9094  poultr9094  fishal9094 regegg9094   nutsal9094 nutleg9094  tdairy9094 ldairy9094 hdairy9094 phy nSES decoff9094 coffee9094  sodium9094,
    numquant=3,
    mscore=T,
    quantname=q,
    cutdat=cutpoints,
    outdat=nhsdata,
    indic=T);
 

data nhs1data;
set nhsdata end=_end_;

agecon=intage;

/*  total red meat   */
%indic3(vbl=qtrmeat9094, prefix=qtrmeat9094, min=1, max=2, reflev=0, missing=., usemiss=0); 

*collapse MA and DR*;
if hiedu in (3,4) then hiedu=3;
%indic3(vbl=hiedu, prefix=hiedu, min=2, max=3, reflev=1, usemiss=0,
        label1='RN', 
        label2='BA', 
        label3='MA or DR');

**HUSBAND'S EDUCATION*;

/*** Recategorize husband's education ***/
 if husbe in (1,2,3) then hused=1;
 if husbe eq 4       then hused=2;
 if husbe eq 5       then hused=3;
 if husbe eq 6       then hused=4;   /* passthrough */

%indic3(vbl=hused, prefix=hused, min=2, max=3, reflev=1, missing=4,
      label1='=<H.S. grad',  
      label2='College grad', 
      label3='Grad school');

/*** family history of dementia ***/
if mdem=0 and fdem=0 and sdem=0 then famdem=0; /*no dem hx*/
else if mdem=1 or fdem=1 or sdem=1 then famdem=1;/*dem <65 y*/
else if mdem=2 or fdem=2 or sdem=2 then famdem=2;/*dem 65+ y*/
else famdem=9;/*missing*/
%indic3(vbl=famdem, prefix=famdem, min=1, max=2, reflev=0, missing=9, usemiss=1, 
        label0='family no dem hx',
        label1='family dem <65 y',
        label2='family dem 65+ y');

  %indic3(vbl=qahei_nomt_9094, prefix=qahei,           min=1, max=2, reflev=0, missing=., usemiss=0);  
  %indic3(vbl=qphy,            prefix=qact,            min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qnSES,           prefix=qnSES,           min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qcalor9094,      prefix=qdaykcal,        min=1, max=2, reflev=0, missing=., usemiss=0);             

  %indic3(vbl=qfruits9094, prefix=qfruits, min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qvegeal9094, prefix=qvegets, min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qwhgrns9094, prefix=qwhgrns, min=1, max=2, reflev=0, missing=., usemiss=1); 
  %indic3(vbl=qpoultr9094, prefix=qpoultr, min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qfishal9094, prefix=qfishs,  min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qregegg9094, prefix=qeggs,   min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qnutsal9094 ,prefix=qnuts ,  min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qtdairy9094, prefix=qtdairy, min=1, max=2, reflev=0, missing=., usemiss=0);  
  %indic3(vbl=qldairy9094, prefix=qldairy, min=1, max=2, reflev=0, missing=., usemiss=0);
  %indic3(vbl=qhdairy9094, prefix=qhdairy, min=1, max=2, reflev=0, missing=., usemiss=0);
  %indic3(vbl=qssb9094,    prefix=qssb,    min=1, max=2, reflev=0, missing=., usemiss=0); 
  %indic3(vbl=qnutleg9094 ,prefix=qnutleg ,  min=1, max=2, reflev=0, missing=., usemiss=0);

%indic3(vbl=dprs, prefix=dprs, min=2, max=2, reflev=1, missing=3,usemiss=1,
        label1='high mental health score',
        label2='low mental health score');
 
%indic3(vbl=alc, prefix=alc,   min=2, max=3, reflev=1, missing=4,usemiss=1,
        label1='alc non-drinkers',
        label2='alc 1-14g',
        label3='alc 15+ g');
%indic3(vbl=smk, prefix=smk,   min=2, max=3, reflev=1, usemiss=0,
        label1='smk never',
        label2='smk past',
        label3='smk current');

if nhor=3 then nhor=2; * to be consistent with other 2 analysis;
%indic3(vbl=nhor, prefix=nhor, min=2, max=2, reflev=1, missing=4,usemiss=1,
        label1='pmh never',
        label2='pmh past or current');
%indic3(vbl=bmi, prefix=bmi,   min=2, max=4, reflev=1, missing=5,usemiss=1,
        label1='bmi <22',
        label2='bmi 22-24',
        label3='bmi 25-29',
        label4='bmi 30+')

/*** alcohol ***/
       if alcoho9094=0.0                          then   alcc9094=1;
         else if alcoho9094>0.0   and alcoho9094<5.0   then   alcc9094=2;
           else if alcoho9094>=5.0  and alcoho9094<10.0  then   alcc9094=3;
             else if alcoho9094>=10.0 and alcoho9094<15.0  then   alcc9094=4;
               else if alcoho9094>=15.0 and alcoho9094<30.0  then   alcc9094=5;
                else if alcoho9094>=30.0                       then   alcc9094=6;   
                 else                                                   alcc9094=.;   
 
       %indic3(vbl=alcc9094, prefix=alcc9094, min=1, max=6, reflev=3, missing=., usemiss=1,
                label1='0 g/d',
                label2='0.1-4.9 g/d',
                label3='5.0-14.9 g/d',
                label4='5.0-14.9 g/d',
                label5='15.0-29.9 g/d',
      	        label6='30+ g/d');
 
/* Apoe_2cat */ 
/* from GWAS */
   if APOE_genotype = 'e2e2' then apoe_2cat = 0;
else if APOE_genotype = 'e2e3' then apoe_2cat=0;
else if APOE_genotype = 'e3e2' then apoe_2cat=0;
else if APOE_genotype = 'e3e3' then apoe_2cat=0;
else if APOE_genotype = 'e2e4' then apoe_2cat=1;
else if APOE_genotype = 'e4e2' then apoe_2cat=1;
else if APOE_genotype = 'e4e3' then apoe_2cat=1;
else if APOE_genotype = 'e3e4' then apoe_2cat=1;
else if APOE_genotype = 'e4e4' then apoe_2cat=1;
else apoe_2cat=.;
run;

 
*******************************************************************************************************
*               (6)  - create dataset "nhscogfn_rmeat" for cross-sectional analysis                   *
*******************************************************************************************************;

/*** MAKE Z-SCORES OF EACH TEST ((using MEAN and SD from baseline scores which are computed in the previous macros ***/

/*** 1st Z scores ***/

 %include '/udd/hpypl/review/macro/zscore.sas';
 %zscore(nhs1data, score,  zscre); 
 %zscore(nhs1data, animal, zaniml); 
 %zscore(nhs1data, words,  zwrds); 
 %zscore(nhs1data, words2, zwrds2); 
 %zscore(nhs1data, story1, zstry1); 
 %zscore(nhs1data, story2, zstry2); 
 %zscore(nhs1data, digbac, zdigbc);

proc means mean std data=nhs1data;
var score animal words words2 story1 story2 digbac;
run;
 

/*** 2ND ZSCORES USING 1ST MEAN AND SD ***/

data nhscogfn_rmeat;
set nhs1data end=_end_;

zscreb=((scoreb-33.72)/2.76); 
zanimlb=((animalb-16.87)/4.65);
zwrdsb=((wordsb-4.61)/1.70); 
zwrds2b=((words2b-2.31)/2.01);
zstry1b=((story1b-9.40)/1.73);
zstry2b=((story2b-8.97)/2.05);
zdigbcb=((digbacb-6.71)/2.43);


/*** CREATE 3rd ZSCORES USING 1ST MEAN AND SD ***/

zscrec=((scorec-33.72)/2.76); 
zanimlc=((animalc-16.87)/4.65);
zwrdsc=((wordsc-4.61)/1.70); 
zwrds2c=((words2c-2.31)/2.01);
zstry1c=((story1c-9.40)/1.73);
zstry2c=((story2c-8.97)/2.05);
zdigbcc=((digbacc-6.71)/2.43);


/*** CREATE 4th ZSCORES USING 1ST MEAN AND SD ***/

zscred=((scored-33.72)/2.76); 
zanimld=((animald-16.87)/4.65);
zwrdsd=((wordsd-4.61)/1.70); 
zwrds2d=((words2d-2.31)/2.01);
zstry1d=((story1d-9.40)/1.73);
zstry2d=((story2d-8.97)/2.05);
zdigbcd=((digbacd-6.71)/2.43);

 
******* AVERAGE OF THE 4 REPEATED COGNITIVE Z-SCORES (USING VALUES THAT ARE AVAILABLE)        ********;
 
av_ntotal=mean(ntotal,ntotalb,ntotalc, ntotald); 
av_nverbl=mean(verbn,verbnb,verbnc,verbnd);
av_zscre=mean(zscre,zscreb,zscrec,zscred);  
av_TICS =mean(score,score,scorec,scored);			/*** TICS CONTINUOUS ***/  
av_zaniml=mean(zaniml,zanimlb,zanimlc,zanimld);
av_zstry1=mean(zstry1,zstry1b,zstry1c,zstry1d);
av_zstry2=mean(zstry2,zstry2b,zstry2c,zstry2d); 
av_zwrds=mean(zwrds,zwrdsb,zwrdsc,zwrdsd);
av_zwrds2=mean(zwrds2,zwrds2b,zwrds2c,zwrds2d); 
av_zdigbc=mean(zdigbc,zdigbcb,zdigbcc,zdigbcd); 


********  - COUNTS AND FLAGS RELATED TO THE COGNITIVE OUTCOMES   ************;
  
/*** Flag cognitive files ***/  

if ntotal=. then ntota=0; else ntota=1;
if ntotalb=. then ntotb=0; else ntotb=1;
if ntotalc=. then ntotc=0; else ntotc=1;
if ntotald=. then ntotd=0; else ntotd=1;

if verbn=. then nverba=0; else nverba=1;
if verbnb=. then nverbb=0; else nverbb=1;
if verbnc=. then nverbc=0; else nverbc=1;
if verbnd=. then nverbd=0; else nverbd=1;

if zscre=. then nscrea=0; else nscrea=1;
if zscreb=. then nscreb=0; else nscreb=1;
if zscrec=. then nscrec=0; else nscrec=1;
if zscred=. then nscred=0; else nscred=1;

/*** patterns of cognitive follow-up - to describe data in the text of the ms ***/

        ************************************** 
          1: followed up the 4th visit 
          2: followed up the third visit
          3: up to second visit baseline only 
        **************************************;

          if zscred ne .      then pattern_zscre=1;
         else if zscrec ne .  then pattern_zscre=2;
         else if zscreb ne .  then pattern_zscre=3;
         else                      pattern_zscre=4;           /* subjects have all TICS at baseline */
       
         if ntotald ne .      then pattern_total=1;
         else if ntotalc ne . then pattern_total=2;
         else if ntotalb ne . then pattern_total=3;
         else if ntotal ne .  then pattern_total=4;
         else                      pattern_total=5;          /* not all the subjects have global score at baseline - but should all have at least one
                                                                during follow-up so pattern_total 5 should have n=0 */
         

/**** count the number of repeated measures per subject, for each cog outcome ***/

ntot = ntota + ntotb + ntotc + ntotd;
nverb = nverba + nverbb + nverbc + nverbd;
nscre = nscrea + nscreb + nscrec + nscred;

       %indic3(vbl=qdecoff9094, prefix=qdecoff9094, min=0, max=2, reflev=0, missing=., usemiss=1);
       %indic3(vbl=qcoffee9094, prefix=qcoffee9094, min=0, max=2, reflev=0, missing=., usemiss=1);
       %indic3(vbl=qsodium9094, prefix=qsodium9094, min=0, max=2, reflev=0, missing=., usemiss=0);

if prmeat9094<0.1 then qprmeat9094=0;
else if prmeat9094<0.25 then qprmeat9094=1;
else qprmeat9094=2;

if urmeat9094<0.5 then qurmeat9094=0;
else if urmeat9094<1 then qurmeat9094=1;
else qurmeat9094=2;

if qprmeat9094=0 then medqprmeat9094=0.035;
if qprmeat9094=1 then medqprmeat9094=0.14;
if qprmeat9094=2 then medqprmeat9094=0.425;

if qurmeat9094=0 then medqurmeat9094=0.316;
if qurmeat9094=1 then medqurmeat9094=0.676;
if qurmeat9094=2 then medqurmeat9094=1.23;


/*  processed red meat   */
%indic3(vbl=qprmeat9094, prefix=qprmeat9094, min=1, max=2, reflev=0, missing=., usemiss=0); 

/*  unprocessed red meat   */
%indic3(vbl=qurmeat9094, prefix=qurmeat9094, min=1, max=2, reflev=0, missing=., usemiss=0); 

  covqprmeat=qprmeat9094;
    %indic3(vbl=covqprmeat, prefix=covqprmeat, min=1, max=2, reflev=0, missing=., usemiss=0); 
  covqurmeat=qurmeat9094;
    %indic3(vbl=covqurmeat, prefix=covqurmeat, min=1, max=2, reflev=0, missing=., usemiss=0); 

  medqprmeat9094_2=medqprmeat9094;
  medqprmeat9094=medqprmeat9094*2; * half serving per day;
  
  av_ntotal_yrs=av_ntotal/(-0.045);
  av_nverbl_yrs=av_nverbl/(-0.045);
  av_zscre_yrs=av_zscre/(-0.069);

ahei8094= mean (of nAHEI80a ahei2010_84 ahei2010_86 ahei2010_90 ahei2010_94);
ahei9094= mean (of ahei2010_90 ahei2010_94);

run;

data nhscogfn_metabolites;
set nhscogfn_rmeat;
keep id av_ntotal_yrs av_nverbl_yrs av_zscre_yrs av_ntotal av_nverbl av_zscre;
run;

data nhscogfn_repeated_metabolites;
set nhscogfn_rmeat;
keep id 
ntotal ntotalb ntotalc ntotald 
verbn verbnb verbnc verbnd
zscre zscreb zscrec zscred;  
run; 

libname curr "data/dementia/";

  data curr.NHS_TICS_only;
  set nhscogfn_metabolites; 
  run;
 