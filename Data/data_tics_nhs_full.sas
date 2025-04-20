options nocenter linesize=125 pagesize=78;

*Read in NHS files*;
 filename nhstools '/proj/nhsass/nhsas00/nhstools/sasautos/';
 filename channing '/usr/local/channing/sasautos/';
 libname readfmt '/proj/nhsass/nhsas00/formats/';

 options mautosource sasautos=(channing nhstools);
 options fmtsearch=(readfmt);
 *options obs=5000;


 %include '/udd/hpypl/review/dietscore/nhs1.exercise.sas';

 %include "scripts/food_nhs.sas";

 %amed8020(); run;

 /**** neighborhood SES - nSES****/ 
 libname curr "";
      data nses8616; set curr.anses8616;  
      proc sort nodupkey; by id; run;


*************************  below is code from nhscogfn_covariates01.sas w.o. any editing until genetic data and data merging  **************************** 
* Purpose: To provide SAS code for reading in NHS cognitive file and deriving commonly used covariates in NHS cognitive analyses;
*   Lisa Li checked the code for nhscogfn_covariates01.090109.sas and approved it, 09/08/09;
*   Please see README file for all the changes made since 09/08/09; 
**********************************************************************************************************************************************************;

data cogfn1234;
infile '/proj/nhalzs/nhalz0a/data/allcogfn1234wv.data.DER.090109' lrecl=173 recfm=d;

input
@1              id                  6.
@7              intage              6.
@31             intdt               4.
@35             intdt2              4.
@39             intdt3              4.
@43             intdt4              4.
@107            animal              2.
@109            animalb             2.
@111            animalc             2.
@113            animald             2.
@115            digbac              2.
@117            digbacb             2.
@119            digbacc             2.
@121            digbacd             2.
@123            score               2.
@125            scoreb              2.
@127            scorec              2.
@129            scored              2.
@131            story1              2.
@133            story1b             2.
@135            story1c             2.
@137            story1d             2.
@139            story2              2.
@141            story2b             2.
@143            story2c             2.
@145            story2d             2.
@147            words               2.
@149            wordsb              2.
@151            wordsc              2.
@153            wordsd              2.
@155            words2              2.
@157            words2b             2.
@159            words2c             2.
@161            words2d             2.
@163            hiedu               1.
@164            edumiss             1.
@165            istatus1            2.
@167            fustatus2           2.
@169            fustatus3           2.
@171            fustatus4           2.
;

label
                id              = 'ID NUMBER'
                intage          = '1st wave: Age in yrs @ interview (to two decimals)'
                intdt           = '1st wave: Interview date (# of months since 1900)'
                intdt2          = '2nd wave: Interview date (# of months since 1900)'
                intdt3          = '3rd wave: Interview date (# of months since 1900)'
                intdt4          = '4th wave: Interview date (# of months since 1900)'
                animal          = '1st wave: Composite score for Verbal fluency'
                animalb         = '2nd wave: Composite score for Verbal fluency'
                animalc         = '3rd wave: Composite score for Verbal fluency'
                animald         = '4th wave: Composite score for Verbal fluency'
                digbac          = '1st wave: Composite score for Digits Backwards'
                digbacb         = '2nd wave: Composite score for Digits Backwards'
                digbacc         = '3rd wave: Composite score for Digits Backwards'
                digbacd         = '4th wave: Composite score for Digits Backwards'
                score           = '1st wave: Composite score for TICS'
                scoreb          = '2nd wave: Composite score for TICS'
                scorec          = '3rd wave: Composite score for TICS'
                scored          = '4th wave: Composite score for TICS'
                story1          = '1st wave: East Boston Memory Test immediate recall'
                story1b         = '2nd wave: East Boston Memory Test immediate recall'
                story1c         = '3rd wave: East Boston Memory Test immediate recall'
                story1d         = '4th wave: East Boston Memory Test immediate recall'
                story2          = '1st wave: Same as story1 but delayed recall'
                story2b         = '2nd wave: Same as story1 but delayed recall'
                story2c         = '3rd wave: Same as story1 but delayed recall'
                story2d         = '4th wave: Same as story1 but delayed recall'
                words           = '1st wave: 10-word list part of the TICS'
                wordsb          = '2nd wave: 10-word list part of the TICS'
                wordsc          = '3rd wave: 10-word list part of the TICS'
                wordsd          = '4th wave: 10-word list part of the TICS'
                words2          = '1st wave: Same as words but delayed recall'
                words2b         = '2nd wave: Same as words but delayed recall'
                words2c         = '3rd wave: Same as words but delayed recall'
                words2d         = '4th wave: Same as words but delayed recall'
                hiedu           = 'Highest attained educational status'
                edumiss         = 'Flag for those whom HIEDU was imputed'
                istatus1        = '1st wave: status of response'
                fustatus2       = '2nd wave: status of response'
                fustatus3       = '3rd wave: status of response'
                fustatus4       = '4th wave: status of response'
;
proc sort; by id;

 %n767880(keep=hbp76 chol76 db76 hbp78 chol78 db78 hbp80 chol80 db80);run;
 %nur82 (keep=hbp82 chol82 db82);run;
 %nur84 (keep=hbp84 chol84 db84);run;
 %nur86 (keep=hbp86 chol86 db86);run;
 %nur88 (keep=hbp88 chol88 db88);run;
 %nur90 (keep=hbp90 chol90 db90);run;
 %nur92 (keep=batch92 husbe92 mdem92 fdem92 sdem92 hbp92 chol92 db92 wt92 marry92 alone92); run;
 %nur94 (keep=batch94 aspd94 ibu94 alz94 alzd94 hbp94 chol94 db94 
              thiaz94 lasix94 ccblo94 betab94 bprx94 wt94); run;
 %nur96 (keep=batch96 antid96 vite96 vitc96 mvt96 aspd96 ibu96 alz96 alzd96 hbp96 chol96 db96
              thiaz96   lasix96    ccblo96    betab96    ace96   bprx96 wt96 marry96 alone96); run;
 %nur98 (keep=q98 antid98 aspd98 ibud98 hbp98 chol98 db98  
              thiaz98   lasix98    ccblo98    betab98    ace98   bprx98 wt98); run;
 %nur00 (keep=q00 antid00 przc00 zol00 paxil00 celex00 vite00 vitc00 mvit00 aspd00 ibud00
              pep00  nerv00 down00  calm00  energ00 blue00 worn00 happy00 tired00 alz00 alzd00 hbp00 chol00 db00
              thiaz00   lasix00    ccblo00    betab00    ace00   bprx00 wt00 marry00 alone00); run;
 %nur02 (keep=alz02 alzd02 hbp02 chol02 db02
              bprx02 lasix02 thiaz02 ace02 ccblo02 betab02); run;
 %nur04 (keep=q04 alz04 alzd04 mdem04 mdemd04 fdem04 fdemd04 sdem04 sdemd04 q28pt04 hbp04 chol04 db04
              bprex04 lasix04 thiaz04 ace04 ccblo04 betab04 k04); run; 

 %ses92(keep= pep92  nerv92 down92  calm92  energ92 blue92 worn92 happy92 tired92); run;
 %ses96(keep= pep96  nerv96 down96  calm96  energ96 blue96 worn96 happy96 tired96); run;
 
 %der7620(keep= irt76   irt78   irt80   irt82   irt84   irt86   irt88   irt90   irt92   irt94   irt96   irt98   irt00   irt02   irt04
                hrt76   hrt78   hrt80   hrt82   hrt84   hrt86   hrt88   hrt90   hrt92   hrt94   hrt96   hrt98   hrt00   hrt02   hrt04
                bmi76   bmi78   bmi80   bmi82   bmi84   bmi86   bmi88   bmi90   bmi92   bmi94   bmi96   bmi98   bmi00   bmi02   bmi04
                can76   can78   can80   can82   can84   can86   can88   can90   can92   can94   can7694
                qt76    qt78    qt80    qt82    qt84    qt86    qt88    qt90    qt92    qt94    qt96    qt98    qt00    qt02    qt04
                                                                                        smkdr94 smkdr96 smkdr98 smkdr00 smkdr02 smkdr04
                                                                                        nhor94  nhor96  nhor98  nhor00  nhor02  nhor04 
                                                                                                                mobf    yobf    namnp04); 
              can7694=0; if  can76=1 or can78=1 or can80=1 or can82=1 or can84=1 or can86=1 or can88=1 or can90=1 or can92=1 or can94=1 then can7694=1;
              
 %act8614 (keep= act94m act96m act98m act00m); run;

 %n94_dt (keep= vite94d vitc94d mvt94d); run;
 %n98_dt (keep= vite98d vitc98d mvt98d); run;

 %n94_nts (keep= alco94n);
 %n98_nts (keep= alco98n);

**************************************    genetic data       ****************************** ;

/*** genetic data ***/
proc import datafile='apoe_with_gsa_after_exclusion.csv'
            out=allnhapoe
            dbms=csv
            replace; 
     getnames=yes;  
run;
/*proc freq; tables APOE; run;*/

    data APOEDATA; 
      set allnhapoe; 
      if Study='NHS';
      id=studyID;
      /* Apoe_2cat */ 
           if APOE = 'e2e2' then apoe_2cat=0;
      else if APOE = 'e2e3' then apoe_2cat=0;
      else if APOE = 'e3e2' then apoe_2cat=0;
      else if APOE = 'e3e3' then apoe_2cat=0;
      else if APOE = 'e2e4' then apoe_2cat=1;
      else if APOE = 'e4e2' then apoe_2cat=1;
      else if APOE = 'e4e3' then apoe_2cat=1;
      else if APOE = 'e3e4' then apoe_2cat=1;
      else if APOE = 'e4e4' then apoe_2cat=1;
      else apoe_2cat=.; 


      if APOE in ('e2e2') then apoe4grp=1;
       else if APOE in ('e2e3','e3e2') then apoe4grp=2;
        else if APOE in ('e3e3') then apoe4grp=3;
         else if APOE in ('e2e4','e4e2') then apoe4grp=4;
          else if APOE in ('e3e4','e4e3') then apoe4grp=5;
           else if APOE in ('e4e4') then apoe4grp=6;
      apoe4con=apoe4grp-1;

        if platform="affy" then platforms=1;
      else if platform="gsa" then platforms=2;
      else if platform="huco2" then platforms=3;
      else if platform="illu" then platforms=4;
      else if platform="omni" then platforms=5;
      else if platform="onco" then platforms=6; 

      /*proc freq;
      tables apoe4grp apoe_2cat APOE4 platform platforms;
      proc means n nmiss min mean std max;
      var PC:;*/
      proc sort;
      by id;
      run;
 
proc import datafile='data/genetic/aPGS002280_scaled_PRS_1000G_after_exclusion.csv'
            out=PRS
            dbms=csv
            replace; 
     getnames=yes;  
     run;

data nhsPRS;
     set PRS;
     if Study='NHS';
     id=studyID; 
     proc sort;
     by id;
     /*proc freq;
     tables platform;
     proc means n nmiss min mean std max DATA=PRS;
      var PC:;*/
      RUN;


libname hppdcase '/proj/hppars/hppar0g/case/';
/**** PD case in NHS between 1976-2012 ****/
data nhspdcase; set hppdcase.pd12nur; dtdx_pd=dtdx; drop dtdx status;
run; 
proc sort nodupkey; by id; run;
  


************************************************************************************************************************************************
***********************************************************         data merging        ********************************************************
************************************************************************************************************************************************;

data allnhs;
merge cogfn1234(in=a) APOEDATA nhsPRS amed8020
      n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 
      der7620(in=b) act8614 n94_dt n98_dt n94_nts n98_nts ses92 ses96 foodsnhs exernhs nhsPRS nhspdcase nses8616; 
by id;
incogfn=a;
roster=b;
if a and yobf>20;
run;

proc datasets nolist;
   delete cogfn1234 APOEDATA nhsPRS 
      n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 
      der7620 act8614 n94_dt n98_dt n94_nts n98_nts ses92 ses96 foodsnhs exernhs nhsPRS; 
run;

data all;
     set allnhs end=_end_;

*Set questionnaire dates that are out of range*;
 array irta {15} irt76 irt78 irt80 irt82 irt84 irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04;

 do i=1 to dim(irta);
   if (irta{i} < (894+24*i)) or (irta{i} >= (918+24*i)) then irta{i} = 894+24*i;
 end;

*Identify women who did not complete a long questionnaire*;
 m92=0; m94=0; m96=0; m98=0; m00=0; m04=0; 

 if batch92 ge 800 or batch92 lt 1 then m92=1;
 if batch94 ge 800 or batch94 lt 1 then m94=1;
 if batch96 ge 800 or batch96 lt 1 then m96=1;
 if q98 < 1 or q98 > 2             then m98=1;
 if q00 < 1 or q00 > 2             then m00=1;
 if q04 < 1 or q04 > 2             then m04=1; 

**Anti-depressants - set missings to no**;
**asked only from 1996                 **;

      if antid96=1  then ad96=1;
 else if antid96^=1 then ad96=0;
 
      if antid98=1  then ad98=1;
 else if antid98^=1 then ad98=0;

 if przc00=1 or zol00=1 or paxil00=1 or celex00=1 or antid00=1 then ad00=1;
 else ad00=0;


**************;
*SF 36 scales*;
**************;
*1992*;
******;
*RECODE ANY VALUES OUT SIDE OF THE PLAUSIBLE RANGE TO MISSING;
array sixpt pep92 energ92 worn92 tired92 nerv92 down92 calm92 blue92 happy92 ;
do over sixpt; if sixpt > 6 or sixpt < 1 then sixpt=.; end;

*Mental Health Score - based on 5 questions;
*this scale is positively scored so that
 the higher the score the better the mental health;
*DPRS92:1=lower score(better),2=higher           *;
**************************************************;
*reverse score two items:  calm92 and happy92;
array mhrev calm92 happy92;
do over mhrev;mhrev=7-mhrev;end;

*determine non-response to questions comprising the index;
*then compute the overall mean of questions;
nmhi=nmiss(nerv92,down92,calm92,blue92,happy92);
     if nmhi > 2 then rmhi=.;
else                  rmhi=mean(nerv92,down92,calm92,blue92,happy92);

*standardize the mean;
min=1;max=6;
mhi92=(rmhi-min)*100/(max-min);
label mhi92='SF36 MENTAL HEALTH INDEX 1992';

*dummy variable for MHI 1992 for logistic regression;
*set up a variable for depressed or not for logistic regression;
     if  0 <= mhi92 <=52 then dprs92=2;
else if 52 <  mhi92      then dprs92=1;
else                          dprs92=3;
label dprs92='SF36 MH Index 1992 Dichotomized at 52';

*Score SF36 energy/fatigue index (based on 4 questions);
*scale is positively scored, so that the higher the score
the better the energy/vitality;
*EFT92:1=lower score(better),2=higher                   *;
*55 is the mean US for those over 65                    *;
*********************************************************;
*reverse score two items:  pep92 and energ92;
array eftrev pep92 energ92;
do over eftrev; eftrev=7-eftrev; end;

*determine non-response to questions comprising the index;
*then compute the overall mean of questions;
   neft=nmiss(pep92,energ92,worn92,tired92);
   if neft > 2 then reft=.;
else                reft=mean(pep92,energ92,worn92,tired92);

*standardize the mean;
min=1;max=6;
eft92=(reft-min)*100/(max-min);
label eft92='SF36 VITALITY SCALE INDEX 1992';
     if  0 <=eft92<50 then eftr92=2;
else if 50 <=eft92    then eftr92=1;
else                       eftr92=3;
label eftr92= 'SF36 EF Index 1992 Dichotomized at 50';

*1996;
*****;
*RECODE ANY VALUES OUT SIDE OF THE PLAUSIBLE RANGE TO MISSING;
array sixpt6 pep96 energ96 worn96 tired96 nerv96 down96 calm96 blue96 happy96 ;
do over sixpt6; if sixpt6 > 6 or sixpt6 < 1 then sixpt6=.; end;

*Mental Health Score - based on 5 questions;
*this scale is positively scored so that
 the higher the score the better the mental health;
*DPRS96:1=lower score(better),2=higher           *;
**************************************************;
*reverse;
array mhrev6 calm96 happy96; do over mhrev6;mhrev6=7-mhrev6;end;

*determine non-response to questions comprising the index;
*compute the overall mean of questions;
   nmhi6=nmiss(nerv96,down96,calm96,blue96,happy96);
   if nmhi6 > 2 then rmhi6=.;
else                 rmhi6=mean(nerv96,down96,calm96,blue96,happy96);

*standardize the mean;
min=1;max=6;
mhi96=(rmhi6-min)*100/(max-min);
label mhi96='SF36 MENTAL HEALTH INDEX 1996';

*dummy variable for MHI 96 for logistic regression;
*set up a variable for depressed or not;
*for logistic regression;
     if  0<=mhi96<=52 then dprs96=2;
else if 52< mhi96     then dprs96=1;
else                       dprs96=3;
label dprs96 ='SF36 MH 1996 Index Dichotomized at 52';

*Score SF36 energy/fatigue index (based on 4 questions)**;
*scale is positively scored, so that the higher the score
the better the energy/vitality;
*EFTR96:1=lower score(better),2=higher                   *;
*********************************************************;
*reverse;
array eftrev6 pep96 energ96;
do over eftrev6; eftrev6=7-eftrev6; end;

*determine non-response to questions comprising the index;
*compute the overall mean of questions;
neft6=nmiss(pep96,energ96,worn96,tired96);
if neft6 > 2 then reft6=.;
else              reft6=mean(pep96,energ96,worn96,tired96);

*standardize the mean;
min=1;max=6;
eft96=(reft6-min)*100/(max-min);
label eft96='SF36 VITALITY SCALE INDEX 1996';

     if  0<=eft96<50 then eftr96=2;
else if 50<=eft96    then eftr96=1;
else                      eftr96=3;
label eftr96= 'SF36 EF Index Dichotomized at 50';

*2000;
*****;
*RECODE ANY VALUES OUT SIDE OF THE PLAUSIBLE RANGE TO MISSING;
array sixpt00 pep00 energ00 worn00 tired00 nerv00 down00 calm00 blue00 happy00 ;
do over sixpt00; if sixpt00 > 6 or sixpt00 < 1 then sixpt00=.; end;

*Mental Health Score - based on 5 questions;
*this scale is positively scored so that
 the higher the score the better the mental health;
*DPRS00:1=lower score(better),2=higher           *;
**************************************************;
*reverse;
array mhrev00 calm00 happy00; do over mhrev00;mhrev00=7-mhrev00;end;

*determine non-response to questions comprising the index;
*compute the overall mean of questions;
   nmhi00=nmiss(nerv00,down00,calm00,blue00,happy00);
   if nmhi00 > 2 then rmhi00=.;
else                  rmhi00=mean(nerv00,down00,calm00,blue00,happy00);

*standardize the mean;
min=1;max=6;
mhi00=(rmhi00-min)*100/(max-min);
label mhi00='SF36 MENTAL HEALTH INDEX 2000';

*dummy variable for MHI 00 for logistic regression;
*set up a variable for depressed or not;
*for logistic regression;
     if  0<=mhi00<=52 then dprs00=2;
else if 52< mhi00     then dprs00=1;
else                       dprs00=3;
label dprs00 ='SF36 MH 2000 Index Dichotomized at 52';

*Score SF36 energy/fatigue index (based on 4 questions)**;
*scale is positively scored, so that the higher the score
the better the energy/vitality;
*EFTR00:1=lower score(better),2=higher                   *;
*********************************************************;
*reverse;
array eftrev00 pep00 energ00;
do over eftrev00; eftrev00=7-eftrev00; end;

*determine non-response to questions comprising the index;
*compute the overall mean of questions;
neft00=nmiss(pep00,energ00,worn00,tired00);
if neft00 > 2 then reft00=.;
else               reft00=mean(pep00,energ00,worn00,tired00);

*standardize the mean;
min=1;max=6;
eft00=(reft00-min)*100/(max-min);
label eft00='SF36 VITALITY SCALE INDEX 1900';

     if  0<=eft00<50 then eftr00=2;
else if 50<=eft00    then eftr00=1;
else                      eftr00=3;
label eftr00= 'SF36 EF Index Dichotomized at 50';

**Alcohol - 0, 1-14, 15+ g/d keep missings separate**;
*****************************************************;
      if        alco94n=0      then alccat94=1;
 else if 0   <  alco94n < 15.0 then alccat94=2;
 else if 15.0<= alco94n        then alccat94=3;
 else                               alccat94=4;

      if        alco98n=0      then alccat98=1;
 else if 0   <  alco98n < 15.0 then alccat98=2;
 else if 15.0<= alco98n        then alccat98=3;
 else                               alccat98=4;

**Smoking - put missings in non-smoking category**;
**************************************************;
 array rsmk{4} smkdr94 smkdr96 smkdr98 smkdr00;
 array dsmk{4} smoke94 smoke96 smoke98 smoke00;

 do k=1 to 4;
        if      rsmk{k} in (0,.) then dsmk{k}=1;
   else if      rsmk{k}  =  1    then dsmk{k}=1;
   else if 2 <= rsmk{k} <=  8    then dsmk{k}=2;
   else if 9 <= rsmk{k} <= 15    then dsmk{k}=3;
 end;

**Physical activity - indicate missings**;
*****************************************;
 array actr{4} act94m act96m act98m act00m;
 array actd{4} act94  act96  act98  act00;

 do k=1 to 4;
   if actr{k} >=998 then actd{k}=.;
   else                  actd{k}=actr{k};
 end;

**BMI - condense categories to <22, 22-24, 25-29, 30+ kg/m2, keep missings separate*;
************************************************************************************;
 array rbmi{4} qt94   qt96   qt98   qt00;
 array dbmi{4} dbmi94 dbmi96 dbmi98 dbmi00;

**COLLAPSE VARIABLES IN CATEGORIES         *;
**1=<22,2=22<=bmi<25,3=25<=bmi<30,4=30<=bmi*;

do k=1 to 4;
     if  1<=rbmi{k}<= 3 then dbmi{k}=1;
else if  4<=rbmi{k}<= 6 then dbmi{k}=2;
else if  7<=rbmi{k}<= 9 then dbmi{k}=3;
else if 10<=rbmi{k}<=13 then dbmi{k}=4;
else                         dbmi{k}=5;
end;

**High blood pressure**;
***********************;
hbp94f=0; hbp96f=0; hbp98f=0; hbp00f=0;
 if hbp76=1 or hbp78=1 or hbp80=1 or hbp82=1 or hbp84=2 or hbp86=1 or 
               hbp88=1 or hbp90=1 or hbp92=1 or hbp94=1 then hbp94f=1; 
 if hbp94f=1  or hbp96 =1                               then hbp96f=1;
 if hbp96f=1  or hbp98 =1                               then hbp98f=1;
 if hbp98f=1  or hbp00 =1                               then hbp00f=1;

**High cholesterol**;
********************;
chol94f=0; chol96f=0; chol98f=0; chol00f=0;
 if chol76=1 or chol78=1 or chol80=1 or chol82=1 or chol84=2 or chol86=1 or 
                chol88=1 or chol90=1 or chol92=1 or chol94=1 then chol94f=1; 
 if chol94f=1 or chol96=1                                    then chol96f=1;
 if chol96f=1 or chol98=1                                    then chol98f=1;
 if chol98f=1 or chol00=1                                    then chol00f=1;

**MI**;
******;
mi94=0; mi96=0; mi98=0; mi00=0;
 if hrt76 in (1,3) or hrt78 in (1,3) or hrt80 in (1,3) or hrt82 in (1,3) or hrt84 in (1,3) or hrt86 in (1,3) or hrt88 in (1,3)
                   or hrt90 in (1,3) or hrt92 in (1,3) or hrt94 in (1,3) then mi94=1; 
 if mi94=1 or hrt96 in (1,3)                                             then mi96=1; 
 if mi96=1 or hrt98 in (1,3)                                             then mi98=1; 
 if mi98=1 or hrt00 in (1,3)                                             then mi00=1; 

**Diabetes**;
************;
db94f=0; db96f=0; db98f=0; db00f=0;
 if db76=1 or db78=1 or db80=1 or db82=1 or db84=2 or db86=1 or 
    db88=1 or db90=1 or db92=1 or db94=1 then db94f=1; 
 if db94f=1 or db96=1                    then db96f=1;
 if db96f=1 or db98=1                    then db98f=1;
 if db98f=1 or db00=1                    then db00f=1;

**Alzheimer's disease**;
***********************;
 if alz94^=1 then alz94=0;
 if alz96^=1 then alz96=0;
 if alz00^=1 then alz00=0;
 if alz02^=1 then alz02=0;
 if alz04^=1 then alz04=0;

alz94c=0; alz96c=0; alz98c=0; alz00c=0;
 if (alz94=1) or (alz96=1 & alzd96=1)                                                then alz94c=1;  
 if (alz94c=1) or (alz96=1)                                                          then alz96c=1;  
 if (alz96c=1) or (alz00=1 & alzd00=1)                                               then alz98c=1;  
 if (alz98c=1) or (alz00=1) or (alz02=1 & alzd02=1) or (alz04=1 and alzd04 in (1,2)) then alz00c=1;  

**Family history of Alzheimer's disease -- set pt's and short qq to missing**;
*****************************************************************************;
 *1992: reported cases had age dx <70 years - therefore, these could be combined with 2004 cases with dx age <65 years, if specific
 family relations were collapsed*;
 *coded separately for mom, dad, and sibling as: 0=no, 1=yes-<70 years, 9=missing*;

 if m92=0 then do;
   if mdem92=1 then mdem92f=1; else mdem92f=0;
   if fdem92=1 then fdem92f=1; else fdem92f=0;
   if sdem92=1 then sdem92f=1; else sdem92f=0;
 end;
 else if m92=1 then do;
   mdem92f=9; fdem92f=9; sdem92f=9;
 end;

 *2004: cases were specified as mom, dad, or sibling and dx age given. due to the small number of missing ages at dx, these were put into the largest category - 65+ years*;
 *Xdem04:  1=no, 2=yes, 3=pt*;
 *Xdemd04: 1=<55, 2=55-64, 3=65+, 4=pt*; 
 *coded separately for mom, dad, and sibling as: 0=no, 1=yes-<65 years, 2=yes-65+ years, 9=unknown case status*; 

 array fam041 {3} mdem04  fdem04  sdem04;
 array fam042 {3} mdemd04 fdemd04 sdemd04;
 array fam043 {3} mdem04f fdem04f sdem04f;

 do k=1 to 3; 
                                                 fam043{k}=0;
        if fam041{k}=2 & fam042{k} in (1,2) then fam043{k}=1;
   else if fam041{k}=2 & fam042{k}=3        then fam043{k}=2;
   else if fam041{k}=2 & fam042{k} in (4,.) then fam043{k}=2;
   else if fam041{k} in (3,.)               then fam043{k}=9;
 end;

 *family history of dementia*;
 array fam921 {3} mdem92f fdem92f sdem92f;
 array famy   {3} mdem    fdem    sdem;
 
 do k=1 to 3;
    
   *no history*;
	if fam043{k}=0 and fam921{k} in (0,9)                       then famy{k}=0;

   *65+ years: if fam member had dementia at age 68, 
    then would be =1 in 92 (ie <70) and =2 in 04, 
    (ie >65)  so =2 in 04 takes priority*;
   else if fam043{k}=2                                              then famy{k}=2;

   *<65 years*;
   else if fam043{k}=1 or (fam921{k}=1 and fam043{k} in (., 1, 9))  then famy{k}=1;

   else                                                                  famy{k}=9;
end;

*POSTMENOPAUSAL HORMONE USE*;
****************************;
array rnhor {4} nhor94  nhor96  nhor98  nhor00;
array dnhor {4} dehor94 dehor96 dehor98 dehor00;
do k=1 to 4;
        if rnhor{k}=2 then dnhor{k}=1;
   else if rnhor{k}=4 then dnhor{k}=2;
   else if rnhor{k}=3 then dnhor{k}=3;
   else                    dnhor{k}=4;
end;


**Vitamin E supplements - set pt's and short qq to missing, missing only vite question to no*;
*********************************************************************************************;

                                   vte94=1;
      if vite94d=2 and m94^=1 then vte94=2;
 else if vite94d=3 or  m94=1  then vte94=3;

                                   vte96=1;
      if vite96 =2 and m96^=1 then vte96=2;
 else if vite96 =3 or  m96=1  then vte96=3;

                                   vte98=1;
      if vite98d=2 and m98^=1 then vte98=2;
 else if vite98d=3 or  m98=1  then vte98=3;

                                   vte00=1;
      if vite00 =2 and m00^=1 then vte00=2;
 else if vite00 =3 or  m00=1  then vte00=3;


**Vitamin C supplements - set pt's and short qq to missing, missing only vitc question to no*;
*********************************************************************************************;


                                   vtc94=1;
      if vitc94d=2 and m94^=1 then vtc94=2;
 else if vitc94d=3 and m94^=1 then vtc94=3;
 else if vitc94d=4 or  m94=1  then vtc94=4;

                                   vtc96=1;
      if vitc96 =2 and m96^=1 then vtc96=2;
 else if vitc96 =3 and m96^=1 then vtc96=3;
 else if vitc96 =4 or  m96=1  then vtc96=4;

                                   vtc98=1;
      if vitc98d=2 and m98^=1 then vtc98=2;
 else if vitc98d=3 and m98^=1 then vtc98=3;
 else if vitc98d=4 or  m98=1  then vtc98=4;

                                   vtc00=1;
      if vitc00 =2 and m00^=1 then vtc00=2;
 else if vitc00 =3 and m00^=1 then vtc00=3;
 else if vitc00 =4 or  m00=1  then vtc00=4;

**Multivitamins - set pt's and short qq to missing, missing only multivit question to no*;
*****************************************************************************************;

                                  mt94=1;
      if mvt94d=2 and m94^=1 then mt94=2;
 else if mvt94d=3 or m94=1   then mt94=3;

                                  mt96=1;
      if mvt96=2 and m96^=1  then mt96=2;
 else if mvt96=3 or m96=1    then mt96=3;

                                  mt98=1;
      if mvt98d=2 and m98^=1 then mt98=2;
 else if mvt98d=3 or m98=1   then mt98=3;

                                  mt00=1;
      if mvit00=1 and m00^=1 then mt00=2;
 else if mvit00=3 or m00=1   then mt00=3;

*ASPIRIN*;
*********;
        if aspd94 in (.,7)                    then aspi94=4;
   else if aspd94=1                           then aspi94=1;
   else if (aspd94=2 or aspd94=3)             then aspi94=2;
   else if (aspd94=4 or aspd94=5 or aspd94=6) then aspi94=3;

        if aspd96 in (.,7)                    then aspi96=4;
   else if aspd96=1                           then aspi96=1;
   else if (aspd96=2 or aspd96=3)             then aspi96=2;
   else if (aspd96=4 or aspd96=5 or aspd96=6) then aspi96=3;

        if aspd98 in (.,7)                    then aspi98=4;
   else if aspd98=1                           then aspi98=1;
   else if (aspd98=2 or aspd98=3)             then aspi98=2;
   else if (aspd98=4 or aspd98=5 or aspd98=6) then aspi98=3;

    	if aspd00=5                           then aspi00=4; 
   else if aspd00=.                           then aspi00=1; 
   else if (aspd00=1 or aspd00=2)             then aspi00=2; 
   else if (aspd00=3 or aspd00=4)             then aspi00=3; 

*IBUPROFEN*;
*asked as yes/no in 94, 96 (1=use), in 98 asked as frequency of use*;
********************************************************************;
     if m94=1        then ibui94=3;
else if ibu94=1      then ibui94=2;
else if ibu94 ne 1   then ibui94=1;

     if m96=1        then ibui96=3;
else if ibu96=1      then ibui96=2;
else if ibu96 ne 1   then ibui96=1;


     if m98=1        then ibui98=3;
else if 2<=ibud98<=4 then ibui98=2;
else                      ibui98=1;

     if m00=1        then ibui00=3;
else if 2<=ibud00<=4 then ibui00=2;
else                      ibui00=1;

 /***antihypertensive***/ 
           if (thiaz94=1 or lasix94 =1 or ccblo94 =1 or betab94 =1 or bprx94=1 ) then htnrx94=1;   else htnrx94=0;
           if (thiaz96=1 or lasix96 =1 or ccblo96 =1 or betab96 =1 or ace96=1 or bprx96=1 ) then htnrx96=1;   else htnrx96=0;
           if (thiaz98=1 or lasix98 =1 or ccblo98 =1 or betab98 =1 or ace98=1 or bprx98=1 ) then htnrx98=1;   else htnrx98=0;
           if (thiaz00=1 or lasix00 =1 or ccblo00 =1 or betab00 =1 or ace00=1 or bprx00=1 ) then htnrx00=1;   else htnrx00=0; 
           if (thiaz02=1 or lasix02 =1 or ccblo02 =1 or betab02 =1 or ace02=1 or bprx02=1 ) then htnrx02=1;   else htnrx02=0;
           if (bprex04=1 or lasix04=1 or  thiaz04=1 or ace04=1 or ccblo04=1 or  betab04=1 or  k04=1) then htnrx04=1;   else htnrx04=0;


**For updated covariates, select most updated value before initial cognitive interview**;

      if irt94 < intdt <= irt96 then typed=1;
 else if irt96 < intdt <= irt98 then typed=2;
 else if irt98 < intdt <= irt00 then typed=3;
 else if irt00 < intdt          then typed=4;

amed9094    = mean (of amed_90    amed_94);
amed9498    = mean (of amed_94    amed_98);

amed_frt9094    = mean (of amed_frt90    amed_frt94);
amed_veg9094    = mean (of amed_veg90    amed_veg94);
amed_leg9094    = mean (of amed_leg90    amed_leg94);
amed_nut9094    = mean (of amed_nut90    amed_nut94);
amed_rmt9094    = mean (of amed_rmt90    amed_rmt94);
amed_fish9094    = mean (of amed_fish90    amed_fish94);
amed_whgrn9094    = mean (of amed_whgrn90    amed_whgrn94);
amed_etoh9094    = mean (of amed_etoh90    amed_etoh94);
amed_ms9094    = mean (of amed_ms90    amed_ms94);

amed_frt9098    = mean (of amed_frt90    amed_frt98);
amed_veg9098    = mean (of amed_veg90    amed_veg98);
amed_leg9098    = mean (of amed_leg90    amed_leg98);
amed_nut9098    = mean (of amed_nut90    amed_nut98);
amed_rmt9098    = mean (of amed_rmt90    amed_rmt98);
amed_fish9098    = mean (of amed_fish90    amed_fish98);
amed_whgrn9098    = mean (of amed_whgrn90    amed_whgrn98);
amed_etoh9098    = mean (of amed_etoh90    amed_etoh98);
amed_ms9098    = mean (of amed_ms90    amed_ms98);

wtchg94 = wt94 - wt92;
wtchg96 = wt96 - wt94;
wtchg98 = wt98 - wt96;
wtchg00 = wt00 - wt98;
 

 array aanti  {4}  ad96     ad96     ad98     ad00;
 array htnrx  {4}  htnrx94  htnrx96  htnrx98  htnrx00;
 array ahbp   {4}  hbp94f   hbp96f   hbp98f   hbp00f;
 array ahchol {4}  chol94f  chol96f  chol98f  chol00f;
 array ami    {4}  mi94     mi96     mi98     mi00;
 array adb    {4}  db94f    db96f    db98f    db00f;
 array aalz   {4}  alz94c   alz96c   alz98c   alz00c;
 array adprs  {4}  dprs92   dprs96   dprs96   dprs00;
 array aeftr  {4}  eftr92   eftr96   eftr96   eftr00;
 array aaspi  {4}  aspi94   aspi96   aspi98   aspi00;
 array aibu   {4}  ibui94   ibui96   ibui98   ibui00;
 array aalco  {4}  alccat94 alccat94 alccat98 alccat98;
 array aalcor {4}  alco94n  alco94n  alco98n  alco98n;
 array calor  {4}  calor94n calor94n calor98n calor98n;
 array asmok  {4}  smoke94  smoke96  smoke98  smoke00;
 array adehor {4}  dehor94  dehor96  dehor98  dehor00;
 array aphy   {4}  act94    act96    act98    act00;
 array abmi   {4}  dbmi94   dbmi96   dbmi98   dbmi00;
 array avite  {4}  vte94    vte96    vte98    vte00;
 array avitc  {4}  vtc94    vtc96    vtc98    vtc00;
 array amvit  {4}  mt94     mt96     mt98     mt00;

 array amed   {4}  amed9094 amed9094 amed9498 amed9498; 
 array actm   {4}  th94     th96     th98     th00;       
 array SES    {4}  nSES_94  nSES_96  nSES_98  nSES_00;
 array bmicons{4}  bmi94    bmi96    bmi98    bmi00;
 array wtchg  {4}  wtchg94  wtchg96  wtchg98  wtchg00;
 array married{4}  marry92  marry92  marry96  marry96;
 array alone  {4}  alone92  alone92  alone96  alone96;
 
 array amed_frt   {4}  amed_frt9094 amed_frt9094 amed_frt9498 amed_frt9498; 
 array amed_veg   {4}  amed_veg9094 amed_veg9094 amed_veg9498 amed_veg9498; 
 array amed_leg   {4}  amed_leg9094 amed_leg9094 amed_leg9498 amed_leg9498; 
 array amed_nut   {4}  amed_nut9094 amed_nut9094 amed_nut9498 amed_nut9498; 
 array amed_rmt   {4}  amed_rmt9094 amed_rmt9094 amed_rmt9498 amed_rmt9498; 
 array amed_fish   {4}  amed_fish9094 amed_fish9094 amed_fish9498 amed_fish9498; 
 array amed_whgrn   {4}  amed_whgrn9094 amed_whgrn9094 amed_whgrn9498 amed_whgrn9498; 
 array amed_etoh   {4}  amed_etoh9094 amed_etoh9094 amed_etoh9498 amed_etoh9498; 
 array amed_ms   {4}  amed_ms9094 amed_ms9094 amed_ms9498 amed_ms9498; 


 do i=1 to 4;
   if typed=i then do;
     phy   =  aphy{i};
     anti  =  aanti{i};
     antihp=  htnrx{i}; if antihp ne 1 then antihp=0;
     hbp   =  ahbp{i};
     hchol =  ahchol{i};
     mi    =  ami{i};
     db    =  adb{i};
     alz   =  aalz{i};
     dprs  =  adprs{i};
     eft   =  aeftr{i};
     asp   =  aaspi{i};
     ibu   =  aibu{i};
     alc   =  aalco{i}; 
     alco  =  aalcor{i}; 
     smk   =  asmok{i};
     nhor  =  adehor{i};
     bmi   =  abmi{i};
     vite  =  avite{i};
     vitc  =  avitc{i};
     mvit  =  amvit{i};
     actcon=  actm{i};
     nSES  =  SES{i};
    amedcon=  amed{i};
    calorie=  calor{i};
    bmicon =  bmicons{i}; 
    wtchgcon= wtchg{i};
    marry   = married{i}; if marry ne 1 then marry=0;
    live_alone=alone{i};  if live_alone ne 1 then live_alone=0;
     amed_frtcon=  amed_frt{i};
     amed_vegcon=  amed_veg{i};
     amed_legcon=  amed_leg{i};
     amed_nutcon=  amed_nut{i};
     amed_rmtcon=  amed_rmt{i};
     amed_fishcon=  amed_fish{i};
     amed_whgrncon=  amed_whgrn{i};
     amed_etohcon=  amed_etoh{i};
     amed_mscon=  amed_ms{i};
   end;
 end;


%indic3(vbl=dprs, prefix=dprs, min=2, max=2, reflev=1, missing=3,
        label1='high mental health score',
        label2='low mental health score');
%indic3(vbl=eft, prefix=eft,   min=2, max=2, reflev=1, missing=3,
        label1='high energy score',
        label2='low energy score');
%indic3(vbl=asp, prefix=asp,   min=2, max=3, reflev=1, missing=4,
        label1='asp non-user',
        label2='asp 1-2/wk',
        label3='asp 3-7/wk');
%indic3(vbl=ibu, prefix=ibu,   min=2, max=2, reflev=1, missing=3,
        label1='ibu non-user',
        label2='ibu user');
%indic3(vbl=alc, prefix=alc,   min=2, max=3, reflev=1, missing=4,
        label1='alc non-drinkers',
        label2='alc 1-14g',
        label3='alc 15+ g');
%indic3(vbl=smk, prefix=smk,   min=2, max=3, reflev=1, usemiss=0,
        label1='smk never',
        label2='smk past',
        label3='smk current');
%indic3(vbl=nhor, prefix=nhor, min=2, max=3, reflev=1, missing=4,
        label1='pmh never',
        label2='pmh past',
        label3='pmh current');

if bmi=5 then bmi=.;
%indic3(vbl=bmi, prefix=bmi,   min=2, max=4, reflev=1, missing=.,
        label1='bmi <22',
        label2='bmi 22-24',
        label3='bmi 25-29',
        label4='bmi 30+')

%indic3(vbl=vite, prefix=vite, min=2, max=2, reflev=1, missing=3,
        label1='vitE non-user',
        label2='vitE user');
%indic3(vbl=vitc, prefix=vitc, min=2, max=3, reflev=1, missing=4,
        label1='vitC non-user',
        label2='vitC seasonal user',
        label3='vitC regular user');
%indic3(vbl=mvit, prefix=mvit, min=2, max=2, reflev=1, missing=3,
        label1='mvit non-user',
        label2='mvit user');


*Create continuous, centered age variable*;
age=intage-74;

*AGE*;
*****;
     if 70.00 <= intage < 72.00 then intagec=1;
else if 72.00 <= intage < 74.00 then intagec=2;
else if 74.00 <= intage < 76.00 then intagec=3;
else if 76.00 <= intage < 78.00 then intagec=4;
else if 78.00 <= intage         then intagec=5;

%indic3(vbl=intagec, prefix=age, min=2, max=5, reflev=1, usemiss=0,
        label1='70-71', 
        label2='72-73', 
        label3='74-75',
        label4='76-77', 
        label5='78+');

**EDUCATION*;
************;
*collapse MA and DR*;
if hiedu in (3,4) then hiedu=3;

%indic3(vbl=hiedu, prefix=edu, min=2, max=3, reflev=1, usemiss=0,
        label1='RN', 
        label2='BA', 
        label3='MA or DR');

**HUSBAND'S EDUCATION*;
*husbe92
        $range 1-6
    1=<high school
    2=some H.S.
    3=H.S. grad
    4=College grad
    5=grad school
    6=passthrough
**********************;
husbe=husbe92;
if husbe92 in (.,6) then husbe=6;

%indic3(vbl=husbe, prefix=husbe, min=2, max=5, reflev=1, missing=6,
      label1='<high school', 
      label2='some H.S.', 
      label3='H.S. grad',
      label4='College grad', 
      label5='Grad school');

****** HUSBAND'S EDUCATION ******;  
if husbe in (1,2,3) then husbedu=1;
if husbe=4 then husbedu=2;
if husbe=5 then husbedu=3;
      %indic3(vbl=husbedu, prefix=husbedu, min=2, max=3, reflev=1, missing=., usemiss=1,
            label1='=<H.S. grad',  
            label2='College grad', 
            label3='Grad school');


*AGE AT MENOPAUSE                                                             *;
*age at menopause as of 1996 was used                                         *;
*all became postmenopausal in 1994 but 2004 data is most updated and accurate *;
*******************************************************************************;
     if 20 <= namnp04 < 50 then amnp=1;
else if 50 <= namnp04 < 53 then amnp=2;
else if 53 <= namnp04 < 95 then amnp=3;
else                            amnp=1;

%indic3(vbl=amnp, prefix=amnp, min=2, max=3, reflev=1, usemiss=0,
        label1='<50', 
        label2='50-52', 
        label3='53+');

*FAMILY HISTORY OF DEMENTIA                                   *;
*coded separately for mom, dad, and sibling as:               *; 
*0=no, 1=yes-<65 years, 2=yes-65+ years, 9=unknown case status*; 
***************************************************************;

%indic3(vbl=mdem, prefix=mdem, min=1, max=2, reflev=0, missing=9,
        label0='mom no dem hx',
        label1='mom dem <65 y',
        label2='mom dem 65+ y');

%indic3(vbl=fdem, prefix=fdem, min=1, max=2, reflev=0, missing=9,
        label0='dad no dem hx',
        label1='dad dem <65 y',
        label2='dad dem 65+ y');

%indic3(vbl=sdem, prefix=sdem, min=1, max=2, reflev=0, missing=9,
        label0='sib no dem hx',
        label1='sib dem <65 y',
        label2='sib dem 65+ y');

if 0=<actcon<0.5  	     then actcc=1;
       else if 0.5=<actcon<2 	     then actcc=2;
       else if 2=<actcon<3.5       then actcc=3;
       else if 3.5=<actcon<5.5     then actcc=4;
       else if actcon>=5.5         then actcc=5;
       else actcc=.;

%indic3(vbl=actcc, prefix=actcc, min=2, max=5, reflev=1, missing=., usemiss=1,
	        label1='<0.5 m/v activity hours/week',
	        label2='0.5 to < 2 hours/week',
	        label3='2 to < 3.5 hours/week',
	        label4='3.5 to < 5.5 hours/week',
	        label5='5.5+ hours/week'); 

if antihp=1 or hbp=1 then highhbp=1; else highhbp=0;
if dprs=2 or anti=1 then depre=1; else depre=0;
   
/*** Indicator for BMI ***/ 
if bmicon<23 then bmic=1;
           else if 23=<bmicon<25 then bmic=2;
             else if 25=<bmicon<30 then bmic=3;
               else if 30=<bmicon<35 then bmic=4;
                 else bmic=5;
          if bmicon=. then bmic=.;
         %indic3(vbl=bmic, prefix=bmic, min=1, max=5, reflev=2, missing=., usemiss=0,
                label1='<23',
                label2='23-24.9',
                label3='25-29.9',
                label4='30-34.9',
      	    label5='35+');
 
if wtchgcon=<-5 then wtchgcat=1;
          else if wtchgcon=<-1 then wtchgcat=2;
          else if wtchgcon<1 then wtchgcat=3;
          else if wtchgcon<5 then wtchgcat=4;
          else wtchgcat=5;
          if wtchgcon=. then wtchgcat=.;
         %indic3(vbl=wtchgcat, prefix=wtchgcat, min=1, max=5, reflev=3, missing=., usemiss=1,
                label1='wtchg =<-5 lbs',
                label2='wtchg -4.9 ~ -1 lbs',
                label3='wtchg >-1 ~ <1 lbs',
                label4='wtchg 1 ~ 4.9 lbs',
      	    label5='wtchg >=5 lbs');  

PGScon=PGS002280_scaled;           
run;



************************************************************************************************************************************************
*****************************************************  average cognative function scores   *****************************************************
************************************************************************************************************************************************;

******select the cognitive population (currently, includes all 23565 who were considered eligible,
      but many have never completed a cognitive interview)
      ------ to select just those who have completed a baseline cognitive assessment;
 
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

proc freq data=all;
tables highhbp  db   anti  antihp actcc  smk   nhor  hbp   hchol   mi  alz   dprs  eft  asp  ibu   alc   alco     bmi  vite  vitc   mvit phy wtchgcat ;
tables antihp*hbp;
proc means n nmiss min mean std median max;
var wtchgcon wt94 wt96 wt98 wt00;
proc means median ;
var actcon;
class actcc;
run; 
 


*******************************************************************************************
*                                        Baseline exclusion                               *
*******************************************************************************************;
   
data nhsdata;
      set  all end=_end_   ;  

%beginex();

*Get counts of starting N*;
	%exclude(id gt 0 or id le 0, nodelete=t); 
*Missing ids*;
	%exclude(id le 0);
*Missing first cognitive interview*;
	%exclude(score eq .);
*Incomplete battery - SELECT POPULATION WITH AT LEAST ONE COMPLETE BATTERY DURING FOLLOW-UP (including baseline)*;
	%exclude(cogm eq 1);
*Missing age at baeline *;
        %exclude(intage eq .); 
*Less than 70 (in master file)*;
	%exclude(intage lt 70); 

* implausible FFQ;
   %exclude(calorie lt 500 );
   %exclude(calorie gt 3500 );
   %exclude(calorie le 0);
   %exclude(calorie eq .);

*Missing amed information *;
        %exclude(amedcon eq .);
 
*Missing smoking information *;
        %exclude(smk eq .); 

*Missing bmi information *;
        %exclude(bmi eq .);  

*Get counts of those left*;
	%exclude(id gt 0, nodelete=t); 

%output();

%endex();                    /* MUST do AFTER all exclusions */
run;


*******************************************************************************************
*                                        baseline covariates                              *
*******************************************************************************************;
 
      
%pctl9(data=nhsdata,
    varlist=  nSES PC1_comb PC2_comb PC3_comb PC4_comb PGS002280_scaled,
    numquant=3,
    mscore=T,
    quantname=q,
    cutdat=cutpoints,
    outdat=nhsdata,
    indic=T);
 
%pctl9(data=nhsdata,
    varlist= amedcon PGScon amed_frtcon amed_vegcon amed_legcon amed_nutcon amed_rmtcon amed_fishcon amed_whgrncon amed_etohcon amed_mscon,
    numquant=5,
    mscore=T,
    quantname=q,
    cutdat=cutpoints,
    outdat=nhsdata,
    indic=T);

data nhs1tics;
set nhsdata end=_end_;

agecon=intage;

%indic3(vbl=qamedcon,    prefix=qamedcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_frtcon,    prefix=qamed_frtcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_vegcon,    prefix=qamed_vegcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_legcon,    prefix=qamed_legcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_nutcon,    prefix=qamed_nutcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_rmtcon,    prefix=qamed_rmtcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_fishcon,    prefix=qamed_fishcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_whgrncon,    prefix=qamed_whgrncon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_etohcon,    prefix=qamed_etohcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
%indic3(vbl=qamed_mscon,    prefix=qamed_mscon,   min=1, max=4, reflev=0, missing=., usemiss=0);

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

fhdem=0;
if famdem in (1,2) then fhdem=1;

if qnSES=. then qnSES=1; *nmiss=11, median imputation;
%indic3(vbl=qnSES,           prefix=qnSES,           min=1, max=2, reflev=0, missing=., usemiss=0);  
if nSES=. then medqnSES=-0.805;
 
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

****** P for trend ******;
if actcc=1 then medactcc=0;
if actcc=2 then medactcc=1;
if actcc=3 then medactcc=2.5;
if actcc=4 then medactcc=5;
if actcc=5 then medactcc=8.5;

if actcc=. then medactcc=99999;
actccmiss=0; if actcc=. then actccmiss=1;

      %indic3(vbl=apoe4grp, prefix=apoe4grp,   min=1, max=6, reflev=3, missing=.,usemiss=1,
            label1='e2e2',
            label2='e2e3',
            label3='e3e3',
            label4='e2e4',
            label5='e3e4',
            label6='e4e4');

      %indic3(vbl=qPC1_comb,     prefix=qPC1_comb,     min=1, max=2, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qPC2_comb,     prefix=qPC2_comb,     min=1, max=2, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qPC3_comb,     prefix=qPC3_comb,     min=1, max=2, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qPC4_comb,     prefix=qPC4_comb,     min=1, max=2, reflev=0, missing=., usemiss=1);

      qPGSmed=medqPGS002280_scaled;
      qPGS=qPGS002280_scaled;
      %indic3(vbl=qPGS,     prefix=qPGS,     min=1, max=2, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qPGScon,  prefix=qPGScon,  min=1, max=4, reflev=0, missing=., usemiss=1);
 
      %indic3(vbl=platforms,     prefix=platforms,     min=2, max=6, reflev=1, missing=., usemiss=1,
            label1='affy',
            label2='gsa',
            label3='huco2',
            label4='illu',
            label5='omni',
            label6='onco' );
run;


 
*******************************************************************************************************
*                             final dataset "nhscogfn" for cross-sectional analysis                   *
*******************************************************************************************************;

/*** MAKE Z-SCORES OF EACH TEST ((using MEAN and SD from baseline scores which are computed in the previous macros ***/

/*** 1st Z scores ***/

 %include '/udd/hpypl/review/macro/zscore.sas';
 %zscore(nhs1tics, score,  zscre); 
 %zscore(nhs1tics, animal, zaniml); 
 %zscore(nhs1tics, words,  zwrds); 
 %zscore(nhs1tics, words2, zwrds2); 
 %zscore(nhs1tics, story1, zstry1); 
 %zscore(nhs1tics, story2, zstry2); 
 %zscore(nhs1tics, digbac, zdigbc);

proc means mean std data=nhs1tics;
var score animal words words2 story1 story2 digbac;
run;
 

/*** 2ND ZSCORES USING 1ST MEAN AND SD ***/

data nhscogfn;
set nhs1tics end=_end_;

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
  
  
  av_ntotal_yrs=av_ntotal/(-0.043);
  av_nverbl_yrs=av_nverbl/(-0.045);
  av_zscre_yrs =av_zscre/(-0.066);

/*** Indicator for baseline BMI ***/ 
if bmi80<23 then bmi80c=1;
           else if 23=<bmi80<25 then bmi80c=2;
             else if 25=<bmi80<30 then bmi80c=3;
               else if 30=<bmi80<35 then bmi80c=4;
                 else bmi80c=5;
          if bmi80=. then bmi80c=.;
         %indic3(vbl=bmi80c, prefix=bmi80c, min=1, max=5, reflev=2, missing=., usemiss=1,
                label1='<23',
                label2='23-24.9',
                label3='25-29.9',
                label4='30-34.9',
      	    label5='35+'); 
run;
 
PROC EXPORT DATA=nhscogfn  
            OUTFILE= 'data/dementia/nhs_tics_full.csv' 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

