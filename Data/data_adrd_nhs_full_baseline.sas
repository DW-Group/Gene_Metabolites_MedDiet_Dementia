options linesize=125 pagesize=78; /*set line and page size*/
filename nhstools '/proj/nhsass/nhsas00/nhstools/sasautos/'; /*READ macros*/
filename local '/usr/local/channing/sasautos/';
filename ehmac '/udd/stleh/ehmac/';
libname readfmt '/proj/nhsass/nhsas00/formats/'; /*READ formats*/
options mautosource sasautos=(local nhstools);          *** path to macros  ***;
options fmtsearch=(readfmt);                            *** path to formats ***; 

************************************************************************************************************************************************
***********************************************************      (1) genetic data       ********************************************************
************************************************************************************************************************************************;

/*** genetic data ***/
proc import datafile='data/genetic/apoe_with_gsa_comb_pcs_after_exclusion.csv'
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

proc import datafile='data/genetic/aPGS000334_scaled_PRS_1000G_after_exclusion.csv'
          out=PRS2
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

data nhsPRS2;
     set PRS2;
     if Study='NHS';
     id=studyID; 
     proc sort;
     by id;
     /*proc freq;
     tables platform;
     proc means n nmiss min mean std max DATA=PRS;
      var PC:;*/
      RUN;

************************************************************************************************************************************************
***********************************************************   READ IN OUTCOME DATA   ***********************************************************
************************************************************************************************************************************************;
 
/****** READ IN memory variables */ 

%nur12(keep = id mchng12 mevent12 mitem12 msec12 minstr12 mconv12 msts12 /* memory */
              dem12 demd12 
              lvig12 lmod12 llift12 lclms12 lclm112 lbend12 lwlkm12 lwlks12 lwlkb12 lbath12 /* physival activity */
              feelpt12 slife12 drint12 empty12 bored12 gdspi12 bdhpn12 happy12 helpl12 styhm12 bdmem12 wndrf12 wrthl12 energ12 hopel12 oboff12 /* mood */
              , out=mem12dat);

%nur14(keep = id mchng14 mevent14 mitem14 msec14 minstr14 mconv14 msts14 /* memory */
              feelpt14 slife14 drint14 empty14 bored14 gdspi14 bdhpn14 happy14 helpl14 styhm14 bdmem14 wndrf14 wrthl14 energ14 hopel14 oboff14 /* mood */
              , out=mem14dat);

data memorynhs;
merge mem12dat mem14dat;
by id;

***********************************************************;
/* memory variables assessed in 2012 */
**********************************************************;

array mem12{7}     mchng12  mevent12  mitem12  msec12  minstr12  mconv12  msts12 ;
array mem12_new{7} mchng12n mevent12n mitem12n msec12n minstr12n mconv12n msts12n ;

        do i = 1 to 7;
                if mem12{i} = . then mem12_new{i} = .;
                else if mem12{i} = 1 then mem12_new{i} = 1; /* yes */
                else if mem12{i} = 2 then mem12_new{i} = 0; /* no */
                else if mem12{i} = 3 then mem12_new{i} = .;
        end;


mem12_6nhs = sum (mevent12n, mitem12n, msec12n, minstr12n, mconv12n, msts12n);

mem12_7nhs = sum (mchng12n, mevent12n, mitem12n, msec12n, minstr12n, mconv12n, msts12n);

**************************************************************;
* cognitive domains - suggested by Olivia Okereke ;
**************************************************************;

memory12 = sum (mevent12n, mitem12n);
executive12 = sum (minstr12n, mconv12n);
visual12=msts12n;
attention12=msec12n;

***********************************************************;
/* memory variables assessed in 2014 */
**********************************************************;

array mem14{7} mchng14 mevent14 mitem14 msec14 minstr14 mconv14 msts14 ;
array mem14_new{7} mchng14n mevent14n mitem14n msec14n minstr14n mconv14n msts14n ;

        do i = 1 to 7;
                if mem14{i} = . then mem14_new{i} = .;
                else if mem14{i} = 1 then mem14_new{i} = 1; /* yes */
                else if mem14{i} = 2 then mem14_new{i} = 0; /* no */
                else if mem14{i} = 3 then mem14_new{i} = .;
        end;


mem14_6nhs = sum (mevent14n, mitem14n, msec14n, minstr14n, mconv14n, msts14n);
mem14_7nhs = sum (mchng14n, mevent14n, mitem14n, msec14n, minstr14n, mconv14n, msts14n);


**************************************************************;
* cognitive domains - suggested by Olivia Okereke ;
**************************************************************;

memory14 = sum (mevent14n, mitem14n);
executive14 = sum (minstr14n, mconv14n);
visual14=msts14n;
attention14=msec14n;

**************************************************************;
* dummy variables for missing on either cognition question ;
**************************************************************;

if mem12_6nhs ne . and mem14_6nhs ne . then mem_miss=0;
else if mem12_6nhs eq . and mem14_6nhs ne . then mem_miss=1;
else if mem14_6nhs eq . and mem12_6nhs eq . then mem_miss=2;

/****     Average of two years     ***/

mem1214_6nhs=mean(mem12_6nhs,mem14_6nhs);
mem1214_7nhs=mean(mem12_7nhs,mem14_7nhs);

if mem1214_6nhs=. then mem6_avg3_3cat = .;
else if mem1214_6nhs = 0 then mem6_avg3_3cat=0;
else if 0 < mem1214_6nhs< 3 then mem6_avg3_3cat=1;
else if 3 <= mem1214_6nhs then mem6_avg3_3cat=2;

if mem1214_7nhs=. then mem7_avg3_3cat = .;
else if mem1214_7nhs = 0 then mem7_avg3_3cat=0;
else if 0 < mem1214_7nhs< 3 then mem7_avg3_3cat=1;
else if 3 <= mem1214_7nhs then mem7_avg3_3cat=2;


if mem1214_7nhs=. then mem7_avg3_2cat = .;
else if  mem1214_7nhs< 3 then mem7_avg3_2cat=0;
else if 3 <= mem1214_7nhs then mem7_avg3_2cat=1;

memory1214=mean(memory12,memory14);
executive1214=mean(executive12, executive14);
visual1214=mean(visual12,visual14);
attention1214=mean(attention12, attention14);

zmem1214_7nhs=mem1214_7nhs;
zmem1214_6nhs=mem1214_6nhs;
zmemory1214=memory1214;
zexecutive1214=executive1214;
zvisual1214=visual1214;
zattention1214=attention1214;


/***     Change of two years     ***/
diffmem1214_7nhs=mem14_7nhs-mem12_7nhs;

if diffmem1214_7nhs=. then diffmem1214_cat=.;
else if diffmem1214_7nhs=0 then diffmem1214_cat=1;
else if diffmem1214_7nhs<0 then diffmem1214_cat=2;
else if diffmem1214_7nhs>0 then diffmem1214_cat=3;


if diffmem1214_7nhs=. then diffmem1214_bin=.;
else if diffmem1214_7nhs=0 then diffmem1214_bin=0;
else if diffmem1214_7nhs<0 then diffmem1214_bin=0;
else if diffmem1214_7nhs>0 then diffmem1214_bin=1;

diffmem1214_6nhs=mem14_6nhs-mem12_6nhs;

if diffmem1214_6nhs=. then diffmem1214_6nhs_bin=.;
else if diffmem1214_6nhs=0 then diffmem1214_6nhs_bin=0;
else if diffmem1214_6nhs<0 then diffmem1214_6nhs_bin=0;
else if diffmem1214_6nhs>0 then diffmem1214_6nhs_bin=1;

diffmemory1214 = memory14 - memory12        ;
diffexecutive1214 = executive14 - executive12     ;
diffvisual1214 = visual14 - visual12        ;
diffattention1214 = attention14 - attention12     ;

if diffmemory1214=. then diffmemory1214_bin=.;
else if diffmemory1214=0 then diffmemory1214_bin=0;
else if diffmemory1214<0 then diffmemory1214_bin=0;
else if diffmemory1214>0 then diffmemory1214_bin=1;

if diffexecutive1214=. then diffexecutive1214_bin=.;
else if diffexecutive1214=0 then diffexecutive1214_bin=0;
else if diffexecutive1214<0 then diffexecutive1214_bin=0;
else if diffexecutive1214>0 then diffexecutive1214_bin=1;

if diffvisual1214=. then diffvisual1214_bin=.;
else if diffvisual1214=0 then diffvisual1214_bin=0;
else if diffvisual1214<0 then diffvisual1214_bin=0;
else if diffvisual1214>0 then diffvisual1214_bin=1;

if diffattention1214=. then diffattention1214_bin=.;
else if diffattention1214=0 then diffattention1214_bin=0;
else if diffattention1214<0 then diffattention1214_bin=0;
else if diffattention1214>0 then diffattention1214_bin=1;
 
mem1214_7nhs_sd=mem1214_7nhs/1.3733482;
proc sort;by id;
run;

/* proc means n nmiss min mean std median max;
   var mem1214_7nhs;                           
   N        NMiss         Minimum            Mean         Std Dev          Median         Maximum
   ----------------------------------------------------------------------------------------------
   63795     6218               0       1.0719570       1.3733482       0.5000000       7.0000000
   ---------------------------------------------------------------------------------------------*/

proc standard data=memorynhs mean=0 std=1 out=memorynhs;
     var zmem1214_7nhs zmem1214_6nhs zmemory1214 zexecutive1214 zvisual1214 zattention1214;
run;

************************************************************************************************************************************************
***********************************************************    (1) Dead and diseases    ********************************************************
************************************************************************************************************************************************;
 
/** Read death files **/
 %deadff (keep= id mod yod icda10 qyr state source conf modx yodx nhsicda enhsicda newicda deadmonth dtdth newicda3, file= /proj/nhdats/nh_dat_cdx/deaths/deadff.current.nhs1)
   dtdth=9999;
   if deadmonth> 0 then dtdth=deadmonth;
    newicda=compress(nhsicda, 'E');
    newicda3=substr(newicda,1,3); 
  run;


/*** Diabetes ***/
 data diabetes;
    %include '/proj/nhdats/nh_dat_cdx/endpoints/diabetes/db.cases.input';
    type2db=0;
    if type=2 and probabil=1 then type2db=1; /*type 2 diabetes*/
    type1db=0;
    if type eq 1 or type eq 5 then type1db=1; /* definite or probable type 1 diabetes */

    dtdxdb=9999;
    if dxmonth > 0 then dtdxdb=dxmonth;

    if dtdxdb<918 then db1976=1;

    dtdxdb2=9999;
    if dxmonth > 0 and type2db=1 then dtdxdb2=dxmonth;

    keep id type2db dtdxdb  type1db dtdxdb2 ;
    run; 

/***Read confirmed MI file, no longer split by nfmi and fmi **/
 data tot_mi;
     %include '/proj/nhdats/nh_dat_cdx/endpoints/mi/mi7620.073120.cases.input';
     mi=0;
    if 11<=conf<=19 then mi=1;

    /*for definite or probable cases*/
     mid=0; if mi=1 and conf=11 then mid=1;
     mip=0; if mi=1 and mid ne 1 then mip=1;

   /*account for missing diagnosis month*/
     dtdxmi=9999;
     if dxmonth > 0 and mi=1 then dtdxmi=dxmonth;

    drop icda nhsicda conf;
    run;

/*** read stroke data ***/
 data tot_str;
    %include '/proj/nhdats/nh_dat_cdx/endpoints/stroke/str7620.031621.cases.input';
    str=0;
    if 11<=conf<=19 then str=1;
    dtdxst=9999;
    if dxmonth>0 and 11=<conf=<19 and str=1 then dtdxst=dxmonth;

    /*for defined or probable*/
     strd=0; if str=1 and conf=11 then strd=1;
     strp=0; if str=1 and strd ne 1 then strp=1;

    ***type of stroke cases;
        if str_type in (1,2) then hem_str=1; else hem_str=0;  /*hemorrhagic stroke cases*/
        if str_type in (3,4,5) then isc_str=1; else isc_str=0;/*ischemic stroke cases*/
        if str_type=6 then unk_str=1; else unk_str=0;         /*unknown type of stroke cases*/

    drop icda nhsicda conf;
    run;

/*** read in cancer data ***/
 data cancer;
    %include '/proj/nhdats/nh_dat_cdx/endpoints/allcancer/canc7622.082823.cases.input';
    dtdxca=9999;
    if dxmonth>0 and 11=<conf=<19 then dtdxca=dxmonth;
    keep id qyr mdx1 ydx1 mdx ydx dtdxca;
    Run;
    proc sort data=cancer out=cancer nodupkey; by id; run;


    data disease;
    merge diabetes deadff tot_mi tot_str cancer ;
    by id;

/*make sure death did not occur before date of diagnosis*/
    if dtdth<dtdxdb and dtdth ne . and dtdxdb ne 9999 then dtdth=dtdxdb;
    if dtdth<dtdxmi and dtdth ne . and dtdxmi ne 9999 then dtdth=dtdxmi;
    if dtdth<dtdxst and dtdth ne . and dtdxst ne 9999 then dtdth=dtdxst;
    if dtdth<dtdxca and dtdth ne . and dtdxca ne 9999 then dtdth=dtdxca;

    /*total CVD (mi+stroke)*/
    dtdxcvd=9999;
    if 0<dtdxmi<9999  and 0<dtdxst<9999 then dtdxcvd=min(dtdxmi,dtdxst);
    else if 0<dtdxmi<9999 then dtdxcvd=dtdxmi;     /* mi*/
    else if 0<dtdxst<9999 then dtdxcvd=dtdxst;     /* stroke */

    keep id dtdth type2db dtdxdb type1db dtdxmi dtdxcvd dtdxst dtdxca newicda
         mi hem_str isc_str unk_str dtdxdb2 newicda3;
    run;

proc sort; by id; run;

    PROC DATASETS;
    delete diabetes deadff tot_mi tot_str cancer;
    RUN;


************************************************************************************************************************************************
*********************************************************** (2) READ IN EXPOSURE DATA   ********************************************************
************************************************************************************************************************************************;

%include "scripts/food_nhs.sas"; 

%amed8020(); run;
 
************************************************************************************************************************************************
*********************************************************** (3) READ IN OTHER COVARIATES *******************************************************
************************************************************************************************************************************************;

    %ahei2010_8410(); run; 

    data ahei2010;
    set ahei2010_8410;
    ahei_nomt_84=ahei2010_84-ahei2010_rmtI84;
    ahei_nomt_86=ahei2010_86-ahei2010_rmtI86;
    ahei_nomt_90=ahei2010_90-ahei2010_rmtI90;
    ahei_nomt_94=ahei2010_94-ahei2010_rmtI94; 
    ahei_nomt_98=ahei2010_98-ahei2010_rmtI98; 
    ahei_nomt_02=ahei2010_02-ahei2010_rmtI02; 
    ahei_nomt_06=ahei2010_06-ahei2010_rmtI06; 
    ahei_nomt_10=ahei2010_10-ahei2010_rmtI10;
     
    %include '/udd/hpypl/review/dietscore/ahei80_mc.sas'; 
    data ahei80; set diet80l; ahei_nomt_80=nAHEI80a-meataI80; keep id nAHEI80a ahei_nomt_80 ssb80; proc sort; by id; run;

    %include '/udd/hpypl/review/dietscore/nhs1.exercise.sas'; 

  /*read in physical activity data: include total activity in mets/week and categorical mets*/
  /*categorical mets are 1=<3 total MET-hrs/wk, 2=3 to <9 MET-hours/week, 3=9 to <18 MET-hours/week, 4=18 to <27 MET-hours/week, 5=27 to <42 MET-hours/week, 
    6=42+ MET-hours/week, 7=skipped activity grid, 8=short of page 1-2 or missing*/
  /*passthru or missing is 998 or 999*/
    %act8614 (keep= act86m act88m act92m act94m act96m act98m act00m act04m act08m
	 act12m act14m actc86m actc88m actc92m actc94m actc96m actc98m actc00m actc04m actc08m actc12m actc14m);
	run;

    data nhs_mets; set act8614;
 	array acts{16} act86m act88m act90m act92m act94m act96m act98m act00m act02m act04m act06m act08m act10m act12m act14m act16m;
	do i=1 to 16; if acts{i}>900 then acts{i}=.; end;
	do i=2 to 16; if acts{i}=. then acts{i}=acts{i-1}; end;

	array actc{16} actc86m actc88m actc90m actc92m actc94m actc96m actc98m actc00m actc02m actc04m actc06m actc08m actc10m actc12m actc14m actc16m;
	do i=1 to 16; if actc{i}=7 or actc{i}=8 then actc{i}=.; end;
	do i=2 to 16; if actc{i}=. then actc{i}=actc{i-1}; end;
        proc sort; by id; run;
	run;

    **************************************************
    *         database of BMI   & covariables        *
    **************************************************;

    %der7620 (keep=irt76 irt78 irt80 irt82 irt84 irt86 irt88 irt90 irt92 irt94 irt96 irt98
		irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 irt20 
		yobf mobf  
		bmiage18 bmi76 bmi78 bmi80 bmi82 bmi84 bmi86 bmi88 bmi90 bmi92 bmi94 bmi96 bmi98 bmi00
		bmi02 bmi04 bmi06 bmi08 bmi10 bmi12 bmi14 bmi16 bmi20 
		/*age76 age78 age80 age82 age84 age86 age88 age90 age92 age94 age96 age98 age00
		age02 age04 age06 age08 age10 age12 age14 age16*/ /*age in years at questionnaire return; some are missing so replace with irt and bdt caculated age** in doloop*/
		qt76 qt78 qt80 qt82 qt84 qt86 qt88 qt90 qt92 qt94 qt96 qt98 qt00
		qt02 qt04 qt06 qt08 qt10 qt12 qt14 qt16 qt20/*categories of BMI*/
		can76 can78 can80 can82 can84 can86 can88 can90 can92 can94 can96 can98 can00
		can02 can04 can06 can08 can10 can12 can14 can16 can20/*cancer*/
		hrt76 hrt78 hrt80 hrt82 hrt84 hrt86 hrt88 hrt90 hrt92 hrt94 hrt96 hrt98 hrt00
		hrt02 hrt04 hrt06 hrt08 hrt10 hrt12 hrt14 hrt16 hrt20/*MI or Angina or both*/
		nhor76 nhor78 nhor80 nhor82 nhor84 nhor86 nhor88 nhor90 nhor92 nhor94 nhor96 nhor98
		nhor00 nhor02 nhor04 nhor06 nhor08 nhor10 nhor12 nhor14 nhor16 nhor20/* new PMH use*/
		smkdr76 smkdr78 smkdr80 smkdr82 smkdr84 smkdr86 smkdr88 smkdr90 smkdr92 smkdr94 smkdr96
		smkdr98 smkdr00 smkdr02 smkdr04 smkdr06 smkdr08 smkdr10 smkdr12 smkdr14 smkdr16 smkdr20/*smoking status w amt of cigarettes*/
		race9204 /*race including 1. White 2. Black 3. American Indian 4. Asian 5. Hawaiian*/ white);

	if race9204=1 then white=1;
	else white=0;
        Run;

    %n767880(keep=mmi76 ammi76 fmi76 afmi76
                  ht76 wt76 wt78 wt80 wt18 chol76 chol78 chol80 hbp76  hbp78  hbp80 db76   db78   db80
                  mi76 mi78 mi80 ang76 ang78 ang80 mar80);  
 
  /*read in medication data (only contains aspirin and NSAID use in meds8016*/
    %meds8016 (keep=aspu80 aspu82 aspu84 aspu86 aspu88 aspu90 aspu92 aspu94 aspu96
	       aspu98 aspu00 aspu02 aspu04 aspu06 aspu08 aspu10 aspu12 aspu14 aspu16);
               /*coding is 1 = current user, 0= nonuser, 9=unknown, KD may want to recode this for aspyn=1 else =0*/
               Run;

    %supp8016(keep=mvitu80 mvyn80 viteu80 vite80d mvitu82  mvyn82 viteu82 vite82d mvitu84 mvyn84 viteu84 vite84d
              mvitu86 mvyn86 viteu86 vite86d mvitu88 mvyn88 viteu88 vite88d
              mvitu90 mvyn90 viteu90 vite90d mvitu92 mvyn92 viteu92 vite92d
              mvitu94 mvyn94 viteu94 vite94d mvitu96 mvyn96 viteu96 vite96d 
              mvitu98 mvyn98 viteu98 vite98d mvitu00 mvyn00 viteu00 vite00d
              mvitu02 mvyn02 viteu02 vite02d mvitu04 mvyn04 viteu04 vite04d 
              mvitu06 mvyn06 viteu06 vite06d mvitu08 mvyn08 viteu08 vite08d
              mvitu10 mvyn10 viteu10 vite10d mvitu12 mvyn12 viteu12 vite12d
              mvitu14 mvyn14 viteu14 vite14d mvitu16 mvyn16 viteu16 vite16d);

              if mvitu80=1 then mvyn80=1; else mvyn80=0; if viteu80=1 then vite80d=1; else vite80d=0;
              if mvitu82=1 then mvyn82=1; else if mvitu82=0 then mvyn82=0; else mvyn82=.; if viteu82=1 then vite82d=1; else if viteu82=0 then vite82d=0; else vite80d=.;
              if mvitu84=1 then mvyn84=1; else if mvitu84=0 then mvyn84=0; else mvyn84=.; if viteu84=1 then vite84d=1; else if viteu84=0 then vite84d=0; else vite84d=.;
              if mvitu86=1 then mvyn86=1; else if mvitu86=0 then mvyn86=0; else mvyn86=.; if viteu86=1 then vite86d=1; else if viteu86=0 then vite86d=0; else vite86d=.;
              if mvitu88=1 then mvyn88=1; else if mvitu88=0 then mvyn88=0; else mvyn88=.; if viteu88=1 then vite88d=1; else if viteu88=0 then vite88d=0; else vite88d=.;

              if mvitu90=1 then mvyn90=1; else if mvitu90=0 then mvyn90=0; else mvyn90=.; if viteu90=1 then vite90d=1; else if viteu90=0 then vite90d=0; else vite90d=.;
              if mvitu92=1 then mvyn92=1; else if mvitu92=0 then mvyn92=0; else mvyn92=.; if viteu92=1 then vite92d=1; else if viteu92=0 then vite92d=0; else vite92d=.;
              if mvitu94=1 then mvyn94=1; else if mvitu94=0 then mvyn94=0; else mvyn94=.; if viteu94=1 then vite94d=1; else if viteu94=0 then vite94d=0; else vite94d=.;
              if mvitu96=1 then mvyn96=1; else if mvitu96=0 then mvyn96=0; else mvyn96=.; if viteu96=1 then vite96d=1; else if viteu96=0 then vite96d=0; else vite96d=.;
              if mvitu98=1 then mvyn98=1; else if mvitu98=0 then mvyn98=0; else mvyn98=.; if viteu98=1 then vite98d=1; else if viteu98=0 then vite98d=0; else vite98d=.;
              if mvitu00=1 then mvyn00=1; else if mvitu00=0 then mvyn00=0; else mvyn00=.; if viteu00=1 then vite00d=1; else if viteu00=0 then vite00d=0; else vite00d=.;
              if mvitu02=1 then mvyn02=1; else if mvitu02=0 then mvyn02=0; else mvyn02=.; if viteu02=1 then vite02d=1; else if viteu02=0 then vite02d=0; else vite02d=.;
              if mvitu04=1 then mvyn04=1; else if mvitu04=0 then mvyn04=0; else mvyn04=.; if viteu04=1 then vite04d=1; else if viteu04=0 then vite04d=0; else vite04d=.;
              if mvitu06=1 then mvyn06=1; else if mvitu06=0 then mvyn06=0; else mvyn06=.; if viteu06=1 then vite06d=1; else if viteu06=0 then vite06d=0; else vite06d=.;
              if mvitu08=1 then mvyn08=1; else if mvitu08=0 then mvyn08=0; else mvyn08=.; if viteu08=1 then vite08d=1; else if viteu08=0 then vite08d=0; else vite08d=.;
              if mvitu10=1 then mvyn10=1; else if mvitu10=0 then mvyn10=0; else mvyn10=.; if viteu10=1 then vite10d=1; else if viteu10=0 then vite10d=0; else vite10d=.;
              if mvitu12=1 then mvyn12=1; else if mvitu12=0 then mvyn12=0; else mvyn12=.; if viteu12=1 then vite12d=1; else if viteu12=0 then vite12d=0; else vite12d=.;
              if mvitu14=1 then mvyn14=1; else if mvitu14=0 then mvyn14=0; else mvyn14=.; if viteu14=1 then vite14d=1; else if viteu14=0 then vite14d=0; else vite14d=.;
              if mvitu16=1 then mvyn16=1; else if mvitu16=0 then mvyn16=0; else mvyn16=.; if viteu16=1 then vite16d=1; else if viteu16=0 then vite16d=0; else vite16d=.;
              proc sort; by id; run;

    %nur82(keep=wt82 str82 stry82 ra82 chol82 hbp82 db82 str80 mi82 ang82
                fdb82 sdb82 mdb82 bdb82 dbfh82 mclc82 fclc82 sclc82 bclc82 mbrcn82 sbrcn82 cafh); 
           if str82=1 and stry82 in (1,2,3,4,5,6,7,8) then str80=1; else str80=0; 
           /* stry82   i2     179-180 YEAR STROKE (L15 S6)
               $range 0 - 11 1 = before 1965  2 = 1965-1969  3 = 1970-1975  4 = 1976   5 = 1977 6 = 1978 7 = 1979   8 = 1980
                            9 = 1981  10 = 1982 11 = 1983 */
       if fdb82=1 or mdb82=1 or sdb82=1 or bdb82=1 then dbfh82=1; else dbfh82=0;
       cafh=0; if mclc82=1 or fclc82=1 or sclc82=1 or bclc82=1 or mbrcn82=1 or sbrcn82=1 then cafh=1;
       run;

    %nur84(keep=wt84  aspd84  db84  hbp84  chol84  mi84  ang84  str84  cabg84 cabgy84 ra84   pe84   ucol84);
       /* cabgy84     i2     247-248 year coronary artery surgery (L10a)
                          $range 0 - 11
                           0 = blank         
                           1 = before 1976  
                           2 = 1976        
                           3 = 1977       
                           4 = 1978      
                           5 = 1979     
                           6 = 1980    
                           7 = 1981   
                           8 = 1982 
                           9 = 1983
                          10 = 1984        
                          11 = 1985 */

     data nur84s;
          set nur84;
                   array dis84 {*} db84  hbp84  chol84  mi84  ang84  str84  cabg84  ra84   pe84   ucol84;
                   array dis84s{*} db84s hbp84s chol84s mi84s ang84s str84s cabg84s ra84s  pe84s  ucol84s;
                   do i=1 to dim(dis84);
                      if dis84{i}=2 then dis84s{i}=1;
                      else dis84s{i}=.;
                   end;
                   drop db84  hbp84  chol84  mi84  ang84  str84  cabg84  ra84   pe84   ucol84;
                   run;
              
      proc datasets; delete nur84; run;
      data nur84;
           set nur84s;
                   array dis84 {*} db84  hbp84  chol84  mi84  ang84  str84  cabg84   ra84   pe84   ucol84;
                   array dis84s{*} db84s hbp84s chol84s mi84s ang84s str84s cabg84s  ra84s  pe84s  ucol84s;
                   do i=1 to dim(dis84);
                      dis84{i}=dis84s{i}; 
                   end;
                   drop db84s  hbp84s  chol84s  mi84s  ang84s  str84s  cabg84s   ra84s   pe84s   ucol84s;
                   run;         


    %nur86(keep=wt86 str86 ra86 cabg86 chol86 hbp86 db86 mi86 ang86); 

    %nur88(keep=wt88 physx88 str88 ra88 cabg88 thiaz88 ccblo88 betab88 ace88 bprx88 clrx88 pct588 pct1088
                chol88 hbp88 db88 mi88  ang88
                fdb88 sdb88 mdb88 bdb88 dbfh88);
                if fdb88=1 or mdb88=1 or sdb88=1 or bdb88=1 then dbfh88=1; else dbfh88=0; 

    %nur90(keep=wt90 physx90 str90 ra90 cabg90 chol90 hbp90 db90 mi90 ang90); run;
 
    %nur92(keep=wt92 physx92 str92 ra92 cabg92 alone92 chol92 hbp92 db92 mi92 ang92
                fdb92 sdb92 mdb92 dbfh92 fhbp92 shbp92 mhbp92 hbpfh92
                husbe92 /* Husbands Educational Level; 1.<high school 2.some HS 3.HS grad 
                           4.College grad 5.grad school 6.PT */ husbedu
                rn92 ba92 ma92 dr92 hiedu mdem92 fdem92 sdem92);
                if fdb92=1 or mdb92=1 or sdb92=1 then dbfh92=1;
                if fhbp92=1 or mhbp92=1 or shbp92=1 then hbpfh92=1; else hbpfh92=0;

                if husbe92=6 or husbe92=. then husbedu=.;
                else if husbe92<=3 then husbedu=1;/*high school and lower*/
                else if husbe92=4 then husbedu=2; /*college*/
                else if husbe92=5 then husbedu=3; /*graduate school*/

                hiedu=.;
                 if dr92=1 then hiedu=4;
                else if ma92=1 then hiedu=3;
                else if ba92=1 then hiedu=2;
                else if rn92=1 then hiedu=1;


    %nur94(keep=wt94 physx94 str94 cabg94 thiaz94 lasix94 ccblo94 betab94 bprx94 couma94  digox94 clrx94 chol94 hbp94 db94 mi94 ang94); run;

    %nur96(keep=wt96 physc96 physy96 str96 ra96  cabg96 thiaz96   lasix96    ccblo96    betab96    ace96   bprx96    couma96   digox96   antia96
                antid96 marry96 alone96  chrx96 chol96 hbp96 db96 mi96 ang96
                brmi96 brmid96 smi96 smid96 mifh96 bstr96 bstrd96 sstr96 sstrd96 strfh96);
                if (brmi96=1 and brmid96 in (1,2)) or (smi96=1 and smid96 in (1,2)) then mifh96=1; else mifh96=0;
                if (bstr96=1 and bstrd96 in (1,2)) or (sstr96=1 and sstrd96 in (1,2)) then strfh96=1; else strfh96=0;

    %nur98(keep=wt98 physc98 physy98 str98 cabg98 thiaz98   lasix98    ccblo98    betab98    ace98   bprx98 couma98   digox98   antia98
                   antid98 chrx98 chrx98 chol98 hbp98 db98 mi98 ang98);  

    %nur00(keep=wt00 physc00 physy00 str00 ra00 cabg00 thiaz00   lasix00    ccblo00    betab00    ace00   bprx00   couma00   digox00   antia00
                   marry00 alone00  chrx00 chrxd00 przc00 zol00 paxil00 celex00 antid00
                   chol00 hbp00 db00 mi00 ang00 depr00 deprd00);   

    %nur02(keep=wt02 physc02 physy02 ra02 sleep02 sleep02c snore02 snore02c str02 renf02  renfd02 iron02d
                cabg02 przc02 zol02 paxil02 celex02 antid02
                bprx02 lasix02 thiaz02 ace02 ccblo02 betab02
                couma02   digox02   antia02 chrx02 oclrx02 chol02 hbp02 db02 mi02 ang02 depr02 ohypo02 hypo02 hypod02 hypo00
                alz02 alzd02 q02);
                                                         /*renfd02 Chronic Renal Failure Dt Dx (L18) $label 1.96 or before, 2.97-99, 3.2000, 4.2001,*/
            if iron02d ne 1 then iron02d=0;
            if sleep02 in (1,2) then sleep02c=5;
               else if sleep02=3 then sleep02c=6;
               else if sleep02=4 then sleep02c=7;
               else if sleep02=5 then sleep02c=8;
               else if sleep02 in (6,7) then sleep02c=9;
               else sleep02c=.;
             if snore02 in (1,2,3,4,5) then snore02c=5-snore02; else snore02c=.;
             if hypo02=1 and hypod02 in (1,2,3) then hypo00=1; else hypo00=0;
             proc freq; tables alz02;
             run;
 
    %nur04(keep=wt04 physc04 physy04 str04 strd04 ra04 rad04 renf04 renfd04 ang04 angd04
                cabg04 cabgd04 antid04 ssri04 bprex04 lasix04 thiaz04 ace04 ccblo04 betab04 k04 iron04 couma04 digox04 antia04
                marry04 alone04 mev04 zoc04 crest04 prav04 lip04 lesc04 oclrx04 chol04 hbp04 db04 mi04
                alz04 alzd04 depr04 ohypo04 mdem04 mdemd04 fdem04 fdemd04 sdem04 sdemd04 q04); 

    %nur06(keep=wt06 physc06 physy06 str06 strd06 ra06 rad06 ang06 angd06 cabg06
                cabgd06  antid06 ssri06 bprx06 lasix06 thiaz06 ace06 ccblo06 betab06 k06 iron06d couma06 digox06 antia06
                mev06 zoc06 crest06 prav06 lip06 ostat06 oclrx06 chol06 hbp06 db06 mi06
                alz06 alzd06 depr06 ohypo06 q06);

    %nur08(keep=physc08 physy08 str08 strd08 ra08 rad08 yr08 mo08 wt08
                db08 chol08 hbp08 mi08 mid08 ang08 angd08 str08 strd08 alone08 marry08
                cabg08 cabgd08 antid08 ssri08 bprx08 lasix08 thiaz08 ace08 ccblo08 couma08 digox08 antia08 betab08 k08
                stat08 mev08 zoc08 crest08 prav08 lip08 ostat08 oclrx08
                brmi08 brmid08 smi08 smid08 omi08 omid08 mifh08 mstr08 mstrd08 fstr08 fstrd08 sstr08 sstrd08 strfh08
                alz08 alzd08 depr08 ohypo08 q08);
                if (brmi08=1 and brmid08 in (1,2)) or (smi08=1 and smid08 in (1,2)) or (omi08=1 and omid08 in (1,2)) then mifh08=1; else mifh08=0;
                if (mstr08=1 and mstrd08 in (1,2)) or (fstr08=1 and fstrd08 in (1,2)) or (sstr08=1 and sstrd08 in (1,2)) then strfh08=1; else strfh08=0;
                if mo08=. then mo08=6;
 
    %nur10(keep=wt10 ra10 rad10 chol10 hbp10 db10 mi10	mid10	mih10
	        ang10	angd10	angc10	cabg10	cabgd10	chf10	chfd10	str10	strd10 antid10 ssri10
                nasp10	ibu10	ibud10	nibu10	celeb10	celed10	ibut10	thiaz10	lasix10
                k10	ccblo10	betab10	ace10	angio10	bprx10	couma10	plavi10	digox10
                antia10	stat10	mev10	prav10	zoc10	lip10	crest10	ostat10	statinpt10	oclrx10
                physx10	physc10	physy10	physpt10 mo10 yr10 
                alz10 alzd10 depr10 ohypo10 q10);
 
   %nur12(keep= wt12 ra12 rad12 chol12 hbp12 db12 mi12	mid12	mih12
	        ang12	angd12	angc12	cabg12	cabgd12	chf12	chfd12	str12	strd12 antid12 ssri12 marry12 alone12
                nasp12	ibu12	ibud12	nibu12	celeb12	celed12	ibut12	thiaz12	lasix12
                k12	ccblo12	betab12	ace12	angio12	bprx12	couma12	plavi12	digox12
                antia12	mev12	prav12	zoc12	lip12	crest12	ostat12	oclrx12
                physx12	physc12	physy12	physpt12 mo12 yr12
                alz12 alzd12 dem12  demd12 depr12 hypo12 hypod12 ohypo12 q12); /*no vitamin data available this cycle*/

   %nur14(keep= wt14 ra14 rad14 chol14 hbp14 db14 mi14	mid14	mih14
	        ang14	angd14	angc14	cabg14	cabgd14	chf14	chfd14	str14	strd14 antid14 ssri14
                nasp14	ibu14	ibud14	nibu14	celeb14	celed14	ibut14	thiaz14	lasix14
                k14	ccblo14	betab14	ace14	angio14	bprx14	couma14	plavi14	digox14
                antia14	mev14	prav14	zoc14	lip14	crest14	ostat14	oclrx14
                physx14	physc14	physy14	physpt14 mo14 yr14
                alz14 alzd14 depr14 ohypo14 park14 q14);

  %nur16(keep= wt16 ra16 rad16 chol16 hbp16 db16 mi16	mid16	mih16 marry16 alone16
	        ang16	angd16	angc16	cabg16	cabgd16	chf16	chfd16	str16	strd16 antid16 ssri16
                nasp16	ibu16	ibud16	nibu16	celeb16	celed16	ibut16	thiaz16	lasix16
                k16	ccblo16	betab16	ace16	angio16	bprx16	couma16	plavi16	digox16
                antia16	mev16	prav16	zoc16	lip16	crest16	ostat16	oclrx16
                physx16	physc16	physy16	physpt16 mo16 yr16
                alz16 alzd16 depr16 ohypo16 park16 q16);
     run; 


   %nur20(keep= alzdem20 alzdemd20 
                wt20    chol20      hbp20       db20 mi20	mid20	mih20 marry20 alone20
	    ang20	angd20	angc20	cabg20	cabgd20	chf20	chfd20	str20	strd20 antida20 
                nasp20	ibu20	ibud20	nibu20	ibut20		 
                thiaz20	ccblo20	betab20	ace20	angio20	bprx20	couma20	pradx20     plavi20  prasu20	digox20
                antia20	stat20      oclrx20
                physx20	physc20	physy20	physpt20  
                depr20  
                insul20 ninsu20     metf20	jard20      invk20 farx20 sita20 ohypo20    
                park20      rasle20 q20);  

      /**** neighborhood SES - nSES****/ 
      libname curr "data/covar/";
      data nses8616; set curr.anses8616;  
      proc sort nodupkey; by id; run;

      /**** PD case in NHS between 1976-2012 ****/
      libname hppdcase '/proj/hppars/hppar0g/case/';
      data nhspdcase; set hppdcase.pd12nur; confpd=conf; s1pd=s1; drop conf s1; 
      proc sort nodupkey; by id; run;

 
 
                              ********************
                              *      database    *
                              ********************;

    data alladrd;
    merge der7620 (in=der)  apoedata foodsnhs ahei2010 ahei80 amed8020
          n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 nur06 nur08 nur10 nur12 nur14 nur16 nur20
          meds8016 disease nhs_mets  exernhs supp8016 nses8616 nhspdcase nhsPRS nhsPRS2
          memorynhs;
     by id;
     exrec=1;
     if der and first.id then exrec=0;  

           height=ht76*25.4/1000;
           if ht76=0 then height=.; 

           if mobf<=0 or mobf>12 then mobf=6;   
               ** birthday in months;   
               bdt=12*yobf+mobf;
               birthday=bdt;

        /* family history */

           if (mmi76=1 and ammi76<=60)  then mothmi=1;
           if (fmi76=1 and afmi76<=60)       then fathmi=1;
           mifh76=0;
           if mothmi=1 or fathmi=1 then mifh76=1;
           if mifh76=1 or mifh96 or mifh08=1 then mifh=1; else mifh=0; /*** note: family history of MI before 60 years, ref: Fung T AJCN 2009paper***/

           if dbfh82=1 or dbfh88=1 or dbfh92=1 then dbfh=1; else dbfh=0;
           if fdb82=1 or mdb82=1 or fdb88=1 or mdb88=1 or fdb92=1 or mdb92=1 then parentaldm=1; else parentaldm=0;
           if mdb82=1 or mdb88=1 or mdb92=1 then mdm=1; else  mdm=0;
           if fdb82=1 or fdb88=1 or fdb92=1 then fdm=1; else fdm=0;

           if hbpfh92=1 then hbpfh=1; else hbpfh=0;

           if strfh96=1 or strfh08=1 then strfh=1; else strfh=0;

          /*** antidepressant medication ***/
           antidep76=0;
           if antid96=1 then antidep96=1; else antidep96=0;
           if antid98=1 then antidep98=1; else antidep98=0;
           antidep00=0;if przc00=1 or zol00=1 or paxil00=1 or celex00=1 or antid00=1 then antidep00=1;
           antidep02=0; if przc02=1 or zol02=1 or paxil02=1 or celex02=1 or antid02=1 then antidep02=1;
           antidep04=0; if antid04=1 or ssri04=1 then antidep04=1;
           antidep06=0; if antid06=1 or ssri06=1 then antidep06=1;
           antidep08=0; if antid08=1 or ssri08=1 then antidep08=1;
           antidep10=0; if antid10=1 or ssri10=1 then antidep10=1;
           antidep12=0; if antid12=1 or ssri12=1 then antidep12=1;
           antidep14=0; if antid14=1 or ssri14=1 then antidep14=1;
           antidep16=0; if antid16=1 or ssri16=1 then antidep16=1;
           antidep20=0; if antida20=1 then antidep20=1;

          /***antihypertensive***/
           htnrx76=0;
           if (thiaz88=1 or ccblo88=1 or betab88=1 or ace88=1 or bprx88=1) then htnrx88=1;   else htnrx88=0;
           if (thiaz94=1 or lasix94 =1 or ccblo94 =1 or betab94 =1 or bprx94=1 ) then htnrx94=1;   else htnrx94=0;
           if (thiaz96=1 or lasix96 =1 or ccblo96 =1 or betab96 =1 or ace96=1 or bprx96=1 ) then htnrx96=1;   else htnrx96=0;
           if (thiaz98=1 or lasix98 =1 or ccblo98 =1 or betab98 =1 or ace98=1 or bprx98=1 ) then htnrx98=1;   else htnrx98=0;
           if (thiaz00=1 or lasix00 =1 or ccblo00 =1 or betab00 =1 or ace00=1 or bprx00=1 ) then htnrx00=1;   else htnrx00=0;
           if (thiaz02=1 or lasix02 =1 or ccblo02 =1 or betab02 =1 or ace02=1 or bprx02=1 ) then htnrx02=1;   else htnrx02=0;
           if (bprex04=1 or lasix04=1 or  thiaz04=1 or ace04=1 or ccblo04=1 or  betab04=1 or  k04=1) then htnrx04=1;   else htnrx04=0;
           if (bprx06=1 or lasix06=1 or thiaz06=1 or ace06=1 or ccblo06=1 or betab06=1 or k06=1) then htnrx06=1;   else htnrx06=0;
           if (bprx08=1 or lasix08=1 or thiaz08=1 or ace08=1 or ccblo08=1 or betab08=1 or k08=1) then htnrx08=1;   else htnrx08=0;
           if (bprx10=1 or lasix10=1 or thiaz10=1 or ace10=1 or ccblo10=1 or betab10=1 or k10=1) then htnrx10=1;   else htnrx10=0;
           if (bprx12=1 or lasix12=1 or thiaz12=1 or ace12=1 or ccblo12=1 or betab12=1 or k12=1) then htnrx12=1;   else htnrx12=0;
           if (bprx14=1 or lasix14=1 or thiaz14=1 or ace14=1 or ccblo14=1 or betab14=1 or k14=1) then htnrx14=1;   else htnrx14=0;
           if (bprx16=1 or lasix16=1 or thiaz16=1 or ace16=1 or ccblo16=1 or betab16=1 or k16=1) then htnrx16=1;   else htnrx16=0;
           if (thiaz20=1 or ace20=1 or ccblo20=1 or betab20=1 or angio20=1 or bprx20=1) then htnrx20=1;   else htnrx20=0;

          /***heart diseases medications ***/
           if (couma94=1 or digox94=1 or betab94 =1 or ccblo94 =1) then chdrx94=1; else chdrx94=0;/***digox94 inclde digoxin or antiarrhythmic***/
           if (couma96=1 or digox96=1 or antia96=1 or betab96 =1 or ccblo96 =1) then chdrx96=1; else chdrx96=0;
           if (couma98=1 or digox98=1 or antia98=1 or betab98 =1 or ccblo98 =1) then chdrx98=1; else chdrx98=0;
           if (couma00=1 or digox00=1 or antia00=1 or betab00 =1 or ccblo00 =1) then chdrx00=1; else chdrx00=0;
           if (couma02=1 or digox02=1 or antia02=1 or betab02 =1 or ccblo02 =1) then chdrx02=1; else chdrx02=0;

	  if (couma04=1 or digox04=1 or antia04=1 or betab04 =1 or ccblo04 =1) then chdrx04=1; else chdrx04=0;
  	  if (couma06=1 or digox06=1 or antia06=1 or betab06 =1 or ccblo06 =1) then chdrx06=1; else chdrx06=0;
           if (couma08=1 or digox08=1 or antia08=1 or betab08 =1 or ccblo08 =1) then chdrx08=1; else chdrx08=0;
	  if (couma10=1 or digox10=1 or antia10=1 or betab10 =1 or ccblo10 =1) then chdrx10=1; else chdrx10=0;
           if (couma12=1 or digox12=1 or antia12=1 or betab12 =1 or ccblo12 =1) then chdrx12=1; else chdrx12=0;
	  if (couma14=1 or digox14=1 or antia14=1 or betab14 =1 or ccblo14 =1) then chdrx14=1; else chdrx14=0;
           if (couma16=1 or digox16=1 or antia16=1 or betab16 =1 or ccblo16 =1) then chdrx16=1; else chdrx16=0;

          /***antihigh cholesterol drug***/
           hchtx76=0;
           if clrx88=1 then hchtx88=1; else hchtx88=0;
           if clrx94=1 then hchtx94=1; else hchtx94=0;
           if chrx96=1 then hchtx96=1; else hchtx96=0;
           if chrx98=1 then hchtx98=1; else hchtx98=0;
           if chrx00=1 or chrxd00=1 then hchtx00=1; else hchtx00=0;
           if chrx02=1 or oclrx02=1 then hchtx02=1; else hchtx02=0;
           if  mev04=1 or zoc04=1 or crest04=1 or prav04=1 or lip04=1 or lesc04=1 or oclrx04=1 then hchtx04=1; else hchtx04=0;
           if  mev06=1 or zoc06=1 or crest06=1 or prav06=1 or lip06=1 or ostat06=1 or oclrx06=1 then hchtx06=1; else hchtx06=0;
           if  mev08=1 or zoc08=1 or crest08=1 or prav08=1 or lip08=1 or ostat08=1 or oclrx08=1 then hchtx08=1; else hchtx08=0;
           if  mev10=1 or zoc10=1 or crest10=1 or prav10=1 or lip10=1 or ostat10=1 or oclrx10=1 then hchtx10=1; else hchtx10=0;
           if  mev12=1 or zoc12=1 or crest12=1 or prav12=1 or lip12=1 or ostat12=1 or oclrx12=1 then hchtx12=1; else hchtx12=0;
           if  mev14=1 or zoc14=1 or crest14=1 or prav14=1 or lip14=1 or ostat14=1 or oclrx14=1 then hchtx14=1; else hchtx14=0;
           if  mev16=1 or zoc16=1 or crest16=1 or prav16=1 or lip16=1 or ostat16=1 or oclrx16=1 then hchtx16=1; else hchtx16=0;
           if  stat20=1 or oclrx20=1 then hchtx20=1; else hchtx20=0;

	 /***marry and living alone***/
           if alone92 ne 1 then alone92=0;
           if alone96 ne 1 then alone96=0;
           if alone00 ne 1 then alone00=0;
           if alone04 ne 1 then alone04=0;
           if alone08 ne 1 then alone08=0;
           if alone12 ne 1 then alone12=0;
           if alone16 ne 1 then alone16=0;
           if alone20 ne 1 then alone20=0;

           if mar80=1   then mar80=1; else mar80=0;
           if marry96=1 then mar96=1; else mar96=0;
           if marry00=1 then mar00=1; else mar00=0;
           if marry04=1 then mar04=1; else mar04=0;
           if marry08=1 then mar08=1; else mar08=0;
           if marry12=1 then mar12=1; else mar12=0;
           if marry16=1 then mar16=1; else mar16=0;
           if marry20=1 then mar20=1; else mar20=0;

        /** physical exam **/
	   pexam88=0; if physx88=2 or physx88=3 then pexam88=1;
	   pexam90=0; if physx90=2 or physx90=3 then pexam90=1;
	   pexam92=0; if physx92=2 or physx92=3 then pexam92=1;
	   pexam94=0; if physx94=2 then pexam94=1;
	   pexam96=0; if physc96=1 or physy96=1 then pexam96=1;
	   pexam98=0; if physc98=1 or physy98=1 then pexam98=1;
	   pexam00=0; if physc00=1 or physy00=1 then pexam00=1;
	   pexam02=0; if physc02=1 or physy02=1 then pexam02=1;
	   pexam04=0; if physc04=1 or physy04=1 then pexam04=1;
	   pexam06=0; if physc06=1 or physy06=1 then pexam06=1;
	   pexam08=0; if physc08=1 or physy08=1 then pexam08=1;
	   pexam10=0; if physc10=1 or physy10=1 then pexam10=1;
	   pexam12=0; if physc12=1 or physy12=1 then pexam12=1;
	   pexam14=0; if physc14=1 or physy14=1 then pexam14=1;
	   pexam16=0; if physc16=1 or physy16=1 then pexam16=1;

       /** physical exam for screening purposes **/
	   psexam88=0; if physx88=3 then psexam88=1;
	   psexam90=0; if physx90=2 then psexam90=1;
	   psexam92=0; if physx92=2 then psexam92=1;
	   psexam96=0; if physc96=1 then psexam96=1;
	   psexam98=0; if physc98=1 then psexam98=1;
	   psexam00=0; if physc00=1 then psexam00=1;
	   psexam02=0; if physc02=1 then psexam02=1;
	   psexam04=0; if physc04=1 then psexam04=1;
	   psexam06=0; if physc06=1 then psexam06=1;
	   psexam08=0; if physc08=1 then psexam08=1;
	   psexam10=0; if physc10=1 then psexam10=1;
	   psexam12=0; if physc12=1 then psexam12=1;
	   psexam14=0; if physc14=1 then psexam14=1;
	   psexam16=0; if physc16=1 then psexam16=1;

      /****** diagnosis time of AD ******/
                                                                  /* alzdemd20 $label 1.<Jun2017;	2.Jun2017-May2019;3.>=Jun2019;	4.pt */
           if alzdem20=1 then do;
                                   if alzdemd20=1 then dtdx_alzd=1410;
                              else if alzdemd20=2 then dtdx_alzd=1422;
                              else if alzdemd20=3 then dtdx_alzd=1434;
                              else                     dtdx_alzd=1422;
                              end;

                                                                  /* alzd16 $label 1.<Jun2014;	2.Jun2014-May2016;3.>=Jun2016;	4.pt */
           if alz16=1 then do;
                                   if alzd16=1 then dtdx_alzd=1374;
                              else if alzd16=2 then dtdx_alzd=1386;
                              else if alzd16=3 then dtdx_alzd=1398;
                              else                  dtdx_alzd=1386;
                              end;
                                                                  /* alzd14 $label 1.<6/12; 2.6/12-5/14; 3.>6/14; 4.pt */                  
           if alz14=1 then do;
                                   if alzd14=1 then dtdx_alzd=1350;
                              else if alzd14=2 then dtdx_alzd=1362;
                              else if alzd14=3 then dtdx_alzd=1374;
                              else                  dtdx_alzd=1362;
                            end;
                                                                   /* alz12 $label 1.2001 or before; 2.2002-2005; 3.2006-2009;4.2010-2011; 5.2012+; 6.pt */
           if dem12=1 then do;                             
                                   if demd12=1 then dtdx_alzd=1212;  
                              else if demd12=2 then dtdx_alzd=1248;
                              else if demd12=3 then dtdx_alzd=1296;
                              else if demd12=4 then dtdx_alzd=1332; 
                              else if demd12=5 then dtdx_alzd=1350;
                              else                  dtdx_alzd=1338;
                            end;
                                                                   /* alz12 $label 1.2001 or before; 2.2002-2005; 3.2006-2009;4.2010-2011; 5.2012+; 6.pt */
           if alz12=1 then do;                             
                                   if alzd12=1 then dtdx_alzd=1212;  
                              else if alzd12=2 then dtdx_alzd=1248;
                              else if alzd12=3 then dtdx_alzd=1296;
                              else if alzd12=4 then dtdx_alzd=1332; 
                              else if alzd12=5 then dtdx_alzd=1350;
                              else                  dtdx_alzd=1338;
                            end;
                                                                   /* alzd10	$label 1.<6/08; 2.6/08-5/10; 3.>6/10; 4.pt */
           if alz10=1 then do;                             
                                   if alzd10=1 then dtdx_alzd=1302;
                              else if alzd10=2 then dtdx_alzd=1314;
                              else if alzd10=3 then dtdx_alzd=1326;
                              else                  dtdx_alzd=1314;
                            end;
                                                                   /* alzd08	$label 1.<6/06; 2.6/06-5/08; 3.>6/08; 4.pt */
           if alz08=1 then do;                             
                                   if alzd08=1 then dtdx_alzd=1278;
                              else if alzd08=2 then dtdx_alzd=1290;
                              else if alzd08=3 then dtdx_alzd=1302;
                              else                  dtdx_alzd=1290;
                            end;
                                                                   /* alzd06	$label 1.<6/04; 2.6/04-5/06; 3.>6/06; 4.pt  */
           if alz06=1 then do;                             
                                   if alzd06=1 then dtdx_alzd=1254;
                              else if alzd06=2 then dtdx_alzd=1266;
                              else if alzd06=3 then dtdx_alzd=1378;
                              else                  dtdx_alzd=1266;
                            end;
                                                                   /* alzd04	$label 1.96 or before; 2.97-01; 3.2002; 4.2003;\	5.2004; 6.pt   */
           if alz04=1 then do;                             
                                   if alzd04=1 then dtdx_alzd=1152;  
                              else if alzd04=2 then dtdx_alzd=1194;
                              else if alzd04=3 then dtdx_alzd=1231; * setup at July of the year for 3,4,5 whole year answers;
                              else if alzd04=4 then dtdx_alzd=1243; 
                              else if alzd04=5 then dtdx_alzd=1255;
                              else                  dtdx_alzd=1243;
                            end;
                                                                   /* alzd02	$label 1.before 6/1/00; 2.6/00-5/02;\3.after 6/1/02; 4.pt date unknown    */
           if alz02=1 then do;                             
                                   if alzd02=1 then dtdx_alzd=1205;  
                              else if alzd02=2 then dtdx_alzd=1218;
                              else if alzd02=3 then dtdx_alzd=1243; 
                              else                  dtdx_alzd=1218;
                            end;

           if alz02=1 or alz04=1 or alz06=1 or alz08=1 or alz10=1 or alz12=1 or dem12=1 or alz14=1 or alz16=1 or alzdem20=1 then alz0220=1; else alz0220=0;

      /*** mortality due to dementia ***/     
          if 0<dtdth<9999 and (newicda3 = 290 OR newicda3 = 331) then dementiadeath=1;  else dementiadeath=0;
      
          if alz0220=1 and dementiadeath=1 and 0<dtdth<dtdx_alzd<9999 then dtdth=dtdx_alzd; /*make sure death did not occur before date of diagnosis*/

      /*** diagnosis time of AD ***/
          if alz0220=1 and dementiadeath=1 then dtdx_ad=min(dtdth,dtdx_alzd);
          else if alz0220=1 then dtdx_ad=dtdx_alzd;
          else if dementiadeath=1 then dtdx_ad=dtdth;
          else dtdx_ad=.;

          if dtdx_ad ne .  and dtdx_ad=dtdx_alzd then nfad=1; else nfad=0;

          if dtdx_ad ne . then prevalent_ad=1; else prevalent_ad=0;

       /*** age of AD diagnosis ***/
          if dtdx_ad ne . and bdt ne . then ADdiagage=int((dtdx_ad-bdt)/12);      
            if ADdiagage =< 0 then ADdiagage=.;


       /*** family history of dementia ***/
          if mdem92=1 or fdem92=1 or sdem92=1 or  mdem04=2 or fdem04=2 or sdem04=2 then fhdem=1; else fhdem=0;

       /*** Parkinson disease - confirmed 1976-2012  or self-reported 2014 2016 ***/ 
          if confpd=11 & s1pd in (1,2) or park14=1 or park16=1 then pdcase=1; else pdcase=0;


/***** Calculate cumulative averages of aMed *****/

array amedorg{10} amed_80   amed_84   amed_86   amed_90   amed_94   amed_98   amed_02   amed_06   amed_10   amed_20; 
array emedcum{10} emed80v  emed84v  emed86v  emed90v  emed94v  emed98v  emed02v  emed06v  emed10v  emed20v; 
   sumvar=0;    n=0;
   do i=1 to 10; 
      if amedorg{i} ne . then do;
         n=n+1; sumvar=sumvar+amedorg{i};
         end; 
      if n ne 0 then emedcum{i}=sumvar/n; else emedcum{i}=.;
      end;

      * calculate weight changes before missing carry forward and readi for missing weight changes carry forward;
      if wt78=. or wt78=0 then wt78=wt76; 
      array wt2      {20} wt80      wt82      wt84      wt86       wt88      wt90      wt92      wt94      wt96      wt98      wt00      wt02      wt04      wt06      wt08      wt10      wt12      wt14      wt16      wt20     ;
      array wt1      {20} wt78      wt80      wt82      wt84       wt86      wt88      wt90      wt92      wt94      wt96      wt98      wt00      wt02      wt04      wt06      wt08      wt10      wt12      wt14      wt16     ;
      array wt21     {20} wtchg80   wtchg82   wtchg84   wtchg86    wtchg88   wtchg90   wtchg92   wtchg94   wtchg96   wtchg98   wtchg00   wtchg02   wtchg04   wtchg06   wtchg08   wtchg10   wtchg12   wtchg14   wtchg16   wtchg20s ;

      do i=1 to 20; 
         if wt1{i}=0 then wt1{i}=.; 
         if wt2{i}=0 then wt2{i}=.; 
         wt21{i}=wt2{i}-wt1{i};
      end;

      wtchg20=wtchg20s/2; *2 years weight changes in lbs;

      array hbp22     {22} hbp76 hbp78 hbp80  hbp82     hbp84      hbp86     hbp88     hbp90     hbp92     hbp94     hbp96     hbp98     hbp00     hbp02     hbp04     hbp06     hbp08     hbp10     hbp12     hbp14     hbp16  hbp20;
      do i=1 to 22; if hbp22{i}=. then hbp22{i}=0; end;

proc means n nmiss min mean std P25 median P75 max;
var wtchg80   wtchg82   wtchg84   wtchg86    wtchg88   wtchg90   wtchg92   wtchg94   wtchg96   wtchg98   wtchg00   wtchg02   wtchg04   wtchg06   wtchg08   wtchg10   wtchg12   wtchg14   wtchg16   wtchg20;      
run; 

proc freq;
tables hbp76 hbp78 hbp80  hbp82     hbp84      hbp86     hbp88     hbp90     hbp92     hbp94     hbp96     hbp98     hbp00     hbp02     hbp04     hbp06     hbp08     hbp10     hbp12     hbp14     hbp16 hbp20;
run;



   /******* get rid of files that are no longer needed ********************************************/
   proc datasets;
        delete der7620 apoedata foodsnhs ahei2010 ahei80   
               n767880 nur82 nur84 nur86 nur88 nur90 nur92 nur94 nur96 nur98 nur00 nur02 nur04 nur06 nur08 nur10 nur12 nur14 nur16 nur20
               meds8016 disease nhs_mets  exernhs supp8016 nses02 nhspdcase nhsPRS nhsPRS2;
   run;

 /*proc means n nmiss min mean std max nolabels;  run;*/

                            ***
                        ***     ***
                    ***            ***
                ***    data-NHS       ***
                    ***            ***
                        ***     ***
                            ***;

data pre_pm;
     set alladrd  end=_end_;

     *diseases at baseline 1980;
      if db76=1 or db78=1 or db80=1 then db_7680=1; else db_7680=0;      
      if hbp76=1 or hbp78=1 or hbp80=1 then hbp_7680=1; else hbp_7680=0;
      if chol76=1 or chol78=1 or chol80=1 then chol_7680=1;else chol_7680=0; 
      if  mi76=1 or mi78=1 or mi80=1 then mi7680=1; else mi7680=0;
      if  can76=1 or can78=1 or can80=1 then can7680=1; else can7680=0; 
      depr80=0; antidep76=0; htnrx76=0; hypo80=0;alone80=0;
      
      irt18=1422; *add 2018 cycle in order to match hp for pooling data;
      
      /* Baseline covariates */
        
      array irt     {22} irt80     irt82     irt84      irt86     irt88     irt90     irt92     irt94     irt96     irt98     irt00     irt02     irt04     irt06     irt08     irt10     irt12     irt14     irt16     irt18     irt20   cutoff;
      array tvar    {21} t80       t82       t84        t86       t88       t90       t92       t94       t96       t98       t00       t02       t04       t06       t08       t10       t12       t14       t16       t18       t20      ;
      array age     {21} age80     age82     age84      age86     age88     age90     age92     age94     age96     age98     age00     age02     age04     age06     age08     age10     age12     age14     age16     age18     age20    ;
      array wt      {21} wt80      wt82      wt84       wt86      wt88      wt90      wt92      wt94      wt96      wt98      wt00      wt02      wt04      wt06      wt08      wt10      wt12      wt14      wt16      wt16      wt20     ;
      array bmi     {21} bmi80     bmi80     bmi80      bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80     bmi80    ; 
      array actm    {21} th80      th80      th80       th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80      th80     ;
      array nsmk    {21} smkdr80   smkdr80   smkdr80    smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80   smkdr80  ;
      array smm     {21} smm80     smm80     smm80      smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80     smm80    ;
      array hor     {21} nhor80    nhor80    nhor80     nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80    nhor80   ; 
      array mvyn    {21} mvyn80    mvyn82    mvyn84     mvyn86    mvyn88    mvyn90    mvyn92    mvyn94    mvyn96    mvyn98    mvyn00    mvyn02    mvyn04    mvyn06    mvyn08    mvyn10    mvyn12    mvyn14    mvyn16    mvyn16    mvyn16   ;
      array aspu    {21} aspu80    aspu82    aspu84     aspu84    aspu88    aspu90    aspu92    aspu94    aspu96    aspu98    aspu00    aspu02    aspu04    aspu06    aspu08    aspu10    aspu12    aspu14    aspu16    aspu16    aspu16   ; 
      array alone   {21} alone80   alone80   alone80    alone80   alone80   alone80   alone92   alone92   alone96   alone96   alone00   alone00   alone04   alone04   alone08   alone08   alone12   alone12   alone16   alone16   alone20  ;       
      array mar     {21} mar80     mar80     mar80      mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80     mar80    ;

      array diab    {21} db_7680   db_7680      db_7680       db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680      db_7680     ;
      array chol    {21} chol_7680 chol82    chol84     chol86    chol88    chol90    chol92    chol94    chol96    chol98    chol00    chol02    chol04    chol06    chol08    chol10    chol12    chol14    chol16    chol16    chol20   ;
      array rphbp   {21} hbp_7680  hbp_7680     hbp_7680      hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680     hbp_7680    ;
      array antid   {21} antidep76 antidep76 antidep76  antidep76 antidep76 antidep76 antidep76 antidep76 antidep96 antidep98 antidep00 antidep02 antidep04 antidep06 antidep08 antidep10 antidep12 antidep14 antidep16 antidep16 antidep20; 
      array antihp  {21} htnrx76   htnrx76   htnrx76    htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76   htnrx76  ; 
      array antitc  {21} hchtx76   hchtx76   hchtx76    hchtx76   hchtx88   hchtx88   hchtx88   hchtx94   hchtx96   hchtx98   hchtx00   hchtx02   hchtx04   hchtx06   hchtx08   hchtx10   hchtx12   hchtx14   hchtx16   hchtx16   hchtx20  ;  
      array repmi   {21} mi7680    mi82      mi84       mi86      mi88      mi90      mi92      mi94      mi96      mi98      mi00      mi02      mi04      mi06      mi08      mi10      mi12      mi14      mi16      mi16      mi20     ;
      array repstrk {21} str80     str82     str84      str86     str88     str90     str92     str94     str96     str98     str00     str02     str04     str06     str08     str10     str12     str14     str16     str16     str20    ;          
      array cance   {21} can7680   can82     can84      can86     can88     can90     can92     can94     can96     can98     can00     can02     can04     can06     can08     can10     can12     can14     can16     can16     can20    ;
      array depr    {21} depr80    depr80    depr80     depr80    depr80    depr80    depr80    depr80    depr80    depr80    depr00    depr02    depr04    depr06    depr08    depr10    depr12    depr14    depr16    depr16    depr20   ;
      array hypo    {21} hypo80    hypo80    hypo80     hypo80    hypo80    hypo80    hypo80    hypo80    hypo80    hypo80    hypo00    hypo02    hypo02    hypo02    hypo02    hypo02    hypo12    hypo12    hypo12    hypo12    hypo12   ;

      array trmeatcon{21} trmeat80d  trmeat80d  trmeat84d  trmeat86d  trmeat86d   trmeat90d  trmeat90d   trmeat94d  trmeat94d  trmeat98d  trmeat98d  trmeat02d  trmeat02d  trmeat06d  trmeat06d  trmeat10d  trmeat10d  trmeat10d  trmeat10d trmeat10d trmeat20d; 
      array prmeatcon{21} prmeat80d  prmeat80d  prmeat84d  prmeat86d  prmeat86d   prmeat90d  prmeat90d   prmeat94d  prmeat94d  prmeat98d  prmeat98d  prmeat02d  prmeat02d  prmeat06d  prmeat06d  prmeat10d  prmeat10d  prmeat10d  prmeat10d prmeat10d prmeat20d; 
      array urmeatcon{21} urmeat80d  urmeat80d  urmeat84d  urmeat86d  urmeat86d   urmeat90d  urmeat90d   urmeat94d  urmeat94d  urmeat98d  urmeat98d  urmeat02d  urmeat02d  urmeat06d  urmeat06d  urmeat10d  urmeat10d  urmeat10d  urmeat10d urmeat10d urmeat20d; 
      array calor    {21} calor80n   calor80n   calor84n   calor86n   calor86n    calor90n   calor90n    calor94n   calor94n   calor98n   calor98n   calor02n   calor02n   calor06n   calor06n   calor10n   calor10n   calor10n   calor10n  calor10n  calor20n;
      array alco     {21} alco80n    alco80n    alco84n    alco86n    alco86n     alco90n    alco90n     alco94n    alco94n    alco98n    alco98n    alco02n    alco02n    alco06n    alco06n    alco10n    alco10n    alco10n    alco10n   alco10n   alco20n;

      array fruitscon{21} fruits80d  fruits80d  fruits84d  fruits86d  fruits86d   fruits90d  fruits90d   fruits94d  fruits94d  fruits98d  fruits98d  fruits02d  fruits02d  fruits06d  fruits06d  fruits10d  fruits10d  fruits10d  fruits10d fruits10d fruits20d; 
      array vegealcon{21} vegeal80d  vegeal80d  vegeal84d  vegeal86d  vegeal86d   vegeal90d  vegeal90d   vegeal94d  vegeal94d  vegeal98d  vegeal98d  vegeal02d  vegeal02d  vegeal06d  vegeal06d  vegeal10d  vegeal10d  vegeal10d  vegeal10d vegeal10d vegeal20d; 
      array whgrnscon{21} whgrns80d  whgrns80d  whgrns84d  whgrns86d  whgrns86d   whgrns90d  whgrns90d   whgrns94d  whgrns94d  whgrns98d  whgrns98d  whgrns02d  whgrns02d  whgrns06d  whgrns06d  whgrns10d  whgrns10d  whgrns10d  whgrns10d whgrns10d whgrns20d; 
      array poultrcon{21} poultr80d  poultr80d  poultr84d  poultr86d  poultr86d   poultr90d  poultr90d   poultr94d  poultr94d  poultr98d  poultr98d  poultr02d  poultr02d  poultr06d  poultr06d  poultr10d  poultr10d  poultr10d  poultr10d poultr10d poultr20d; 
      array fishalcon{21} fishal80d  fishal80d  fishal84d  fishal86d  fishal86d   fishal90d  fishal90d   fishal94d  fishal94d  fishal98d  fishal98d  fishal02d  fishal02d  fishal06d  fishal06d  fishal10d  fishal10d  fishal10d  fishal10d fishal10d fishal20d; 
      array regeggcon{21} regegg80d  regegg80d  regegg84d  regegg86d  regegg86d   regegg90d  regegg90d   regegg94d  regegg94d  regegg98d  regegg98d  regegg02d  regegg02d  regegg06d  regegg06d  regegg10d  regegg10d  regegg10d  regegg10d regegg10d regegg20d; 
      array nutsalcon{21} nutsal80d  nutsal80d  nutsal84d  nutsal86d  nutsal86d   nutsal90d  nutsal90d   nutsal94d  nutsal94d  nutsal98d  nutsal98d  nutsal02d  nutsal02d  nutsal06d  nutsal06d  nutsal10d  nutsal10d  nutsal10d  nutsal10d nutsal10d nutsal20d; 
      array nutlegcon{21} nutleg80d  nutleg80d  nutleg84d  nutleg86d  nutleg86d   nutleg90d  nutleg90d   nutleg94d  nutleg94d  nutleg98d  nutleg98d  nutleg02d  nutleg02d  nutleg06d  nutleg06d  nutleg10d  nutleg10d  nutleg10d  nutleg10d nutleg10d nutleg20d; 
      array tdairycon{21} tdairy80d  tdairy80d  tdairy84d  tdairy86d  tdairy86d   tdairy90d  tdairy90d   tdairy94d  tdairy94d  tdairy98d  tdairy98d  tdairy02d  tdairy02d  tdairy06d  tdairy06d  tdairy10d  tdairy10d  tdairy10d  tdairy10d tdairy10d tdairy20d; 
      array ldairycon{21} ldairy80d  ldairy80d  ldairy84d  ldairy86d  ldairy86d   ldairy90d  ldairy90d   ldairy94d  ldairy94d  ldairy98d  ldairy98d  ldairy02d  ldairy02d  ldairy06d  ldairy06d  ldairy10d  ldairy10d  ldairy10d  ldairy10d ldairy10d ldairy20d; 
      array hdairycon{21} hdairy80d  hdairy80d  hdairy84d  hdairy86d  hdairy86d   hdairy90d  hdairy90d   hdairy94d  hdairy94d  hdairy98d  hdairy98d  hdairy02d  hdairy02d  hdairy06d  hdairy06d  hdairy10d  hdairy10d  hdairy10d  hdairy10d hdairy10d hdairy20d; 
      array legumecon{21} legume80d  legume80d  legume84d  legume86d  legume86d   legume90d  legume90d   legume94d  legume94d  legume98d  legume98d  legume02d  legume02d  legume06d  legume06d  legume10d  legume10d  legume10d  legume10d legume10d legume20d; 
      array soyprocon{21} soypro80d  soypro80d  soypro84d  soypro86d  soypro86d   soypro90d  soypro90d   soypro94d  soypro94d  soypro98d  soypro98d  soypro02d  soypro02d  soypro06d  soypro06d  soypro10d  soypro10d  soypro10d  soypro10d soypro10d soypro20d; 
      array dcoffcon {21} dcaf_s84   dcaf_s84   dcaf_s84   dcaf_s86   dcaf_s86    dcaf_s90   dcaf_s90    dcaf_s94   dcaf_s94   dcaf_s98   dcaf_s98   dcaf_s02   dcaf_s02   dcaf_s06   dcaf_s06   dcaf_s10   dcaf_s10   dcaf_s10   dcaf_s10  dcaf_s10  dcaf_s20 ;
      array coffcon  {21} coff_s80   coff_s84   coff_s84   coff_s86   coff_s86    coff_s90   coff_s90    coff_s94   coff_s94   coff_s98   coff_s98   coff_s02   coff_s02   coff_s06   coff_s06   coff_s10   coff_s10   coff_s10   coff_s10  coff_s10  coff_s20 ;
      array nacon    {21} na80a      na84a      na84a      na86a      na86a       na90a      na90a       na94a      na94a      na98a      na98a      na02a      na02a      na06a      na06a      na10a      na10a      na10a      na10a     na10a     na20a    ;
 

      array ssbcon   {21} ssb80      ssb80      ahei2010_ssb84        ahei2010_ssb86         ahei2010_ssb86         ahei2010_ssb90        ahei2010_ssb90        ahei2010_ssb94        ahei2010_ssb94        ahei2010_ssb98        ahei2010_ssb98 
                                                ahei2010_ssb02        ahei2010_ssb02         ahei2010_ssb06         ahei2010_ssb06        ahei2010_ssb10        ahei2010_ssb10        ahei2010_ssb10        ahei2010_ssb10        ahei2010_ssb10        ssb20;
      array aheis    {21} ahei_nomt_80          ahei_nomt_80          ahei_nomt_84           ahei_nomt_86           ahei_nomt_86          ahei_nomt_90          ahei_nomt_90          ahei_nomt_94          ahei_nomt_94          ahei_nomt_98
                          ahei_nomt_98          ahei_nomt_02          ahei_nomt_02           ahei_nomt_06           ahei_nomt_06          ahei_nomt_10          ahei_nomt_10          ahei_nomt_10          ahei_nomt_10          ahei_nomt_10          ahei_nomt_10;
      array aheia    {21} nAHEI80a              nAHEI80a              ahei2010_84            ahei2010_86            ahei2010_86           ahei2010_90           ahei2010_90           ahei2010_94           ahei2010_94           ahei2010_98
                          ahei2010_98           ahei2010_02           ahei2010_02            ahei2010_06            ahei2010_06           ahei2010_10           ahei2010_10           ahei2010_14           ahei2010_14           ahei2010_10           ahei2010_10;

      array ameds    {21} amed_80    amed_80      amed_84     amed_86     amed_86      amed_90     amed_90       amed_94    amed_94      amed_98    amed_98     amed_02     amed_02     amed_06     amed_06     amed_10     amed_10      amed_10    amed_10    amed_10      amed_20; 
      array amedv    {21} emed80v   emed80v     emed84v    emed86v    emed86v     emed90v    emed90v      emed94v   emed94v     emed98v   emed98v    emed02v    emed02v    emed06v    emed06v    emed10v    emed10v     emed10v   emed10v   emed10v     emed20v; 
      array wtchg    {21} wtchg80   wtchg82     wtchg84    wtchg86    wtchg88     wtchg90    wtchg92      wtchg94   wtchg96     wtchg98   wtchg00    wtchg02    wtchg04    wtchg06    wtchg08    wtchg10    wtchg12     wtchg14   wtchg16   wtchg16     wtchg20 ;
      array SES      {21} nSES_86   nSES_86     nSES_86    nSES_86    nSES_86     nSES_86    nSES_86      nSES_86   nSES_86     nSES_86   nSES_86    nSES_86    nSES_86    nSES_86    nSES_86    nSES_86    nSES_86     nSES_86   nSES_86   nSES_86     nSES_86 ;

      array amed_frt    {21} amed_frt80    amed_frt80      amed_frt84     amed_frt86     amed_frt86      amed_frt90     amed_frt90       amed_frt94    amed_frt94      amed_frt98    amed_frt98     amed_frt02     amed_frt02     amed_frt06     amed_frt06     amed_frt10     amed_frt10      amed_frt10    amed_frt10    amed_frt10      amed_frt20; 
      array amed_veg    {21} amed_veg80    amed_veg80      amed_veg84     amed_veg86     amed_veg86      amed_veg90     amed_veg90       amed_veg94    amed_veg94      amed_veg98    amed_veg98     amed_veg02     amed_veg02     amed_veg06     amed_veg06     amed_veg10     amed_veg10      amed_veg10    amed_veg10    amed_veg10      amed_veg20; 
      array amed_leg    {21} amed_leg80    amed_leg80      amed_leg84     amed_leg86     amed_leg86      amed_leg90     amed_leg90       amed_leg94    amed_leg94      amed_leg98    amed_leg98     amed_leg02     amed_leg02     amed_leg06     amed_leg06     amed_leg10     amed_leg10      amed_leg10    amed_leg10    amed_leg10      amed_leg20; 
      array amed_nut    {21} amed_nut80    amed_nut80      amed_nut84     amed_nut86     amed_nut86      amed_nut90     amed_nut90       amed_nut94    amed_nut94      amed_nut98    amed_nut98     amed_nut02     amed_nut02     amed_nut06     amed_nut06     amed_nut10     amed_nut10      amed_nut10    amed_nut10    amed_nut10      amed_nut20; 
      array amed_rmt    {21} amed_rmt80    amed_rmt80      amed_rmt84     amed_rmt86     amed_rmt86      amed_rmt90     amed_rmt90       amed_rmt94    amed_rmt94      amed_rmt98    amed_rmt98     amed_rmt02     amed_rmt02     amed_rmt06     amed_rmt06     amed_rmt10     amed_rmt10      amed_rmt10    amed_rmt10    amed_rmt10      amed_rmt20; 
      array amed_fish   {21} amed_fish80    amed_fish80      amed_fish84     amed_fish86     amed_fish86      amed_fish90     amed_fish90       amed_fish94    amed_fish94      amed_fish98    amed_fish98     amed_fish02     amed_fish02     amed_fish06     amed_fish06     amed_fish10     amed_fish10      amed_fish10    amed_fish10    amed_fish10      amed_fish20; 
      array amed_whgrn    {21} amed_whgrn80    amed_whgrn80      amed_whgrn84     amed_whgrn86     amed_whgrn86      amed_whgrn90     amed_whgrn90       amed_whgrn94    amed_whgrn94      amed_whgrn98    amed_whgrn98     amed_whgrn02     amed_whgrn02     amed_whgrn06     amed_whgrn06     amed_whgrn10     amed_whgrn10      amed_whgrn10    amed_whgrn10    amed_whgrn10      amed_whgrn20; 
      array amed_etoh    {21} amed_etoh80    amed_etoh80      amed_etoh84     amed_etoh86     amed_etoh86      amed_etoh90     amed_etoh90       amed_etoh94    amed_etoh94      amed_etoh98    amed_etoh98     amed_etoh02     amed_etoh02     amed_etoh06     amed_etoh06     amed_etoh10     amed_etoh10      amed_etoh10    amed_etoh10    amed_etoh10      amed_etoh20; 
      array amed_ms    {21} amed_ms80    amed_ms80      amed_ms84     amed_ms86     amed_ms86      amed_ms90     amed_ms90       amed_ms94    amed_ms94      amed_ms98    amed_ms98     amed_ms02     amed_ms02     amed_ms06     amed_ms06     amed_ms10     amed_ms10      amed_ms10    amed_ms10    amed_ms10      amed_ms20; 

       
   do i=2 to 21;
      if trmeatcon{i}=.           then trmeatcon{i} = trmeatcon{i-1};
      if urmeatcon{i}=.           then urmeatcon{i} = urmeatcon{i-1};
      if prmeatcon{i}=.           then prmeatcon{i} = prmeatcon{i-1};
      if fruitscon{i}=.           then fruitscon{i} = fruitscon{i-1};
      if vegealcon{i}=.           then vegealcon{i} = vegealcon{i-1};
      if whgrnscon{i}=.           then whgrnscon{i} = whgrnscon{i-1};
      if poultrcon{i}=.           then poultrcon{i} = poultrcon{i-1};
      if fishalcon{i}=.           then fishalcon{i} = fishalcon{i-1};
      if regeggcon{i}=.           then regeggcon{i} = regeggcon{i-1};
      if nutsalcon{i}=.           then nutsalcon{i} = nutsalcon{i-1};
      if nutlegcon{i}=.           then nutlegcon{i} = nutlegcon{i-1};
      if tdairycon{i}=.           then tdairycon{i} = tdairycon{i-1};
      if ssbcon{i}=.              then ssbcon{i}    = ssbcon{i-1};
      if ldairycon{i}=.           then ldairycon{i} = ldairycon{i-1};
      if hdairycon{i}=.           then hdairycon{i} = hdairycon{i-1};
      if legumecon{i}=.           then legumecon{i} = legumecon{i-1};
      if coffcon{i}=.             then coffcon{i}   = coffcon{i-1};
      if dcoffcon{i}=.            then dcoffcon{i}  = dcoffcon{i-1};
      if nacon{i}=.               then nacon{i}     = nacon{i-1};
      if ssbcon{i}=.              then ssbcon{i}    = ssbcon{i-1};
      if aheis{i}=.               then aheis{i}     = aheis{i-1}; 
      if aheia{i}=.               then aheia{i}     = aheia{i-1};
      if ameds{i}=.               then ameds{i}     = ameds{i-1};
      if amedv{i}=.               then amedv{i}     = amedv{i-1};
      if wt{i}=0 or wt{i}=.       then wt{i}     = wt{i-1};
      if wtchg{i}=.               then wtchg{i}  = wtchg{i-1};
      if nsmk{i}=0 or nsmk{i}=.   then nsmk{i}   = nsmk{i-1}; 
      if actm{i}=.                then actm{i} = actm{i-1};
      if aheis{i}=.               then aheis{i}= aheis{i-1};
      if hor{i}=0 or hor{i}=.     then hor{i}    = hor{i-1};  
      if mvyn{i}=.                then mvyn{i}   = mvyn{i-1};
      if aspu{i}=.                then aspu{i}   = aspu{i-1};
      if calor{i}=.               then calor{i}  = calor{i-1};
      if alco{i}=.                then alco{i}   = alco{i-1};
      if SES{i}=.                 then SES{i}    = SES{i-1};
      if cance{i-1}=1     then cance{i}=1; 
      if repmi{i-1}=1     then repmi{i}=1;
      if repstrk{i-1}=1   then repstrk{i}=1;
      if diab{i-1}=1      then diab{i}=1;
      if chol{i-1}=1      then chol{i}=1;
      if rphbp{i-1}=1     then rphbp{i}=1;  
      if depr{i-1}=1      then depr{i}=1;
      /*if antid{i-1}=1     then antid{i}=1;
      if antihp{i-1}=1    then antihp{i}=1;  
      if antitc{i-1}=1    then antitc{i}=1;*/

      if amed_frt{i}=.               then amed_frt{i}     = amed_frt{i-1};
      if amed_veg{i}=.               then amed_veg{i}     = amed_veg{i-1};
      if amed_leg{i}=.               then amed_leg{i}     = amed_leg{i-1};
      if amed_nut{i}=.               then amed_nut{i}     = amed_nut{i-1};
      if amed_rmt{i}=.               then amed_rmt{i}     = amed_rmt{i-1};
      if amed_fish{i}=.               then amed_fish{i}     = amed_fish{i-1};
      if amed_whgrn{i}=.               then amed_whgrn{i}     = amed_whgrn{i-1};
      if amed_etoh{i}=.               then amed_etoh{i}     = amed_etoh{i-1};
      if amed_ms{i}=.               then amed_ms{i}     = amed_ms{i-1};

    end;
     
/*** Set cutoff at 2023 Jan 31 ***/
/*************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.     ****
*************************************************************************/

 cutoff=1477;    
   do i=1 to 19;
      if (irt{i}<(942+24*i) | irt{i}>=(966+24*i)) then irt{i}=942+24*i;
  end;

  if irt20<1439 or irt20>=1457 then irt20=1439;
  /* irt20:  	$range 1439-1457 	$label 1439.nov 19;\ \   1457.may 21   */

%beginex();

/********* DO LOOP OVER TIME PERIOD **************/

   do i=1 to DIM(irt)-1;     

    period=i;

    do j=1 to DIM(tvar);
    tvar{j}=0;
    end;
   
    tvar{i}=1; 


    /*** total AD ***/
      adcase=0;  tad=irt{i+1}-irt{i};    
      if irt{i}<dtdx_ad<=irt{i+1} then do; adcase=1; tad=dtdx_ad-irt{i};  end;
      if irt{i} lt dtdth le irt{i+1} then tad=min(tad, dtdth-irt{i});

    /*** confirmed AD ***/
      confadcase=0;      
      if irt{i}<dtdx_ad<=irt{i+1} and alz0220=1 and dementiadeath=1 then confadcase=1; 

    /* incident non-fatal AD */
      nfadcase=0;     
      if irt{i}<dtdx_ad<=irt{i+1} and alz0220=1 then  nfadcase=1;  

    /* incident fatal dementia */
      fadcase=0;   
      if irt{i}<dtdx_ad<=irt{i+1} and dementiadeath=1 then fadcase=1;   
  
      if dtdth eq 9999 then dtdth=.;
 
  ****** main exposure ******; 
  trmeat = trmeatcon{i}; 
  prmeat = prmeatcon{i};
  urmeat = urmeatcon{i}; 

  fruits = fruitscon{i};
  vegets = vegealcon{i};
  whgrns = whgrnscon{i};
  poultr = poultrcon{i};
  fishs  = fishalcon{i};
  eggs   = regeggcon{i};
  nuts   = nutsalcon{i};
  nutleg = nutlegcon{i};
  tdairy = tdairycon{i};
  daykcal= calor{i}; 
  alcocon= alco{i}; 
  ssb    = ssbcon{i};
  ldairy = ldairycon{i};
  hdairy = hdairycon{i};
  legume = legumecon{i};
  soypro = soyprocon{i};
  legsoy = sum(0,soypro,legume);
  coffee = coffcon{i};
  dcoffee= dcoffcon{i};
  sodium = nacon{i};
  ahei8010=aheia{i};
  amedcon =ameds{i};
  amedcov =amedv{i};

  nSES     =SES{i};
  scfmem   = mem1214_7nhs;  
  scfmem_z = zmem1214_7nhs;
  scfmem_sd= mem1214_7nhs_sd;
  PGScon   = PGS002280_scaled;
  PGScon2   = PGS000334_scaled;

  amed_frtcon =amed_frt{i};
  amed_vegcon =amed_veg{i};
  amed_legcon =amed_leg{i};
  amed_nutcon =amed_nut{i};
  amed_rmtcon =amed_rmt{i};
  amed_fishcon =amed_fish{i};
  amed_whgrncon =amed_whgrn{i};
  amed_etohcon =amed_etoh{i};
  amed_mscon =amed_ms{i};


  ****** covariates ******;
      ****** Nurse's education level ******;
      if hiedu in (3,4) then hiedu=3;
      %indic3(vbl=hiedu, prefix=hiedu, min=2, max=3, reflev=1, missing=., usemiss=1,
        label1='RN', 
        label2='BA', 
        label3='MA or DR');

      ****** HUSBAND'S EDUCATION ******;  
      %indic3(vbl=husbedu, prefix=husbedu, min=2, max=3, reflev=1, missing=., usemiss=1,
            label1='=<H.S. grad',  
            label2='College grad', 
            label3='Grad school');

      /***smoking status***/
      if nsmk{i}=1 then smm{i}=1;  *never;
        else if nsmk{i} in(2,3,4,5,6,7,8) then smm{i}=2; *past smoker;
        else if nsmk{i} in (9,10,15) then smm{i}=3; *current, 1-14 cigs ;
        else if nsmk{i}=11 then smm{i}=4; *current, 15-24 cigs ;
        else if nsmk{i} in(12,13,14) then smm{i}=5; *current, 25+ cigs ;
        else smm{i}=.; *missing;
        smkk=smm{i};
        %indic3(vbl=smkk, prefix=smkk, min=2, max=5, reflev=1, missing=., usemiss=1,
	        label1='Never smoke',
	        label2='Past smoker',
	        label3='Current smoker 1-14 cigs',
	        label4='Current smoker 15-24 cigs',
	        label5='Current smoker 25+ cigs');

        if smkk in (1,2) then smkcat=smkk;
        else if smkk in (3,4,5) then smkcat=3;
        else smkcat=.;
        %indic3(vbl=smkcat, prefix=smkcat, min=2, max=3, reflev=1, missing=., usemiss=0,
	        label1='Never smoke',
	        label2='Past smoker',
	        label3='Current smoker');

      /***exercise***/
       actcon=actm{i};
       if 0=<actcon<2  	     then actcc=1;
       else if 2=<actcon<3.5       then actcc=2;
       else if actcon>=3.5     then actcc=3;
       else actcc=2;
       %indic3(vbl=actcc, prefix=actcc, min=2, max=3, reflev=1, missing=., usemiss=0,
	        label1='<2 m/v activity hours/week',
	        label2='2 to < 3.5 hours/week',
	        label3='3.5+ hours/week');                                                               
       /* actcon=actm{i};
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
	        label5='5.5+ hours/week');  */
/*** alcohol ***/
       if alcocon=0.0                          then   alcc=1;
         else if alcocon>0.0   and alcocon<5.0   then   alcc=2;
           else if alcocon>=5.0  and alcocon<10.0  then   alcc=3;
             else if alcocon>=10.0 and alcocon<15.0  then   alcc=4;
               else if alcocon>=15.0 and alcocon<30.0  then   alcc=5;
                else if alcocon>=30.0                    then   alcc=6;   
                 else                                              alcc=.;   
 
       %indic3(vbl=alcc, prefix=alcc, min=1, max=6, reflev=3, missing=., usemiss=0,
                label1='0 g/d',
                label2='0.1-4.9 g/d',
                label3='5.0-14.9 g/d',
                label4='5.0-14.9 g/d',
                label5='15.0-29.9 g/d',
      	        label6='30+ g/d');

       /*** dietary pattern***/
        ahei=aheis{i}; 

      /****** DEFINE BMI ******/ 
         bmi{i}=(wt{i}*0.45359237)/(height*height);
         if bmi{i}=0 then bmi{i}=.;
   
      /*** Indicator for BMI ***/ 
             bmicon=bmi{i}; 
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

         basebmicon=bmi80; 
         if bmi80<23 then basebmic=1;
           else if 23=<bmi80<25 then basebmic=2;
             else if 25=<bmi80<30 then basebmic=3;
               else if 30=<bmi80<35 then basebmic=4;
                 else basebmic=5;
          if bmi80=. then basebmic=.;
         %indic3(vbl=basebmic, prefix=basebmic, min=1, max=5, reflev=2, missing=., usemiss=0,
                label1='<23',
                label2='23-24.9',
                label3='25-29.9',
                label4='30-34.9',
      	    label5='35+');

          wtchgcon=wtchg{i};
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

 
  ******covariables******;
        /****** DEFINE AGE ******/ 
         age{i}=int((irt{i}-bdt)/12);      
            if age{i} =< 0 then age{i}=.;

            agecon=age{i};

            agegp=int((age{i}-30)/5);     **** Define the agegp in the i-th period;
            if agegp<=0 then agegp=1;
            if agegp>8 then agegp=8;
           %indic3(vbl=agegp, prefix=agegp, reflev=1, min=1, max=8, usemiss=0,
                     label1='age<40',
                     label2='age40-44',
                     label3='age45-49',
                     label4='age50-54',
                     label5='age55-59',
                     label6='age60-64',
                     label7='age65-69',
                     label8='age >= 70');

       /*** phms ***/ 
         if hor{i} =1 then phmsstatus=1;      
         else if hor{i} in (2,5) then phmsstatus=2;      
         else if hor{i} =3  then phmsstatus=3;      
         else if hor{i} =4  then phmsstatus=4;     
         else phmsstatus=.; 
     %indic3(vbl=phmsstatus, prefix=phmsstatus, reflev=1, missing=., usemiss=1, min=2, max=4,      
                     label1='pre mnp',      
                     label2='never/unknown pmhuser, post mnp',      
                     label3='curr pmh user and post mnp',      
                     label4='past pmhuser and post mnp'); 

       ****** living alone ******;
        live_alone=alone{i};
        if live_alone ne 1 then live_alone=0;   

       ****** current married or not ******;
        marry=mar{i};
        if marry ne 1 then marry=0; 

      /*** family history of diabetes***/
         IF dbfh NE 1 THEN dbfh=0;

      /*** family history of MI***/
         IF mifh NE 1 THEN mifh=0;

      /*** family history of cancer***/
         IF cafh NE 1 THEN cafh=0;
 
     /*** aspirin ***/
            select(aspu{i});
                when (1)     aspirin=1; 
      	         otherwise    aspirin=0;
            end;  

      /*** multiple vitamin ***/
             if mvyn{i}=1 then mvit=1;
             else mvit=0;
 
     /******************sef-reported chronic conditions ******************/
 
      /****** Indicator for History of High Blood Pressure *******/
            select(rphbp{i});
                when (1)     htn=1;
      	         otherwise    htn=0;
            end;
            

      /****** Indicator for History of diabetes *******/
            select(diab{i});
                when (1)     diabetes=1;
      	         otherwise    diabetes=0;
            end;
          
      /****** Indicator for History of High TC ******/
            select(chol{i});
                when (1)     hchol=1;
      	         otherwise    hchol=0;
            end;

            if chol{i}=1 or antitc{i}=1 then hightc=1; else hightc=0; 

      /****** Indicator for History of selfreport heart disease MI & Angina ******/
            select(repmi{i});
                when (1)     chd=1;
      	         otherwise    chd=0;
            end;
 
      /****** Indicator for History of selfreport stroke ******/
            select(repstrk{i});
                when (1)     stroke=1;
      	         otherwise    stroke=0;   
            end;
 

      /****** Indicator for History of Cancer ******/
            select(cance{i});
                when (1)     cancer=1;
      	         otherwise    cancer=0;
            end;

      /****** Indicator for History of hypothyroid ******/
            select(hypo{i});
                when (1)     hypothyroid=1;
      	         otherwise    hypothyroid=0;
            end;

      /****** Indicator for History of depression ******/
            select(depr{i});
                when (1)     depression=1;
      	         otherwise    depression=0;
            end;
 
     /******************regular use of medications ******************/
          
     /****** Indicator for antidepressant user *******/
            if antid{i}=1 then antidep=1;
               else   antidep=0;
                      
      /****** Indicator for antihypertensive user in last two years, no carry forward *******/
            select(antihp{i});
                when (1)     antihbp=1;
      	         otherwise    antihbp=0;
            end;
      
       if htn=1 or antihbp=1 then highhbp=1; else highhbp=0;
 
      /****** Indicator for TC lowering *******/
            select(antitc{i});
                when (1)     antihtc=1;
      	         otherwise    antihtc=0;
            end;

       if dtdth in (0,.)then dtdth=.;
       if dtdth eq 9999 then dtdth=.;

  
    
      /****************  BASELINE EXCLUSIONS ********************************************/
      if i=1 then do;

   %exclude(exrec eq 1);                 *multiple records and not in master file nur02;

   %exclude(bdt eq .);
   %exclude(age80 eq .);

   %exclude(pdcase eq 1);
   %exclude(str80 eq 1);
   %exclude(can7680 eq 1);

   %exclude(calor80n lt 500 );
   %exclude(calor80n gt 3500 );
   %exclude(calor80n le 0);
   %exclude(calor80n eq .);  

   %exclude(bmi80 eq .);
   %exclude(emed80 eq .);
   %exclude(smkdr80 eq .);
   %exclude(smkdr80 eq 0);

   %exclude(0 lt dtdth    le irt{i} );    
   %exclude(0 lt dtdx_ad  le irt{i});     

   %output();
  end;

  else if i> 1 then do;
   %exclude(irt{i-1} le dtdth lt irt{i});
   %exclude(irt{i-1} lt dtdx_ad le irt{i});  
   %output();
  end;

 end;

          /* END OF DO-LOOP OVER TIME PERIODs */
      
      %endex();
keep id period period agecon cutoff confadcase adcase tad dtdx_ad dtdth nfadcase  fadcase 
     t80 t82 t84 t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18 t20
     irt80 irt82 irt84 irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 irt18 irt20
     trmeat prmeat urmeat 
     &hiedu_ &husbedu_ fhdem marry live_alone nSES aspirin mvit &smkk_ actcc &actcc_ &phmsstatus_ &bmic_  ahei
     bmi80 hbp_7680 db_7680 chol_7680 mi7680 can7680  
     htn highhbp diabetes hchol hightc hypothyroid depression antidep antihtc antihbp chd stroke cancer
     fruits  vegets  whgrns  poultr  fishs  eggs   nuts nutleg tdairy daykcal ssb smkk  bmic 
     bmicon actcon phmsstatus  hiedu  husbedu &alcc_ alcocon alcc
     ldairy  hdairy  legume  soypro  legsoy coffee dcoffee sodium ahei8010
     PGS002280_scaled PGS000334_scaled apoe4grp apoe_2cat APOE4 platform platforms PC1_comb PC2_comb PC3_comb PC4_comb apoe4con amedcon amedcov ADdiagage bdt
     amed_frtcon amed_vegcon amed_legcon amed_nutcon amed_rmtcon amed_fishcon amed_whgrncon amed_etohcon amed_mscon
     basebmicon basebmic &basebmic_ smkcat &smkcat_ wtchgcon wtchgcat &wtchgcat_
     scfmem scfmem_z  scfmem_sd PGScon PGScon2 dementiadeath alz0220 q02 q04 q06 q08 q10 q12 q14 q16 q20;
run; 
 
proc means n nmiss min mean std median max;
var nSES actcon hiedu smkk;
class period;
run;



%pctl9(data=pre_pm,
    varlist= nSES PC1_comb PC2_comb PC3_comb PC4_comb PGS002280_scaled PGS000334_scaled,
    numquant=3,
    mscore=T,
    quantname=q,
    cutdat=cutpoints,
    outdat=pre_pm,
    indic=T);

%pctl9(data=pre_pm,
    varlist= ahei8010 amedcon amedcov PGScon PGScon2 amed_frtcon amed_vegcon amed_legcon amed_nutcon amed_rmtcon amed_fishcon amed_whgrncon amed_etohcon amed_mscon,
    numquant=5,
    mscore=T,
    quantname=q,
    cutdat=cutpoints,
    outdat=pre_pm,
    indic=T);

data pre_pm;
  set pre_pm end=_end_; 

  if bmi80<23 then bmi80c=1;
           else if 23=<bmi80<25 then bmi80c=2;
             else if 25=<bmi80<30 then bmi80c=3;
               else if 30=<bmi80<35 then bmi80c=4;
                 else bmi80c=5;
          if bmi80=. then bmi80c=.;
         %indic3(vbl=bmi80c, prefix=bmi80c, min=1, max=5, reflev=2, missing=., usemiss=1,
                label1='18.5-22.9',
                label2='23-24.9',
                label3='25-29.9',
                label4='30-34.9',
      	        label5='35+');

      ****** neiborghood SES index   ******;
      %indic3(vbl=qnSES, prefix=qnSES, min=1, max=2, reflev=0, missing=., usemiss=1); 

      if antidep=1 or  depression=1 then depre=1; else depre=0;

      %indic3(vbl=qahei8010,   prefix=qahei8010,   min=1, max=4, reflev=0, missing=., usemiss=1); 
      %indic3(vbl=qamedcon,    prefix=qamedcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamedcov,    prefix=qamedcov,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_frtcon,    prefix=qamed_frtcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_vegcon,    prefix=qamed_vegcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_legcon,    prefix=qamed_legcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_nutcon,    prefix=qamed_nutcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_rmtcon,    prefix=qamed_rmtcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_fishcon,    prefix=qamed_fishcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_whgrncon,    prefix=qamed_whgrncon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_etohcon,    prefix=qamed_etohcon,    min=1, max=4, reflev=0, missing=., usemiss=0); 
      %indic3(vbl=qamed_mscon,    prefix=qamed_mscon,    min=1, max=4, reflev=0, missing=., usemiss=0); 


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
      %indic3(vbl=platforms,     prefix=platforms,     min=2, max=6, reflev=1, missing=., usemiss=0,
            label1='affy',
            label2='gsa',
            label3='huco2',
            label4='illu',
            label5='omni',
            label6='onco' ); 
run;

proc means median; var nSES;    class qnSES; run;
proc means median; var amedcon; class qamedcon; run;
proc means median; var actcon;  class actcc; run;

proc means median; var amed_etohcon; class qamed_etohcon; run;
proc means median; var amed_mscon; class qamed_mscon; run;
proc means median; var amed_fishcon; class qamed_fishcon; run;
proc means median; var amed_rmtcon; class qamed_rmtcon; run;
proc means median; var amed_whgrncon; class qamed_whgrncon; run;
proc means median; var amed_legcon; class qamed_legcon; run;
proc means median; var amed_frtcon; class qamed_frtcon; run;
proc means median; var amed_nutcon; class qamed_nutcon; run;
proc means median; var amed_vegcon; class qamed_vegcon; run;

/* proc means nolabels n nmiss min mean std max;run; */

data nhs1data;
set pre_pm;
cohort=1;
sex=1;
id=100000000+id;
interval=period;
birthday=bdt;
proff=1;
FORMAT _all_;
INFORMAT _all_;

****** P for trend: beta of per unit is missingless as 99999 is counted, just P for trend ******;
 
if wtchgcon=. then wtchgcon=99999;
wtchgmiss=0; if wtchgcon=. then wtchgmiss=1;

****** P for trend ******;
if actcc=1 then medactcc=0;
if actcc=2 then medactcc=1;
if actcc=3 then medactcc=2.5;
if actcc=4 then medactcc=4.79;
if actcc=5 then medactcc=8;
if actcc=. then medactcc=99999;
actccmiss=0; if actcc=. then actccmiss=1;

if nSES=. then medqnSES=99999;
nSESmiss=0; if nSES=. then nSESmiss=1;

run;

proc datasets;
  delete pre_pm;
  run;
