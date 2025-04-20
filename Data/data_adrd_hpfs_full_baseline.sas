options nocenter ls=150 ps=100;
filename hpstools '/proj/hpsass/hpsas00/hpstools/sasautos';  /*** These are for read macros ***/
filename channing '/usr/local/channing/sasautos';
libname readfmt '/proj/hpsass/hpsas00/formats';
options mautosource sasautos=(channing hpstools);
options fmtsearch=(readfmt) nofmterr;
  
**************************************************
*                    Outcome                     *
**************************************************;

%hp_dead (keep= mmdth yydth icda newicda dtdth );
   
   if 0=<yydth<50 then yydth=yydth+100;
   
   newicda=compress(icda, 'E');
   if newicda='V812' then newicda=812;
   
   if mmdth<0 and yydth>0 then mod=6;
   dtdth=(yydth*12)+mmdth;
   if dtdth=. then delete ;
   /*proc freq; tables dtdth newicda yydth;*/
   proc sort;
      by id;
      run;
 
 
  /***********read in CVD**********/
   filename cardio10 '/proj/hpchds/hpchd0q/CARDIO2016/cardio2016.02062020';                  

   data cvd10;
        infile cardio10 lrecl=700;
        input
        id 1-6
        mi_nf 7  /*0=noncase, 1=definite, 2=probable, 3=selfreport, 4=rejected*/
        mnyr_nmi 8-11
        fat_chd 27 /* fatal mi; 0=noncase, 1=def, 2=presumed, 3=sudden*/
        mnyr_fat 28-31
        cabg 12 /*same as nf mi*/
        mnyr_ca 13-16
        angina 17  /*same*/
        mnyr_an 18-21
        tot_chd 32 /*nf, fatal mi + cabg*/
        mnyr_chd 33-36
        tot_cvd 37
        mnyr_cvd 38-41
        nf_str 22 /*0=noncase, 1=definite, 2=probable, 3=selfreport, 4=rejected**/
        fat_str 42 /* fatal mi; 0=noncase, 1=def, 2=presumed, 3=sudden*/
        mnyr_nfst 23-26
        mnyr_tst 47-50 /** p-t total stroke **/
        s1 $ 44
         /*
         Stroke Subcode 1, type of stroke
         $label 1.430: subarachnoid hemorrhage;\
         2.icda 431:intracerebral hemorrhage;\
         3.432-4,6:thrombotic stroke;\
         4.432-4,6:embolic stroke;\
         5.432-4,6:thrombotic/embolic stroke;\
         6.436-8:-unknown type;\
         */
        ;

        /** MI **/
        nfmi=0;
        if mi_nf in (1,2) then nfmi=1;

           nfmid=0; if mi_nf=1 then nfmid=1; /***for method part***/
           nfmip=0; if mi_nf=2 and nfmid ne 1 then nfmip=1;

        fatchd=0;
        if fat_chd in (1,2,3) then fatchd=1;

           fatchdd=0;   if fat_chd=1 then fatchdd=1;/***for method part***/
           fatchdp=0;   if fat_chd=2 and fatchdd ne 1 then fatchdp=1;
           fatchds=0;   if fat_chd=3 and fatchdd ne 1 and fatchds ne 1 then fatchds=1;

        tot_mi = sum (nfmi,fatchd);

        totmi=0;
        if tot_mi in (1,2) then totmi=1;


        /** stroke **/

        if s1 in (1,2) then hem_str=1; else hem_str=0;
        if s1 in (3,4,5) then isc_str=1; else isc_str=0;
        if s1=6 then unk_str=1; else unk_str=0;


        nfstr=0;
        if nf_str in (1,2) then nfstr=1;

        fatstr=0;
        if fat_str in (1,2) then fatstr=1;

        tot_str=sum(fatstr, nfstr);

        totstr=0;
        if tot_str in (1,2) then totstr=1;

        /* cabg */
        totcabg=0;
        if cabg in (1,2) then totcabg=1;
 
        if totmi=1 or totstr=1 or totcabg=1
           then do;totcvd=1;dtdxcvd=mnyr_cvd; end; 
           else do;totcvd=0;dtdxcvd=.;end;

        /*label
        nfstr='non-fatal stroke'
        fatstr='fatal stroke'
        totstr='nonfatal + fatal stroke'
        hem_str='hemorrhgic stroke'
        isc_str='ischemic stroke'
        unk_str='unknown type'
        totcvd='mi+stk+cabg';*/ 
    proc sort; by id;
    run;



/******diabetes cases ******/
     data diabetes;
         infile '/proj/hpdbxs/hpdbx00/diab_cases_8620/diab8620.05MAY23' lrecl=49 recfm=d;
             input
             @1		hpfsid		$   8.
             @1		id		    6.
             @7		cd		    2.
             @10		yrdx		    2.
             @13		dbconf		    2.
             @16		diag_mn		    2.
             @19		diag_yr		    2.
             @22		ketoac		    1.
             @23		coma		    1.
             @24		loss		    1.
             @25		hunger		    1.
             @26		thirst		    1.
             @27		urin		    1.
             @28		visual		    1.
             @30		type		    1.
             @32		prob		    1.
             @34		dbcase		    1.
             @36		dtdx		    4.
             @41		symptoms	    1.
             @43		fu		    4.
             @47		a1cmid		    1.
             @48		a1chigh		    1.
             ;
                  dtdxdb=9999;
                  if dtdx>0 and dbcase=1 then dtdxdb=dtdx;

                  type2db=0;
                  if type=2 and prob=1 then type2db=1; *type 2 diabees;

                  type1db=0;
                  if type eq 1 or type eq 5 then type1db=1;

                  dtdxdb2=9999;
                  if dtdx > 0 and type2db=1 then dtdxdb2=dtdx;  

                  keep  id type2db dtdxdb type1db dtdxdb2;
                  run;
                  proc sort; by id; run; 

/***********read in total cancer**********/

             data totalcancer;
                  infile '/proj/hpchds/hpchd0y/Cancer_2018_data/Cancer2016_Apr2021_withduplicate.dat';
                  input
                  id 1-6
                  dateca 8-11  /* date of cancer */
                  dod 12-15  /* date of death for all deaths */
                  cancer   16 /* all cancers 1,0 */
                  confca 47-48
                  cancerdth 49
                  can_icdx      52-54
                  ;

                  dtdxcan=9999;
                  if dateca>0 and 11=<confca=<19 then dtdxcan=dateca;

                  run;

******************************************
*    Food and dietary patterns           *
******************************************;

    %include "food_hpfs.sas"; *including ahei14 ahei18 emed14 emed18;

    %include "mind_hpfs.sas"; 

    /****** AHEI ******/
    libname AHEIdat "/proj/hpalcs/hpalc0b/DIETSCORES/HPFS/"; 
    data ahei86; set AHEIdat.hnahei86l; ahei86=nAHEI86a; ahei86_noal=nAHEI86_noal; ahei_nomt_86=ahei86-meataI86; keep id ahei86 ahei86_noal ssb86 ahei_nomt_86; proc sort; by id; run;
    data ahei90; set AHEIdat.hnahei90l; ahei90=nAHEI90a; ahei90_noal=nAHEI90_noal; ahei_nomt_90=ahei90-meataI90; keep id ahei90 ahei90_noal ssb90 ahei_nomt_90; proc sort; by id; run;
    data ahei94; set AHEIdat.hnahei94l; ahei94=nAHEI94a; ahei94_noal=nAHEI94_noal; ahei_nomt_94=ahei94-meataI94; keep id ahei94 ahei94_noal ssb94 ahei_nomt_94; proc sort; by id; run;
    data ahei98; set AHEIdat.hnahei98l; ahei98=nAHEI98a; ahei98_noal=nAHEI98_noal; ahei_nomt_98=ahei98-meataI98; keep id ahei98 ahei98_noal ssb98 ahei_nomt_98; proc sort; by id; run;
    data ahei02; set AHEIdat.hnahei02l; ahei02=nAHEI02a; ahei02_noal=nAHEI02_noal; ahei_nomt_02=ahei02-meataI02; keep id ahei02 ahei02_noal ssb02 ahei_nomt_02; proc sort; by id; run;
    data ahei06; set AHEIdat.hnahei06l; ahei06=nAHEI06a; ahei06_noal=nAHEI06_noal; ahei_nomt_06=ahei06-meataI06; keep id ahei06 ahei06_noal ssb06 ahei_nomt_06; proc sort; by id; run;
    libname dietsc '/udd/hpypl/review/dietscore/';
    data ahei10; set dietsc.hnahei10l; ahei10=nAHEI10a; ahei10_noal=nAHEI10_noal; ahei_nomt_10=ahei10-meataI10; keep id ahei10 ssb10 ahei10_noal ahei_nomt_10; proc sort; by id; run;  

    /****** PDIs ******/ 
    data p86hpfs; infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX86HPFSMB.dat' lrecl=400;
         input id 1-6
            @8 pdi86 6.2
           @18 hpdi86 6.2
           @28 updi86 6.2; proc sort nodupkey; by id; run;

    data p90hpfs; infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX90HPFSMB.dat' lrecl=400;
         input id 1-6
            @8 pdi90 6.2
           @18 hpdi90 6.2       
           @28 updi90 6.2; proc sort nodupkey; by id; run;

    data p94hpfs;  infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX94HPFSMB.dat' lrecl=400;
         input id 1-6
            @8 pdi94 6.2
           @18 hpdi94 6.2
           @28 updi94 6.2; proc sort nodupkey; by id; run;

    data p98hpfs;   infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX98HPFSMB.dat' lrecl=400;
         input id 1-6
            @8 pdi98 6.2
           @18 hpdi98 6.2
           @28 updi98 6.2; proc sort nodupkey; by id; run;

    data p02hpfs;   infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX02HPFSMB.dat' lrecl=400;
         input id 1-6
            @8 pdi02 6.2
           @18 hpdi02 6.2
           @28 updi02 6.2; proc sort nodupkey; by id; run;

    data p06hpfs;   infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX06HPFSMB.dat' lrecl=400;
        input id 1-6
            @8 pdi06 6.2
           @18 hpdi06 6.2
           @28 updi06 6.2; proc sort nodupkey; by id; run;

    data p10hpfs;   infile '/udd/hpbme/PDIs/HPFS/PLANTINDEX10HPFSMB.dat' lrecl=400;
        input id 1-6
            @8 pdi10 6.2
           @18 hpdi10 6.2
           @28 updi10 6.2; proc sort nodupkey; by id; run;
 
    /**AMED**/ 
    data amed86hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED86HPFS' lrecl=200;
    input id 1-6  
      @8 emed86h 6.2  
     @18 alcoh286 6.2 
     @28 emsrat86h 6.2 
     @38 efsh86h 6.2 
     @48 emt86h 6.2 
     @58 ewgr86h 6.2 
     @68 eleg86h 6.2 
     @78 efru86h 6.2 
     @88 enut86h 6.2 
     @98 evegg86h 6.2;
     proc sort nodupkey; by id; run; 

    data amed90hpfs; 
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED90HPFS' lrecl=200;
    input id 1-6  
      @8 emed90h 6.2  
     @18 alcoh290 6.2 
     @28 emsrat90h 6.2 
     @38 efsh90h 6.2 
     @48 emt90h 6.2 
     @58 ewgr90h 6.2 
     @68 eleg90h 6.2 
     @78 efru90h 6.2 
     @88 enut90h 6.2 
     @98 evegg90h 6.2 ; 
    proc sort nodupkey; by id; run; 

    data amed94hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED94HPFS' lrecl=200;
    input id 1-6  
      @8 emed94h 6.2  
     @18 alcoh294 6.2 
     @28 emsrat94h 6.2 
     @38 efsh94h 6.2 
     @48 emt94h 6.2 
     @58 ewgr94h 6.2 
     @68 eleg94h 6.2 
     @78 efru94h 6.2 
     @88 enut94h 6.2 
     @98 evegg94h 6.2;  
     proc sort nodupkey; by id; run; 

    data amed98hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED98HPFS' lrecl=200;
    input id 1-6  
      @8 emed98h 6.2  
     @18 alcoh298 6.2 
     @28 emsrat98h 6.2 
     @38 efsh98h 6.2 
     @48 emt98h 6.2 
     @58 ewgr98h 6.2 
     @68 eleg98h 6.2 
     @78 efru98h 6.2 
     @88 enut98h 6.2 
     @98 evegg98h 6.2;
     proc sort nodupkey; by id; run; 

    data amed02hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED02HPFS' lrecl=200;
    input id 1-6  
      @8 emed02h 6.2  
     @18 alcoh202 6.2 
     @28 emsrat02h 6.2 
     @38 efsh02h 6.2 
     @48 emt02h 6.2 
     @58 ewgr02h 6.2 
     @68 eleg02h 6.2 
     @78 efru02h 6.2 
     @88 enut02h 6.2 
     @98 evegg02h 6.2;
     proc sort nodupkey; by id; run; 

    data amed06hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED06HPFS' lrecl=200;
    input id 1-6  
      @8 emed06h 6.2  
     @18 alcoh206 6.2 
     @28 emsrat06h 6.2 
     @38 efsh06h 6.2 
     @48 emt06h 6.2 
     @58 ewgr06h 6.2 
     @68 eleg06h 6.2 
     @78 efru06h 6.2 
     @88 enut06h 6.2 
     @98 evegg06h 6.2; 
     proc sort nodupkey; by id; run; 

    data amed10hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/EMED10HPFS' lrecl=200;
    input id 1-6
      @8 emed10h 6.2
     @18 alcoh210 6.2
     @28 emsrat10h 6.2
     @38 efsh10h 6.2
     @48 emt10h 6.2
     @58 ewgr10h 6.2
     @68 eleg10h 6.2
     @78 efru10h 6.2
     @88 enut10h 6.2
     @98 evegg10h 6.2;
     proc sort nodupkey; by id; run; 

    /**DASH**/ 
    data dash86hpfs;
       infile '/proj/hpalcs/hpalc0b/stroke/paper/2004_programreview/HPFS/DIETS/diets86.dat' lrecl=250;
    input id 1-6 
       @8 dash86q 6.2;
    proc sort nodupkey; by id; run; 

    data dash90hpfs;
       infile '/proj/hpalcs/hpalc0b/stroke/paper/2004_programreview/HPFS/DIETS/diets90.dat' lrecl=250;
    input id 1-6 
       @8 dash90q 6.2;
    proc sort nodupkey; by id; run; 

    data dash94hpfs;
       infile '/proj/hpalcs/hpalc0b/stroke/paper/2004_programreview/HPFS/DIETS/diets94.dat' lrecl=250;
    input id 1-6 
       @8 dash94q 6.2;
    proc sort nodupkey; by id; run; 

    data dash98hpfs;
       infile '/proj/hpalcs/hpalc0b/stroke/paper/2004_programreview/HPFS/DIETS/diets98.dat' lrecl=250;
    input id 1-6 
       @8 dash98q 6.2;
    proc sort nodupkey; by id; run; 

    data dash02hpfs;
       infile '/proj/hpalcs/hpalc0b/stroke/paper/2004_programreview/HPFS/DIETS/diets02.dat' lrecl=250;
    input id 1-6 
       @8 dash02q 6.2;
    proc sort nodupkey; by id; run; 

    data dash06hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/DASHHPFS06' lrecl=150; 
    input id 1-6 
       @8 dash06q 6.2;
    proc sort nodupkey; by id; run; 

    data dash10hpfs;
       infile '/proj/hpbios/hpbio00/hpfs.dietqual/DASHHPFS10' lrecl=150; 
    input id 1-6 
       @8 dash10q 6.2;
    proc sort nodupkey; by id; run;  



   /****** Moderate to vigorous Exercise (to match with NHS where no act80m, only th80 ******/
   %include '/udd/hpypl/review/dietscore/hpfs.exercise.sas'; proc sort; by id; run;

                ******************************************
                *                Covariates              *
                ******************************************;


   /****read in self-reported diabetes, cancer, CVD and hypertension data***/
   %hp_der (keep= dbmy86 dbmy09 rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06 
                  height bmi2186 wt86 wt88 wt90 wt92 wt94 wt96 wt98 wt00 wt02 wt04 wt06 
                  smoke86 smoke88 smoke90 smoke92 smoke94 smoke96 smoke98 smoke00 smoke02 smoke04 smoke06 
                  cgnm86 cgnm88 cgnm90 cgnm92 cgnm94 cgnm96 cgnm98 cgnm00 cgnm02 cgnm04 cgnm06
                  pckyr86 pckyr88 pckyr90 pckyr92 pckyr94 pckyr96 pckyr98 pckyr00 pckyr02 pckyr04 pckyr06
                  fdb3087 mdb3087 sdb3087 bdb3087 fdb2987 mdb2987 sdb2987 bdb2987 dbfh87
                  fclc86 afclc86 mclc86 amclc86 cafh86
                  fclc90 fpro90  mclc90 sclc90 spro90 cafh90
                  profssn);
                  /* Profession: 1.dentist 2.hosp pharmacist 3.optometrist 4.osteopath 5.pharmacist 6.podiatrist 7.vet*/
                  /** family history of db **/
                  dbfh87=0;    * onset after age 30 presumed to be NIDDM;                  
                  if fdb3087=1 or mdb3087=1 or sdb3087=1 or bdb3087=1 then dbfh87=1;
                  /**family history of Cancer***/
                  cafh86=0; if fclc86=1 or mclc86=1 then cafh86=1;
                  cafh90=0;
                  if fclc90=1 or fpro90=1 or  mclc90=1 or sclc90=1 or spro90=1 then cafh90=1;
                  proc freq; tables profssn; run;

   %hp_der_2 (keep= rtmnyr08 wt08 smoke08 cgnm08 rtmnyr10 wt10 smoke10 cgnm10 rtmnyr12 wt12 smoke12 cgnm12 rtmnyr14 wt14 smoke14 cgnm14
                    rtmnyr16 wt16 smoke16 cgnm18 rtmnyr18 wt18 smoke18 cgnm20 rtmnyr20 wt20 smoke20 cgnm20);

   %hp86(keep=hbp86 db86 dbd86 colc86  pros86  lymp86  ocan86  mel86 chol86 trig86 mi86 mid86 ang86 angd86 str86  strd86 cabg86 cabgd86
              seuro86 scand86 ocauc86 afric86 asian86 oanc86 mmi86 fmi86 
              asp86 mar86  livng86
              betab86 lasix86 diur86 calcb86 ald86
              hbpd86 renf86 renfd86 antch86  renf86 pe86 ucol86 asth86 
              prosd86 lympd86  ocand86 colcd86 meld86 ); 
              if db86 ne 1 then db86=0; if mi86 ne 1 then mi86=0; if str86 ne 1 then str86=0;
              if colc86 ne 1 then colc86=0; if mel86 ne 1 then mel86=0; if pros86 ne 1 then pros86=0;
              if lymp86 ne 1 then lymp86=0; if ocan86 ne 1 then ocan86=0;



   %h86_dt(keep=mvt86d mvt86 vite86d vite86);  if mvt86d=0 then mvt86=1; else if mvt86d in (1,2) then mvt86=2; else mvt86=3;
                                               vite86=9; if vite86d=0 or vite86d=1 then vite86=1; else if vite86d=2 then vite86=2;run;  
              /***vite86d             i1             228            Use of Vitamin E (35)
	            $label 0.never;\
	             1.past only;\
	             2.yes current;\
	             9.passthru***/

   %hp87(keep=slp87 snore87);

   %hp88(keep=hbp88 db88 colc88  pros88  lymp88  ocan88  mel88 chol88 trig88 mi88 ang88 str88 cabg88   
              mvt88 asp88 mar88  livng88
              betab88 lasix88 diur88 ccblo88 ald88
              hbpd88  renf88 renfd88 
              pct588 pct1088 pct2088 pct3088 pct4088 pct88 physc88 chrx88 renf88 pe88 ucol88 park88); mvt88=3-mvt88;
 
   %hp90(keep=hbp90 db90 colc90  pros90  lymp90  ocan90  mel90 chol90 trig90 mi90 ang90 str90 cabg90   
              asp90 mar90 livng90 antid90 betab90 lasix90 diur90 ccblo90 ald90
              renfd90  hbpd90
              mdb90 fdb90 sdb90 dbfh90 physc90  chrx90 
              renf90 pe90 ucol90 park90 mdem90 fdem90 sdem90); if fdb90=1 or mdb90=1 or sdb90=1 then dbfh90=1; else dbfh90=0;
 
   %h90_dt(keep=mvt90d mvt90 vite90d vite90);  mvt90=mvt90d;  vite90=9; if vite90d=2 then vite90=2; else if vite90d=1 then vite90=1;run;
                                                              /****vite90d	$range 1-3 ;1 = No 	2 = Yes 	3 = Missing, PASSTHU*****/ 
   %hp92(keep=hbp92 db92 colc92  pros92  lymp92  ocan92  mel92 chol92 trig92 mi92 ang92 str92 cabg92 
              renf92 pe92 ucol92 asth92 pneu92 park92 
              mvt92 asp92 mar92 livng92
              antid92 betab92 lasix92 thiaz92 ccblo92 ald92
              fdb92 mdb92 sdb92 dbfh92
              renf92 renfd92 hbpd92 vite92 
              mlng92 mclc92 flng92 fclc92 fpro92 slng92 sclc92 spro92 cafh92 physx92 chrx92 
              mdem92 fdem92 sdem92
              tbpos92 tbneg92 tbunk92 tbpt92 pre8792 aft8792  tbdpt92 );
              if fdb92=1 or mdb92=1 or sdb92=1 then dbfh92=1; else dbfh92=0; vite92=3-vite92;
              if mlng92=1 or mclc92=1 or flng92=1 or fclc92=1 or fpro92=1 or slng92=1 or sclc92=1 or spro92=1 then cafh92=1; else cafh92=0;

   %hp94(keep=hbp94 db94 colc94  pros94  lymp94  ocan94  mel94 chol94 trig94 mi94 ang94 str94 cabg94 
              pe94 ucol94 pneu94 park94 ms94 ms94s 
              asp94 mar94 livng94
              antid94 betab94 lasix94 thiaz94 ccblo94 ald94 hbpd94 physx94 chrx94 ); if ms94=2 then ms94s=1;


   %h94_dt(keep=mvt94d mvt94 vite94d vite94);mvt94=mvt94d;vite94=9; if vite94d=2 then vite94=2; else if vite94d=1 then vite94=1;

   %hp96(keep=hbp96 db96 colc96  pros96  lymp96  ocan96  mel96 chol96 trig96 mi96 ang96 str96 cabg96 renf96
              pe96 ucol96 asth96 pneu96  emph96 copd96   
              asp96 aspd96 mvt96 mar96  livar96 przc96 tcyc96 antid96 betab96 lasix96 thiaz96 ccblo96 ald96
              hbpd96 renf96 renfd96 vite96
              pclc96 sclc196 sclc296 fpro96 bpro196 bpro296 sbrc196 mbrcn96 cafh96 physx96 chrx96 
              pe96 ucol96 asth96 pneu96  emph96  copd96   renf96 park96 ms96
              brmi96  smi96  bstr96 sstr96); 
              copd96=emph96; if copd96 ne 1 then copd96=0; 
              vite96=3-vite96;
              copd96=emph96; if copd96 ne 1 then copd96=0;
              if aspd96 in (2,3,4,5) then asp96=1; else asp96=0;
              if pclc96=1 or sclc196=1 or sclc296=1 or fpro96=1 or bpro196=1 or bpro296=1 or sbrc196=1 or mbrcn96=1 then cafh96=1; else cafh96=0;

   %hp98(keep=hbp98 db98 colc98  pros98  lymp98  ocan98  mel98 chol98 trig98 mi98 ang98 str98 cabg98 renf98
              pe98 ucol98 asth98 pneu98  emph98 copd98  hbpd98 
              mar98  livng98  przc98 tcyc98 antid98 betab98 lasix98 thiaz98 ccblo98 ald98                
              pe98 ucol98 asth98 pneu98  emph98 copd98 renf98 park98 ms98   chrx98 physx98
              flag98); copd98=emph98; if copd98 ne 1 then copd98=0;

   %h98_dt(keep=mvt98d mvt98 vite98d vite98 melat98d melat98);
              mvt98=mvt98d;vite98=9; if vite98d=2 then vite98=2; else if vite98d=1 then vite98=1; 
              melat98=melat98d; if melat98 ne 1 then melat98=0; /*proc freq; tables melat98; run;*/

   %hp00(keep=hbp00 db00 colc00  pros00  lymp00  ocan00  mel00 chol00 trig00 mi00 ang00 strk00 cabg00
              slp00 snore00 bus00 indrs00 space00 ill00  htsc00 panic00 worry00 out00 
              slp00 snore00 renf00
              pe00 ucol00 asth00 pneu00  emph00 copd00  renf00 park00 als00 ms00 
              asp00 mvt00 mar00 livng00 przc00 tcyc00 antid00 betab00 lasix00 thiaz00 ccblo00 oanth00
              hbpd00 renfd00 vite00 physx00 chrx00 ochrx00 melat00 flag00
              );
              copd00=emph00; if copd00 ne 1 then copd00=0;
              if melat00 ne 1 then melat00=0; /*proc freq; tables melat00; run;*/

   %hp02(keep=hbp02 db02 colc02  pros02  lymp02  ocan02  mel02 chol02 trig02 mi02 ang02 strk02 cabg02 przc02 tcyc02 antid02 ra02 oa02 smk02 ncig02 renal02
              betab02 lasix02 thiaz02 calcb02 anthp02    iron02 asp02 mvt02
              legpn02 rest02 night02
              pe02 ucol02 asth02 pneu02 pern02 copd02 renal02 park02 als02 ms02  
              mar02 lalon02  depr02 
              renld02  hbpd02 vite02 physc02 stat02 ochrx02  melat02 flag02); 
              if copd02 ne 1 then copd02=0;
              if melat02 ne 1 then melat02=0; /*proc freq; tables melat02; run;*/

   %hp04(keep=hbp04 db04 colc04  pros04  lymp04  ocan04  mel04 chol04 trig04 mi04 mi1d04 ang04 angd04 strk04 cabg04 cabgd04 przc04 tcyc04 antid04 ra04 oa04 smk04 ncig04 
              dfslp04 wake04 early04 nap04 restd04 renal04 betab04 ace04 lasix04 thiaz04 calcb04 anthp04 iron04 asp04 mvt04
              pe04 ucol04 asth04 pneu04 pern04  copd04 renal04 park04 als04 ms04 
              hbpd04 renld04 vite04 physc04 depr04  mar04 lalon04
              stat04 mev04 zoc04 crest04 prav04 lip04 ostat04 ochrx04  melat04 flag04
              dfslp04 wake04 early04 nap04 restd04 shing04 shngd04);if copd04 ne 1 then copd04=0;
              if melat04 ne 1 then melat04=0; /*proc freq; tables melat04; run;*/

   %hp06(keep=hbp06 hbpd06 db06 colc06  pros06  lymp06  ocan06  mel06 chol06 trig06 mi06 ang06 angd06 strk06 cabg06 cabgd06 mib06 angd06 przc06 tcyc06 antid06 ra06 oa06 smk06 ncig06 renal06
              betab06 ace06 lasix06 thiaz06 calcb06 anthp06 iron06 asp06 mvt06
              pe06 ucol06 asth06 pneu06 pern06  copd06  renal06 als06 park06 
              hbpd06 renld06 vite06 physc06 mar06 lalon06
              stat06 mev06 zoc06 crest06 prav06 lip06 ostat06 ochrx06  melat06 flag06
              shing06 shingd06);if copd06 ne 1 then copd06=0;
              if melat06 ne 1 then melat06=0; /*proc freq; tables melat06; run;*/

   %hp08(keep=hbp08 hbpd08 db08 colc08  pros08  lymp08  ocan08  mel08 chol08 trig08 mi08 ang08 strk08 cabg08 cabgd08 mi1d08 angd08 ssri08 snri08 tcyc08 maoi08 antid08 
              ra08 oa08 smk08 ncig08 renal08  betab08 ace08 lasix08 thiaz08 calcb08 arb08 anthp08 iron08 asp08 mvt08 wt08
              pe08 ucol08 asth08 pneu08 pern08 copd08 renal08 park08 als08 
              group08 social mar08 lalon08 depr08 ncig08
              hbpd08      renld08
              drt08 vite08 physc08 stat08 mev08 zoc08 crest08 prav08 lip08 ostat08 ochrx08  melat08 flag08
              sleep08 slpmed08 svacc08 shingd08  shing08 );
              if copd08 ne 1 then copd08=0; 
              if group08 in (1,2,3,4,5,6) then social=group08; else social=.;
              if melat08 ne 1 then melat08=0; /*proc freq; tables melat08; run;*/
              /*group08 # hrs/wk in groups (social,work,volunteer,etc)? (L42) $noblank $label 1.None;\ 2.1 - 2 hours;\ 3.3 - 5 hours;\
                        4.6 - 10 hours;\ 5.11 - 15 hours;\  6.16+ hours;\ 7.Passthru */
            

   %hp10(keep=hbp10 hbpd10 db10 chol10 trig10 mi10 ang10 strk10 cabg10 cabgd10  
              physc10 asp10 mvt10 mar10 lalon10
              pros10  lymp10 ocan10 colc10 mel10
              andep10 
              betab10 ace10 lasix10 thiaz10  calcb10  arb10 anthp10
              stat10 mev10 zoc10 crest10 prav10 lip10 ostat10 ochrx10
              pe10 ucol10 asth10 pern10  copd10 renal10 park10 als10  melat10 flag10
              slpmed10);
              if melat10 ne 1 then melat10=0; /*proc freq; tables melat10; run;*/ 

   %hp12(keep=hbp12 hbpd12 db12 chol12 trig12 mi12 ang12 strk12 cabg12 cabgd12  
              physc12   mvt12 mar12 lalon12
              pros12  lymp12 ocan12 colc12 mel12
              andep12 ssri12 
              betab12 ace12 lasix12 thiaz12  calcb12  arb12 anthp12
              stat12 mev12 zoc12 crest12 prav12 lip12 ostat12 ochrx12
              flag12 slphrs12 snore12 slpfal12 slpwak12 slperl12 slpnap12 slprst12 slpmed12 slpmed12yn
              demmo12 demfad12 demsi12);
              if slpmed12 in (1,2,3,4) then slpmed12yn=1; else slpmed12yn=0; 


   %hp14(keep=hbp14 hbpd14 db14 chol14  mi14 ang14 strk14 cabg14 cabgd14 
              physc14 asp14 mvt14 mar14 lalon14
              pros14  lymp14 ocan14 colc14 mel14
              andep14 ssri14
              betab14 ace14 lasix14 thiaz14  calcb14  arb14 anthp14
              stat14 mev14 zoc14 crest14 prav14 lip14 ostat14 ochrx14 
              melat14 flag14 park14);
              if flag14=2 then flag14=1; *only for this study;
              if melat14 ne 1 then melat14=0; /*proc freq; tables melat14; run;*/



   %hp16(keep=cabg16 cabgd16 alz16 alzd16 park16 mar16 asp16 lalon16
              hbp16 hbpd16 db16 chol16  mi16 ang16 stk16
              pros16 lymp16 ocan16 colc16 mel16 kidca16 blad16 depres16 depresd16
              ssris16  tric16  antidep16  thiaz16   lasix16  potas16   cblock16   bblock16  aceinhb16  angio16  antihy16  
              stat16  zocor16   crestor16   prava16  lipit16 othmed16  ochrx16 flag16) ; 
              /* alz16           i1              Alzheimer's or other type of dementia */

   %hp18(keep=cabg18 cabgd18 alz18 alzd18 park18
              mvt18 mar18 asp18 lalon18
              hbp18 hbpd18 db18 chol18  mi18 ang18 strk18
              pros18 lymp18 ocan18 colc18 mel18 kidca18 blad18 depres18 depresd18
              ssris18  tric18  antidep18  thiaz18   lasix18  potas18   cblock18   bblock18  aceinhb18  angio18  antihy18  
              stat18  zocor18   crestor18   prava18  lipit18 othmed18  ochrx18 flag18); 
              /* alz18 Alzheimer's or other type of dementia */

   %hp20(keep=cabg20 cabgd20 alz20 alzd20 park20
              mar20 asp20 lalon20
              hbp20 hbpd20 db20 chol20  mi20 ang20 strk20
              pros20 lymp20 ocan20 colc20 mel20 kidca20 blad20 depres20 depresd20
              ssris20  tric20  antidep20  thiaz20   cblock20   bblock20  aceinhb20  angio20  antihy20  
              stat20  zocor20   crestor20   prava20  lipit20 othmed20  ochrx20 flag20); /* Alzheimer's or other type of dementia */

   %hmet8616(keep=act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16);
           array actm{*} act86 act88 act90 act92 act94 act96 act98 act00 act02 act04 act06 act08 act10 act12 act14 act16;
           do i=1 to dim(actm);
           if actm{i}>250 then actm{i}=.;
           end;



                ******************************************
                *              genetic data              *
                ******************************************; 

proc import datafile='data/genetic/apoe_with_gsa_comb_pcs_after_exclusion.csv'
            out=allnhapoe
            dbms=csv
            replace; 
     getnames=yes;  
run;
/*proc freq; tables APOE; run;*/

    data hpAPOEDATA; 
      set allnhapoe; 
      if Study='HPFS';
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
      proc sort;
      by id;

      /*proc freq;
      tables apoe4grp apoe_2cat APOE4 platform platforms;
      proc means n nmiss min mean std max;
      var PC:;*/
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

data hpPRS;
     set PRS;
     if Study='HPFS';
     id=studyID;  
     proc sort;
     by id; 
      RUN;

data hpPRS2;
     set PRS2;
     if Study='HPFS';
     id=studyID;  
     proc sort;
     by id; 
      RUN;

libname hppdcase '/proj/hppars/hppar0g/case/';
/*****READ IN Parkson Disease before 2012***/

data hppdcase; set hppdcase.hppd12;
keep id conf s1;run;

data hppdcase; set hppdcase;
rename conf=pdconf;
rename s1=pds1;run;
  
/**** neighborhood SES - nSES****/ 
libname curr "";
data hpses; set curr.ahpses8818;  
      proc sort nodupkey; by id; run;

         
                ******************************************
                *               merge data               *
                ******************************************;   

data base;
   merge hp_der(in=mst) diabetes cvd10 hp_der_2 hp_dead 
                        hp86 hp87 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp10 hp12 hp14 hp16 hp18 hp20
                        foodshpfs h86_dt h90_dt h94_dt h98_dt
                        ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10 
                        p86hpfs p90hpfs p94hpfs p98hpfs p02hpfs p06hpfs p10hpfs
                        amed86hpfs amed90hpfs amed94hpfs amed98hpfs amed02hpfs amed06hpfs amed10hpfs
                        dash86hpfs dash90hpfs dash94hpfs dash98hpfs dash02hpfs dash06hpfs dash10hpfs
                        totalcancer hpAPOEDATA  hpPRS hpPRS2 
                        hppdcase hmet8616 pal  hpses 
                        mindhpfs
                        end=_end_;
    by id;
    exrec=1;
    if first.id and mst then exrec=0; 

    birthday=dbmy86;

   can86=0;   if pros86=1 or lymp86=1 or ocan86=1 or colc86=1 or mel86=1 then can86=1;
   can88=0;   if pros88=1 or lymp88=1 or ocan88=1 or colc88=1 or mel88=1 then can88=1;
   can90=0;   if pros90=1 or lymp90=1 or ocan90=1 or colc90=1 or mel90=1 then can90=1;
   can92=0;   if pros92=1 or lymp92=1 or ocan92=1 or colc92=1 or mel92=1 then can92=1;
   can94=0;   if pros94=1 or lymp94=1 or ocan94=1 or colc94=1 or mel94=1 then can94=1;
   can96=0;   if pros96=1 or lymp96=1 or ocan96=1 or colc96=1 or mel96=1 then can96=1;
   can98=0;   if pros98=1 or lymp98=1 or ocan98=1 or colc98=1 or mel98=1 then can98=1;
   can00=0;   if pros00=1 or lymp00=1 or ocan00=1 or colc00=1 or mel00=1 then can00=1;
   can02=0;   if pros02=1 or lymp02=1 or ocan02=1 or colc02=1 or mel02=1 then can02=1;
   can04=0;   if pros04=1 or lymp04=1 or ocan04=1 or colc04=1 or mel04=1 then can04=1;
   can06=0;   if pros06=1 or lymp06=1 or ocan06=1 or colc06=1 or mel06=1 then can06=1;
   can08=0;   if pros08=1 or lymp08=1 or ocan08=1 or colc08=1 or mel08=1 then can08=1;
   can10=0;   if pros10=1 or lymp10=1 or ocan10=1 or colc10=1 or mel10=1 then can10=1;
   can12=0;   if pros12=1 or lymp12=1 or ocan12=1 or colc12=1 or mel12=1 then can12=1;
   can14=0;   if pros14=1 or lymp14=1 or ocan14=1 or colc14=1 or mel14=1 then can14=1;
   can16=0;   if pros16=1 or lymp16=1 or ocan16=1 or colc16=1 or mel16=1 or kidca16=1 or blad16=1 then can16=1;
   can18=0;   if pros18=1 or lymp18=1 or ocan18=1 or colc18=1 or mel18=1 or kidca18=1 or blad18=1 then can18=1;
   can20=0;   if pros20=1 or lymp20=1 or ocan20=1 or colc20=1 or mel20=1 or kidca20=1 or blad20=1 then can20=1;

  /***ethnicity data****/
     if seuro86=1 then ethnic =1;
     else if scand86=1 then ethnic =1;
     else if ocauc86=1 then ethnic =1;
     else if afric86=1 then ethnic =2;
     else if asian86=1 then ethnic =3;
     else if oanc86=1  then ethnic =4;
     else ethnic=.;
     %indic3(vbl=ethnic, prefix=ethnic, reflev=1, min=2, max=4, missing=., usemiss=1,
           label1='causcasian',
           label2='african',
           label3='asian',
           label4='other');
     if ethnic>3 then eth3g=3;
     else             eth3g=ethnic;
     %indic3(vbl=eth3g, prefix=eth3g, reflev=1, min=2, max=3, missing=., usemiss=1,
           label1='causcasian',
           label2='african',
           label3='asian & others');

     if eth3g=1 then white=1; else white=0; if ethnic=. then white=.;

/***feel depressed***/
   if depr02=2 then dep02=1;
   else if depr02=1 then dep02=0;
   else dep02=.;

   if depr04=2 then dep04=1;
   else if depr04=1 then dep04=0;
   else dep04=.;

   if depr08=2 then dep08=1;
   else if depr08=1 then dep08=0;
   else dep08=.;

   if depres16=1 then dep16=1; else dep16=0; /*Depression, clinician-diagnosed (ever) */
   if depres18=1 then dep18=1; else dep18=0;
   if depres20=1 then dep20=1; else dep20=0;

/***antidepressant***/
   if antid90=1 then antd90=1; else antd90=0;
   if antid92=1 then antd92=1; else antd92=0;
   if antid94=1 then antd94=1; else antd94=0;
   if przc96=1 or tcyc96=1 or antid96=1 then antd96=1; else antd96=0;
   if przc98=1 or tcyc98=1 or antid98=1 then antd98=1; else antd98=0;
   if przc00=1 or tcyc00=1 or antid00=1 then antd00=1; else antd00=0;
   if przc02=1 or tcyc02=1 or antid02=1 then antd02=1; else antd02=0;
   if przc04=1 or tcyc04=1 or antid04=1 then antd04=1; else antd04=0; 
   if przc06=1 or tcyc06=1 or antid06=1 then antd06=1; else antd06=0; 
   if ssri08=1 or snri08=1 or tcyc08=1 or maoi08=1 or antid08=1 then antd08=1; else antd08=0;
   if andep10=1 then antd10=1; else antd10=0;
   if ssri12=1 or andep12=1 then antd12=1; else antd12=0;
   if ssri14=1 or andep14=1 then antd14=1; else antd14=0;
   if ssris16=1 or tric16=1 or antidep16=1 then antd16=1; else antd16=0; 
   if ssris18=1 or tric18=1 or antidep18=1 then antd18=1; else antd18=0; 
   if ssris20=1 or tric20=1 or antidep20=1 then antd20=1; else antd20=0; 

/***antihypertensive***/
   if (betab86 =1 or lasix86=1 or diur86 =1 or calcb86=1 or ald86=1) then htnrx86=1;       else htnrx86=0;
   if (betab88 =1 or lasix88=1 or diur88 =1 or ccblo88=1 or ald88=1) then htnrx88=1;       else htnrx88=0;
   if (betab90 =1 or lasix90=1 or diur90 =1 or ccblo90=1 or ald90=1) then htnrx90=1;       else htnrx90=0;
   if (betab92 =1 or lasix92=1 or thiaz92=1 or ccblo92=1 or ald92=1) then htnrx92=1;       else htnrx92=0;
   if (betab94 =1 or lasix94=1 or thiaz94=1 or ccblo94=1 or ald94=1) then htnrx94=1;       else htnrx94=0;
   if (betab96 =1 or lasix96=1 or thiaz96=1 or ccblo96=1 or ald96=1) then htnrx96=1;       else htnrx96=0;
   if (betab98 =1 or lasix98=1 or thiaz98=1 or ccblo98=1 or ald98=1) then htnrx98=1;       else htnrx98=0;
   if (betab00 =1 or lasix00=1 or thiaz00=1 or ccblo00=1 or oanth00=1) then htnrx00=1;       else htnrx00=0;
   if (betab02 =1 or lasix02=1 or thiaz02 =1 or calcb02=1 or anthp02=1) then htnrx02=1;       else htnrx02=0;
   if (betab04 =1 or ace04=1 or lasix04=1 or thiaz04 =1 or calcb04=1 or anthp04=1) then htnrx04=1;       else htnrx04=0; 
   if (betab06 =1 or ace06=1 or lasix06=1 or thiaz06 =1 or calcb06=1 or anthp06=1) then htnrx06=1;       else htnrx06=0; 
   if (betab08 =1 or ace08=1 or lasix08=1 or thiaz08 =1 or calcb08=1 or arb08=1 or anthp08=1) then htnrx08=1;       else htnrx08=0; 
   if (betab10 =1 or ace10=1 or lasix10=1 or thiaz10 =1 or calcb10=1 or arb10=1 or anthp10=1) then htnrx10=1;       else htnrx10=0;
   if (betab12 =1 or ace12=1 or lasix12=1 or thiaz12 =1 or calcb12=1 or arb12=1 or anthp12=1) then htnrx12=1;       else htnrx12=0;
   if (betab14 =1 or ace14=1 or lasix14=1 or thiaz14 =1 or calcb14=1 or arb14=1 or anthp14=1) then htnrx14=1;       else htnrx14=0;
   if (thiaz16 =1 or lasix16=1 or potas16=1 or cblock16=1 or bblock16 =1 or aceinhb16=1 or angio16=1 or antihy16=1) then htnrx16=1;       else htnrx16=0;
   if (thiaz18 =1 or lasix18=1 or potas18=1 or cblock18=1 or bblock18 =1 or aceinhb18=1 or angio18=1 or antihy18=1) then htnrx18=1;       else htnrx18=0;
   if (thiaz20 =1 or cblock20=1 or bblock20 =1 or aceinhb20=1 or angio20=1 or antihy20=1) then htnrx20=1;       else htnrx20=0;

/*antihbc*/
	if antch86=1 then hbcrx86=1; else hbcrx86=0;
	if chrx88=1  then hbcrx88=1; else hbcrx88=0;
	if chrx90=1  then hbcrx90=1; else hbcrx90=0;
	if chrx92=1  then hbcrx92=1; else hbcrx92=0;
	if chrx94=1  then hbcrx94=1; else hbcrx94=0;
	if chrx96=1  then hbcrx96=1; else hbcrx96=0;
	if chrx98=1  then hbcrx98=1; else hbcrx98=0;
	if chrx00=1 or ochrx00=1 then hbcrx00=1; else hbcrx00=0;
	if stat02=1 or ochrx02=1 then hbcrx02=1; else hbcrx02=0; 
	if stat04=1 or mev04=1 or zoc04=1 or crest04=1 or prav04=1 or lip04=1 or ostat04=1 or ochrx04=1 then hbcrx04=1; else hbcrx04=0;
	if stat06=1 or mev06=1 or zoc06=1 or crest06=1 or prav06=1 or lip06=1 or ostat06=1 or ochrx06=1 then hbcrx06=1; else hbcrx06=0;
	if stat08=1 or mev08=1 or zoc08=1 or crest08=1 or prav08=1 or lip08=1 or ostat08=1 or ochrx08=1 then hbcrx08=1; else hbcrx08=0;
	if stat10=1 or mev10=1 or zoc10=1 or crest10=1 or prav10=1 or lip10=1 or ostat10=1 or ochrx10=1 then hbcrx10=1; else hbcrx10=0;
 	if stat12=1 or mev12=1 or zoc12=1 or crest12=1 or prav12=1 or lip12=1 or ostat12=1 or ochrx12=1 then hbcrx12=1; else hbcrx12=0;
 	if stat14=1 or mev14=1 or zoc14=1 or crest14=1 or prav14=1 or lip14=1 or ostat14=1 or ochrx14=1 then hbcrx14=1; else hbcrx14=0;
 	if stat16=1 or zocor16=1 or crestor16=1 or prava16=1 or lipit16=1 or othmed16=1 or ochrx16=1 then hbcrx16=1; else hbcrx16=0;
 	if stat18=1 or zocor18=1 or crestor18=1 or prava18=1 or lipit18=1 or othmed18=1 or ochrx18=1 then hbcrx18=1; else hbcrx18=0;
 	if stat20=1 or zocor20=1 or crestor20=1 or prava20=1 or lipit20=1 or othmed20=1 or ochrx20=1 then hbcrx20=1; else hbcrx20=0;


/* Set up family hx of MI */
   mifh=0;
   if mmi86=1 or fmi86=1 or brmi96=1 or  smi96=1  then mifh=1;

   strfh=0;
   if bstr96 or sstr96=1 then strfh=1;

   cvdfh=0; 
   if mifh=1 or strfh=1 then cvdfh=1;

/* Set up family hx of diabetes */
   famhxdb=0;
   if dbfh87=1 or dbfh90=1 or dbfh92=1 then famhxdb=1;
   parentaldb=0;
   if fdb3087=1 or mdb3087=1 or fdb90=1 or mdb90=1 or fdb92=1 or mdb92=1 then parentaldb=1;
   mdm=0;
   if mdb3087=1 or mdb90=1 or mdb92=1 then mdm=1;
   fdm=0;
   if fdb3087=1 or fdb90=1 or fdb92=1 then fdm=1;  

  /*** Set up family hx of cancer ***/
   if cafh86=1 or cafh90=1 or cafh92=1 or cafh96=1 then famhxca=1; else famhxca=0;
 

  /***sleep duration***/
     if slp00='1' then sleep00c=4.5;
     if slp00='2' then sleep00c=6;
     if slp00='3' then sleep00c=7;
     if slp00='4' then sleep00c=8;
     if slp00='5' then sleep00c=9;
     if slp00='6' then sleep00c=10;
     if slp00='7' then sleep00c=11.5;
     if slp00='8' then sleep00c=.;
     if slp00='A' then sleep00c=5.5;
     if slp00='B' then sleep00c=6.5;
     if slp00='C' then sleep00c=7.5;
     if slp00='D' then sleep00c=8.5;
     if slp00='E' then sleep00c=9.5;
     if slp00='F' then sleep00c=10.5;

     if slp87='1' then sleep87c=4.5;
     if slp87='2' then sleep87c=6;
     if slp87='3' then sleep87c=7;
     if slp87='4' then sleep87c=8;
     if slp87='5' then sleep87c=9;
     if slp87='6' then sleep87c=10;
     if slp87='7' then sleep87c=11.5;
     if slp87='8' then sleep87c=.;
     if slp87='A' then sleep87c=5.5;
     if slp87='B' then sleep87c=6.5;
     if slp87='C' then sleep87c=7.5;
     if slp87='D' then sleep87c=8.5;
     if slp87='E' then sleep87c=9.5;
     if slp87='F' then sleep87c=10.5;

     if sleep08='1' then sleep08c=4.5;
     if sleep08='2' then sleep08c=6;
     if sleep08='3' then sleep08c=7;
     if sleep08='4' then sleep08c=8;
     if sleep08='5' then sleep08c=9;
     if sleep08='6' then sleep08c=10;
     if sleep08='7' then sleep08c=11.5;
     if sleep08='8' then sleep08c=.;
     if sleep08='A' then sleep08c=5.5;
     if sleep08='B' then sleep08c=6.5;
     if sleep08='C' then sleep08c=7.5;
     if sleep08='D' then sleep08c=8.5;
     if sleep08='E' then sleep08c=9.5;
     if sleep08='F' then sleep08c=10.5;

     if slphrs12='1' then sleep12c=4.5;
     if slphrs12='2' then sleep12c=6;
     if slphrs12='3' then sleep12c=7;
     if slphrs12='4' then sleep12c=8;
     if slphrs12='5' then sleep12c=9;
     if slphrs12='6' then sleep12c=10;
     if slphrs12='7' then sleep12c=11.5;
     if slphrs12='8' then sleep12c=.;
     if slphrs12='a' then sleep12c=5.5;
     if slphrs12='b' then sleep12c=6.5;
     if slphrs12='c' then sleep12c=7.5;
     if slphrs12='d' then sleep12c=8.5;
     if slphrs12='e' then sleep12c=9.5;
     if slphrs12='f' then sleep12c=10.5;

     if dfslp04 in (1,2,3) then dfslp04c=3-dfslp04; else dfslp04c=.;
     dfslpbase=dfslp04c; * carry backward as no data available before 04;
     if slpfal12 in (1,2,3) then slpfal12c=3-slpfal12; else slpfal12c=.;

     if snore87=5 then snore87c=0; else if snore87 in (3,4) then snore87c=1; else if snore87 in (1,2) then snore87c=2; else snore87c=.;
     if snore00=5 then snore00c=0; else if snore00 in (3,4) then snore00c=1; else if snore00 in (1,2) then snore00c=2; else snore00c=.;
     if snore12=5 then snore12c=0; else if snore12 in (3,4) then snore12c=1; else if snore12 in (1,2) then snore12c=2; else snore12c=.;

     if melat98 ne 1 then melat98=0;
     if melat00 ne 1 then melat00=0;
     if melat02 ne 1 then melat02=0;
     if melat04 ne 1 then melat04=0;
     if melat06 ne 1 then melat06=0;
     if melat08 ne 1 then melat08=0;
     if melat10 ne 1 then melat10=0;
     if melat14 ne 1 then melat14=0;

     if melat00=1 then melat9800=melat98+melat00;                                                 else melat9800=0; 
     if melat02=1 then melat9802=melat98+melat00+melat02;                                         else melat9802=0; if melat9802>2 then melat9802=2;
     if melat04=1 then melat9804=melat98+melat00+melat02+melat04;                                 else melat9804=0; if melat9804>2 then melat9804=2;     
     if melat06=1 then melat9806=melat98+melat00+melat02+melat04+melat06;                         else melat9806=0; if melat9806>2 then melat9806=2;     
     if melat08=1 then melat9808=melat98+melat00+melat02+melat04+melat06+melat08;                 else melat9808=0; if melat9808>2 then melat9808=2;  
     if melat10=1 then melat9810=melat98+melat00+melat02+melat04+melat06+melat08+melat10;         else melat9810=0; if melat9810>2 then melat9810=2;  
     if melat14=1 then melat9814=melat98+melat00+melat02+melat04+melat06+melat08+melat10+melat14; else melat9814=0; if melat9814>2 then melat9814=2;
        
     if cabg86=1 or cabg88=1 or cabg90=1 or cabg92=1 or cabg94=1 or cabg96=1 or cabg98=1 then basecabg=1; else basecabg=0;


       /* cabg update */
        if dtdxcvd=. then do;
                              if cabg20=1 then do;
                                                 if cabgd20=1 then dtdxcvd=1410;
                                            else if cabgd20=2 then dtdxcvd=1422;
                                            else if cabgd20=3 then dtdxcvd=1434;
                                            else if cabgd20=4 then dtdxcvd=1441;
                                            else                   dtdxcvd=1434;
                                                end;

                              if cabg18=1 then do;
                                                 if cabgd18=1 then dtdxcvd=1386;
                                            else if cabgd18=2 then dtdxcvd=1398;
                                            else if cabgd18=3 then dtdxcvd=1410;
                                            else if cabgd18=4 then dtdxcvd=1417;
                                            else                   dtdxcvd=1410;
                                                end;

                              if cabg16=1 then do;
                                                 if cabgd16=1 then dtdxcvd=1362;
                                            else if cabgd16=2 then dtdxcvd=1374;
                                            else if cabgd16=3 then dtdxcvd=1386;
                                            else if cabgd16=4 then dtdxcvd=1393;
                                            else                   dtdxcvd=1386;
                                                end;
                            end;

       /* diagnosis time of ADRD */

                              if alz20=1 then do; /* 1.Before 2018;\2.2018;\3.2019;\4.2020;\5.Passthru */
                                                 if alzd20=1 then dtdx_alzd=1410;
                                            else if alzd20=2 then dtdx_alzd=1422;
                                            else if alzd20=3 then dtdx_alzd=1434;
                                            else if alzd20=4 then dtdx_alzd=1441;
                                            else                  dtdx_alzd=1434;
                                                end;

                              if alz18=1 then do; /* 1.Before 2016;\	2.2016;\	3.2017;\	4.2018;\ 	5.Passthru*/

                                                 if alzd18=1 then dtdx_alzd=1386;
                                            else if alzd18=2 then dtdx_alzd=1398;
                                            else if alzd18=3 then dtdx_alzd=1410;
                                            else if alzd18=4 then dtdx_alzd=1417;
                                            else                  dtdx_alzd=1410;
                                                end;

                              if alz16=1 then do; /*$label 1.Before 2014;\2.2014;\3.2015;\  4.2016;\  5.Passthru */
                                                 if alzd16=1 then dtdx_alzd=1362;
                                            else if alzd16=2 then dtdx_alzd=1374;
                                            else if alzd16=3 then dtdx_alzd=1386;
                                            else if alzd16=4 then dtdx_alzd=1393;
                                            else                  dtdx_alzd=1386;
                                                end;


     if  alz16=1  or alz18=1  or alz20=1  then alz1620=1;  else alz1620=0; 

      /*** mortality due to dementia or AD***/     
          if 0<dtdth<9999 and (newicda = 290 OR newicda = 331) then dementiadeath=1;  else dementiadeath=0;

          if alz1620=1 or dementiadeath=1 then adrd=1; else adrd=0;

          if alz1620=1 and dementiadeath=1 and 0<dtdth<dtdx_alzd<9999 then dtdth=dtdx_alzd; /*make sure death did not occur before date of diagnosis*/

      /*** diagnosis time of AD ***/
          if alz1620=1 and dementiadeath=1 then dtdx_ad=min(dtdth,dtdx_alzd);
          else if alz1620=1 then dtdx_ad=dtdx_alzd;
          else if dementiadeath=1 then dtdx_ad=dtdth;
          else dtdx_ad=.;

          if dtdx_ad ne .  and dtdx_ad=dtdx_alzd then nfad=1; else nfad=0;

          if dtdx_ad ne . then prevalent_ad=1; else prevalent_ad=0;

       /*** age of AD diagnosis ***/
          if dtdx_ad ne . and birthday ne . then ADdiagage=int((dtdx_ad-birthday)/12);      
            if ADdiagage =< 0 then ADdiagage=.;


      pd_self=0;
      if park14=1 or park16=1  or park18=1 or park20=1 then pd_self=1; 

      if pdconf=11 & pds1 in (1,2) then pd=1;
            else if pd_self=1 then pd=1;
            else pd=0;

  /***** Calculate cumulative averages of aMed *****/

     array emedorg{9} emed86h emed90h emed94h emed98h emed02h emed06h emed10h emed14h emed18h; 
     array emedcum{9} emed86v emed90v emed94v emed98v emed02v emed06v emed10v emed14v emed18v; 
      sumvar=0;    n=0;
      do i=1 to 9; 
         if emedorg{i} ne . then do;
         n=n+1; sumvar=sumvar+emedorg{i};
         end; 
      if n ne 0 then emedcum{i}=sumvar/n; else emedcum{i}=.;
      end;

      * calculate weight changes before missing carry forward and readi for missing weight changes carry forward;

      array wt2      {18}    wt86     wt88     wt90     wt92     wt94     wt96     wt98     wt00     wt02     wt04     wt06      wt08      wt10      wt12      wt14      wt16      wt18       wt20       ;
      array wt1      {18}    wt86     wt86     wt88     wt90     wt92     wt94     wt96     wt98     wt00     wt02     wt04      wt06      wt08      wt10      wt12      wt14      wt16       wt18       ;
      array wt21     {18}    wtchg86  wtchg88  wtchg90  wtchg92  wtchg94  wtchg96  wtchg98  wtchg00  wtchg02  wtchg04  wtchg06   wtchg08   wtchg10   wtchg12   wtchg14   wtchg16   wtchg18    wtchg20    ;

      array hbp     {18}    hbp86    hbp88    hbp90    hbp92    hbp94    hbp96    hbp98    hbp00    hbp02    hbp04    hbp06     hbp08     hbp10    hbp12     hbp14     hbp16     hbp18      hbp20      ;

      do i=1 to 18; 
         if wt1{i}=0 then wt1{i}=.; 
         if wt2{i}=0 then wt2{i}=.; 
         wt21{i}=wt2{i}-wt1{i};
         if hbp{i} ne 1 then hbp{i}=0;
      end; 

      if tbpos92=1 then tbcat=3; else if tbneg92=1 then  tbcat=1;  else tbcat=2;  
      %indic3(vbl=tbcat, prefix=tbcat, min=1, max=3, reflev=1, missing=., usemiss=0,
                label1='TB negative',
                label2='non-response',
      	    label3='TB positive');

      if shing08=1 or shing06=1 or shing04=1 then shingle=1; else shingle=0;
      if shingle=1 then shing2vac=1; 
              else if shingle=0 and svacc08=2 then shing2vac=2;
              else if shingle=0 and svacc08=1 then shing2vac=3;
              else shing2vac=.; 

      if shingle=1 and svacc08=2 then shingesvac=1; 
              else if shingle=1 and svacc08=1 then shingesvac=2;
              else if shingle=0 and svacc08=2 then shingesvac=3;
              else if shingle=0 and svacc08=1 then shingesvac=4;
              else shingesvac=.; 


proc freq;
tables tbcat shing2vac;
run; 


  /*   proc freq;
   tables alz1620 dementiadeath adrd nfad prevalent_ad pd park14 park16 park18 park20;
   run;  
 
  proc means nolabels;    run; 

  * time trend for quality control of the emed score coding;
  proc freq;
  tables emed:  alcoh2:  emsrat:  efsh:  emt:  ewgr:  eleg:  efru: enut: evegg:;
  RUN;

  proc means n nmiss min mean std median max;
  var emed86h emed90h emed94h emed98h emed02h emed06h emed10h emed14h emed18h ;
  RUN;*/  

proc datasets nolist;
   delete hp_der diabetes cvd10 hp_der_2 hp_dead 
                        hp86 hp87 hp88 hp90 hp92 hp94 hp96 hp98 hp00 hp02 hp04 hp06 hp08 hp10 hp12 hp14 hp16 hp18 hp20
                        foodshpfs h86_dt h90_dt h94_dt h98_dt
                        ahei86 ahei90 ahei94 ahei98 ahei02 ahei06 ahei10 
                        p86hpfs p90hpfs p94hpfs p98hpfs p02hpfs p06hpfs p10hpfs
                        amed86hpfs amed90hpfs amed94hpfs amed98hpfs amed02hpfs amed06hpfs amed10hpfs
                        dash86hpfs dash90hpfs dash94hpfs dash98hpfs dash02hpfs dash06hpfs dash10hpfs
                        totalcancer hpAPOEDATA  hpPRS hpPRS2 hppdcase hmet8616 pal  hpses;
   run;

data pre_pm;
     set base  end=_end_;
 
   if chol86 ne 1 then chol86=0;
   if hbp86  ne 1 then hbp86=0;
   if white  ne 1 then white=0;

   %indic3(vbl=eth3g, prefix=eth3g, reflev=1, min=2, max=3, missing=., usemiss=1,
           label1='causcasian',
           label2='african',
           label3='asian & others');

   depr86=0;

   /* current dietary patterns are just carry forwaed by 2010, will update later */

   array rtmnyr  {19}    rtmnyr86 rtmnyr88 rtmnyr90 rtmnyr92 rtmnyr94 rtmnyr96 rtmnyr98 rtmnyr00 rtmnyr02 rtmnyr04 rtmnyr06  rtmnyr08  rtmnyr10 rtmnyr12  rtmnyr14  rtmnyr16  rtmnyr18   rtmnyr20   cutoff; 
   array irt     {19}    irt86    irt88    irt90    irt92    irt94    irt96    irt98    irt00    irt02    irt04    irt06     irt08     irt10    irt12     irt14     irt16     irt18      irt20      cutoff;
   array tvar    {18}    t86      t88      t90      t92      t94      t96      t98      t00      t02      t04      t06       t08       t10      t12       t14       t16       t18        t20        ;
   array age     {18}    age86    age88    age90    age92    age94    age96    age98    age00    age02    age04    age06     age08     age10    age12     age14     age16     age18      age20      ;
   array mvyn    {18}    mvt86    mvt88    mvt90    mvt92    mvt94    mvt96    mvt98    mvt00    mvt02    mvt04    mvt06     mvt08     mvt10    mvt12     mvt14     mvt14     mvt18      mvt18      ;
   array marriage{18}   mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86 mar86;
   array aspi    {18}    asp86    asp88    asp90    asp92    asp94    asp96    asp96    asp00    asp02    asp04    asp06     asp08     asp10    asp12     asp14     asp16     asp18      asp20      ;
   array alone    {18} livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86 livng86;
   array hbp     {18} hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86 hbp86;
   array chol    {18}    chol86   chol88   chol90   chol92   chol94   chol96   chol98   chol00   chol02   chol04   chol06    chol08    chol10   chol10    chol14    chol16    chol18     chol20     ;
   array antitc  {18}    hbcrx86  hbcrx88  hbcrx90  hbcrx92  hbcrx94  hbcrx96  hbcrx98  hbcrx00  hbcrx02  hbcrx04  hbcrx06   hbcrx08   hbcrx10  hbcrx12   hbcrx14   hbcrx16   hbcrx18    hbcrx20    ;
   array tg      {18}    trig86   trig88   trig90   trig92   trig94   trig96   trig98   trig00   trig02   trig04   trig06    trig08    trig10   trig12    trig12    trig12    trig12     trig12     ;
   array diab    {18} db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86 db86;
   array repmi   {18}    mi86     mi88     mi90     mi92     mi94     mi96     mi98     mi00     mi02     mi04     mi06      mi08      mi10     mi12      mi14      mi16      mi18       mi20       ;
   array repstrk {18}    str86    str88    str90    str92    str94    str96    str98    strk00   strk02   strk04   strk06    strk08    strk10   strk12    strk14    stk16     strk18     strk20     ;
   array rpang   {18}    ang86    ang88    ang90    ang92    ang94    ang96    ang98    ang00    ang02    ang04    ang06     ang08     ang10    ang12     ang14     ang16     ang18      ang20      ;
   array rpcabg  {18}    cabg86   cabg88   cabg90   cabg92   cabg94   cabg96   cabg98   cabg00   cabg02   cabg04   cabg06    cabg08    cabg10   cabg12    cabg14    cabg16    cabg18     cabg20     ;
   array canc    {18}    can86    can88    can90    can92    can94    can96    can98    can00    can02    can04    can06     can08     can10    can12     can14     can16     can18      can20      ;
   array antd    {18}    antd90   antd90   antd90   antd92   antd94   antd96   antd98   antd00   antd02   antd04   antd06    antd08    antd10   antd12    antd14    antd16    antd18     antd20     ;
   array antihp {18} htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86 htnr86;
   array depr    {18} depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86 depr86;
   array weight  {18}    wt86     wt88     wt90     wt92     wt94     wt96     wt98     wt00     wt02     wt04     wt06      wt08      wt10     wt12      wt14      wt16      wt18       wt20       ;
   array bmi     {18}    bmi86    bmi86    bmi86    bmi86    bmi86    bmi86    bmi86    bmi86    bmi86    bmi86    bmi86     bmi86     bmi86    bmi86     bmi86     bmi86     bmi86      bmi86      ;  
   array smok    {18}    smoke86  smoke86  smoke86  smoke86  smoke86  smoke86  smoke86  smoke86  smoke86  smoke86  smoke86   smoke86   smoke86  smoke86   smoke86   smoke86   smoke86    smoke86    ;
   array cignm   {18}    cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86 cgnm86;
   array smkm    {18}    smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86 smkm86;
   array smm     {18}    smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86  smm86;  
   array alcon   {18}    alco86n  alco86n  alco90n  alco90n  alco94n  alco94n  alco98n  alco98n  alco02n  alco02n  alco06n   alco06n   alco10n  alco10n   alco14n   alco14n   alco18n    alco18n    ;
   array actm    {18}    th86     th86     th86     th86     th86     th86     th86     th86     th86     th86     th86      th86      th86     th86      th86      th86      th86       th86       ;  
   array aheim   {18}    ahei_nomt_86   ahei_nomt_86   ahei_nomt_90   ahei_nomt_90   ahei_nomt_94   ahei_nomt_94   ahei_nomt_98   ahei_nomt_98   ahei_nomt_02   ahei_nomt_02   ahei_nomt_06    ahei_nomt_06    ahei_nomt_10   
                         ahei_nomt_10   ahei_nomt_14   ahei_nomt_14   ahei_nomt_14   ahei_nomt_14     ;
   array dashm   {18}    dash86q  dash86q  dash90q  dash90q  dash94q  dash94q  dash98q  dash98q  dash02q  dash02q  dash06q   dash06q   dash10q  dash10q   dash10q   dash10q   dash10q    dash10q    ;
   array pdim    {18}    pdi86    pdi86    pdi90    pdi90    pdi94    pdi94    pdi98    pdi98    pdi02    pdi02    pdi06     pdi06     pdi10    pdi10     pdi10     pdi10     pdi10      pdi10      ;
   array hpdim   {18}    hpdi86   hpdi86   hpdi90   hpdi90   hpdi94   hpdi94   hpdi98   hpdi98   hpdi02   hpdi02   hpdi06    hpdi06    hpdi10   hpdi10    hpdi10    hpdi10    hpdi10     hpdi10     ;
   array updim   {18}    updi86   updi86   updi90   updi90   updi94   updi94   updi98   updi98   updi02   updi02   updi06    updi06    updi10   updi10    updi10    updi10    updi10     updi10     ;

   array calor    {18}  calor86n  calor86n    calor90n   calor90n    calor94n   calor94n   calor98n   calor98n   calor02n   calor02n   calor06n   calor06n   calor10n   calor10n   calor14n   calor14n  calor18n  calor18n;
   array trmeatcon{18} trmeat86d  trmeat86d   trmeat90d  trmeat90d   trmeat94d  trmeat94d  trmeat98d  trmeat98d  trmeat02d  trmeat02d  trmeat06d  trmeat06d  trmeat10d  trmeat10d  trmeat14d  trmeat14d trmeat18d trmeat18d; 
   array prmeatcon{18} prmeat86d  prmeat86d   prmeat90d  prmeat90d   prmeat94d  prmeat94d  prmeat98d  prmeat98d  prmeat02d  prmeat02d  prmeat06d  prmeat06d  prmeat10d  prmeat10d  prmeat14d  prmeat14d prmeat18d prmeat18d; 
   array urmeatcon{18} urmeat86d  urmeat86d   urmeat90d  urmeat90d   urmeat94d  urmeat94d  urmeat98d  urmeat98d  urmeat02d  urmeat02d  urmeat06d  urmeat06d  urmeat10d  urmeat10d  urmeat14d  urmeat14d urmeat18d urmeat18d; 

   array fruitscon{18} fruits86d  fruits86d   fruits90d  fruits90d   fruits94d  fruits94d  fruits98d  fruits98d  fruits02d  fruits02d  fruits06d  fruits06d  fruits10d  fruits10d  fruits14d  fruits14d fruits18d fruits18d; 
   array vegealcon{18} vegeal86d  vegeal86d   vegeal90d  vegeal90d   vegeal94d  vegeal94d  vegeal98d  vegeal98d  vegeal02d  vegeal02d  vegeal06d  vegeal06d  vegeal10d  vegeal10d  vegeal14d  vegeal14d vegeal18d vegeal18d; 
   array whgrnscon{18} whgrns86d  whgrns86d   whgrns90d  whgrns90d   whgrns94d  whgrns94d  whgrns98d  whgrns98d  whgrns02d  whgrns02d  whgrns06d  whgrns06d  whgrns10d  whgrns10d  whgrns14d  whgrns14d whgrns18d whgrns18d; 
   array poultrcon{18} poultr86d  poultr86d   poultr90d  poultr90d   poultr94d  poultr94d  poultr98d  poultr98d  poultr02d  poultr02d  poultr06d  poultr06d  poultr10d  poultr10d  poultr14d  poultr14d poultr18d poultr18d; 
   array fishalcon{18} fishal86d  fishal86d   fishal90d  fishal90d   fishal94d  fishal94d  fishal98d  fishal98d  fishal02d  fishal02d  fishal06d  fishal06d  fishal10d  fishal10d  fishal14d  fishal14d fishal18d fishal18d; 
   array regeggcon{18} regegg86d  regegg86d   regegg90d  regegg90d   regegg94d  regegg94d  regegg98d  regegg98d  regegg02d  regegg02d  regegg06d  regegg06d  regegg10d  regegg10d  regegg14d  regegg14d regegg18d regegg18d; 
   array nutsalcon{18} nutsal86d  nutsal86d   nutsal90d  nutsal90d   nutsal94d  nutsal94d  nutsal98d  nutsal98d  nutsal02d  nutsal02d  nutsal06d  nutsal06d  nutsal10d  nutsal10d  nutsal14d  nutsal14d nutsal18d nutsal18d; 
   array nutlegcon{18} nutleg86d  nutleg86d   nutleg90d  nutleg90d   nutleg94d  nutleg94d  nutleg98d  nutleg98d  nutleg02d  nutleg02d  nutleg06d  nutleg06d  nutleg10d  nutleg10d  nutleg14d  nutleg14d nutleg18d nutleg18d; 
   array tdairycon{18} tdairy86d  tdairy86d   tdairy90d  tdairy90d   tdairy94d  tdairy94d  tdairy98d  tdairy98d  tdairy02d  tdairy02d  tdairy06d  tdairy06d  tdairy10d  tdairy10d  tdairy14d  tdairy14d tdairy18d tdairy18d; 
   array ldairycon{18} ldairy86d  ldairy86d   ldairy90d  ldairy90d   ldairy94d  ldairy94d  ldairy98d  ldairy98d  ldairy02d  ldairy02d  ldairy06d  ldairy06d  ldairy10d  ldairy10d  ldairy14d  ldairy14d ldairy18d ldairy18d; 
   array hdairycon{18} hdairy86d  hdairy86d   hdairy90d  hdairy90d   hdairy94d  hdairy94d  hdairy98d  hdairy98d  hdairy02d  hdairy02d  hdairy06d  hdairy06d  hdairy10d  hdairy10d  hdairy14d  hdairy14d hdairy18d hdairy18d; 
   array legumecon{18} legume86d  legume86d   legume90d  legume90d   legume94d  legume94d  legume98d  legume98d  legume02d  legume02d  legume06d  legume06d  legume10d  legume10d  legume14d  legume14d legume18d legume18d; 
   array soyprocon{18} soypro86d  soypro86d   soypro90d  soypro90d   soypro94d  soypro94d  soypro98d  soypro98d  soypro02d  soypro02d  soypro06d  soypro06d  soypro10d  soypro10d  soypro14d  soypro14d soypro18d soypro18d; 
   array ssbcon   {18} ssb86      ssb86       ssb90      ssb90       ssb94      ssb94      ssb98      ssb98      ssb02      ssb02      ssb06      ssb06      ssb10      ssb10      ssb14      ssb14     ssb18     ssb18;

   array emedm   {18}    emed86h  emed86h  emed90h  emed90h  emed94h  emed94h  emed98h  emed98h  emed02h  emed02h  emed06h   emed06h   emed10h  emed10h   emed14h   emed14h   emed20h    emed20h    ;
   array emedv   {18}    emed86v  emed86v  emed90v  emed90v  emed94v  emed94v  emed98v  emed98v  emed02v  emed02v  emed06v   emed06v   emed10v  emed10v   emed14v   emed14v   emed20v    emed20v    ;
   array mindm   {18}    MIND86  MIND86  MIND90  MIND90  MIND94  MIND94  MIND98  MIND98  MIND02  MIND02  MIND06   MIND06   MIND10  MIND10   MIND14   MIND14   MIND18    MIND18    ;
   array wtchg   {18}    wtchg86  wtchg88  wtchg90  wtchg92  wtchg94  wtchg96  wtchg98  wtchg00  wtchg02  wtchg04  wtchg06   wtchg08   wtchg10  wtchg12   wtchg14   wtchg16   wtchg18    wtchg20    ;
   array SES {18} nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88 nSES_88;
  array emedm_alco   {18}    alcoh286  alcoh286  alcoh290  alcoh290  alcoh294  alcoh294  alcoh298  alcoh298  alcoh202  alcoh202  alcoh206   alcoh206   alcoh210  alcoh210   alcoh214   alcoh214   alcoh220    alcoh220    ;
  array emedm_emsrat {18}    emsrat86h  emsrat86h  emsrat90h  emsrat90h  emsrat94h  emsrat94h  emsrat98h  emsrat98h  emsrat02h  emsrat02h  emsrat06h   emsrat06h   emsrat10h  emsrat10h   emsrat14h   emsrat14h   emsrat20h    emsrat20h    ;
  array emedm_efsh   {18}    efsh86h  efsh86h  efsh90h  efsh90h  efsh94h  efsh94h  efsh98h  efsh98h  efsh02h  efsh02h  efsh06h   efsh06h   efsh10h  efsh10h   efsh14h   efsh14h   efsh20h    efsh20h    ;
  array emedm_emt    {18}    emt86h  emt86h  emt90h  emt90h  emt94h  emt94h  emt98h  emt98h  emt02h  emt02h  emt06h   emt06h   emt10h  emt10h   emt14h   emt14h   emt20h    emt20h    ;
  array emedm_ewgr   {18}    ewgr86h  ewgr86h  ewgr90h  ewgr90h  ewgr94h  ewgr94h  ewgr98h  ewgr98h  ewgr02h  ewgr02h  ewgr06h   ewgr06h   ewgr10h  ewgr10h   ewgr14h   ewgr14h   ewgr20h    ewgr20h    ;
  array emedm_eleg   {18}    eleg86h  eleg86h  eleg90h  eleg90h  eleg94h  eleg94h  eleg98h  eleg98h  eleg02h  eleg02h  eleg06h   eleg06h   eleg10h  eleg10h   eleg14h   eleg14h   eleg20h    eleg20h    ;
  array emedm_efru   {18}    efru86h  efru86h  efru90h  efru90h  efru94h  efru94h  efru98h  efru98h  efru02h  efru02h  efru06h   efru06h   efru10h  efru10h   efru14h   efru14h   efru20h    efru20h    ;
  array emedm_enut   {18}    enut86h  enut86h  enut90h  enut90h  enut94h  enut94h  enut98h  enut98h  enut02h  enut02h  enut06h   enut06h   enut10h  enut10h   enut14h   enut14h   enut20h    enut20h    ;
  array emedm_evegg  {18}    evegg86h  evegg86h  evegg90h  evegg90h  evegg94h  evegg94h  evegg98h  evegg98h  evegg02h  evegg02h  evegg06h   evegg06h   evegg10h  evegg10h   evegg14h   evegg14h   evegg20h    evegg20h    ;

  array mind_whsvg   {18}  whsvg86  whsvg86  whsvg90  whsvg90  whsvg94  whsvg94  whsvg98  whsvg98  whsvg02  whsvg02  whsvg06  whsvg06  whsvg10  whsvg10  whsvg14  whsvg14  whsvg18  whsvg18;
  array mind_GVsvg   {18}  GVsvg86  GVsvg86  GVsvg90  GVsvg90  GVsvg94  GVsvg94  GVsvg98  GVsvg98  GVsvg02  GVsvg02  GVsvg06  GVsvg06  GVsvg10  GVsvg10  GVsvg14  GVsvg14  GVsvg18  GVsvg18;
  array mind_vgsvg   {18}  vgsvg86  vgsvg86  vgsvg90  vgsvg90  vgsvg94  vgsvg94  vgsvg98  vgsvg98  vgsvg02  vgsvg02  vgsvg06  vgsvg06  vgsvg10  vgsvg10  vgsvg14  vgsvg14  vgsvg18  vgsvg18;
  array mind_brsvg   {18}  brsvg86  brsvg86  brsvg90  brsvg90  brsvg94  brsvg94  brsvg98  brsvg98  brsvg02  brsvg02  brsvg06  brsvg06  brsvg10  brsvg10  brsvg14  brsvg14  brsvg18  brsvg18;
  array mind_mtsvg   {18}  mtsvg86  mtsvg86  mtsvg90  mtsvg90  mtsvg94  mtsvg94  mtsvg98  mtsvg98  mtsvg02  mtsvg02  mtsvg06  mtsvg06  mtsvg10  mtsvg10  mtsvg14  mtsvg14  mtsvg18  mtsvg18;
  array mind_fssvg   {18}  fssvg86  fssvg86  fssvg90  fssvg90  fssvg94  fssvg94  fssvg98  fssvg98  fssvg02  fssvg02  fssvg06  fssvg06  fssvg10  fssvg10  fssvg14  fssvg14  fssvg18  fssvg18;
  array mind_ptsvg   {18}  ptsvg86  ptsvg86  ptsvg90  ptsvg90  ptsvg94  ptsvg94  ptsvg98  ptsvg98  ptsvg02  ptsvg02  ptsvg06  ptsvg06  ptsvg10  ptsvg10  ptsvg14  ptsvg14  ptsvg18  ptsvg18;
  array mind_bnsvg   {18}  bnsvg86  bnsvg86  bnsvg90  bnsvg90  bnsvg94  bnsvg94  bnsvg98  bnsvg98  bnsvg02  bnsvg02  bnsvg06  bnsvg06  bnsvg10  bnsvg10  bnsvg14  bnsvg14  bnsvg18  bnsvg18;
  array mind_ntsvg   {18}  ntsvg86  ntsvg86  ntsvg90  ntsvg90  ntsvg94  ntsvg94  ntsvg98  ntsvg98  ntsvg02  ntsvg02  ntsvg06  ntsvg06  ntsvg10  ntsvg10  ntsvg14  ntsvg14  ntsvg18  ntsvg18;
  array mind_Ffsvg   {18}  Ffsvg86  Ffsvg86  Ffsvg90  Ffsvg90  Ffsvg94  Ffsvg94  Ffsvg98  Ffsvg98  Ffsvg02  Ffsvg02  Ffsvg06  Ffsvg06  Ffsvg10  Ffsvg10  Ffsvg14  Ffsvg14  Ffsvg18  Ffsvg18;
  array mind_bu_amt_tot {18}  bu_amt_tot86  bu_amt_tot86  bu_amt_tot90  bu_amt_tot90  bu_amt_tot94  bu_amt_tot94  bu_amt_tot98  bu_amt_tot98  bu_amt_tot02  bu_amt_tot02  bu_amt_tot06  bu_amt_tot06  bu_amt_tot10  bu_amt_tot10  bu_amt_tot14  bu_amt_tot14  bu_amt_tot18  bu_amt_tot18;
  array mind_o_amt_tot  {18}  o_amt_tot86  o_amt_tot86  o_amt_tot90  o_amt_tot90  o_amt_tot94  o_amt_tot94  o_amt_tot98  o_amt_tot98  o_amt_tot02  o_amt_tot02  o_amt_tot06  o_amt_tot06  o_amt_tot10  o_amt_tot10  o_amt_tot14  o_amt_tot14  o_amt_tot18  o_amt_tot18;
  array mind_CHsvg   {18}  CHsvg86  CHsvg86  CHsvg90  CHsvg90  CHsvg94  CHsvg94  CHsvg98  CHsvg98  CHsvg02  CHsvg02  CHsvg06  CHsvg06  CHsvg10  CHsvg10  CHsvg14  CHsvg14  CHsvg18  CHsvg18;
  array mind_SWsvg   {18}  SWsvg86  SWsvg86  SWsvg90  SWsvg90  SWsvg94  SWsvg94  SWsvg98  SWsvg98  SWsvg02  SWsvg02  SWsvg06  SWsvg06  SWsvg10  SWsvg10  SWsvg14  SWsvg14  SWsvg18  SWsvg18;
  array mind_alsvg   {18}  alsvg86  alsvg86  alsvg90  alsvg90  alsvg94  alsvg94  alsvg98  alsvg98  alsvg02  alsvg02  alsvg06  alsvg06  alsvg10  alsvg10  alsvg14  alsvg14  alsvg18  alsvg18;

   /****** rtmnyr to irt ******/
    do c=1 to 18;
     irt{c}=rtmnyr{c};
    end;
    drop c;

   /****** missing replace ******/
    do b=2 to 18;

       if actm{b}=.     then actm{b}=actm{b-1};
       if weight{b}=.   then weight{b}=weight{b-1};
       if wtchg{b}=.    then wtchg{b}=wtchg{b-1};
       if mvyn{b}=.     then mvyn{b}= mvyn{b-1};
       if aspi{b}=.     then aspi{b}= aspi{b-1};
       if calor{b}=.    then calor{b}= calor{b-1};

       if alcon{b}=.    then alcon{b}=alcon{b-1}; 
       if smm{b}=.      then smm{b}=smm{b-1};
       if aheim{b}=.    then aheim{b}=aheim{b-1}; 
       if dashm{b}=.    then dashm{b}=dashm{b-1}; 
       if pdim{b}=.     then pdim{b}=pdim{b-1}; 
       if hpdim{b}=.    then hpdim{b}=hpdim{b-1}; 
       if updim{b}=.    then updim{b}=updim{b-1};  
       if emedm{b}=.    then emedm{b}=emedm{b-1}; 
       if emedv{b}=.    then emedv{b}=emedv{b-1};  
       if mindm{b}=.    then mindm{b}=mindm{b-1}; 

       if emedm_alco{b}=.    then emedm_alco{b}=emedm_alco{b-1}; 
       if emedm_emsrat{b}=.    then emedm_emsrat{b}=emedm_emsrat{b-1}; 
       if emedm_efsh{b}=.    then emedm_efsh{b}=emedm_efsh{b-1}; 
       if emedm_emt{b}=.    then emedm_emt{b}=emedm_emt{b-1}; 
       if emedm_ewgr{b}=.    then emedm_ewgr{b}=emedm_ewgr{b-1}; 
       if emedm_eleg{b}=.    then emedm_eleg{b}=emedm_eleg{b-1}; 
       if emedm_efru{b}=.    then emedm_efru{b}=emedm_efru{b-1}; 
       if emedm_enut{b}=.    then emedm_enut{b}=emedm_enut{b-1}; 
       if emedm_evegg{b}=.    then emedm_evegg{b}=emedm_evegg{b-1}; 

        if mind_whsvg{b} =. then mind_whsvg{b} = mind_whsvg{b-1}; 
        if mind_GVsvg{b} =. then mind_GVsvg{b} = mind_GVsvg{b-1}; 
        if mind_vgsvg{b} =. then mind_vgsvg{b} = mind_vgsvg{b-1}; 
        if mind_brsvg{b} =. then mind_brsvg{b} = mind_brsvg{b-1}; 
        if mind_mtsvg{b} =. then mind_mtsvg{b} = mind_mtsvg{b-1}; 
        if mind_fssvg{b} =. then mind_fssvg{b} = mind_fssvg{b-1}; 
        if mind_ptsvg{b} =. then mind_ptsvg{b} = mind_ptsvg{b-1}; 
        if mind_bnsvg{b} =. then mind_bnsvg{b} = mind_bnsvg{b-1}; 
        if mind_ntsvg{b} =. then mind_ntsvg{b} = mind_ntsvg{b-1}; 
        if mind_Ffsvg{b} =. then mind_Ffsvg{b} = mind_Ffsvg{b-1}; 
        if mind_bu_amt_tot{b} =. then mind_bu_amt_tot{b} = mind_bu_amt_tot{b-1}; 
        if mind_o_amt_tot{b} =. then mind_o_amt_tot{b} = mind_o_amt_tot{b-1}; 
        if mind_CHsvg{b} =. then mind_CHsvg{b} = mind_CHsvg{b-1}; 
        if mind_SWsvg{b} =. then mind_SWsvg{b} = mind_SWsvg{b-1}; 
        if mind_alsvg{b} =. then mind_alsvg{b} = mind_alsvg{b-1}; 

       if smok{b}=. or smok{b}=5   then smok{b}=smok{b-1};
       if cignm{b}=. or cignm{b}=0 or cignm{b}=7  then cignm{b}=cignm{b-1}; 

       if trmeatcon{b}=.           then trmeatcon{b} = trmeatcon{b-1};
       if urmeatcon{b}=.           then urmeatcon{b} = urmeatcon{b-1};
       if prmeatcon{b}=.           then prmeatcon{b} = prmeatcon{b-1};
       if fruitscon{b}=.           then fruitscon{b} = fruitscon{b-1};
       if vegealcon{b}=.           then vegealcon{b} = vegealcon{b-1};
       if whgrnscon{b}=.           then whgrnscon{b} = whgrnscon{b-1};
       if poultrcon{b}=.           then poultrcon{b} = poultrcon{b-1};
       if fishalcon{b}=.           then fishalcon{b} = fishalcon{b-1};
       if regeggcon{b}=.           then regeggcon{b} = regeggcon{b-1};
       if nutsalcon{b}=.           then nutsalcon{b} = nutsalcon{b-1};
       if nutlegcon{b}=.           then nutlegcon{b} = nutlegcon{b-1};
       if tdairycon{b}=.           then tdairycon{b} = tdairycon{b-1};
       if ldairycon{b}=.           then ldairycon{b} = ldairycon{b-1};
       if hdairycon{b}=.           then hdairycon{b} = hdairycon{b-1};
       if legumecon{b}=.           then legumecon{b} = legumecon{b-1}; 
       if soyprocon{b}=.           then soyprocon{b} = soyprocon{b-1};
       if ssbcon{b}=.              then ssbcon{b}    = ssbcon{b-1};
       if SES{b}=.                 then SES{b}       = SES{b-1};

       if diab{b-1}=1   then diab{b}=1;
       if chol{b-1}=1   then chol{b}=1;
       if hbp{b-1}=1    then hbp{b}=1;
       if depr{b-1}=1   then depr{b}=1;   
       /*if antitc{b-1}=1 then antitc{b}=1;
       if antihp{b-1}=1 then antihp{b}=1;
       if antd{b-1}=1   then antd{b}=1;*/
    end;
 
/*** Set cutoff at 2023 Jan 31 ***/

/*************************************************************************
***** If an irt date is before June of that qq year or after or equal ****
***** to the next qq year it is incorrect and should be defaulted to  ****
***** June of that qq year.    Make time period indicator tvar=0.     ****
*************************************************************************/
 cutoff=1477;    
   do i=1 to DIM(irt)-1;
      if (irt{i}<(1009+24*i) | irt{i}>=(1033+24*i)) then irt{i}=1009+24*i;
   end;  

%beginex();

   *****************Do-Loop over time periods*****************************;
   
   do i=1 to dim(irt)-1;
      period=i;
      do j=1 to dim(tvar);
         tvar{j}=0;
         end;
      tvar{i}=1;

/*** outcomes ***/

    /*** total AD ***/
      adcase=0;  tad=irt{i+1}-irt{i};    
      if irt{i}<dtdx_ad<=irt{i+1} then do; adcase=1; tad=dtdx_ad-irt{i};  end;
      if irt{i} lt dtdth le irt{i+1} then tad=min(tad, dtdth-irt{i});

    /*** confirmed AD ***/
      confadcase=0;      
      if irt{i}<dtdx_ad<=irt{i+1} and alz1620=1 and dementiadeath=1 then confadcase=1; 

    /* incident non-fatal AD */
      nfadcase=0;     
      if irt{i}<dtdx_ad<=irt{i+1} and alz1620=1 then  nfadcase=1;  

    /* incident fatal dementia */
      fadcase=0;   
      if irt{i}<dtdx_ad<=irt{i+1} and dementiadeath=1 then fadcase=1;   
      /* note, this fetal case may not include all fatal ADRD as who had ADRD before death were censored in doloop by %exclude(irt{i-1} lt dtdx_ad le irt{i}) */
  
      if dtdth eq 9999 then dtdth=.;


 /*Define cvd, cancer and total mortality case variables and time to follow for MPHREG model*/

      dead=0;    tdead=irt{i+1}-irt{i};   
      if irt{i}<dtdth<=irt{i+1} then do; dead=1;tdead=dtdth-irt{i}; end;
      
/******Primary category according to diagnosis coding of ICD 9 ******/
 
      /*** cancer ***/
      dead_can=0; tdead_can=irt{i+1}-irt{i};  
      if irt{i}<dtdth<=irt{i+1} and 140<=newicda<=208 then do; dead_can=1;tdead_can=dtdth-irt{i}; end;

      /*** Diseases of the Circulatory System ***/ 
      dead_cvd=0;  tdead_cvd=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and (390<=newicda<=459) then do; dead_cvd=1;tdead_cvd=dtdth-irt{i}; end;

          /*** Diseases of CHD ***/ 
          dead_chd=0;  tdead_chd=irt{i+1}-irt{i};        
          if irt{i}<dtdth<=irt{i+1} and (390<=newicda<=429 or 440<=newicda<=459) then do; dead_chd=1;tdead_chd=dtdth-irt{i}; end;
          /*** Diseases of stroke ***/ 
          dead_str=0;  tdead_str=irt{i+1}-irt{i};        
          if irt{i}<dtdth<=irt{i+1} and (430<=newicda<=438) then do; dead_str=1;tdead_str=dtdth-irt{i}; end;

      /*** Diseases of the Respiratory System ***/ 
      dead_resp=0;  tdead_resp=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and (460<=newicda<520) then do; dead_resp=1;tdead_resp=dtdth-irt{i}; end;

      /*** Diseases of neurodegenerative disease ***/ 
      dead_neuro=0;  tdead_neuro=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and newicda in (290,332,335,340,342,348) then do; dead_neuro=1;tdead_neuro=dtdth-irt{i}; end;

      /*** Infectious and Parasitic Diseases ***/
      dead_inf=0;  tdead_inf=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and (0<newicda<140) then do; dead_inf=1;tdead_inf=dtdth-irt{i}; end;   
        
      /*** kidney disease ***/ 
      dead_rf=0;  tdead_rf=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and (580<=newicda=<593) then do; dead_rf=1;tdead_rf=dtdth-irt{i}; end;
                                       
      /*** diabetes ***/ 
      dead_dm=0;  tdead_dm=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and (newicda=250) then do; dead_dm=1;tdead_dm=dtdth-irt{i}; end;

      dead_oth=0;tdead_oth=irt{i+1}-irt{i}; 
      if irt{i}<dtdth<=irt{i+1} and (dead_can=0 and dead_cvd=0 and dead_resp=0  and dead_neuro=0 and
                                     dead_inf=0 and dead_rf=0  and dead_dm=0)  then do; 
                                                                        dead_oth=1;   tdead_oth=dtdth-irt{i}; end;


      /*** others besides cvd and cancer ***/ 
      dead_oths=0;  tdead_oths=irt{i+1}-irt{i};        
      if irt{i}<dtdth<=irt{i+1} and (dead_can=0 and dead_cvd=0 )  then do; 
                                                                        dead_oths=1;   tdead_oths=dtdth-irt{i}; end;
    
      if dtdth eq 9999 then dtdth=.;

                    
   /******AGE GROUP******/
      age{i}= int((irt{i}-birthday)/12);

      agecon=age{i};

      if age{i}<60 then agegp=0;     **** Define the agegp in the i-th period;
            else if age{i}<65 then agegp=1;
            else if age{i}<70 then agegp=2;
            else if age{i}<75 then agegp=3;
            else if age{i}<80 then agegp=4;
            else agegp=5;
            if age{i}=. then agegp=.;

      %indic3(vbl=agegp, prefix=agegp, min=1, max=5, reflev=0,missing=., usemiss=0,
                   label1='age<60ys',
                   label2='age60-64ys',
                   label3='age65-69ys',
                   label4='age70-74ys',
                   label5='age75-79ys',
                   label6='age>=80ys');

  ****** main exposure ******; 
        ahei=aheim{i}; 
        dashcon=dashm{i};
        pdicon =pdim{i};
        hpdicon=hpdim{i};
        updicon=updim{i};
        amedcon=emedm{i}; 
        amedcov=emedv{i};  
        mindcon=mindm{i};  

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
        alcocon= alcon{i}; 
        ssb    = ssbcon{i};
        ldairy = ldairycon{i};
        hdairy = hdairycon{i};
        legume = legumecon{i};
        soypro = soyprocon{i};
        legsoy = sum(0,soypro,legume);
        nSES   = SES{i};

        scfmem   = mem0820;  
        scfmem_z = zmem;
        scfmem_sd= mem0820_sd;
        PGScon   =PGS002280_scaled;
        PGScon2   =PGS000334_scaled;

      /****** Indicator for Smoking ******/ 
      /****** Indicator for Smoking ******/ 
      /***smoke  $label 1.Never;  2.Past; 3.No, unknown past history;   4.Current; 5.PASSTHRU;
          cgnm $label 1.1-4/d; 2.5-14; 3.15-24; 4.25-34; 5.35-44;6.45 or more; 7.PASSTHRU; 0.not app or no qx***/ 

      if smok{i}=5 then smok{i}=.;      
      if smok{i}<4 then cignm{i}=0;      
     
      if cignm{i}=7 then cignm{i}=.;      
      if cignm{i}=0 then cignm{i}=.;      
    
      if smok{i}=1 then smkm(i)=1;                    /**   smkm=1, never*/
      if smok{i}=2 then smkm(i)=2;                    /**   smkm=2, past*/
      if smok{i}=3 then smkm(i)=2;      
      if smok{i}=. then smkm(i)=7;      
      if smok{i}=4 then 
         do;	 
         if 1<=cignm{i}<=2 then smkm(i)=3;           /**     smkm=3, 1-14*/
         else if cignm{i}=3 then smkm(i)=4;          /*      smkm=4,15-24*/
         else if 4<=cignm{i}<=6 then smkm(i)=5;	    /** smkm=5: >=25/d */
         else if cignm{i}=. then smkm(i)=6;         /** smkm=6: currnt smker, amount unknown */
         end;  

      select(smkm{i});
	 when (1)     smm{i}=1;
	 when (2)     smm{i}=2;
	 when (3,6)   smm{i}=3;
	 when (4)     smm{i}=4;	
	 when (5)     smm{i}=5;	
	 otherwise    smm{i}=.;
	 end;

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


       /****** Indicator for BMI ******/
     if height>0 and weight{i}>0 then bmi{i}=(weight{i}*0.45359237)/   
        ((height*25.4/1000)*(height*25.4/1000));

      bmicon=bmi{i};if bmi{i}=0 then bmicon=.;
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

         basebmicon=bmi86;
         if bmi86<23 then basebmic=1;
           else if 23=<bmi86<25 then basebmic=2;
             else if 25=<bmi86<30 then basebmic=3;
               else if 30=<bmi86<35 then basebmic=4;
                 else basebmic=5;
          if bmi86=. then basebmic=.;
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

         if i>1 then do; bmichg=bmi{i}-bmi{i-1};end;
            if bmichg=<0 then bmichgcat=0; else bmichgcat=1; if bmichg=. then bmichgcat=.;
         %indic3(vbl=bmichgcat, prefix=bmichgcat, min=1, max=1, reflev=0, missing=., usemiss=1,
                label1='bmichg=<0',
                label2='bmichg>0' );


      /****** Indicator Physical Activity ******/
       actcon=actm{i}; /***continue, hours/wk***/  
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

      /****** Indicator Alcohol consumption ******/
       alcocon=alcon{i}; /***continuous variable, g/day***/                                                               
 
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

****** covariates ******;

       if white ne 1 then white=0;

      /*** marriage status ***/

       if marriage{i}=1 then marry=1; else marry=0; 

      /*** living alone ***/
            select(alone{i});
                when (1)     live_alone=1;
      	         otherwise    live_alone=0;
            end;

      /*** multiple vitamin ***/
             if mvyn{i}=2 then mvit=1; else   mvit=0;   * consistent with NHS;     
      
     /*** aspirin ***/
             aspirin=aspi{i};
             if aspirin ne 1 then aspirin=0;    

     /***family history of MI***/
       if famhxmi ne 1 then famhxmi=0;

     /***family history of cancer***/
       if famhxca ne 1 then famhxca=0;
             
      /****** Indicator for History of antihypertensive user *******/
            select(antihp{i});
                when (1)     antihbp=1;
      	         otherwise    antihbp=0;
            end;
     
      /****** Indicator for History of High Blood Pressure *******/
            select(hbp{i});
                when (1)     htn=1;
      	         otherwise    htn=0;
            end;

            if hbp{i}=1 or antihp{i}=1 then highhbp=1; else highhbp=0;
            if hbp86=1 or htnr86 then basehbp=1; else basehbp=0;       
           
      /****** Indicator for History of High TC ******/
            select(chol{i});
                when (1)     hchol=1;
      	         otherwise    hchol=0;
            end;

            IF antitc{i}=1 then antihtc=1; else antihtc=0;

            if chol{i}=1 or antitc{i}=1 then hightc=1; else hightc=0;
            if chol86=1 or hbcrx86=1 then basechol=1; else basechol=0; 
   
      /****** Indicator for History of diabetes *******/
            select(diab{i});
                when (1)     diabetes=1;
      	         otherwise    diabetes=0;
            end;

            if db86=1 then basedb=1; else basedb=0;

      /****** Indicator for History of selfreport heart disease MI & Angina ******/
            select(repmi{i});
                when (1)     chd=1;
      	         otherwise    chd=0;
            end;

            select(rpang{i});
                when (1)     angina=1;
      	         otherwise    angina=0;
            end;

            select(rpcabg{i});
                when (1)     cabg=1;
      	         otherwise    cabg=0;
            end;

      /****** Indicator for History of selfreport stroke ******/
            select(repstrk{i});
                when (1)     stroke=1;
      	         otherwise    stroke=0;   
            end;


      /****** Indicator for History of Cancer ******/
            select(canc{i});
                when (1)     cancer=1;
      	         otherwise    cancer=0;
            end;

      /****** Indicator for History of depression ******/
            select(depr{i});
                when (1)     depression=1;
      	         otherwise    depression=0;
            end;

     /****** Indicator for History of antidepressant user *******/
            select(antd{i});
                when (1)     antidep=1;
      	         otherwise    antidep=0;
            end;
      
       if dtdth in (0,.)then dtdth=.;
       if dtdth eq 9999 then dtdth=.;


     if white ne 1 then white=0;
     if mvit ne 1 then mvit=0;  

     /*  Profession  */
     /* profssn: 1.dentist 2.hosp pharmacist 3.optometrist 4.osteopath 5.pharmacist 6.podiatrist 7.vet*/

     if profssn=1 then proff=1;
          else if profssn=7 then proff=2;
          else proff=3;
     if profssn=. then proff=.; 
     %indic3(vbl=proff, prefix=proff, min=2, max=3, reflev=1, missing=., usemiss=0,
          label1='Dentist',
          label2='Vet',
          label3='Pharmacist/Optometrist/Osteopath/Podiatrist');


     /* family history of demenia*/
       if mdem90=1 or mdem92=1 or demmo12=2 or fdem90=1 or fdem92=1 or demfad12=2 or sdem90=1 or sdem92=1 or demsi12=2 then fhdem=1; /* dem hx*/
       else  fhdem=0;
  
     /****************  BASELINE EXCLUSIONS ********************************************/
      if i=1 then do;

   %exclude(exrec eq 1);                 *multiple records and not in master file; 

   %exclude(age86 eq .);   
   %exclude(pd eq 1);
   %exclude(str86 eq 1);
   %exclude(can86 eq 1);
  
   %exclude(calor86n gt 4200); 
   %exclude(calor86n lt 800);
   %exclude(calor86n eq .);  

   %exclude(calor86n eq .); 
   %exclude(calor86n eq .); 

   %exclude(bmi86 eq .);
   %exclude(emed86h eq .);
   %exclude(MIND86 eq .);   
   %exclude(smoke86 eq .);
   %exclude(smoke86 eq 5);

   %exclude(0 lt dtdth le irt{i});    
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

keep id period dtdth agecon agegp age86  birthday newicda
     confadcase adcase tad dtdx_ad dtdth nfadcase  fadcase  dead tdead

     irt86 irt88 irt90 irt92 irt94 irt96 irt98 irt00 irt02 irt04 irt06 irt08 irt10 irt12 irt14 irt16 irt18 irt20 cutoff
     t86 t88 t90 t92 t94 t96 t98 t00 t02 t04 t06 t08 t10 t12 t14 t16 t18 t20

     dead  dead_can  dead_cvd  dead_chd  dead_str  dead_resp  dead_neuro  dead_inf   dead_rf   dead_dm  dead_oth  dead_oths 
     tdead tdead_can tdead_cvd tdead_chd tdead_str tdead_resp tdead_neuro tdead_inf  tdead_rf  tdead_dm tdead_oth tdead_oths  

     trmeat prmeat urmeat 
     fruits vegets whgrns  poultr  fishs  eggs   nuts nutleg tdairy ldairy  hdairy  legume  soypro  legsoy daykcal ssb
     ahei dashcon pdicon  hpdicon updicon  amedcon amedcov mindcon

     actcon  actcc &actcc_ alcocon alcc &alcc_ smkk &smkk_ white marry live_alone mvit aspirin  
     bmicon bmic &bmic_  
     antihbp htn highhbp basehbp hchol hightc antihtc  diabetes angina chd cabg stroke cancer depression antidep

     profssn  fhdem nSES proff &proff_
     
     PGS002280_scaled PGS000334_scaled apoe4grp apoe_2cat APOE4 platform platforms PC1_comb PC2_comb PC3_comb PC4_comb apoe4con ADdiagage

     basebmicon basebmic &basebmic_  smkcat &smkcat_ wtchgcon wtchgcat &wtchgcat_

     scfmem scfmem_z  scfmem_sd PGScon PGScon2 dementiadeath alz1620 flag16 flag18 flag20 tbcat shing2vac;
run;

/* check missing */ 
proc means n nmiss min mean std median max;
var nSES actcon proff smkk;
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
    varlist= ahei dashcon pdicon  hpdicon updicon  amedcon amedcov mindcon PGScon PGScon2,
    numquant=5,
    mscore=T,
    quantname=q,
    cutdat=cutpoints,
    outdat=pre_pm,
    indic=T);
 

data pre_pm;
  set pre_pm end=_end_;  

  ****** neiborghood SES index   ******;
      %indic3(vbl=qnSES, prefix=qnSES, min=1, max=2, reflev=0, missing=., usemiss=1);  

      if antidep=1 or  depression=1 then depre=1; else depre=0;

      %indic3(vbl=qahei,       prefix=qahei,   min=1, max=4, reflev=0, missing=., usemiss=1);  
      %indic3(vbl=qamedcon,    prefix=qamedcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
      %indic3(vbl=qamedcov,    prefix=qamedcov,   min=1, max=4, reflev=0, missing=., usemiss=0);  
      %indic3(vbl=qmindcon,    prefix=qmindcon,   min=1, max=4, reflev=0, missing=., usemiss=0);
      %indic3(vbl=qdashcon,    prefix=qdashcon,   min=1, max=4, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qpdicon,     prefix=qpdicon,    min=1, max=4, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qhpdicon,    prefix=qhpdicon,   min=1, max=4, reflev=0, missing=., usemiss=1);
      %indic3(vbl=qupdicon,    prefix=qupdicon,   min=1, max=4, reflev=0, missing=., usemiss=1);

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

/* for trend */
proc means median; var nSES;    class qnSES; run;
proc means median; var amedcon; class qamedcon; run;
proc means median; var mindcon; class qmindcon; run;
proc means median; var actcon;  class actcc; run;

data hpfsdata;
set pre_pm  end=_end_;
cohort=2;
sex=2;
id=200000000+id; 
interval=period+3;
phmsstatus=1;
hiedu=1;
husbedu=1;
FORMAT _all_;
INFORMAT _all_;

****** P for trend ******;
if actcc=1 then medactcc=0;
if actcc=2 then medactcc=1.2;
if actcc=3 then medactcc=2.5;
if actcc=4 then medactcc=5.0;
if actcc=5 then medactcc=9.7;
if actcc=. then medactcc=99999;
actccmiss=0; if actcc=. then actccmiss=1;

if nSES=. then medqnSES=99999;
nSESmiss=0; if nSES=. then nSESmiss=1;
 
if wtchgcon=. then wtchgcon=99999;
wtchgmiss=0; if wtchgcon=. then wtchgmiss=1;

if tbcat=. then tbcat=2;
%indic3(vbl=tbcat, prefix=tbcat, min=1, max=3, reflev=1, missing=., usemiss=0,
                label1='TB negative',
                label2='non-response',
      	    label3='TB positive');

run;

proc datasets;
  delete pre_pm;
  run;