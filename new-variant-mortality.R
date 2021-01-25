
library(tidyverse)
library(patchwork)
library(survival)

## Data loading functions ----
getDeathsLineList = function(path, ...) {
  tmp = readxl::read_excel(path, col_types = "text")
  datecols = c(colnames(tmp) %>% stringr::str_subset("date"),"dod")
    
  for(datecol in datecols) {
    tmp[[datecol]] = suppressWarnings(as.Date(as.numeric(tmp[[datecol]]),"1899-12-30"))
  }
  tmp = tmp %>% 
      dplyr::mutate(
        age = as.integer(age),
        gender = normaliseGender(ifelse(is.na(gender),"unknown",gender))
      )
  return(tmp %>% mutate(FINALID=as.numeric(finalid)) %>% dplyr::ungroup())
}

#' @description Load line list
#' 
#' @return raw line list data set
getSGeneLineList = function(path,...) {
  tmp = readr::read_csv(path)
  return(tmp %>% dplyr::ungroup())
}

#' @description Load line list
#' 
#' @return raw line list data set
getLineList = function(path,...) {
    tmp = readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()))
    tmp = tmp %>% 
      dplyr::mutate(
        Onsetdate = as.Date(Onsetdate,"%d/%m/%Y"),
        specimen_date = as.Date(specimen_date,"%d/%m/%Y"),
        lab_report_date = as.Date(lab_report_date,"%d/%m/%Y")
      ) 
    if(any(is.na(tmp$specimen_date))) stop("Problem parsing dates")
    return(tmp %>% mutate(
      pillar_2_testingkit = tolower(pillar_2_testingkit),
      age = suppressWarnings(as.numeric(age)),
      FINALID = as.numeric(FINALID),
      ethnicity_final = case_when(
        ethnicity_final %in% c("African (Black or Black British)","Any other Black background","Caribbean (Black or Black British)") ~ "Afro-carribbean",
        ethnicity_final %in% c("Any other Asian background","Bangladeshi (Asian or Asian British)","Indian (Asian or Asian British)","Pakistani (Asian or Asian British)") ~ "Asian",
        ethnicity_final %in% c("Any other White background","British (White)","Irish (White)") ~ "White",
        ethnicity_final %in% c("Any other Mixed background","Any other ethnic group","White and Black Caribbean (Mixed)","White and Black African (Mixed)","Chinese (other ethnic group)") ~ "Other",
        TRUE ~ "Unknown"),
      residential_category = case_when(
        cat == 'Residential dwelling (including houses, flats, sheltered accommodation)' ~ "Residential",
        cat == 'Care/Nursing home' ~ "Care home",
        cat == 'Undetermined'~"Other/Unknown",
        cat == 'Medical facilities (including hospitals and hospices, and mental health)'~"Other/Unknown",
        cat == 'Other property classifications'~"Other/Unknown",
        cat == 'House in multiple occupancy (HMO)' ~ "Residential",
        cat == 'Prisons, detention centres, secure units'~"Other/Unknown",
        cat == 'Residential institution (including residential education)'~"Other/Unknown",
        cat == 'No fixed abode'~"Other/Unknown",
        cat == 'Overseas address'~"Other/Unknown",
        TRUE~"Other"
      )
    ) %>% dplyr::ungroup())
}


normaliseGender = function(gender) {
  case_when(
    is.na(gender) ~ NA_character_,
    gender %>% stringr::str_detect("f|F") ~ "female",
    gender %>% stringr::str_detect("m|M") ~ "male",
    gender %>% stringr::str_detect("u|U") ~ "unknown",
    TRUE ~ "unknown")
}

## Load the data
# check if already loaded
# we discard S gene data for patients who have had more than one S gene test result. 
# this will put these patients into the "Unknown / P1" category, and will be excluded from the case control part of study but included in the main analysis.

loadData = function(dir="~/Data/new-variant",date="20210118", censorLength = 28, B117Date = as.Date("2020-10-01"), earliestDate = B117Date, latestDate = NULL, minAge = 30) {
  ll = getLineList(path=paste0(dir,"/Anonymised Combined Line List ",date,".csv"))
  message("Loaded cases line list ",date,":",nrow(ll))
  dll = getDeathsLineList(path = paste0(dir,"/",date," COVID19 Deaths.xlsx"))
  message("Loaded deaths line list ",date,":",nrow(dll))
  sgll = getSGeneLineList(path = paste0(dir,"/SGTF_linelist_",date,".csv"))
  message("Loaded sgene line list ",date,":",nrow(sgll)," entries, ",n_distinct(sgll$FINALID)," patients")
  sgllUniq = sgll %>% filter(!is.na(FINALID)) %>% group_by(FINALID) %>% mutate(sGeneResultCount = n()) %>% filter(sGeneResultCount == 1) %>% ungroup()
  message("Selected patients with single s gene results ",nrow(sgllUniq))
  #sgllUniq = sgll %>% filter(!is.na(FINALID)) %>% group_by(FINALID) %>% arrange(desc(specimen_date)) %>% mutate(sGeneResultCount = n()) %>% filter(row_number() == 1) %>% ungroup()
  
  
  if (identical(latestDate,NULL)) latestDate = max(dll$dod,na.rm=TRUE)
  
  # Is this a potential source of bias?
  excludedSgll = sgll %>% anti_join(sgllUniq, by = "FINALID")
  if(nrow(excludedSgll) > 0) {
    message("we excluded ",excludedSgll %>% pull(FINALID) %>% unique() %>% length()," patients with multiple S Gene test results in pillar 2")
    excludedDeaths = dll %>% inner_join(excludedSgll, by="FINALID", suffix=c("",".sgene")) #%>% filter(specimen_date > "2020-10-01" & age>30)
    message("of these ",excludedDeaths %>% filter(specimen_date.sgene > B117Date) %>% pull(FINALID) %>% unique() %>% length()," died during the study period")
  }
  
  # get the S gene line list and remove all patients who have more than one S Gene result. Otherwise we don;t have a single source in time for when to start the survival clock.
  # interpretSGene assigns sGene "Positive" for S+N+ORF+ with cutOff for S at 40 and N & ORF at 30; Negative is S-N+ORF+; equivocal is anything else (e.g. S-N-ORF+)
  # We exclude anything which does not match the line list.
  
  # Join the line list, deaths line list and uniquified sGene line list by FINALID
  coxData = ll %>% 
    left_join(sgllUniq, by="FINALID", suffix=c("",".sgene")) %>% 
    left_join(dll %>% mutate(FINALID=as.numeric(finalid)), by=c("FINALID"),suffix=c("",".death")) %>% 
    filter(specimen_date >= B117Date & specimen_date <= latestDate) # make sure we are not including tests for which we could not have death info.
  
  message("select cases with specimen_date >",earliestDate," & <=",latestDate," :",nrow(coxData))
  
  coxData %>% filter(is.na(specimen_date.sgene)) %>% nrow() %>% message(" patients with missing S gene result excluded")
  coxData %>% filter(!is.na(specimen_date.sgene)) %>% nrow() %>% message(" patients with known S gene result included")
  
  coxData %>% filter(is.na(dod)) %>% nrow() %>% message(" patients who have no date of death (and are presumed alive)")
  
  coxData2 = coxData %>% 
    select(
      FINALID, 
      specimen_date, 
      lab_report_date, 
      specimen_date, 
      specimen_date.sgene, 
      dateadmission_NHSE, 
      dod,
      death_type28, 
      age,sex,imd_decile,
      P2CH3CQ,P2CH1CQ,P2CH2CQ,P2CH4CQ, # gene copy numbers
      sgtf, sgtf_under30CT, # gene copy number interpretation
      NHSER_name,NHSER_code,
      ethnicity_final, # reduced set of ethnicities
      residential_category, # reduced set of residential categories
      LTLA_name) %>% 
    mutate(
      died = !is.na(dod),
      #diedWithin28days = !is.na(death_type28), # This is the line list definition of death within 28 days
      diedWithin28days = died & as.numeric(dod-specimen_date) <= censorLength, # This is a definition of death within 28 days not relying on line list definition
      dodAt28days = ifelse(diedWithin28days, dod, NA), # this is the date of death (or NA if it is >28 days after sgene specimen date or there is no specimen)
      censoredDate = pmin(latestDate, specimen_date+censorLength,na.rm = TRUE), # earliest of last recorded death or specimen date
    ) %>% 
    mutate(
      reportingDelay = as.numeric(lab_report_date-specimen_date),
      sGeneDelay = as.numeric(specimen_date.sgene-specimen_date),
      admissionDelay = as.numeric(dateadmission_NHSE-specimen_date),
      time = ifelse(diedWithin28days, dodAt28days-as.numeric(specimen_date), censoredDate-as.numeric(specimen_date)), 
      status = ifelse(!diedWithin28days,0,1)
    ) %>%
    mutate(
      ageCat = cut(age, breaks = c(-Inf,30,60,70,80,Inf), labels = c("<30","30-59","60-69","70-79","80+"), right=FALSE),
      LTLA_name =  as.factor(LTLA_name),
      preB117 = specimen_date < B117Date
    )
  
  # TODO: https://cran.r-project.org/web/packages/Gmisc/vignettes/Grid-based_flowcharts.html
  # https://www.r-bloggers.com/2018/05/flow-charts-in-r/
  
  tmpCox = coxData2 %>% interpretSGene()
  tmpCox %>% filter(sGene=="Equivocal") %>% nrow() %>% message(" with equivocal results excluded")
  tmpCox %>% filter(sGene=="Negative") %>% nrow() %>% message(" S- cases")
  tmpCox %>% filter(sGene=="Positive") %>% nrow() %>% message(" S+ controls")
  
  reportExclusions = function(cox,type) {
    type = paste0(" ",type)
    cox %>% filter(is.na(imd_decile)) %>% nrow() %>% message(type," with missing IMD")
    cox %>% filter(is.na(sex) | sex=="Unknown") %>% nrow() %>% message(type," with missing or unknown gender")
    cox %>% filter(is.na(age)) %>% nrow() %>% message(type," with unknown age")
    cox %>% filter(age<minAge) %>% nrow() %>% message(type," with age <",minAge)
    cox %>% filter(is.na(ethnicity_final)) %>% nrow() %>% message(type," with unknown ethnicity")
    cox %>% filter(admissionDelay < 0) %>% nrow() %>% message(type," admitted before pillar 2 test taken")
    cox %>% filter(reportingDelay < 0 | reportingDelay>=20) %>% nrow() %>% message(type," whose tests are reported before specimen taken or >19 days after specimen taken")
    cox %>% filter(sGeneDelay < 0 | sGeneDelay>=20) %>% nrow() %>% message(type," whose sGene tests are >19 days after original positive test")
    cox %>% filter(time<=0) %>% nrow() %>% message(type," whose time from test to event (death) is 0 or less (e.g. postmortem tests)")
  }
  
  tmpCox %>% filter(sGene=="Negative") %>% reportExclusions("S- cases")
  tmpCox %>% filter(sGene=="Positive") %>% reportExclusions("S+ controls")
  
  coxData2 = coxData2 %>%
    filter(!is.na(imd_decile)) %>%
    filter(!is.na(sex) & sex != "Unknown") %>% mutate() %>%
    filter(!is.na(ageCat) & !is.na(age)) %>%
    filter(!is.na(ethnicity_final) & !ethnicity_final %in% c("Unknown")) %>%
    filter(is.na(admissionDelay) | admissionDelay >= 0) %>%
    filter(is.na(reportingDelay) | (reportingDelay >= 0 & reportingDelay <20) ) %>%
    filter(is.na(sGeneDelay) | (sGeneDelay >= 0 & sGeneDelay <20) ) %>%
    filter(time>0) %>% # quite a few specimens take post mortem?
    mutate(
      ethnicity_final = ethnicity_final %>% forcats::as_factor() %>% forcats::fct_relevel("White"),
      sex = sex %>% forcats::as_factor(),
      imd_decile = forcats::fct_relevel(forcats::as_factor(as.numeric(imd_decile)),"5"),
      imd_decile_2 = forcats::as_factor(as.numeric(imd_decile)),
      deathStatus = ifelse(diedWithin28days,"Dead <28 days","Other")
    ) %>% filter(age >= minAge)
  
  tmpCox = coxData2 %>% interpretSGene()
  tmpCox %>% filter(sGene=="Negative") %>% nrow() %>% message(" S- cases eligable for inclusion")
  tmpCox %>% filter(sGene=="Positive") %>% nrow() %>% message(" S+ controls eligable for inclusion")
  
  ## quality control
  if (any(is.na(coxData2$LTLA_name))) stop("unknowns for LTLA")
  
  return(coxData2)
}

## Analysis functions: ----

# function to help interpret S gene N gene and ORF gene in 
interpretSGene = function(sGeneLineList, S_CT = 30, ORF1ab_CT = 30, N_CT = 30, Control_CT = Inf) {
  sGeneLineList %>% mutate(
    ORF1ab_CT_threshold = ORF1ab_CT,
    N_CT_threshold = N_CT,
    S_CT_threshold = S_CT,
    S_pos = P2CH3CQ > 0 & P2CH3CQ <= S_CT,
    S_undetect = P2CH3CQ == 0,
    N_pos = P2CH2CQ > 0 & P2CH2CQ <= N_CT,
    ORF1ab_pos = P2CH1CQ > 0 & P2CH1CQ <= ORF1ab_CT,
    Control_pos = P2CH4CQ > 0 & P2CH4CQ <= Control_CT,
    sGene = case_when(
      is.na(P2CH1CQ) ~ "Unknown",
      S_pos & N_pos & ORF1ab_pos & Control_pos ~ "Positive",
      S_undetect & N_pos & ORF1ab_pos & Control_pos ~ "Negative",
      TRUE ~ "Equivocal"
    ),
    CT_N = ifelse(P2CH2CQ > 0, P2CH2CQ, 40)
  ) %>% mutate(
    result = ifelse(!Control_pos,"No control",paste0(ifelse(S_pos,"S+","S-"),ifelse(N_pos,"N+","N-"),ifelse(ORF1ab_pos,"ORF+","ORF-")))
  ) %>% mutate(
      sGeneEra = case_when(
        is.na(sGene) ~ "Unknown (and P1)",
        sGene == "Negative" & !preB117 ~ "Neg Post B.1.1.7",
        sGene == "Negative" & preB117 ~ "Neg Pre B.1.1.7",
        TRUE ~ as.character(sGene)
      ),
      sGene = sGene %>% forcats::fct_relevel("Positive"),
      sGeneEra = sGeneEra %>% forcats::as_factor() %>% forcats::fct_relevel("Positive"),
      relativeCopyNumber = 2^(median(CT_N,na.rm=TRUE)-CT_N)
    )
}

bootstrappedCoxModel = function(coxFinal, modelName, modelFormula) {
  message("Calculating ",modelName)
  bootSurvModels = coxFinal %>% group_by(boot) %>% group_modify(function(d,g,...) {
    survModel = survival::coxph(modelFormula, d)
    out = tidyCoxmodel(survModel) %>% mutate(model=modelName)
    rm(survModel)
    return(out)
  })
}

tidyCoxmodel = function(survModel) {
  broom::tidy(survModel) %>% bind_cols(as.data.frame(confint(survModel))) %>% 
    mutate(
      HR = exp(estimate), 
      HR.lower=exp(`2.5 %`), 
      HR.upper=exp(`97.5 %`), 
      HR.pValue = p.value)
}

createPairedMatches = function(coxData, ageTolerance, specimenDateTolerance,matchOn, max) {
  
  posSgene = coxData %>% filter(sGene=="Positive")
  negSgene = coxData %>% filter(sGene=="Negative")
  message("S gene positives: ",posSgene %>% nrow())
  message("S gene negatives: ",negSgene %>% nrow())
  
  if(as.numeric(nrow(negSgene)) * as.numeric(nrow(posSgene)) > max) stop("this would try and create a very large cross join of ",nrow(posSgene),"x",nrow(negSgene))
  
  # pairwise matching
  ## case matching: 
  # match cases by "sex","imd_decile","LTLA_name","ethnicity_final","ageCat"
  # then ensure age difference < 2 years & specimen date < 4 days
  # select S negatives and match to S positives
  # randomly select a single match out of possible matches
  matches = 
    negSgene %>% 
    inner_join(
      posSgene, 
      by=matchOn,suffix=c(".neg",".pos")) %>% 
    filter(
      abs(age.neg-age.pos) <= ageTolerance & 
        abs(as.numeric(specimen_date.neg-specimen_date.pos)) <= specimenDateTolerance
    )
  
  return(matches)
}

uniquifyPairs = function(pairedMatches, resolutionStrategy, bootstraps, ratio.neg, ratio.pos) {
  message("Before sampling the matched set contains ",n_distinct(pairedMatches$FINALID.neg)," unique S gene negatives and ",n_distinct(pairedMatches$FINALID.pos)," unique S gene positives")
  
  if (resolutionStrategy == "bootstrap") {
    # TODO: I am not sure if there is a risk here of a systematic bias, or at least a random bias that is systematically present
    # throughout the analysis.
    # What happens here is the majority of matches are 1:1 and will always be picked in each bootstrap.
    # These matches are essentially oversampled
    # Other match ratios eg. 1:2 or 2:1 will be randomly picked for each boot strap - the lines starting group_by(FINALID.neg) and group_by(FINALID.pos) do this.
    # There are not equal amounts of potential matches from both classes before sampling.
    
    # @Leon
    # randomly select a single match out of possible matches for each bootstrap
    set.seed(101)
    matchesSelected = NULL
    message("bootstrapping",appendLF = FALSE)
    for(i in 1:bootstraps) {
      message("..",i,appendLF = FALSE)
      tmp = pairedMatches %>% ungroup() %>%
        mutate(rnd = runif(nrow(pairedMatches)), wt.neg=1, wt.pos=1) %>% 
        arrange(rnd) %>% # put the dataset in a random order
        group_by(FINALID.neg) %>% filter(row_number()==1) %>% # filter out a single matching at random from pos to neg
        group_by(FINALID.pos) %>% filter(row_number()==1) %>%  # filter out a single matching at random from neg to pos
        ungroup() %>% mutate(boot = i)
      matchesSelected = bind_rows(matchesSelected,tmp)
    }
    message("..finished")
    
    # check for some systematic bias as a result of matching process
    # This is biased by censoring - TODO visualise and explain.
    # negExcluded = negSgene %>% left_join(negMatched %>% select("FINALID") %>% mutate(matched=1), by="FINALID")
    # posExcluded = posSgene %>% left_join(posMatched %>% select("FINALID") %>% mutate(matched=1), by="FINALID")
    # inclusions = bind_rows(negExcluded,posExcluded) %>% mutate(matched = ifelse(is.na(matched), 0, matched))
    # inclusions %>% group_by(sGene, died, matched) %>% summarise(count = n()) %>% clipr::write_clip() #%>% mutate(percentOfSGene = count/sum(count)) %>% group_by(died) %>% mutate(percentByDiedStatus = count/sum(count)) %>% arrange(sGene,died)
    # p_died_given_sGene = inclusions %>% group_by(sGene,died) %>% summarise(count = n()) %>% mutate(percent = count/sum(count))
    # p_died_given_sGene_and_matched = inclusions %>% filter(matched) %>% group_by(sGene,died) %>% summarise(count = n()) %>% mutate(percent = count/sum(count))
    # p_died_given_sGene %>% clipr::write_clip()
    # p_died_given_sGene_and_matched %>% clipr::write_clip()
  } else if(resolutionStrategy == "drop") {
    # drop all pairs which have more than one match
    matchesSelected = pairedMatches %>%
      group_by(FINALID.pos) %>% filter(n() == ratio.pos) %>%
      group_by(FINALID.neg) %>% filter(n() == ratio.neg) %>%
      ungroup() %>%
      mutate(boot=1, wt.neg=1, wt.pos=1)
  } else if(resolutionStrategy == "none") {
    # ignore duplicates
    matchesSelected = pairedMatches %>% mutate(boot=1) %>%
      group_by(FINALID.neg) %>% mutate(wt.neg=1/n()) %>%
      group_by(FINALID.pos) %>% mutate(wt.pos=1/n()) %>%
      ungroup()
  } else {
    stop("No such strategy: ",resolutionStrategy)
  }
  message("After sampling the matched set contains ",n_distinct(matchesSelected$FINALID.neg)," unique S gene negatives and ",n_distinct(matchesSelected$FINALID.pos)," unique S gene positives")
  return(matchesSelected)
}

pairedMatchesToCoxData = function(matches, matchOn) {
  if (!"boot" %in% colnames(matches)) matches = matches %>% mutate(boot=1)
  sPos = matches %>% select(boot,ends_with(".pos"),all_of(matchOn)) %>% rename_with(.fn=function(x) stringr::str_remove(x,".pos"),.cols=ends_with(".pos")) %>% mutate(sGene="Positive") %>% distinct()
  sNeg = matches %>% select(boot,ends_with(".neg"),all_of(matchOn)) %>% rename_with(.fn=function(x) stringr::str_remove(x,".neg"),.cols=ends_with(".neg")) %>% mutate(sGene="Negative") %>% distinct()
  return(bind_rows(sPos,sNeg))
}

runAnalysis = function(
    coxData, sGeneCtThreshold = ctThreshold, ctThreshold=30, 
    bootstraps=10, ageTolerance = 3, specimenDateTolerance = 5,
    modelFormula = list(`S gene only`=Surv(time,status) ~ sGene),
    max = 450000*450000, includeRaw=FALSE, includeMatched=FALSE,
    matchOn = c("sex","imd_decile","LTLA_name","ethnicity_final","ageCat"), #"residential_category"),
    resolutionStrategy = "bootstrap",ratio.neg=1,ratio.pos=1
  ) {
  result = list(params=tibble(
    sGeneCtThreshold = sGeneCtThreshold,
    ctThreshold = ctThreshold,
    bootstraps=bootstraps,
    ageTolerance = ageTolerance,
    specimenDateTolerance = specimenDateTolerance
    ))
  
  message("Interpreting data with CT threshold of S<",sGeneCtThreshold,":N<",ctThreshold,":ORF<",ctThreshold,
          "; tolerances - age=",ageTolerance,":specimen_date=",specimenDateTolerance)
  
  # Interpret S Gene cutoffs
  coxData3 = coxData %>% interpretSGene(S_CT = sGeneCtThreshold,N_CT = ctThreshold,ORF1ab_CT = ctThreshold)
  if (includeRaw) result$rawData = coxData3
  
  # make a set of matched pairs
  matches = coxData3 %>% createPairedMatches(ageTolerance = ageTolerance, specimenDateTolerance = specimenDateTolerance, matchOn = matchOn, max=max)
  if (includeMatched) result$pairedMatchedUnfilteredData = matches
  
  # resolve duplicates in matched pairs
  matchesSelected = matches %>% uniquifyPairs(resolutionStrategy=resolutionStrategy, bootstraps = bootstraps, ratio.neg=ratio.neg, ratio.pos=ratio.pos)
  if (includeMatched) result$pairedMatchedData = matchesSelected
  
  # prepare data for cox PH models
  coxFinal = matchesSelected %>% pairedMatchesToCoxData(matchOn = matchOn) %>% 
    mutate(sGene = sGene %>% forcats::as_factor() %>% forcats::fct_relevel("Positive")) # equivocal category is dropped
  
  negMatched = coxFinal %>% filter(sGene=="Negative")
  posMatched = coxFinal %>% filter(sGene=="Positive")
  if (includeMatched) result$coxMatchedData = coxFinal
  
  ## execute survival models ----
  message("Execute survival models")
  result$table4Data = NULL
  for (modelName in names(modelFormula)) {
    result$table4Data = result$table4Data %>% bind_rows(
      bootstrappedCoxModel(coxFinal, modelName, modelFormula[[modelName]])
    )
  }
  
  # https://ete-online.biomedcentral.com/articles/10.1186/s12982-017-0060-8
  # this implements a paired data direct measure of hazard rate
  message("Paired data hazard rate")
  result$table5Data = matchesSelected %>% group_by(boot) %>% group_modify(function(d,g,...) {
    tibble(
      model = c("Unweighted PLME HR","Weighted PLME HR"),
      G= c(
        d %>% filter(time.neg < time.pos & status.neg==1) %>% nrow(),
        d %>% filter(time.neg < time.pos & status.neg==1) %>% summarise(G = sum(wt.neg)) %>% pull(G)
      ),
      H=c(
        d %>% filter(time.pos < time.neg & status.pos==1) %>% nrow(),
        d %>% filter(time.pos < time.neg & status.pos==1) %>% summarise(H = sum(wt.pos)) %>% pull(H)
      )
    ) %>% mutate(
      HR = G/H,
      estimate = log(G/H),
      std.error = sqrt(1/G + 1/H),
      HR.lower = exp(estimate-1.96*std.error),
      HR.upper = exp(estimate+1.96*std.error),
      HR.pValue = pnorm(0,estimate,std.error),
      term = "sGeneNegative"
    )
  })
  
  ## Pull results
  ## high level stats about number of deaths in each arm of control study where 
  combinedSummary1 = coxFinal %>% group_by(sGene,boot) %>% summarise(pairs = n()) %>% left_join(
    coxFinal %>% filter(diedWithin28days) %>% group_by(sGene,boot) %>% summarise(deaths = n()),
    by = c("sGene","boot")
  )
  result$table1Data = combinedSummary1

  ## summary of pairwise matching data set ----
  message("Summarising results")
  combinedSummary2 = summaryTable(posMatched, bootstraps) %>% rename_with(~paste0("S pos ",.x),.cols=c(N,`%age`,`mean (SD)`)) %>% #,by=c("category","value")) %>%
    full_join(summaryTable(negMatched, bootstraps) %>% rename_with(~paste0("S neg ",.x),.cols=c(N,`%age`,`mean (SD)`)),by=c("category","value")) %>% 
    full_join(summaryTable(coxFinal %>% filter(diedWithin28days), bootstraps) %>% rename_with(~paste0("Died ",.x),.cols=c(N,`%age`,`mean (SD)`)),by=c("category","value"))
  result$table2Data = combinedSummary2 %>% select(category,everything()) %>% arrange(category)

  # detect residual biases from matching with tolerances
  deltas = matchesSelected %>% group_by(boot) %>% summarise(
    ageDifference = mean(age.neg-age.pos),
    specimenDateDifference = mean(as.numeric(specimen_date.neg-specimen_date.pos))
  )
  result$table3Data = deltas
  
  return(result)
}





runAnalysisQuick = function(
  coxData, sGeneCtThreshold = ctThreshold, ctThreshold=30, 
  ageTolerance = 3, specimenDateTolerance = 5,
  matchOn = c("sex","imd_decile","LTLA_name","ethnicity_final","ageCat"), #"residential_category"),
  max = 450000*450000, includeRaw=FALSE, includeMatched=FALSE
) {
  result = list(params=tibble(
    sGeneCtThreshold = sGeneCtThreshold,
    ctThreshold = ctThreshold,
    ageTolerance = ageTolerance,
    specimenDateTolerance = specimenDateTolerance
  ))
  
  message("Interpreting data with CT threshold of S<",sGeneCtThreshold,":N<",ctThreshold,":ORF<",ctThreshold,
          "; tolerances - age=",ageTolerance,":specimen_date=",specimenDateTolerance)
  
  coxData3 = coxData %>% interpretSGene(S_CT = sGeneCtThreshold,N_CT = ctThreshold,ORF1ab_CT = ctThreshold)
  if (includeRaw) result$rawData = coxData3
  
  matches = coxData3 %>% createPairedMatches(ageTolerance = ageTolerance, specimenDateTolerance = specimenDateTolerance, matchOn = matchOn, max=max)
  
  rm(coxData3)
  
  # calculate a weighing for each side of each match based on frequecy
  matches = matches %>% 
    group_by(FINALID.pos) %>% mutate(wt.pos = 1/n()) %>%
    group_by(FINALID.neg) %>% mutate(wt.neg = 1/n()) %>%
    ungroup()
  
  if (includeMatched) {
    result$weightedMatchedData = matched
  }
  
  # https://ete-online.biomedcentral.com/articles/10.1186/s12982-017-0060-8
  # TODO: explain weighing modification to pairwise matched HR calculation to account for over-matching
  PLME_HR = matches %>% 
    mutate(
      weightedG = ifelse(time.neg < time.pos & status.neg==1, wt.neg, 0), 
      weightedH = ifelse(time.pos < time.neg & status.pos==1, wt.pos, 0),
    ) %>%
    summarise(
      G=sum(weightedG),
      H=sum(weightedH)
    ) %>% mutate(
      mean.HR = G/H,
      mean.estimate = log(G/H),
      mean.std.error = sqrt(1/G + 1/H),
      mean.HR.lower = exp(mean.estimate-1.96*mean.std.error),
      mean.HR.upper = exp(mean.estimate+1.96*mean.std.error),
      mean.HR.pValue = pnorm(0,mean.estimate,mean.std.error),
      model = "Weighted PLME HR",
      term = "sGeneNegative"
    )
  
  return(PLME_HR)
}




## Formatting functions ----
summaryTable = function(tmp,bootstraps=1) {
  groupPercent = function(df) df %>% summarise(N=as.integer(round(n()/bootstraps))) %>% mutate(`%age`=sprintf("%1.1f%%",N/sum(N)*100))
  groupMean = function(df,col) df %>% filter(is.finite({{col}})) %>% summarise(N=as.integer(round(n()/bootstraps)),  `mean (SD)`=sprintf("%1.1f (\u00B1%1.1f)",mean({{col}},na.rm=TRUE),sd({{col}},na.rm=TRUE)))
  bind_rows(
    tmp %>% group_by(value = NHSER_name) %>% groupPercent() %>% mutate(category="Region"),
    tmp %>% group_by(value = ethnicity_final) %>% groupPercent() %>% mutate(category="Ethnicity"),
    tmp %>% group_by(value = sex) %>% groupPercent() %>% mutate(category="Gender"),
    tmp %>% group_by(value = sGeneEra) %>% groupPercent() %>% mutate(category="S gene"),
    tmp %>% group_by(value = deathStatus) %>% groupPercent() %>% mutate(category="Status"),
    tmp %>% group_by(value = ageCat) %>% groupPercent() %>% mutate(category="Age by category"),
    tmp %>% groupMean(age) %>% mutate(category="Age"),
    tmp %>% groupMean(CT_N) %>% mutate(category="N gene CT"),
    tmp %>% group_by(value=imd_decile_2) %>% groupPercent() %>% mutate(category="IMD")
  ) %>% ungroup()
}


summariseAnalysisRun = function(run) {
  if(is.data.frame(run)) {
    tmp = run
  } else {
    tmp = run$table1Data
  }
  tmp = tmp %>% pivot_wider(names_from = sGene, values_from = deaths, names_prefix="deaths.")
  tmp = tmp %>% summarise(mean.matchedPairs = mean(pairs), mean.matchedDeaths.Positive = mean(deaths.Positive), mean.matchedDeaths.Negative = mean(deaths.Negative))
  tmp = tmp %>% bind_cols(
    run$table3Data %>% summarise(mean.ageDifference = mean(ageDifference), mean.specimenDateDifference = mean(specimenDateDifference)),
  )
  # combining bootstrap results
  tmp2 = run$table4Data %>% combineBootstraps()
  
  tmp2 = tmp2 %>% left_join(tmp, by=character()) %>% left_join(run$params, by=character())
  return(tmp2 %>% select(model,everything()))
}

combineBootstraps = function(analysisResult) {
  analysisResult %>% group_by(model,term) %>% summarise(
    across(
      .cols = c(estimate,std.error,HR.pValue), 
      .fns = c(mean=mean,sd=sd), .names="{.fn}.{.col}" # calculate mean and SD for all bootstrapped parameters
    )
  ) %>% mutate(
    mean.HR = exp(mean.estimate),
    mean.HR.lower = exp(mean.estimate-1.96*mean.std.error),
    mean.HR.upper = exp(mean.estimate+1.96*mean.std.error)
  )
}

prettyPrintSummary = function(summaryDf) {
  out = summaryDf %>% group_by(model,term) %>% summarise(
    `Hazard rate (95% CI)`=sprintf("%1.1f (%1.1f \u2013 %1.1f)",mean.HR,mean.HR.lower,mean.HR.upper),
    `p value`=scales::pvalue(mean.HR.pValue,0.001)) %>%
    mutate(
      Value = case_when(
        term %>% stringr::str_starts("sGene") ~ term %>% stringr::str_remove("sGene"),
        term %>% stringr::str_starts("sex") ~ term %>% stringr::str_remove("sex"),
        term %>% stringr::str_starts("imd_decile") ~ term %>% stringr::str_remove("imd_decile"),
        term %>% stringr::str_starts("ethnicity_final") ~ term %>% stringr::str_remove("ethnicity_final"),
        TRUE~NA_character_),
      Predictor = case_when(
        term %>% stringr::str_starts("sGene")~"S gene status",
        term == "CT_N"~"N gene CT value",
        term == "age"~"Age",
        term %>% stringr::str_starts("sex")~"Gender",
        term %>% stringr::str_starts("imd_decile")~"IMD",
        term %>% stringr::str_starts("ethnicity_final")~"Ethnicity",
        TRUE~NA_character_
      ),
      Model = model) %>%
    ungroup() %>%
    select(-model,-term)
  out = out %>% bind_rows(
    out %>% select(Model) %>% distinct() %>% mutate(Value = "Positive (ref)",Predictor="S gene status",
                                                    `Hazard rate (95% CI)`="\u2014",`p value`="\u2014")
  ) %>%
    mutate(
      Predictor = Predictor %>% forcats::as_factor() %>% forcats::fct_relevel("S gene status"),
      Value = Value %>% forcats::as_factor() %>% forcats::fct_relevel("Positive (ref)"),
      Model = Model %>% forcats::as_factor() %>% forcats::fct_relevel("S gene only")
    )
  return(out %>% select(Model,Predictor,Value,`Hazard rate (95% CI)`,`p value`) %>% group_by(Model,Predictor) %>% arrange(Model,Predictor,Value))
}

