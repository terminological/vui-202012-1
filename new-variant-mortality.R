
library(tidyverse)
library(patchwork)
library(survival)

setwd("~/Git/vui-202012-1/")
source("./provenance.R")

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
  return(tmp %>% mutate(FINALID=as.numeric(finalid)) %>% dplyr::ungroup() %>% p_comment(glue::glue("Deaths from: {path}",path=path)) %>% p_status())
}

#' @description Load line list
#' 
#' @return raw line list data set
getSGeneLineList = function(path,...) {
  tmp = readr::read_csv(path)
  return(tmp %>% dplyr::ungroup() %>% p_comment(glue::glue("S gene from: {path}",path=path)) %>% p_status())
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
        TRUE ~ "Other/Unknown"
      )
    ) %>% dplyr::ungroup() %>% p_comment(glue::glue("Cases from: {path}",path=path)) %>% p_status())
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

loadData = function(dir="~/Data/new-variant",date, censorLength = 28, B117Date = as.Date("2020-10-01"), earliestDate = B117Date, latestDate = NULL, minAge = 30) {
  
  ll = getLineList(path=paste0(dir,"/Anonymised Combined Line List ",date,".csv"))
  dll = getDeathsLineList(path = paste0(dir,"/",date," COVID19 Deaths.xlsx"))
  sgll = getSGeneLineList(path = paste0(dir,"/SGTF_linelist_",date,".csv"))
  
  sgllUniq = sgll %>% 
    p_clear() %>%
    p_status(TRUE ~ "{count} S gene test results") %>%
    p_exclude(is.na(FINALID) ~ "{count} tests with no id") %>%
    p_modify(
      function(d) {
        d %>% group_by(FINALID) %>% 
          mutate(sGeneResultCount = n()) %>% 
          filter(sGeneResultCount == 1) %>% 
          ungroup()
      },
      .message = "removing {count.in-count.out} multiple tests\n{count.out} patients with single results"
    )
  
  if (identical(latestDate,NULL)) latestDate = max(dll$dod,na.rm=TRUE)
  
  # Is this a potential source of bias?
  # excludedSgll = sgll %>% anti_join(sgllUniq, by = "FINALID")
  # if(nrow(excludedSgll) > 0) {
  #   message("we excluded ",excludedSgll %>% pull(FINALID) %>% unique() %>% length()," patients with multiple S Gene test results in pillar 2")
  #   excludedDeaths = dll %>% inner_join(excludedSgll, by="FINALID", suffix=c("",".sgene")) #%>% filter(specimen_date > "2020-10-01" & age>30)
  #   message("of these ",excludedDeaths %>% filter(specimen_date.sgene > B117Date) %>% pull(FINALID) %>% unique() %>% length()," died during the study period")
  # }
  
  # get the S gene line list and remove all patients who have more than one S Gene result. Otherwise we don;t have a single source in time for when to start the survival clock.
  # interpretSGene assigns sGene "Positive" for S+N+ORF+ with cutOff for S at 40 and N & ORF at 30; Negative is S-N+ORF+; equivocal is anything else (e.g. S-N-ORF+)
  # We exclude anything which does not match the line list.
  
  # Join the line list, deaths line list and uniquified sGene line list by FINALID
  coxData = ll %>% 
    left_join(sgllUniq, by="FINALID", suffix=c("",".sgene")) %>% 
    left_join(dll %>% mutate(FINALID=as.numeric(finalid)), by=c("FINALID"),suffix=c("",".death")) %>% 
    p_copy(sgllUniq) %>%
    p_exclude(
      specimen_date < earliestDate ~ "excluded {count} samples taken before {earliestDate}",
      specimen_date > latestDate ~ "excluded {count} samples taken after {latestDate}",
    ) %>%
    p_status(
      !is.na(specimen_date.sgene) ~ "{count} patients with known S gene result",
      is.na(specimen_date.sgene) ~ "{count} patients with unknown S gene result",
      !is.na(dod) ~ "{count} who died",
      is.na(dod) ~ "{count} who remain alive"
    )
  
  coxData2 = coxData %>% 
    p_select(
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
    p_mutate(
      died = !is.na(dod),
      diedWithin28days = died & as.numeric(dod-specimen_date) <= censorLength, # This is a definition of death within 28 days not relying on line list definition
      dodAt28days = ifelse(diedWithin28days, dod, NA), # this is the date of death (or NA if it is >28 days after sgene specimen date or there is no specimen)
      censoredDate = pmin(latestDate, specimen_date+censorLength,na.rm = TRUE) # earliest of last recorded death or specimen date
    ) %>% 
    p_mutate(
      reportingDelay = as.numeric(lab_report_date-specimen_date),
      sGeneDelay = as.numeric(specimen_date.sgene-specimen_date),
      admissionDelay = as.numeric(dateadmission_NHSE-specimen_date),
      time = ifelse(diedWithin28days, dodAt28days-as.numeric(specimen_date), censoredDate-as.numeric(specimen_date)), 
      status = ifelse(!diedWithin28days,0,1)
    ) %>%
    p_mutate(
      ageCat = cut(age, breaks = c(-Inf,30,60,70,80,Inf), labels = c("<30","30-59","60-69","70-79","80+"), right=FALSE),
      LTLA_name =  as.factor(LTLA_name),
      preB117 = specimen_date < B117Date
    )
  

  
  # tmpCox = coxData2 %>% interpretSGene()
  # tmpCox %>% filter(sGene=="Equivocal") %>% nrow() %>% message(" with equivocal results excluded")
  # tmpCox %>% filter(sGene=="Negative") %>% nrow() %>% message(" S- cases")
  # tmpCox %>% filter(sGene=="Positive") %>% nrow() %>% message(" S+ controls")
  
  # reportExclusions = function(cox,type) {
  #   type = paste0(" ",type)
  #   cox %>% filter(is.na(imd_decile)) %>% nrow() %>% message(type," with missing IMD")
  #   cox %>% filter(is.na(sex) | sex=="Unknown") %>% nrow() %>% message(type," with missing or unknown gender")
  #   cox %>% filter(is.na(age)) %>% nrow() %>% message(type," with unknown age")
  #   cox %>% filter(age<minAge) %>% nrow() %>% message(type," with age <",minAge)
  #   cox %>% filter(is.na(ethnicity_final)) %>% nrow() %>% message(type," with unknown ethnicity")
  #   cox %>% filter(admissionDelay < 0) %>% nrow() %>% message(type," admitted before pillar 2 test taken")
  #   cox %>% filter(reportingDelay < 0 | reportingDelay>=20) %>% nrow() %>% message(type," whose tests are reported before specimen taken or >19 days after specimen taken")
  #   cox %>% filter(sGeneDelay < 0 | sGeneDelay>=20) %>% nrow() %>% message(type," whose sGene tests are >19 days after original positive test")
  #   cox %>% filter(time<=0) %>% nrow() %>% message(type," whose time from test to event (death) is 0 or less (e.g. postmortem tests)")
  # }
  # 
  # tmpCox %>% filter(sGene=="Negative") %>% reportExclusions("S- cases")
  # tmpCox %>% filter(sGene=="Positive") %>% reportExclusions("S+ controls")
  # 
  coxData3 = coxData2 %>%
    p_exclude(
      is.na(imd_decile) ~ "{count} with missing IMD",
      is.na(sex) | sex == "Unknown" ~ "{count} with missing or unknown gender",
      is.na(age) ~ "{count} with unknown age",
      age<minAge ~ "{count} with age<{minAge}",
      is.na(ethnicity_final) ~ "{count} with unknown ethnicity",
      admissionDelay < 0 ~ "{count} admitted before pillar 2 test taken",
      reportingDelay < 0 | reportingDelay>=20 ~ "{count} with reporting delay <0 or >19 days",
      sGeneDelay < 0 | sGeneDelay>=20 ~ "{count} whose sGene tests are <0 or >19 days after original positive test",
      time<=0 ~ "{count} whose time from test to event (death) is <0 (e.g. postmortem tests)",
      is.na(LTLA_name) ~ "{count} with unknown location"
      ) %>%
    p_mutate(
      ethnicity_final = ethnicity_final %>% forcats::as_factor() %>% forcats::fct_relevel("White"),
      sex = sex %>% forcats::as_factor(),
      imd_decile = forcats::fct_relevel(forcats::as_factor(as.numeric(imd_decile)),"5"),
      imd_decile_2 = forcats::as_factor(as.numeric(imd_decile)),
      deathStatus = ifelse(diedWithin28days,"Dead <28 days","Other")
    )
  
  # tmpCox = coxData2 %>% interpretSGene()
  # tmpCox %>% filter(sGene=="Negative") %>% nrow() %>% message(" S- cases eligible for inclusion")
  # tmpCox %>% filter(sGene=="Positive") %>% nrow() %>% message(" S+ controls eligible for inclusion")
  
  ## quality control
  
  
  return(coxData3)
}

## Analysis functions: ----

# function to help interpret S gene CT values in context of N gene and ORF gene to give S gene status. 
# With the defaults this produces the same result as the sgtf_30 column in the sounrge SGTF line list
interpretSGene = function(sGeneLineList, S_CT = 30, ORF1ab_CT = 30, N_CT = 30, Control_CT = Inf) {
  params = glue::glue("CT thresholds S<={S_CT}:N<={N_CT}:ORF<={ORF1ab_CT}")
  sGeneLineList %>% 
    p_mutate(
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
    ) %>% 
    p_mutate(
      result = ifelse(!Control_pos,"No control",paste0(ifelse(S_pos,"S+","S-"),ifelse(N_pos,"N+","N-"),ifelse(ORF1ab_pos,"ORF+","ORF-")))
    ) %>% 
    p_mutate(
      sGeneEra = case_when(
        is.na(sGene) ~ "Unknown",
        sGene == "Negative" & !preB117 ~ "Neg Post B.1.1.7",
        sGene == "Negative" & preB117 ~ "Neg Pre B.1.1.7",
        TRUE ~ as.character(sGene)
      ),
      sGene = sGene %>% forcats::fct_relevel("Positive"),
      sGeneEra = sGeneEra %>% forcats::as_factor() %>% forcats::fct_relevel("Positive"),
      relativeCopyNumber = 2^(median(CT_N,na.rm=TRUE)-CT_N)
    ) %>%
    p_status(
      .headline = params,
      sGene == "Unknown" ~ "{count} with unknown S status",
      sGene == "Positive" ~ "{count} with S+ status",
      sGene == "Negative" ~ "{count} with S- status",
      sGene == "Equivocal" ~ "{count} with equivocal S status",
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
  m = paste0(matchOn,collapse=', ')
  out = coxData %>% 
    p_exclude(
      preB117 ~ "{count} predating B.1.1.7 emergence",
      sGene == "Unknown" ~ "{count} with unknown S status",
      sGene == "Equivocal" ~ "{count} with equivocal S status"
    ) %>%
    p_modify(
      .headline = "Matching by:\n{m}\nwith tolerances:\nage<={ageTolerance} & specimen<={specimenDateTolerance}",
      .message = "{count.in} patients \n resulting in {count.out} pairs",
      function(d) {
        posSgene = d %>% filter(sGene=="Positive")
        negSgene = d %>% filter(sGene=="Negative")
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
          ) %>% 
          group_by(FINALID.neg) %>% mutate(order.pos = n()) %>%
          group_by(FINALID.pos) %>% mutate(order.neg = n()) %>%
          ungroup() %>% mutate(
            mapType.neg = paste0(order.neg,":",order.pos),
            mapType.pos = paste0(order.pos,":",order.neg)
          )
        
        return(matches)
      }
    )
  return(out)
}

uniquifyPairs = function(pairedMatches, resolutionStrategy, bootstraps, ratio.neg=NULL, ratio.pos=NULL) {
  message = glue::glue("Selecting unique pairs by {resolutionStrategy} with {bootstraps} replicates.\n",
                      "Before sampling the matched set contains:\n",
                      "{n_distinct(pairedMatches$FINALID.neg)} unique S-\n",
                      "and {n_distinct(pairedMatches$FINALID.pos)} unique S+")
  
  if (!is.null(ratio.neg) & ! is.null(ratio.pos)) {
    pairedMatched = pairedMatches %>%
        p_exclude(
          !(order.pos <= ratio.pos) ~ "{count} pairs where S+ is connected to more than {ratio.pos} S-",
          !(order.neg <= ratio.neg) ~ "{count} pairs where S- is connected to more than {ratio.neg} S+"
        )
  }
  
  out = pairedMatches %>% p_modify(
    .headline = message,
    .message = "{count.in} pairs input:\n{count.out} pairs over all replicates",
    function(d) {
  
  
      if (resolutionStrategy == "bootstrap") {
        warning("This option is here for legacy purposes. It is superceded by the edge sample option which does what this should have done")
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
          tmp = d %>% ungroup() %>%
            mutate(rnd = runif(nrow(d))) %>% 
            arrange(rnd) %>% # put the dataset in a random order
            group_by(FINALID.neg) %>% filter(row_number()==1) %>% # filter out a single matching at random from pos to neg
            group_by(FINALID.pos) %>% filter(row_number()==1) %>%  # filter out a single matching at random from neg to pos
            # N.B. I'm not sure this is as random as it looks
            # It selects random matches for LHS and then selects random reverse matches that are paired with
            # LHS random match. however this leads to situations where there is less chance of
            # a RHS appearing if it matches a LHS that is over-represented.
            # it's not symmetrical.
            ungroup() %>% mutate(boot = i)
          matchesSelected = bind_rows(matchesSelected,tmp)
        }
        message("..finished")
        
      } else if(resolutionStrategy == "none") {
        # ignore duplicates
        matchesSelected = d %>% mutate(boot=1) %>% ungroup()
        
      } else if(resolutionStrategy == "edge sample") {
        
        # Get the edge list
        # randomly select a single match out of possible matches for each bootstrap
        set.seed(101)
        matchesSelected = NULL
        message("bootstrapping",appendLF = FALSE)
        for(i in 1:bootstraps) {
          message("..",i,appendLF = FALSE)
          
          tmp = d %>% ungroup() %>%
            # randomly order the pairs
            mutate(ordering = sample.int(n=nrow(d))) %>% 
            # include only the first of each of the S negs. excluding all others, essentially removing duplicates.
            group_by(FINALID.neg) %>% mutate(
              chosen.neg = ifelse(ordering==min(ordering), "inc", "exc")
            ) %>% 
            # include only the first of each of the S pos. excluding all others, essentially removing duplicates.
            group_by(FINALID.pos) %>% mutate( 
              chosen.pos = ifelse(ordering==min(ordering), "inc", "exc")
            ) %>% ungroup() %>% 
            # pick the first matches of both. This will pick out a pair if both ends are ranked low.
            filter((chosen.neg=="inc" & chosen.pos=="inc")) %>%
            mutate(boot = i)
          matchesSelected = bind_rows(matchesSelected,tmp)
        }
        message("..finished")
        
      } else if(resolutionStrategy == "node sample") {
    
        # Get a list of unique nodes S+ and S- combined
        nodes = bind_rows(
          d %>% select(FINALID = FINALID.neg),# %>% mutate(side="neg"),
          d %>% select(FINALID = FINALID.pos)# %>% mutate(side="pos"),
        ) %>% distinct()
        
        # initialise bootstrapping loop
        matchesSelected = NULL
        set.seed(101)
        message("bootstrapping",appendLF = FALSE)
        for (i in 1:bootstraps) {
          message("..",i,appendLF = FALSE)
          #browser()
          
          # We are picking nodes at random from the whole set of possible nodes regardless of Sneg or Spos status
          # all we are doing here is assigning a random order to the node list
          nodes = nodes %>% mutate(ordering = sample.int(n=nrow(nodes)))
          
          # create a directed edge list for this ordered set of nodes
          # by making an undirected edgelist and picking only edges that go from lower to higher in rank order
          # This duplicates the edges by adding in the reverse direction
          edges = bind_rows( 
              d %>% select(FINALID.source = FINALID.neg, FINALID.target = FINALID.pos) %>% mutate(direction="neg->pos"),
              d %>% select(FINALID.source = FINALID.pos, FINALID.target = FINALID.neg) %>% mutate(direction="pos->neg")
            )
          
          # Assign the node rank to both sides of the edge
          edges = edges %>%
            inner_join(nodes %>% rename(ordering.source=ordering), by=c("FINALID.source"="FINALID")) %>%
            inner_join(nodes %>% rename(ordering.target=ordering), by=c("FINALID.target"="FINALID")) 
          
          # This filter de-duplicates the edges and creates an directed graph. from source to target (source can be S- or S+)
          edges = edges %>%
            filter(ordering.source < ordering.target) 
          
          # Iteratively select edges for removal from edgelist
          output = NULL
          while( nrow(edges) > 0) {
            # for every FINALID.source(==ordering.source) pick out the first edge target combination
            # i.e. for every source pick edges where the target is the lowest rank entry for a given source (i.e. the first one)
            selected = edges %>% group_by(ordering.source) %>% filter(ordering.target == min(ordering.target))
            # if an edge is selected its target cannot be a source of another edge, so deselect any where this happens
            # this also takes care of the reverse situation where the source of one edge is the target of another edge.
            # but this is also prevented by the filter which removes edges with source > target
            selected = selected %>% anti_join(selected, by=c("ordering.source" = "ordering.target"))
            # there are also situation where we have selected the same target for multiple sources.
            # in this case we want to select only the highest ranked source
            selected = selected %>% group_by(ordering.target) %>% filter(ordering.source == min(ordering.source))
            # the selection is guaranteed to have distinct sources and targets
            if(any(duplicated(selected$FINALID.source)) | any(duplicated(selected$FINALID.target))) stop("Aargh...")
            # take all edges containing any selected source or target nodes away from edges.
            edges = edges %>% 
              # remove selected source nodes from edge list
              anti_join(selected, by=c("ordering.source")) %>% 
              # remove selected target nodes that are sources in the edge list
              anti_join(selected, by=c("ordering.source"="ordering.target")) %>%
              # remove selected target nodes from edge list
              anti_join(selected, by=c("ordering.target")) %>% 
              # remove selected source nodes that are targets in the edge list
              anti_join(selected, by=c("ordering.target"="ordering.source"))
              
            if (nrow(selected)==0) break #lets hope this doesn't happen if it does we just probably have to abort
            output = output %>% bind_rows(selected)
          }
          #output now has edgelist of selected edges in it but the source and targets need to be mapped back to original order
          output = output %>% mutate(
            FINALID.neg = ifelse(direction == "neg->pos", FINALID.source, FINALID.target),
            FINALID.pos = ifelse(direction == "pos->neg", FINALID.source, FINALID.target)
          )
          # selected output rows from d
          thisIteration = d %>% semi_join(output,by=c("FINALID.neg","FINALID.pos")) %>% mutate(boot=i)
          matchesSelected = matchesSelected %>% bind_rows(thisIteration)
        }
        rm(nodes,edges,thisIteration,output,selected)
        message("..finished")
        
        
      } else {
        stop("No such strategy: ",resolutionStrategy)
      }
      return(matchesSelected)
    })
  message("After sampling the matched set contains ",n_distinct(out$FINALID.neg)," unique S- and ",n_distinct(out$FINALID.pos)," unique S+")
  return(out)
}

pairedMatchesToCoxData = function(matches, matchOn) {
  out = matches %>% p_modify(
    .message = "{count.in} pairs converted to {count.out} individuals",
    function (d) {
      if (!"boot" %in% colnames(d)) d = d %>% mutate(boot=1)
      sPos = d %>% select(boot,ends_with(".pos"),all_of(matchOn)) %>% 
        rename_with(.fn=function(x) stringr::str_remove(x,".pos"),.cols=ends_with(".pos")) %>% 
        mutate(sGene="Positive") %>% distinct()
      sNeg = d %>% select(boot,ends_with(".neg"),all_of(matchOn)) %>% 
        rename_with(.fn=function(x) stringr::str_remove(x,".neg"),.cols=ends_with(".neg")) %>% 
        mutate(sGene="Negative") %>% distinct()
      return(bind_rows(sPos,sNeg))
    }
  ) %>% p_group_by(sGene) %>% p_status(
    TRUE ~ "{count} in all replicates",
    diedWithin28days ~ "of whom {count} died within 28 days",
    !diedWithin28days ~ "and {count} still alive at 28 days",
    boot==1 ~ "{count} in 1st replicate",
    diedWithin28days & boot==1 ~ "of whom {count} died within 28 days",
    !diedWithin28days & boot==1 ~ "and {count} were still alive at 28 days"
  ) %>% p_ungroup()
  return(out)
}

runAnalysis = function(
    coxData, sGeneCtThreshold = ctThreshold, ctThreshold=30, 
    bootstraps=10, ageTolerance = 5, specimenDateTolerance = 1,
    modelFormula = list(`S gene + age`=Surv(time,status) ~ sGene+age),
    max = 450000*450000, includeRaw=FALSE, includeMatched=FALSE,
    matchOn = c("sex","imd_decile","LTLA_name","ethnicity_final","ageCat"), #,"residential_category"),
    resolutionStrategy = "none",ratio.neg=1,ratio.pos=1
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
  if (includeRaw) result$pairedMatchedUnfilteredData = matches
  
  # resolve duplicates in matched pairs
  matchesSelected = matches %>% uniquifyPairs(resolutionStrategy=resolutionStrategy, bootstraps = bootstraps, ratio.neg=ratio.neg, ratio.pos=ratio.pos)
  if (includeRaw) result$pairedMatchedData = matchesSelected
  
  # prepare data for cox PH models
  coxFinal = matchesSelected %>% pairedMatchesToCoxData(matchOn = matchOn) %>% 
    p_mutate(sGene = sGene %>% forcats::as_factor() %>% forcats::fct_relevel("Positive")) # equivocal category is dropped
  
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
  
  negMatched = coxFinal %>% filter(sGene=="Negative")
  posMatched = coxFinal %>% filter(sGene=="Positive")
  if(nrow(posMatched) == nrow(negMatched)) {
    message("Paired data hazard rate can be directly calculated as data is 1:1")
    result$table5Data = matchesSelected %>% group_by(boot) %>% group_modify(function(d,g,...) {
      tibble(
        model = c("PLME HR"),
        G= c(
          d %>% filter(time.neg < time.pos & status.neg==1) %>% nrow()
        ),
        H=c(
          d %>% filter(time.pos < time.neg & status.pos==1) %>% nrow()
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
  }
  
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





# runAnalysisQuick = function(
#   coxData, sGeneCtThreshold = ctThreshold, ctThreshold=30, 
#   ageTolerance = 3, specimenDateTolerance = 5,
#   matchOn = c("sex","imd_decile","LTLA_name","ethnicity_final","ageCat"), #"residential_category"),
#   max = 450000*450000, includeRaw=FALSE, includeMatched=FALSE
# ) {
#   result = list(params=tibble(
#     sGeneCtThreshold = sGeneCtThreshold,
#     ctThreshold = ctThreshold,
#     ageTolerance = ageTolerance,
#     specimenDateTolerance = specimenDateTolerance
#   ))
#   
#   message("Interpreting data with CT threshold of S<",sGeneCtThreshold,":N<",ctThreshold,":ORF<",ctThreshold,
#           "; tolerances - age=",ageTolerance,":specimen_date=",specimenDateTolerance)
#   
#   coxData3 = coxData %>% interpretSGene(S_CT = sGeneCtThreshold,N_CT = ctThreshold,ORF1ab_CT = ctThreshold)
#   if (includeRaw) result$rawData = coxData3
#   
#   matches = coxData3 %>% createPairedMatches(ageTolerance = ageTolerance, specimenDateTolerance = specimenDateTolerance, matchOn = matchOn, max=max)
#   
#   rm(coxData3)
#   
#   # calculate a weighing for each side of each match based on frequecy
#   matches = matches %>% 
#     group_by(FINALID.pos) %>% mutate(wt.pos = 1/n()) %>%
#     group_by(FINALID.neg) %>% mutate(wt.neg = 1/n()) %>%
#     ungroup()
#   
#   if (includeMatched) {
#     result$weightedMatchedData = matched
#   }
#   
#   # https://ete-online.biomedcentral.com/articles/10.1186/s12982-017-0060-8
#   # TODO: explain weighing modification to pairwise matched HR calculation to account for over-matching
#   PLME_HR = matches %>% 
#     mutate(
#       weightedG = ifelse(time.neg < time.pos & status.neg==1, wt.neg, 0), 
#       weightedH = ifelse(time.pos < time.neg & status.pos==1, wt.pos, 0),
#     ) %>%
#     summarise(
#       G=sum(weightedG),
#       H=sum(weightedH)
#     ) %>% mutate(
#       mean.HR = G/H,
#       mean.estimate = log(G/H),
#       mean.std.error = sqrt(1/G + 1/H),
#       mean.HR.lower = exp(mean.estimate-1.96*mean.std.error),
#       mean.HR.upper = exp(mean.estimate+1.96*mean.std.error),
#       mean.HR.pValue = pnorm(0,mean.estimate,mean.std.error),
#       model = "Weighted PLME HR",
#       term = "sGeneNegative"
#     )
#   
#   return(PLME_HR)
# }




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
    `Hazard rate (95% CI)`=sprintf("%1.1f (%1.1f \u2013 %1.1f)",
                                   ifelse(mean.HR<10,mean.HR,Inf),
                                   ifelse(mean.HR.lower<10,mean.HR.lower,Inf),
                                   ifelse(mean.HR.upper<10,mean.HR.upper,Inf)),
    `p value`=scales::pvalue(mean.HR.pValue,0.001)) %>%
    mutate(
      Value = case_when(
        term %>% stringr::str_starts("sGene") ~ term %>% stringr::str_remove("sGene"),
        term %>% stringr::str_starts("sex") ~ term %>% stringr::str_remove("sex"),
        term %>% stringr::str_starts("imd_decile") ~ term %>% stringr::str_remove("imd_decile"),
        term %>% stringr::str_starts("ethnicity_final") ~ term %>% stringr::str_remove("ethnicity_final"),
        term %>% stringr::str_starts("residential_category") ~ term %>% stringr::str_remove("residential_category"),
        TRUE~NA_character_),
      Predictor = case_when(
        term %>% stringr::str_starts("sGene")~"S gene status",
        term == "CT_N"~"N gene CT value",
        term == "age"~"Age",
        term %>% stringr::str_starts("sex")~"Gender",
        term %>% stringr::str_starts("imd_decile")~"IMD",
        term %>% stringr::str_starts("ethnicity_final")~"Ethnicity",
        term %>% stringr::str_starts("residential_category")~"Residence",
        TRUE~term
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
      Model = Model %>% forcats::as_factor() %>% forcats::fct_relevel("S gene + age")
    )
  return(out %>% select(Model,Predictor,Value,`Hazard rate (95% CI)`,`p value`) %>% group_by(Model,Predictor) %>% arrange(Model,Predictor,Value))
}

