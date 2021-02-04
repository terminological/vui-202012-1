
library(tidyverse)
library(patchwork)
library(survival)

setwd("~/Git/vui-202012-1/")
source("./provenance.R")

## Data loading functions ----

#' Loads the line list of deaths in the format provided by PHE
#'
#' @param path - where on the local filesystem
#' @param ... - not used
#' @return a dataframe of death data
#' @examples dll = getDeathsLineList("~/Data/new-variant/20210129 COVID19 Deaths.xlsx")
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

#' Loads the raw line list of s Gene results in the format provided by PHE
#'
#' @param path - where on the local filesystem
#' @param ... - not used
#' @return a dataframe of s gene data
#' @examples sgll = getSGeneLineList("~/Data/new-variant/SGTF_linelist_20210129.csv")
getSGeneLineList = function(path,...) {
  tmp = readr::read_csv(path)
  return(tmp %>% dplyr::ungroup() %>% p_comment(glue::glue("S gene from: {path}",path=path)) %>% p_status())
}

#' Loads the raw line list of cases data in the format provided by PHE
#' Normalise ethnicity to 4 categories, residential category to 3 categories
#'
#' @param path - where on the local filesystem
#' @param ... - not used
#' @return a dataframe of case data
#' @examples sgll = getLineList("~/Data/new-variant/SGTF_linelist_20210129.csv")#' /home/terminological/Data/new-variant/Anonymised Combined Line List 20210129.csv
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

#' Normalise a vector of strings specifying genders strings to a minimal format
#'
#' @param gender - string vector
#' @return a string vector
#' @examples tmp %>% mutate(gender = normaliseGender(ifelse(is.na(gender),"unknown",gender)))
normaliseGender = function(gender) {
  case_when(
    is.na(gender) ~ NA_character_,
    gender %>% stringr::str_detect("f|F") ~ "female",
    gender %>% stringr::str_detect("m|M") ~ "male",
    gender %>% stringr::str_detect("u|U") ~ "unknown",
    TRUE ~ "unknown")
}

#' Load data and assemble the data used in the analysis.
#' The resulting data includes the raw CT values.
#'
#' @param dir - the data directory
#' @param date - the data date as a %Y%m%d string which is how it is represented in the PHE files
#' @param censorLength - how long before you censor observations of death - default is 28
#' @param B117Date - the date before which B.1.1.7 was not circulating - default "2020-10-01"
#' @param earliestDate - the earliest date to consider in the analysis - defaults to B117Date
#' @param latestDate - the latest date to consider in the analysis - defaults to latest date in death line list
#' @param minAge - the smallest age to consider
#'
#' @return a data set including a time, and status column, for cox modelling, and raw CT values. The dataframe retains provenance (available through df %>% p_flowchart() )
#'
#' @examples coxData = loadData("~/Data/new-variant","20200129")
loadData = function(dir="~/Data/new-variant",date, censorLength = 28, B117Date = as.Date("2020-10-01"), earliestDate = B117Date, latestDate = NULL, minAge = 30) {
  
  ll = getLineList(path=paste0(dir,"/Anonymised Combined Line List ",date,".csv"))
  dll = getDeathsLineList(path = paste0(dir,"/",date," COVID19 Deaths.xlsx"))
  sgll = getSGeneLineList(path = paste0(dir,"/SGTF_linelist_",date,".csv"))
  
  # exclude any non unique results
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
  
  # Join the line list, deaths line list and uniquified sGene line list by FINALID
  coxData = ll %>% 
    left_join(sgllUniq, by="FINALID", suffix=c("",".sgene")) %>% 
    left_join(dll %>% mutate(FINALID=as.numeric(finalid)), by=c("FINALID"),suffix=c("",".death")) %>% 
    p_clear() %>%
    p_status(
      TRUE ~ "{count} linked patients"
    ) %>% 
    p_exclude(
      is.na(specimen_date.sgene) ~ "{count} patients with unknown S gene status",
      specimen_date < earliestDate ~ "{count} diagnosed before {earliestDate}",
      specimen_date > latestDate ~ "{count} diagnosed after {latestDate}",
    ) %>%
    p_status(
      TRUE ~ "{count} included patients",
      !is.na(dod) ~ "of whom {count} died",
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
      sgtf, sgtf_under30CT, # gene copy number interpretation from sGene line list
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
  
  # TODO: We may need to consider the data stratified by sGene status from this point in the future, but this would need to assume a CT threshold
  # although we could use official (sgtf, sgtf_under30CT) flags for this:
  # coxData2 = coxData2 %>% interpretSGene()
  # coxData2 = coxData2 %>% p_stratify(sgtf_under30CT) %>% p_status_summary(count = n(), .glue="{count} with s gene status {sgtf_under30CT}")
  
  coxData3 = coxData2 %>%
    p_exclude(
      is.na(imd_decile) ~ "{count} with missing IMD",
      is.na(sex) | sex == "Unknown" ~ "{count} with missing or unknown gender",
      is.na(age) ~ "{count} with unknown age",
      age<minAge ~ "{count} with age<{minAge}",
      is.na(ethnicity_final) ~ "{count} with unknown ethnicity",
      admissionDelay < 0 ~ "{count} admitted before pillar 2 test taken",
      reportingDelay < 0 | reportingDelay>=20 ~ "{count} with reporting delay <0 or >19 days",
      sGeneDelay < 0 | sGeneDelay>=20 ~ "{count} with sGene <0 or >19 days after first test",
      time<=0 ~ "{count} whose outcome occurred before sGene test",
      is.na(LTLA_name) ~ "{count} with unknown location"
      ) %>%
    p_mutate(
      ethnicity_final = ethnicity_final %>% forcats::as_factor() %>% forcats::fct_relevel("White"),
      sex = sex %>% forcats::as_factor(),
      imd_decile = forcats::fct_relevel(forcats::as_factor(as.numeric(imd_decile)),"5"),
      imd_decile_2 = forcats::as_factor(as.numeric(imd_decile)),
      deathStatus = ifelse(diedWithin28days,"Dead <28 days","Other")
    ) %>%
    p_status(
      TRUE ~ "{count} passing data quality checks",
      !is.na(dod) ~ "of whom {count} died",
    )
  
  return(coxData3)
}

## Analysis functions: ----

#' Interpret S gene status according to various cut off values
#' function to help interpret S gene CT values in context of N gene and ORF gene to give S gene status. 
#' With the defaults this produces the same result as the sgtf_30 column in the source SGTF line list
#' Defaults are S:30,ORF:30,N:30,Control:Inf
#'
#' @param sGeneLineList - a dataframe includeing 
#' @param S_CT - S gene detected when P2CH3CQ <= this value
#' @param ORF1ab_CT - ORF1ab gene detected when P2CH1CQ <= this value
#' @param N_CT - N gene detected when P2CH2CQ <= this value
#' @param Control_CT - control sample is positive when P2CH4CQ <= this value
#'
#' @return - the same dataframe with additional columns including "sGene" and "result"
#'
#' @examples coxData = coxData %>% interpretSGene()
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
      # sGene == "Unknown" ~ "{count} with unknown S status",
      sGene == "Positive" ~ "{count} with S+ status",
      sGene == "Negative" ~ "{count} with S- status",
      sGene == "Equivocal" ~ "{count} with equivocal S status",
    )
}



#' Matches coxData entries with each other to form pairs with similar properties but differing by S gene status (Positives and Negatives only)
#'
#' @param coxData - data frame to match containing all the colums specified in matchOn and age, specimen_date, sGene, status (died or censored), time and FINALID as a minimum
#' @param ageTolerance - maximum delta in year of age allowed in matching pairs. an age 75 and 80 will match if ageTolerance is 5 (80-75)<=5
#' @param specimenDateTolerance - maximum delta in specimen date allowed
#' @param matchOn - a set of columns to perform exact matching on as a character vector, e.g. c("sex","imd_decide","LTLA_name","ethnicity_final",...)
#' @param max - a safety net for max numbers of items to attempt to join.
#'
#' @return - essentially an edgelist or pair list consisting of the matchOn columns and additional columns mentioned above.
#'
#' @examples
createPairedMatches = function(coxData, ageTolerance, specimenDateTolerance,matchOn, max) {
  m = paste0(matchOn,collapse=', ')
  out = coxData %>% 
    p_exclude(
      # preB117 ~ "{count} predating VOC emergence",
      # sGene == "Unknown" ~ "{count} with unknown S status",
      sGene == "Equivocal" ~ "{count} with equivocal S status"
    ) %>%
    p_modify(
      .headline = "Matching by:\n{m}\nwith tolerances: age \u00b1{ageTolerance} years & sample date \u00b1{specimenDateTolerance} days",
      .message = "from {count.in} individuals, we found {count.out} pairs",
      function(d) {
        
        d = d %>% select(all_of(matchOn),age,specimen_date,FINALID,sGene,status,time)
        
        posSgene = d %>% filter(sGene=="Positive") %>% 
          rename_with(.cols=!all_of(matchOn), .fn = function(x) paste0(x,".pos")) %>% 
          data.table::as.data.table()
        posSgene %>% data.table::setkeyv(matchOn)
        negSgene = d %>% filter(sGene=="Negative") %>% 
          rename_with(.cols=!all_of(matchOn), .fn = function(x) paste0(x,".neg")) %>% 
          data.table::as.data.table()
        negSgene %>% data.table::setkeyv(matchOn)
        
        if(as.numeric(nrow(negSgene)) * as.numeric(nrow(posSgene)) > max) {
          stop("this would try and create a very large cross join of ",nrow(posSgene),"x",nrow(negSgene))
        }
        message("Joining: ",nrow(posSgene),"x",nrow(negSgene))
        
        # pairwise matching
        ## case matching: 
        # match cases by "sex","imd_decile","LTLA_name","ethnicity_final","ageCat"
        # then ensure age difference < 2 years & specimen date < 4 days
        # select S negatives and match to S positives
        # randomly select a single match out of possible matches
        matches = negSgene[posSgene, on=matchOn, nomatch=0, allow.cartesian=TRUE]
        matches = matches[abs(age.neg-age.pos) <= ageTolerance][abs(as.numeric(specimen_date.neg-specimen_date.pos)) <= specimenDateTolerance]
        # this is done with data.table for performance reasons. the original dplyr code was:
        # matches = negSgene %>% 
        # inner_join(
        #   posSgene, 
        #   by=matchOn,suffix=c(".neg",".pos")) %>% 
        # filter(
        #   abs(age.neg-age.pos) <= ageTolerance & 
        #     abs(as.numeric(specimen_date.neg-specimen_date.pos)) <= specimenDateTolerance
        # ) 
        
        matches = matches %>% 
          group_by(FINALID.neg) %>% mutate(order.pos = n()) %>%
          group_by(FINALID.pos) %>% mutate(order.neg = n()) %>%
          ungroup() %>% mutate(
            mapType.neg = paste0(order.neg,":",order.pos),
            mapType.pos = paste0(order.pos,":",order.neg)
          )
        
        rm(posSgene,negSgene)
        gc()
        
        return(matches)
      }
    )
  return(out)
}


#' Resolution of multiple matches in the dataset.
#' We find that especially in the early and late stages of the transition to S Negative, there are some patients who closely match large numbers of other patients.
#' This is doubly true the more we relax our matching tolerances.
#' Ideally we don't individuals patients in the data set multiple times (although one school of thought is that given they are matched multiple times it doesn;t matter)
#' There are 3 strategies to remove duplicates:
#' "none" - don't remove them
#' "node sample" - randomly sample the unique cases and controls who have not yet been selected, pick at random a match associated with that patient. select that pair, ensuring both case and control are unique
#' "edge sample" - randomly sample the pairs not yet selected, select a pair and make sure that any other pair involving either case or control is removed.
#' In testing there is not much difference between the last 2 options in analysis.
#'
#' @param pairedMatches a dataframe of paired matches
#' @param resolutionStrategy a name of a resolution strategy to adopt - currently supported are "none","edge sample","node sample", and a legacy - "bootstrap" option which was out first attempt at this.
#' @param bootstraps - for edge or node sampling the number of replicates to produce - for none there will only ever be 1 replicate
#' @param ratio.neg - a filter to exclude S Negatives that are matched to a large number of S Positive - these can be over-represented in the early stage of the transition to S Negative
#' @param ratio.pos - a filter to exclude S Positives that are matched to a large number of S Negatives - these can be over-represented in the late stage of the transition to S Negative
#'
#' @return an paired list that (depending on options) includes unique cases and controls. This has a "boot" column which defines which replicate it is. For very many bootstraps this might get unmanageably large.
#'
#' @examples
uniquifyPairs = function(pairedMatches, resolutionStrategy, bootstraps, ratio.neg=NULL, ratio.pos=NULL) {
  message = glue::glue("Selecting unique pairs by {resolutionStrategy} with {bootstraps} replicates.\n",
                      "Before sampling the matched set contains:\n",
                      "{n_distinct(pairedMatches$FINALID.neg)} unique S- cases and {n_distinct(pairedMatches$FINALID.pos)} unique S+ controls")
  
  if (!is.null(ratio.neg) & ! is.null(ratio.pos)) {
    pairedMatched = pairedMatches %>%
        p_exclude(
          !(order.pos <= ratio.pos) ~ "{count} pairs where S+ is connected to more than {ratio.pos} S-",
          !(order.neg <= ratio.neg) ~ "{count} pairs where S- is connected to more than {ratio.neg} S+"
        )
  }
  
  out = pairedMatches %>% p_modify(
    .message = message,
    function(d) {
  
  
      if (resolutionStrategy == "bootstrap") {
        warning("This option is here for legacy purposes. It is superceded by the edge sample option which does what this should have done")
        # NB: I am not sure if there is a risk here of a systematic bias, or at least a random bias that is systematically present
        # throughout the analysis. What happens here is the majority of matches are 1:1 and will always be picked in each bootstrap.
        # These matches are essentially oversampled
        # Other match ratios eg. 1:2 or 2:1 will be randomly picked for each boot strap - the lines starting group_by(FINALID.neg) and group_by(FINALID.pos) do this.
        # There are not equal amounts of potential matches from both classes before sampling.
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
            # it's not symmetrical. It is one source of our initially high estimates.
            ungroup() %>% mutate(boot = i)
          matchesSelected = bind_rows(matchesSelected,tmp)
        }
        message("..finished")
        
      } else if(resolutionStrategy == "none") {
        # ignore duplicates altogether
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
  return(out)
}

#' Pivots a wide pairlist into a long list suitable for default cox model analysis, and restores any data stripped out when matching pairs.
#'
#' @param matches - the paired data
#' @param coxData - the original unpaired data
#'
#' @return a cox model friendly dataframe
#'
#' @examples
pairedMatchesToCoxData = function(matches, coxData) {
  out = matches %>% p_modify(
    function (d) {
      if (!"boot" %in% colnames(d)) d = d %>% mutate(boot=1)
      return(bind_rows(
        coxData %>% inner_join(d %>% select(boot,FINALID=FINALID.pos), by=c("FINALID")) %>% distinct(),
        coxData %>% inner_join(d %>% select(boot,FINALID=FINALID.neg), by=c("FINALID")) %>% distinct()
      ))
    }
  ) %>% 
    p_group_by(sGene) %>% 
    p_status_summary(
      count = n(),
      boots = n_distinct(boot),
      died = sum(ifelse(diedWithin28days,1,0)),
      .glue = c(
        "{sprintf('%1.0f',count/boots)} patients",
        "of whom {sprintf('%1.0f',died/boots)} died within 28 days",
        "and {sprintf('%1.0f',(count-died)/boots)} survived for 28 days"
      )
    ) %>% p_ungroup(
      count = n(),
      boots = n_distinct(boot),
      .glue = c("{sprintf('%1.0f',count/boots)} matched patients.")
    )
  return(out)
}



#' The main entry point for the analysis and all control options.
#'
#' @param coxData - the loaded data from loadData() function
#' @param sGeneCtThreshold - defaults to the same value as ctThreshold
#' @param ctThreshold - CT values defining positivity CT <= ctThreshold
#' @param bootstraps - Number of replicates if needed - test with a small number
#' @param ageTolerance - delta age <= age tolerance (years)
#' @param specimenDateTolerance  - delta specimen_date <= specimen date tolerance (days)
#' @param modelFormula - a list of cox model formulae to try. Columns allowed are any from the coxData input, should be specified as a named list of formulae e.g. list(`S gene + age`=Surv(time,status) ~ sGene+age),
#' @param max - maximum number of sGene positive x sGene negative you are allowing to be matched before the algorithm aborts. With loose mathcing criteria very large joins can cause crashes.
#' @param includeRaw - include as output the raw unpaired cox Data - for further analysis
#' @param includeMatched  - include as output the paired cox Data - for further analysis 
#' @param matchOn - columns to use to generate pairs (see createPairedMatches())
#' @param resolutionStrategy - strategy to resolve pair duplicates (see uniquifyPairs())
#' @param ratio.neg - filters for excluding cases with many matches to controls (see uniquifyPairs())
#' @param ratio.pos - filters for excluding controls with many matches to cases (see uniquifyPairs())
#'
#' @return - a summary output of all the analysis
#' @export
#'
#' @examples
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
  if (includeMatched) result$pairedMatchedData = matchesSelected
  
  # prepare data for cox PH models
  coxFinal = matchesSelected %>% pairedMatchesToCoxData(coxData = coxData3) %>% 
    p_mutate(sGene = sGene %>% forcats::as_factor() %>% forcats::fct_drop() %>% forcats::fct_relevel("Positive")) # equivocal category is dropped
  
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
  
  ## high level stats about number of deaths in each arm of control study where 
  combinedSummary1 = coxFinal %>% group_by(sGene,boot) %>% summarise(pairs = n()) %>% left_join(
    coxFinal %>% filter(diedWithin28days) %>% group_by(sGene,boot) %>% summarise(deaths = n()),
    by = c("sGene","boot")
  )
  result$table1Data = combinedSummary1

  ## summary of pairwise matching data set
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
  
  rm(matchesSelected, matches, coxFinal, coxData3, posMatched, negMatched)
  gc()
  
  return(result)
}


#' Conduct a survival analysis multiple times over bootstrapped data
#'
#' @param coxFinal - multiple replicates of survival data containing a "boot" column specifying the replicate nymber
#' @param modelName - a description of the model
#' @param modelFormula - the specification of the model
#'
#' @return - a tidy summary table of 1 row per HR from the model
#'
#' @examples
bootstrappedCoxModel = function(coxFinal, modelName, modelFormula) {
  message("Calculating ",modelName)
  bootSurvModels = coxFinal %>% group_by(boot) %>% group_modify(function(d,g,...) {
    survModel = survival::coxph(modelFormula, d)
    out = tidyCoxmodel(survModel) %>% mutate(model=modelName, N=d %>% nrow())
    rm(survModel)
    return(out)
  })
}

## Formatting functions ----

#' Extract the headline values from the cox model result
#' @param survModel - the cox model
#'
#' @return a dataframe of HR and beta values for each component of the model
#'
#' @examples
tidyCoxmodel = function(survModel) {
  # TODO: Dig N out from here
  broom::tidy(survModel) %>% bind_cols(as.data.frame(confint(survModel))) %>% 
    mutate(
      HR = exp(estimate), 
      HR.lower=exp(`2.5 %`), 
      HR.upper=exp(`97.5 %`), 
      HR.pValue = p.value)
}

# Internal function to generate a printable summary table of information from raw cox data.
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

# Internal function to generate a printable summary table of the analysis.
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

# combine cox model results from bootstraps combined Beta is the mean or mean, and mean of std.error, and pValue as the mean of pValue. HR is then recalculated from Betas 
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

# pretty print the hazard table, re-labelling colums as required
# TODO: if we change the cox models we use we need to make sure new predictors are added here, 
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

