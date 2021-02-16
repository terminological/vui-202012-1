## A set of functions to perform tidy manipulations on data and retain information about them for the purposes of documentation
# This is a proof of concept
# requires the following packages
# install.packages("DiagrammeR")
# install.packages("DiagrammeRsvg")

# df = tibble(a=c(1,1,1,2,2,2), b=c(1,2,3,1,2,3), c=c(1,2,3,4,5,6))
# df = df %>% group_by(a)
# df %>% p_clear() %>% p_generate(nrow) %>% provenance()
# df %>% p_clear() %>% prov("test") %>% prov("test2") %>% provenance()
# df %>% p_clear() %>% p_comment(df %>% summarise(message = n())) %>% provenance()
p_comment = function(.data, .messageDf=NULL, type="info", .headline=NULL, terminal = (type=="exclusion")) {
  
  if (identical(.messageDf,NULL) & identical(.headline,NULL)) stop("one of .messageDf or .headline nust be defined")
  env = rlang::caller_env()
  if (!identical(.headline,NULL)) {.headline = glue::glue(.headline,.envir = env)}

  grps = .data %>% groups()
  tmp = .messageDf
  env$total = nrow(.data)
  # convert .message into a data frame if not already
  if (!is.data.frame(tmp)) {
    tmp = tibble(.message=as.character(tmp))
  }
  
  if (!".message" %in% colnames(tmp)) stop(".messageDf dataframes must contain at least the .message column")
  
  if (grps %>% length() > 0) {
    if(!(".strata" %in% colnames(tmp))) {
      # merge groups and values into a single column.
      strata = .data %>% group_data() %>% 
        select(-.rows)
      if (all(grps %in% (tmp %>% colnames()))) {
        # .messageDf was grouped
        tmp = strata %>% full_join(tmp,by=grps %>% sapply(as_label) %>% as.character())
      } else if ((tmp %>% groups() %>% length())==0) {
        # assume .message will either be 1 or exactly the right length
        # TODO: if .messageDf is NULL this fails.
        tmp = strata %>% mutate(.message = tmp$.message)
      }
      tmp = tmp %>% 
        mutate(across(.cols = !!!grps,.fns=~paste0(cur_column(),":",.x), .names=".strata.{.col}")) %>% 
        unite(col = ".strata",starts_with(".strata."),sep="; ")
    }
  } else {
    tmp = tmp %>% mutate(.strata = NA_character_)
  }
  
  tmp = tmp %>% ungroup() %>% select(.strata,.message) 
  tmp = tmp %>% group_by(.strata) %>% summarise(.message = paste0(.message,collapse = "\n"))
  # add in .headline if present
  if (!identical(.headline,NULL)) {tmp = tmp %>% mutate(.message = paste0(.headline,"\n",.message))}
  
  # add the message to the .data@prov attribute which is a graph structure of nodes, edges and a head 
  if (identical(attr(.data,"prov"),NULL)) {
    # as a new message
    nodes=tmp %>% mutate(rank=1,type=type,id=row_number())
    edges=tibble(id=integer(),to=integer(),from=integer(),rel=character(), .strata=character()) # empty edges
    head=nodes %>% select(from = id,.strata) %>% distinct()
    attr(.data,"prov") = list(nodes=nodes,edges=edges,head=head)
  } else {
    current = attr(.data,"prov")
    currentRank = max(current$nodes$rank)
    currentMaxId = max(current$nodes$id)
    nodes=tmp %>% mutate(rank=currentRank+1,type=type,id=row_number()+currentMaxId)
    new = list()
    new$nodes = current$nodes %>% bind_rows(nodes)
    #TODO: if the new node is terminal put in an invisible extra node for the branch point.
    newEdges = nodes %>% 
      select(to = id, rel = type, .strata) %>% 
      inner_join(current$head, by=character(),suffix=c("",".head")) %>% 
      filter(is.na(.strata) | is.na(.strata.head) | .strata == .strata.head)
    new$edges = current$edges %>% bind_rows(newEdges)
    if (!terminal) {
      # this is 
      new$head = nodes %>% select(from = id, .strata)
    } else {
      # this rank is a terminal part of the flowchart (i.e. an exclusion)
      new$head = current$head
    }
    attr(.data,"prov") = new
  }
  return(.data)
}

provenance = function(.data) {
  return(attr(.data,"prov"))
}

p_clear = function(.data) {
  attr(.data,"prov")=NULL
  return(.data)
}

p_generate = function(.data, messageFn, type="info") {
  out = .data %>% group_modify(function(d,g,...) {
    tibble(.message=messageFn(d))
  }) 
  return(p_comment(.data,out,type))
}

# e.g. df %>% p_clear() %>% p_glue(nrow,"has {.nrow} rows") %>% provenance()
# TODO: make this into p_status_summary where messageFn is a ... named list and glue is a list.
p_glue = function(.data, messageFn, .glue = paste0("{",fnName,"}"), type="info") {
  default_env = rlang::caller_env()
  fnName = paste0(".",as.character(substitute(messageFn)))
  out = .data %>% group_modify(function(d,g,...) {
    tibble(!!fnName:=messageFn(d))
  }) 
  messages = glue::glue_data(.x=out, .glue)
  return(p_comment(.data,messages,type))
}

# e.g. iris %>% p_group_by(Species) %>% 
#   p_status_summary(
#     count = n(), 
#     meanLength = mean(Sepal.Length), 
#     .glue = c("{meanLength} average length","{count} items")
#   ) %>% provenance()
# TODO: make this into p_status_summary where messageFn is a ... named list and glue is a list.
p_status_summary = function(.data, ..., .glue, type="info") {
  default_env = rlang::caller_env()
  grps = .data %>% groups()
  out = .data %>% summarise(...) %>% group_by(!!!grps)
  out = out %>% group_by_all() %>% group_modify(function(d,g,..) {
    tibble(.message = sapply(.glue, function(gs) {glue::glue_data(.x=g,gs)}))
  })
  return(p_comment(.data,out,type))
}


p_copy = function(.data, from) {
  attr(.data,"prov") = attr(from,"prov")
  return(.data)
}

# df %>% p_clear() %>% p_comment("test") %>% p_exclude(c%%2==0 ~ `excluding {count} even items`) %>% provenance()
p_exclude = function(.data, ..., na.rm=FALSE, .headline="Exclude:") {
  messages=NULL
  grps = .data %>% groups()
  default_env = rlang::caller_env()
  out = .data %>% mutate(.retain = TRUE)
  filters = rlang::list2(...)
  default_env$total = nrow(.data)
  for(filter in filters) {
    glueSpec = rlang::f_rhs(filter)
    filt = rlang::f_lhs(filter)
    out = out %>% group_modify(function(d,g,...) {
      d %>% 
        mutate(.excl = rlang::eval_tidy(filt,data = d, env=default_env)) %>% 
        mutate(.retain = .retain & !ifelse(is.na(.excl),na.rm,.excl))
    })
    tmp = out %>% 
      filter(.excl) %>% 
      summarise(count = n()) %>% 
      ungroup() %>%
      tidyr::complete(!!!grps,fill=list(count=0))
    tmp$.message = rlang::eval_tidy(glue::glue_data(tmp,glueSpec,.envir = default_env), data=tmp, env = default_env)
    messages = messages %>% bind_rows(tmp)
  }
  out = out %>% filter(.retain) %>% select(-.retain,-.excl) %>% p_copy(.data) %>% p_comment(messages, .headline = .headline, type="exclusion")
  return(out)
}

# x = function() {filterCond = 1; df %>% p_clear() %>% p_status(c%%2==filterCond ~ `consisting of {count} even items+{filterCond}`,c%%2!=0 ~ `and {count} odd items+{filterCond}`) %>% provenance()}
# x()
# df %>% p_clear() %>% p_status(c%%2==0 ~ `consisting of {count} even items`,c%%2!=0 ~ `and {count} odd items`) %>% provenance()
p_status = function(.data, ..., na.count=FALSE, .headline=NULL) {
  messages = NULL
  grps = .data %>% groups()
  default_env = rlang::caller_env()
  out = .data
  default_env$total = nrow(.data)
  filters = rlang::list2(...)
  if(length(filters)==0) filters = list(TRUE ~ "{count} items")
  for(filter in filters) {
    glueSpec = rlang::f_rhs(filter)
    filt = rlang::f_lhs(filter)
    tmp = out %>% group_modify(function(d,g,...) {
      z = d %>% 
        mutate(.excl = rlang::eval_tidy(filt, data=d, env=default_env)) %>% 
        filter(.excl) %>% 
        summarise(count = n()) 
      return(z)
    }) 
    tmp = tmp %>% ungroup() %>% 
      tidyr::complete(!!!grps,fill=list(count=0)) 
    
    tmp$.message = rlang::eval_tidy(glue::glue_data(tmp,glueSpec,.envir = default_env), data=tmp, env = default_env)
    
    
    messages = messages %>% bind_rows(tmp)
    # filtBool = rlang::eval_tidy(filt, data=.data, env = default_env)
    # filtBool = ifelse(is.na(filtBool),!na.count,filtBool)
    # tmp = .data %>% ungroup() %>% filter(filtBool) %>% group_by(!!!grps) %>% summarise(count=n())
    # tmp$.message = rlang::eval_tidy(glue::glue_data(tmp,glueSpec,.envir = default_env), data=tmp, env = default_env)
    # messages = messages %>% bind_rows(tmp)
  }
  out = .data %>% p_comment(.messageDf = messages,.headline = .headline, type="summary")
  return(out)
}

# df %>% p_clear() %>% p_status(c%%2==0 ~ "consisting of {count} even items",c%%2!=0 ~ "and {count} odd items") %>% p_ungroup() %>% provenance()
p_ungroup = function(.data, ..., .glue = "totalling {count} items") {
  filters = enexprs(...)
  if(length(filters)==0) filters = list(count=expr(n()))
  out = .data %>% ungroup() %>% p_copy(.data) 
  out = p_status_summary(out, !!!filters, .glue=.glue)
  return(out)
}

# df %>% p_clear() %>% p_comment("test") %>% p_summarise(.message = "combining as mean",mean.a = mean(a)) %>% provenance()
p_summarise = function(.data, ..., .message = NULL) {
  out = .data %>% summarise(...)
  out = out %>% p_copy(.data)
  if (!identical(.message,NULL)) out = out %>% p_comment(.message, "summarise")
  return(out)
}

# df %>% p_clear() %>% p_comment("test") %>% p_mutate(.message = "add date",d = 123) %>% p_comment("test 2") %>% provenance()
p_mutate = function(.data, ..., .message = NULL) {
  out = .data %>% mutate(...)
  out = out %>% p_copy(.data) 
  if (!identical(.message,NULL)) out = out %>% p_comment(.message, "modify")
  return(out)
}

# df %>% ungroup() %>% p_clear() %>% p_status() %>% p_group_by(a) %>% p_comment("test") %>% p_group_by(b) %>% p_comment("test 2") %>% p_status() %>% p_flowchart()
# TODO: multiple columns
p_group_by = function(.data, col, .message = NULL) {
  col = ensym(col)
  total = nrow(.data)
  if(identical(.message,NULL)) .message = glue::glue("stratifying by {col}",col=col)
  tmp = .data %>% ungroup() %>% p_copy(.data)
  tmp = tmp %>% p_comment(.message, "stratify") %>% group_by(!!col)
  return(tmp)
}

p_select = function(.data, ...) {
  out = .data %>% select(...) %>% p_copy(.data)
  return(out)
}

# df %>% p_clear() %>% p_modify(function(d) { d %>% filter(c==2) }, .message="was {count.in}, now {count.out}") %>% provenance()
# TODO: FAILS: df %>% p_clear() %>% p_modify(function(d) { d %>% filter(c==2) }, .headline="was {nrow(df)}") %>% provenance()
# df %>% p_clear() %>% p_modify(function(d) { d %>% filter(c==2) }) %>% provenance() # NULL
p_modify = function(.data, dplyrFn, .message=NULL, .headline=NULL) {
  env = rlang::caller_env()
  if (!identical(.headline,NULL)) {.headline = glue::glue(.headline,.envir = env)}
  out = dplyrFn(.data)
  out = out %>% p_copy(.data)
  grps = .data %>% groups()
  tmp = NULL
  if (!identical(.message,NULL)) {
    countIn = .data %>% summarise(count.in=n()) %>% ungroup() %>% tidyr::complete(!!!grps,fill=list(count.in=0))
    countOut = out %>% summarise(count.out = n()) %>% ungroup() %>% tidyr::complete(!!!grps,fill=list(count.out=0))
    tmp = countIn %>% left_join(countOut, by=(grps %>% sapply(as_label) %>% as.character())) %>% mutate(count.out = ifelse(is.na(count.out),0,count.out)) 
    tmp = tmp %>% mutate(.message = glue::glue(.message))
  }
  if (!identical(.message,NULL) | !identical(.headline,NULL)) out = out %>% p_comment(.messageDf = tmp,.headline = .headline, type = "modify")
  return(out)
}

# TODO: https://cran.r-project.org/web/packages/Gmisc/vignettes/Grid-based_flowcharts.html
# https://www.r-bloggers.com/2018/05/flow-charts-in-r/


# df %>% ungroup() %>% p_clear() %>% p_status() %>% p_group_by(a) %>% p_comment("test") %>% p_group_by(b) %>% p_comment("test 2") %>% p_status() %>% p_flowchart()
p_flowchart = function(.data, filename = NULL, fill="lightgrey", fontsize="8", colour="black", ...) {
  
  nodesDf = (.data %>% provenance())$nodes %>% rename(label = .message) %>% mutate(
    label = ifelse(is.na(.strata),label,paste0(.strata,"\n",label)),
    fillcolor= case_when(type=="summary"~"grey90",type=="exclusion"~"grey80",TRUE~"white")
  )
  edgesDf = (.data %>% provenance())$edges %>% mutate(
    headport=ifelse(rel=="exclusion","w","n"),
    weight=ifelse(rel=="exclusion","1","10"),
    tailport="s",
    colour="black",
    id = row_number()
  )
  
  outNode = nodesDf %>% group_by(desc(id)) %>% mutate(
    #html = label %>% htmltools::htmlEscape() %>% stringr::str_replace_all("\\n","</TD></TR><TR><TD ALIGN='LEFT'>")
  #) %>% mutate(nodeSpec = glue::glue("'{id}' [label=<<TABLE><TR><TD ALIGN='LEFT'>{html}</TD></TR></TABLE>>,group='{.strata}',fillcolor='{fillcolor}'];")) %>%
    html = label %>% htmltools::htmlEscape() %>% stringr::str_replace_all("\\n","<BR ALIGN='LEFT'/>")
  ) %>% mutate(nodeSpec = glue::glue("'{id}' [label=<{html}<BR ALIGN='LEFT'/>>,group='{.strata}',fillcolor='{fillcolor}'];")) %>%
    group_by(rank) %>% 
    summarise(rankSpec = paste0(nodeSpec,collapse = "\n")) %>% 
    mutate(rankSpec = paste(
      "{ rank='same';",
      rankSpec,
      "}",sep="\n")
    ) %>% arrange(desc(rank)) %>% pull(rankSpec) %>% paste0(collapse="\n")
  
  outEdge = edgesDf %>% arrange(desc(id)) %>% mutate(edgeSpec = glue::glue("'{from}' -> '{to}' [tailport='{tailport}',weight='{weight}']")) %>%
    pull(edgeSpec) %>% paste0(collapse="\n")
  
  outGraph = paste(
    "digraph {
     graph [layout = 'dot',
        splines='ortho',
        rankdir = 'TB',
        outputorder = 'edgesfirst',
        bgcolor = 'white',
        ranksep = '0.25',
        nodesep = '0.2',
        newrank='true']

    node [fontname = 'Helvetica',
        fontsize = '8',
        shape='box',
        fixedsize = 'false',
        margin = '0.1,0.1',
        width = '0',
        height = '0',
        style = 'filled',
        color = 'black',
        fontcolor = 'black',
        labeljust='l']

    edge [fontname = 'Helvetica',
        fontsize = '8',
        len = '0.5',
        color = 'black',
        arrowsize = '0.5']
    ",outNode,"\n",outEdge,sep="\n","}")
  
  #cat(outGraph)
  
  if (!identical(filename,NULL)) {
    filename = filename %>% stringr::str_remove("\\..*$")
    
    rsvg::rsvg_png(charToRaw(DiagrammeRsvg::export_svg(DiagrammeR::grViz(outGraph))),file = paste0(filename,".png"),...)
    rsvg::rsvg_pdf(charToRaw(DiagrammeRsvg::export_svg(DiagrammeR::grViz(outGraph))),file = paste0(filename,".pdf"),...)
    
    unlink(paste0(filename,".dot"))
    outGraph %>% writeChar(paste0(filename,".dot"))
    unlink(paste0(filename,".svg"))
    DiagrammeRsvg::export_svg(DiagrammeR::grViz(outGraph)) %>% writeChar(paste0(filename,".svg"))
    
  } else {
    filename = tempfile()
    rsvg::rsvg_png(charToRaw(DiagrammeRsvg::export_svg(DiagrammeR::grViz(outGraph))),file = paste0(filename,".png"),...)
    rsvg::rsvg_pdf(charToRaw(DiagrammeRsvg::export_svg(DiagrammeR::grViz(outGraph))),file = paste0(filename,".pdf"),...)
  }
  
  if (isTRUE(getOption("knitr.in.progress"))) {
    fmt <- rmarkdown::default_output_format(knitr::current_input())$name
    return(knitr::include_graphics(path = normalizePath(paste0(filename,".png"),mustWork = FALSE),auto_pdf = TRUE))
  } else {
    return(DiagrammeR::grViz(outGraph))
  }
  
  #tmp = DiagrammeR::generate_dot(graph)
  #tmp2 = tmp %>% stringr::str_replace("graph \\[","graph [splines='ortho',")
  #tmp2 %>% DiagrammeR::grViz()
  #DiagrammeR::export_graph(graph,file_name = "graph.svg")
}

demo = function() {
  cutOff = 3
  iris %>% ungroup() %>% p_clear() %>%
    p_status() %>%
    p_group_by(Species) %>%
    p_status(
      Sepal.Width<cutOff ~ "consisting of {count} short sepal <{cutOff}",
      Sepal.Width>=cutOff ~ "and {count} long sepal >={cutOff}"
    )  %>%
    p_exclude(
      Petal.Width<0.3 ~ "excluding {count} with narrow petals",
      Petal.Width == max(Petal.Width) ~ "and {count} outlier"
    ) %>%
    p_comment("test message") %>%
    p_status(TRUE ~ "{count} of type {Species}") %>%
    p_ungroup() %>%
    p_status(TRUE ~ "{count} together with cutOff {cutOff}") %>%
    p_flowchart()
}