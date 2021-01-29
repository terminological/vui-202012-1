## A set of functions to perform tidy maniulations on data and retain information about them for the purposes of documentation
# This is a proof of concept
# requires the following packages
# install.packages("DiagrammeR")
# install.packages("DiagrammeRsvg")

# df = tibble(a=c(1,1,1,2,2,2), b=c(1,2,3,1,2,3), c=c(1,2,3,4,5,6))
# df = df %>% group_by(a)
# df %>% p_clear() %>% p_generate(nrow) %>% provenance()
# df %>% p_clear() %>% prov("test") %>% prov("test2") %>% provenance()
# df %>% p_clear() %>% p_comment(df %>% summarise(message = n())) %>% provenance()
p_comment = function(.data, message, type="info") {
  if (is.data.frame(message)) {
    
    if (!".message" %in% colnames(message)) stop("dataframes as messages must contain at least the .message column")
    tmp = message %>% ungroup() %>% select(c(!!!.data %>% groups(),.message)) 
    tmp = tmp %>% group_by(!!!.data %>% groups()) %>% summarise(.message = paste0(.message,collapse = "\n")) %>% ungroup()
    tmpNames = colnames(tmp %>% select(-.message))
    if (length(tmpNames) > 0) {
      tmpPrefix = paste0(tmpNames,":",collapse = "")
      message = tmp %>% pivot_wider(names_from = -.message, values_from = .message, names_prefix = tmpPrefix) %>% as.list()
    } else {
      message = tmp %>% pull(.message) %>% as.list()
    }
    message = list(message=message,type=type)
    
  } else if (.data %>% groups() %>% length() == 0 & !is.list(message)) {
    
    message=list(message=message,type=type)
  } else if (.data %>% groups() %>% length() == 1) {
    
    tmp = .data %>% group_data() %>% select(-.rows) %>% mutate(message = message)
    message = tmp %>% pivot_wider(names_from = -message, values_from = message, names_prefix = paste0(.data %>% groups(),":",collapse = "")) %>% as.list()
    message = list(message=message,type=type)
  } else {
    stop("multiply grouped dataframes are not supported")
  }
  # add the message to the .data
  if (identical(attr(.data,"prov"),NULL)) {
    # as a new message
    attr(.data,"prov") = list(message)
  } else {
    # to existing messages
    current = attr(.data,"prov")
    attr(.data,"prov") = c(current,list(message))
  }
  return(.data)
}

provenance = function(.data) return(attr(.data,"prov"))

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
p_glue = function(.data, messageFn, glue = paste0("{",fnName,"}"), type="info") {
  fnName = paste0(".",as.character(substitute(messageFn)))
  out = .data %>% group_modify(function(d,g,...) {
    tibble(!!fnName:=messageFn(d))
  }) 
  messages = glue::glue_data(.x=out, glue)
  return(p_comment(.data,messages,type))
}

p_copy = function(.data, from) {
  attr(.data,"prov") = attr(from,"prov")
  return(.data)
}

# df %>% p_clear() %>% p_comment("test") %>% p_exclude(c%%2==0 ~ `excluding {count} even items`) %>% provenance()
p_exclude = function(.data, ...) {
  messages = NULL
  filterList = NULL
  filters = rlang::list2(...)
  for(filter in filters) {
    glueSpec = as.character(attr(terms(filter),"variables")[3])
    filt = as.character(attr(terms(filter),"variables")[2])
    message = .data %>% filter(eval(rlang::parse_expr(filt))) %>% summarise(count=n()) %>% mutate(.message = glue::glue(glueSpec))
    messages = messages %>% bind_rows(message)
    filterList = c(filterList,filt)
  }
  filterList2 = unname(sapply(filterList, function(x) paste0("!(",x,")")))
  out = .data %>% filter(!!!rlang::parse_exprs(filterList2)) %>% p_copy(.data) %>% p_comment(messages, type="exclusion")
  return(out)
}


# df %>% p_clear() %>% p_status(c%%2==0 ~ `consisting of {count} even items`,c%%2!=0 ~ `and {count} odd items`) %>% provenance()
p_status = function(.data, ...) {
  messages = NULL
  filters = rlang::list2(...)
  if(length(filters)==0) filters = list(TRUE ~ `{count} items`)
  for(filter in filters) {
    glueSpec = as.character(attr(terms(filter),"variables")[3])
    filt = as.character(attr(terms(filter),"variables")[2])
    message = .data %>% filter(eval(rlang::parse_expr(filt))) %>% summarise(count=n()) %>% mutate(.message = glue::glue(glueSpec))
    messages = messages %>% bind_rows(message)
  }
  out = .data %>% p_comment(messages, type="summary")
  return(out)
}

# df %>% p_clear() %>% p_status(c%%2==0 ~ `consisting of {count} even items`,c%%2!=0 ~ `and {count} odd items`) %>% p_ungroup() %>% provenance()
p_ungroup = function(.data, ...) {
  messages = NULL
  filters = rlang::list2(...)
  if(length(filters)==0) filters = list(TRUE ~ `combining {count} items`)
  out = .data %>% ungroup()
  for(filter in filters) {
    glueSpec = as.character(attr(terms(filter),"variables")[3])
    filt = as.character(attr(terms(filter),"variables")[2])
    message = out %>% filter(eval(rlang::parse_expr(filt))) %>% summarise(count=n()) %>% mutate(.message = glue::glue(glueSpec))
    messages = messages %>% bind_rows(message)
  }
  out = out %>% p_copy(.data) %>% p_comment(messages, type="summary")
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
  dots = rlang::list2(...)
  out = .data %>% mutate(!!!dots)
  out = out %>% p_copy(.data) 
  if (!identical(.message,NULL)) out = out %>% p_comment(.message, "modify")
  return(out)
}

# df %>% ungroup() %>% p_clear() %>% p_status() %>% p_group_by(a) %>% p_comment("test") %>% p_group_by(b) %>% p_comment("test 2") %>% p_status() %>% p_flowchart()
p_group_by = function(.data, col, .message = NULL) {
  col = ensym(col)
  if(identical(.message,NULL)) .message = glue::glue("stratifying by {col}",col=col)
  tmp = .data %>% ungroup() %>% p_copy(.data)
  tmp = tmp %>% p_comment(.message, "stratify") %>% group_by(!!col)
  return(tmp)
}

p_modify = function(.data, dplyrFn, .message=NULL) {
  out = dplyrFn(.data)
  out = out %>% p_copy(.data)
  if (!identical(.message,NULL)) out = out %>% p_comment(.message, "modify")
}

# TODO: https://cran.r-project.org/web/packages/Gmisc/vignettes/Grid-based_flowcharts.html
# https://www.r-bloggers.com/2018/05/flow-charts-in-r/

# df %>% ungroup() %>% p_clear() %>% p_status() %>% p_group_by(a) %>% p_comment("test") %>% p_group_by(b) %>% p_comment("test 2") %>% p_status() %>% p_flowchart()
# p_flowchart = function(.data, fill="lightgrey", fontsize="8", colour="black", ...) {
#   bgp = grid::gpar(fill=fill, col=colour,...)
#   tgp = grid::gpar(fontsize=fontsize,...)
#   grobList = NULL
#   lists = .data %>% provenance()
#   prevLevel = NULL
#   allLevels = NULL
#   for (item in lists) {
#     type = item$type
#     thisLevel = NULL
#     for (i in 1:length(item$message)) {
#       m = item$message[[i]]
#       n = names(item$message)[[i]]
#       if (!identical(n,NULL)) m = paste0(n,"\n",m)
#       (tmp = Gmisc::boxGrob(label = m,txt_gp = tgp,box_gp = bgp))
#       if (length(prevLevel) == length(item$message)) {
#         linkTo = prevLevel[[i]]
#         tmp2 = Gmisc::connectGrob(linkTo,tmp,"N")
#         grobList = c(grobList,list(tmp2))
#       } else if (length(prevLevel) == 1) {
#         linkTo = prevLevel[[1]]
#         tmp2 = Gmisc::connectGrob(linkTo,tmp,"N")
#         grobList = c(grobList,list(tmp2))
#       } else if (length(item$message) == 1 & length(prevLevel) > 1) {
#         for (j in 1:length(prevLevel)) {
#           linkTo = prevLevel[[j]]
#           tmp2 = Gmisc::connectGrob(linkTo,tmp,"N")
#           grobList = c(grobList,list(tmp2))
#         }
#       }
#       thisLevel = c(thisLevel,list(tmp))
#     }
#     if(length(thisLevel) >1) {
#       tmp3 = Gmisc::spreadHorizontal(thisLevel)
#     } else {
#       tmp3=thisLevel[[1]]
#     }
#     allLevels = c(allLevels,list(tmp3))
#     prevLevel = thisLevel
#   }
#   
#   if(length(allLevels) >1) {
#     tmp4 = Gmisc::spreadVertical(allLevels)
#   } else {
#     tmp4=allLevels[[1]]
#   }
#   
#   grobList = c(grobList,list(tmp4))
#   grid::grid.newpage()
#   browser()
#   sapply(grobList,plot)
# }


# df %>% ungroup() %>% p_clear() %>% p_status() %>% p_group_by(a) %>% p_comment("test") %>% p_group_by(b) %>% p_comment("test 2") %>% p_status() %>% p_flowchart()
p_flowchart = function(.data, filename = NULL, fill="lightgrey", fontsize="8", colour="black", ...) {
  
  nodesDf = tibble()
  edgesDf = tibble()
  .createNode = function(label, type) {
    id = nrow(nodesDf)+1
    nodesDf <<- nodesDf %>% bind_rows(tibble(id = id, label=label, type=type))
    return(id)
  }
  .createEdge = function(from,to,rel=NA) {
    id = nrow(edgesDf)+1
    edgesDf <<- edgesDf %>% bind_rows(tibble(id=id,to=to,from=from,rel=rel))
    return(id)
  }
  
  lists = .data %>% provenance()
  prevLevel = NULL
  for (item in lists) {
    type = item$type
    thisLevel = NULL
    for (i in 1:length(item$message)) {
      m = item$message[[i]]
      n = names(item$message)[[i]]
      if (!identical(n,NULL)) m = paste0(n,"\n",m)
      tmp = .createNode(m,type)
      if (length(prevLevel) == length(item$message)) {
        linkTo = prevLevel[[i]]
        .createEdge(from = linkTo, to = tmp, rel=type)
      } else if (length(prevLevel) == 1) {
        linkTo = prevLevel[[1]]
        .createEdge(from = linkTo, to = tmp, rel=type)
      } else if (length(item$message) == 1 & length(prevLevel) > 1) {
        for (j in 1:length(prevLevel)) {
          linkTo = prevLevel[[j]]
          .createEdge(from = linkTo, to = tmp, rel=type)
        }
      }
      thisLevel = c(thisLevel,tmp)
    }
    if (type != "exclusion") prevLevel = thisLevel
  }
  browser()
  graph = DiagrammeR::create_graph(
    #graph [splines=ortho]
    nodes_df = nodesDf %>% mutate(
      shape="rectangle",
      color="black",
      fontcolor="black",
      fixedsize=FALSE, 
      fillcolor= case_when(type=="summary"~"grey90",type=="exclusion"~"grey80",TRUE~"white")
    ),
    edges_df = edgesDf %>% mutate(
      headport=ifelse(rel=="exclusion","w","n"),
      tailport="s",
      colour="black"
    ), attr_theme = "tb")
  if (!identical(filename,NULL)) {
    DiagrammeR::export_graph(graph,file_name = filename,...)
  } else {
    return(graph %>% DiagrammeR::render_graph())
  }
  #tmp = DiagrammeR::generate_dot(graph)
  #tmp2 = tmp %>% stringr::str_replace("graph \\[","graph [splines='ortho',")
  #tmp2 %>% DiagrammeR::grViz()
  #DiagrammeR::export_graph(graph,file_name = "graph.svg")
}


# df %>% ungroup() %>% p_clear() %>%
#   p_status() %>%
#   p_group_by(a) %>%
#   p_status(
#     c%%2==0 ~ `consisting of {count} even items`,
#     c%%2!=0 ~ `and {count} odd items`
#   )  %>%
#   p_exclude(
#     c%%2==0 ~ `excluding {count} even items`
#   ) %>%
#   p_comment("test") %>%
#   p_group_by(b) %>%
#   p_comment("test 2") %>%
#   p_status() %>%
#   p_ungroup() %>%
#   p_flowchart()
