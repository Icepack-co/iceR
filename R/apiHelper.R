library(RProtoBuf)
library(dplyr)
library(purrr)
library(httr)
library(sf)
library(ggplot2)
library(leaflet)
library(magrittr)

#' @export
#' @rawNamespace useDynLib(iceR)
apiHelper <- setRefClass("apiHelper",
           fields = list(
                        activeProblems = 'list',
                        modelType = "character",
                        apiToken = "character",
                        endpoint = "character",
                        route = "character"),
           methods = list(
             initialize = function(configFile = 'config.json', ...) {
               if(!file.exists(configFile)){
                 writeLines(con = configFile,
                            text = jsonlite::toJSON(list(apiToken = "", endpoint = 'https://api.icepack.ai/'),
                                                    auto_unbox = T))
                 stop(paste0("file not found: '",configFile,"'. Creating a blank file. Please provide a valid apiToken in file."))
               }
               x <- jsonlite::fromJSON(configFile)
               if(!('apiToken' %in% names(x))){
                 stop(paste0("'apiToken' definition not found in config file: ", configFile))
               }
               if(!('endpoint' %in% names(x))){
                 stop(paste0("'endpoint' definition not found in config file: ", configFile))
               }
               if(x$apiToken == ''){
                 stop(paste0("'apiToken' definition in ", configFile," cannot be blank"))
               }
               if(x$endpoint == ''){
                 stop(paste0("'endpoint' definition in ", configFile," cannot be blank"))
               } else {
                 # check if the endpoint needs to be padded.
                 if(substr(x$endpoint, nchar(x$endpoint), nchar(x$endpoint))[[1]] != "/"){
                   x$endpoint <- paste0(x$endpoint, "/")
                 }
               }
               readProtoFiles(system.file("data", 'problem.proto', package = "iceR", mustWork = TRUE))
               .self$apiToken = x$apiToken
               .self$endpoint = x$endpoint
               .self$initFields(...)
               models<- c("tsp-mcvfz472gty6",
                          'tsptw-kcxbievqo879',
                          'cvrp-jkfdoctmp51n',
                          'cvrptw-acyas3nzweqb',
                          'ivr7-kt461v8eoaif',
                          'ivr8-yni1c9k2swof',
                          'nvd-hap0j2y4zlm1',
                          'ndd-cmibu6krtqja',
                          'ns3-tbfvuwtge2iq',
                          'isr-z4foi53qznrv',
                          'matrix-vyv95n7wchpl',
                          'ivrdata-o43e0dvs78zq')

               if(length(modelType) == 0){
                 stop(paste0("No model type provided. Should be one of the following:\n", paste0(models, collapse =',')))
               }
               if(!(modelType %in% models)){
                 stop(paste0("Invalid model type provided. Should be one of the following:\n", paste0(models, collapse =',')))
               } else {
                 readProtoFiles(system.file("data", paste0(modelType, '.proto'), package = "iceR", mustWork = TRUE))
               }
               # grab the route associated with the model type.
               modelroutes <- rep(x = 'vehicle-router/solve/', times = length(models))
                 modelroutes[models == 'ns3-tbfvuwtge2iq'] <- 'network-sourcing/solve/'
               modelroutes[models == 'matrix-vyv95n7wchpl'] <- 'matrix/'
               modelroutes[models == 'ivrdata-o43e0dvs78zq'] <- 'vehicle-router/data/'
               .self$route <- modelroutes[match(modelType, models)]
             }
           ))

#' @export
#' @rawNamespace useDynLib(iceR)
postSolveRequest <- function(apiHelper, solveRequest){
  if(solveRequest$isInitialized()) {
    p <- new (problem.ProblemEnvelope)
    p$type <- apiHelper$modelType
    p$subType <- 0
    p$content <- solveRequest$serialize(NULL)
    body<- p$serialize(NULL)
    postResp <- POST(url = paste0(apiHelper$endpoint, apiHelper$route),
                     body = body,
                     add_headers(`Content-Type`="application/protobuf",
                                 Authorization=paste0("Apitoken ",apiHelper$apiToken)))
    if(postResp$status_code == 200) {
      requestID <- content(postResp)$requestid
      apiHelper$activeProblems[[length(apiHelper$activeProblems) + 1]] <- requestID
      return(requestID)
    } else {
      if(postResp$status_code == 405 && length(content(postResp)) == 0){
        stop('Unexpected http error code: ', postResp$status_code, "\nDid you specify the end-point correctly?\n")
      }else{
        stop('Unexpected http error code: ', postResp$status_code, "\n", content(postResp))
      }
    }
  } else {
    stop("Problem payload not initialised. Have you set all the fields on the problem object?")
  }
}

#' @export
#' @rawNamespace useDynLib(iceR)
getSolutionInstance <- function(apiHelper, solveRequest, paretoResponse, solutionIndex){
  if(solutionIndex > paretoResponse$frontier %>% length) {
    stop("A solution was selected outside of the defined frontier")
  }
  # this method was previously

  nsr <- solveRequest
  if('visitSequence' %in% names(nsr$model)) {
    nsr$model$visitSequence <- paretoResponse$frontier[[solutionIndex]]$visitSequence
  } else { # then it's ndd not nvd
    nsr$model$taskSequence <- paretoResponse$frontier[[solutionIndex]]$taskSequence
  }
  nsr$solveType <- 1                                      # pop this in evaluate mode.
  requestID <- postSolveRequest(apiHelper, nsr)           # submit the model to the api
  resp <- getResponse(apiHelper, requestID)               # retrieve the model response
  return(resp)
}

#' @export
#' @rawNamespace useDynLib(iceR)
getResponse <- function(apiHelper, requestID, getSolverTime = FALSE){
	solResLogs <- c()
  repeat {
    getResp <- GET(url = paste0(apiHelper$endpoint,apiHelper$route, requestID),
                   add_headers(`Content-Type` = "application/protobuf",
                               Authorization = paste0("Apitoken ", apiHelper$apiToken)))
    if (getResp$status_code != 200) {
      break;
    }
    resp_payload <- content(getResp)

    result <- read(descriptor = problem.ProblemEnvelope, resp_payload)

    # unmarshal solver response object
    solRes <- read(problem.SolverResponse, result$content)

    solResStr <- map(solRes$logs, function(i) { i$toString() }) %>% unlist 
	solResLogs <- append(solResLogs, solResStr)
    cat(solResStr)

    if (solRes$state != 1 &&
        solRes$logs[[length(solRes$logs)]]$type != 2) {
      Sys.sleep(0.05) # exponential back-off is implemented on the api, this isn't required.
    } else {
      apiHelper$activeProblems <- apiHelper$activeProblems[apiHelper$activeProblems != requestID]
      if(length(solRes$solution) == 0){
        return(NULL)
      } else {
        solutionResponse <- (switch(apiHelper$modelType,
                       "tsp-mcvfz472gty6"={read(TSP.SolutionResponse,solRes$solution)},
                       "cvrp-jkfdoctmp51n" ={read(CVRP.SolutionResponse,solRes$solution)},
                       'tsptw-kcxbievqo879'={read(TSPTW.SolutionResponse, solRes$solution)},
                       "cvrptw-acyas3nzweqb" ={read(CVRPTW.SolutionResponse,solRes$solution)},
                       'ivr8-yni1c9k2swof' = {read(IVR8.SolutionResponse, solRes$solution)},
                       'ivr7-kt461v8eoaif' = {read(IVR7.SolutionResponse, solRes$solution)},
                       'nvd-hap0j2y4zlm1' = {read(NVD.SolutionResponse, solRes$solution)},
                       'ndd-cmibu6krtqja' = {read(NDD.SolutionResponse, solRes$solution)},
                       'isr-z4foi53qznrv' = {read(ISR.SolutionResponse, solRes$solution)},
                       'ns3-tbfvuwtge2iq' = {read(NS3.SolutionResponse, solRes$solution)},
                       'matrix-vyv95n7wchpl' = {read(Matrix.MatrixResponse, solRes$solution)},
        )
        )
		if(getSolverTime){
			log_line <- solResLogs[grepl("Solved in: ", solResLogs)][1]
		    duration <- sub("s.*", "", sub(".*Solved in: ", "", log_line)) %>% as.numeric()
			return(list(resp=solutionResponse, solverTime=duration))
		} else {
			return(solutionResponse)
		}
      }
      break;
    }
  }
  return(solRes)
}

#' @export
#' @rawNamespace useDynLib(iceR)
plotResponseLeaflet <- function(solutionResponse, solveRequest){
  l<- tabulate(solutionResponse, solveRequest)
  m <- leaflet() %>%
    addTiles()
  if('id' %in% names(l$nodes)){
    m <- m %>%
      addCircleMarkers(data = l$nodes, lng = ~x, lat = ~y, color = 'black', radius = 5,  popup = ~id)
  }else{
    # then it's a vehicle-centric plot
    if(all(c('vehicleId', 'sequence', 'taskId', 'jobId') %in% names(l$nodes))){
      l$nodes <- l$nodes %>%
        mutate(label = paste0("Vehicle: ",
                              vehicleId,'<br>Stop: ',
                              sequence,"<br>Task: ",
                              taskId, "<br>Job: ", jobId))
    }else{
      l$nodes <- l$nodes %>%
        mutate(label = paste0("Vehicle: ",vehicleId,'<br>Stop: ',sequence))
    }
    m <- m %>%
      addCircleMarkers(data = l$nodes, lng = ~x, lat = ~y, color = 'black', radius = 5,  popup = ~label)
  }
  if(solveRequest@type == 'NS3.SolveRequest'){
    m <- m %>% addPolylines(data = st_as_sf(l$routes %>% select(geometry)))
  }else{
    if(!is.null(l$edges$geometry)){
      if('vehicleId' %in% names(l$edges)){
        uv <- unique(l$edges$vehicleId)
        cols <- colorRampPalette(c("blue", "red"))(length(uv))
        for(i in 1:length(uv)){
          m <- m %>% addPolylines(data = l$edges %>%
                                    filter(vehicleId == uv[i]) %>% select(geometry), color = ~cols[i])
        }
      }else{
        m <- m %>% addPolylines(data = l$edges$geometry)
      }
    }else{
      cat("Seems like there are no geometries defined for this model. Plotting points only\n")
    }
  }
  return(m)
}

#' @export
#' @rawNamespace useDynLib(iceR)
plotResponse <- function(solutionResponse, solveRequest,
                         addNodes = TRUE, addNodeLabels = FALSE,
                         addGeometry = TRUE, addSegments = TRUE){
  l<- tabulate(solutionResponse, solveRequest)
  if(length(names(l)) == 1 & 'frontier' %in% names(l)){
    objNames = names(l$frontier)[2:ncol(l$frontier)]
    p <- ggplot() +
      theme_bw() +
      geom_point(data = l$frontier, aes_string(x = objNames[1], y = objNames[2])) +
      scale_x_continuous(labels = comma) +
      scale_y_continuous(labels = comma)
    return (p)
  }else{
    p <- ggplot() + theme_bw() + xlab('') + ylab('')
    if(addSegments){
      if(solveRequest@type == 'NS3.SolveRequest'){
        ss<- l$assignments %>% left_join(l$nodes %>%
                                        select(id, FX = x, FY = y),
                                        by = c('source' = "id")) %>%
                                left_join(l$nodes %>%
                                          select(id, TX = x, TY = y),
                                          by = c('destination' = "id"))
        p <- p + geom_segment(data = ss, aes(x = FX, xend = TX, y = FY, yend = TY), col = 'grey')
        if(is.null(l$routes$geometry) || !addGeometry){
          p <- p + coord_equal()
        }
      }else{
        p <- p + geom_segment(data = l$edges, aes(x = fx, xend = tx, y = fy, yend = ty), col = 'grey')
        if(is.null(l$edges$geometry) || !addGeometry){
          p <- p + coord_equal()
        }
      }
    }
    if(addGeometry){
      if(solveRequest@type == 'NS3.SolveRequest'){
        if(!is.null(l$routes$geometry)){
          p <- p + geom_sf(data = l$routes, aes(geometry = geometry))
        }
      }else{
        if(!is.null(l$edges$geometry)){
          p <- p + geom_sf(data = l$edges, aes(geometry = geometry))
        }
      }
    }
    if(addNodes){
      p <- p + geom_point(data = l$nodes, aes(x = x, y = y))
    }
    if(addNodeLabels){
      if(solveRequest@type == 'NS3.SolveRequest'){
        p <- p + geom_label(data = l$nodes, aes(x = x, y = y, label = id))
      }else{
        if(all(c('vehicleId', 'sequence', 'taskId', 'jobId') %in% names(l$nodes))){
          l$nodes %<>% mutate(label = paste0("V: ",vehicleId,'\nStop:',sequence,"\nT:",taskId, "\nJ:", jobId))
        }else{
          l$nodes %<>% mutate(label = paste0("V: ",vehicleId,'\nStop:',sequence))
        }
        p <- p + geom_label(data = l$nodes, aes(x = x, y = y, label = label))
      }
    }
    if('day' %in% names(l$nodes)){
      # this is then a ndd plot at we can facet by the day
      p <- p + facet_wrap(day ~ ., ncol = 5)
    }
    return (p)
  }
}

#' @export
#' @rawNamespace useDynLib(iceR)
toBytes <- function(m){ return(m$serialize(NULL)) }

#' @export
#' @rawNamespace useDynLib(iceR)
tabulatecvrp <- function(sr, resp){
  pts<-sr$model$points %>% lapply(function(i){return(data.frame(id = i$id, x = i$x, y = i$y))})
  pts[[length(pts) + 1]] <- data.frame(id = sr$model$depot$id, x = sr$model$depot$x, y = sr$model$depot$y)
  pts<- do.call(rbind, pts)
  pts$id %<>% as.character()

  edges<- lapply(resp$routes, function(i){
    do.call(rbind, lapply(i$edges, function(j){
      if(length(j$geometry) == 0){
        return(NULL)
      }
      xs <- lapply(j$geometry, function(k) { return(k$x)}) %>% unlist
      ys <- lapply(j$geometry, function(k) { return(k$y)}) %>% unlist
      geomdf <- data.frame(xs, ys) %>%  st_as_sf(coords = c(1,2)) %>% summarise(do_union = FALSE) %>%  st_cast("LINESTRING")
      geomdf$fromId <- j$from
      geomdf$toId <- j$to
      geomdf$distance <- j$distance
      return(geomdf)
    }))
  })

  for(i in 1:length(edges)){
    if(!is.null(edges[[i]])){
      edges[[i]]$vehicleId <- i
    }
  }
  edges<- do.call(rbind, edges)
  if(!is.null(edges)){
    edges<- edges %>% left_join(pts %>% select(fromId = id, fx = x, fy = y), by = 'fromId') %>%
      left_join(pts %>% select(toId = id, tx = x, ty = y), by = 'toId')
    edges$vehicleId  %<>% as.character()
  }
  nodes <- lapply(resp$routes, function(i){
    if("arrivalTimes" %in% names(i)){
      return(data.frame(locationId =  i$sequence,
                        quantity = cumsum(i$visitCapacities),
                        arrivalTime = i$arrivalTimes))
    }else{
      return(data.frame(locationId =  i$sequence,
                        quantity = cumsum(i$visitCapacities)))
    }
  })

  for(i in 1:length(nodes)){
    if(nrow(nodes[[i]]) > 0){
      nodes[[i]]$vehicleId <- i
      nodes[[i]]$sequence <- 1:nrow(nodes[[i]])
    }
  }
  nodes <- do.call(rbind, nodes)
  nodes$locationId %<>% as.character()
  nodes$vehicleId  %<>% as.character()
  nodes %<>% left_join(pts, by = c('locationId' = 'id'))
  return(list(nodes = nodes, edges = edges))
}

#' @export
#' @rawNamespace useDynLib(iceR)
infeasibilitiesToTable <- function(resp){
  infeasibilityRow <- function(item){
    ctasks = paste0(item$constrainingTaskIds, collapse = ",")
    return( data.frame(dimId = item$dimId, message = item$message,
                       limit = item$limit, value =item$value, count = item$count,
                       constrainingTask = ctasks ) )
  }
  infeasibilityForTask <- function(task){
    lst<- lapply(task$infeasibilityInfo, infeasibilityRow)
    lst<- do.call(rbind, lst)
    lst$taskId <- task$taskId
    lst<- lst[,c(1:5,7,6)]
    return (lst)
  }
  return(do.call(rbind, lapply(resp$infeasibilities, infeasibilityForTask)))
}

#' @export
#' @rawNamespace useDynLib(iceR)
transitRulesToTable <- function(resp, sr){
  do.call(rbind, lapply(resp$routes, function (i){
    do.call(rbind, lapply(i$transitRuleAttributes, function(j){
      return(data.frame(vehicleId = i$vehicleId,
                        ruleId = j$ruleId,
                        dimId = j$dimId,
                        fromStopId = j$fromStopId,
                        toStopId = j$toStopId,
                        startValue = j$startValue,
                        endValue = j$endValue,
                        cost = j$cost))
    }))
  }))
}

#' @export
#' @rawNamespace useDynLib(iceR)
compartmentsToTable <- function(nodes, solResp, solveReq){
  # table the stop information.
  # then for each vehicle, extract the compartments on the vehicle.
  # once we have those compartments, we can
  vids <- nodes$vehicleId %>% unique()
  m <- solveReq$model
  cdims <- m$dimensions$capacityDimensions %>% lapply(function(i){i$id}) %>% unlist
  res <- list()
  if(m$compartments %>% length == 0 & m$compartmentSets %>% length == 0){
    return (res)
  }

  for(v in vids){
    # find the compartments associated with this vehicle.
    vindex <- which(m$vehicles %>% lapply(function(i){i$id == v}) %>% unlist())
    veh<- m$vehicles[[vindex]]
    cset <- NULL
    if(veh$compartmentSetId != ''){
      cset <- m$compartmentSets[[m$compartmentSets %>% lapply(function(i){i$id == veh$compartmentSetId}) %>% unlist]]
    }else{
      vclass <- m$vehicleClasses[[m$vehicleClasses %>% lapply(function(i){i$id == veh$classId }) %>% unlist]]
      if(vclass$compartmentSetId != ''){
        cset <- m$compartmentSets[[m$compartmentSets %>% lapply(function(i){i$id == vclass$compartmentSetId}) %>% unlist]]
      }
    }
    if(is.null(cset)){
      next # there is no compartment set assigned to this vehicle, nothing to do here.
    }
    for(d in cdims){
      # we have to loop through all the capacitated dimensions here as well.
      l <- list()
      for(i in 1:length(m$compartments)){
        if(m$compartments[[i]]$id %in% cset$compartmentIds){
          for(j in 1:length(m$compartments[[i]]$capacities)){
            if(m$compartments[[i]]$capacities[[j]]$dimensionId == d){
              l[[length(l) + 1]] <- data.frame(compartmentId = m$compartments[[i]]$id,
                                               capacity = m$compartments[[i]]$capacities[[j]]$capacity,
                                               stringsAsFactors = F)
            }
          }
        }
      }
      compCaps <- do.call(rbind, l)
      nodes[["tmp_delta"]] <- nodes[[paste0(d, "_end")]] - nodes[[paste0(d, "_start")]]

      stab <- nodes %>% filter(vehicleId == v) %>% arrange(stopId)
      if(nrow(stab) > 2 ){
        for(i in 2:(nrow(stab) - 1)){ # normally do this for each vehicle, but we know we only
          # have one vehicle in this example
          stage = paste0('stop.',stab$stopId[i])
          compCaps<- compCaps %>%
            left_join(stab[2:i, ] %>% # 2 because we skip the vehicle start.
                        select(compartmentId,tmp_delta) %>%
                        group_by(compartmentId) %>%
                        summarise(!!stage := sum(tmp_delta)), #not sure dplyr are proud of this syntax but it works?
                      by = 'compartmentId')
        }
        cnames <- compCaps$compartmentId
        compCaps<- t(compCaps[,2:ncol(compCaps)])
        colnames(compCaps) <- cnames
        compCaps[is.na(compCaps)] <- 0
        compCaps %<>% as.data.frame()
        compCaps$taskId <- c(NA, stab$taskId[2:(nrow(stab)- 1)])
        compCaps$vehicleId <- v
        compCaps$dimension <- d
        compCaps$vehicleId[1] <- compCaps$dimension[1] <- NA
        # so for each vehicle, for each capacitated dimension
        res[[paste0(v, ' dimension:', d)]]<-compCaps
      }
    }
  }
  return(res)
}
