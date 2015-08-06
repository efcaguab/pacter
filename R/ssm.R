#' ssm
#'
#' @return jj
#' @export
#'
#' @examples
#' #
ssm <- function(detections,
                acoustic_array,
                range_test,
                sentinel,
                delta_t,
                lon_grid = list(from = NA, to = NA, extra = 0.2),
                lat_grid = list(from = NA, to = NA, extra = 0.2), 
                cell_width = NULL, 
                n_cells = 2500,
                start_time = NULL,
                end_time = NULL,
                params = list(logbx = -9, logby = -9, logitcor = 0, logsigma = 0.4),
                safe=FALSE,
                verbose=FALSE,
                admb_errors=c("stop","warn","ignore")){
  
  # Ensure that the working directory is the same than at the beginning
  cwd <- getwd()
  on.exit(setwd(cwd))
  
  # prepare the grid
  state_space_grid <- ssm_grid(acoustic_array, 
                               lon_grid,
                               lat_grid,
                               cell_width,
                               n_cells)
  write_ssm_grid(state_space_grid)

  # prepare the data
  state_space_data <- ssm_data(detections, 
                               acoustic_array,
                               range_test,
                               sentinel,
                               delta_t,
                               ssm_grid = state_space_grid,
                               start_time = NULL, end_time = NULL)
  write_ssm_data(state_space_data)
  
  default_params <- identical(params, list(logbx = -9, logby = -9, logitcor = 0, logsigma = 0.4))
  
  if (default_params){
    
    # cf it has never been compiled, do so
    if(system.file("onssm", package= "pacter") == ""){
      setwd(system.file(package = "pacter"))
      compile_admb("onssm", 
                   safe = safe, 
                   re = FALSE,
                   verbose = verbose,
                   admb_errors = admb_errors)
    } 
    
    # copy the compiled software into the temporary directory
    file.copy(system.file("onssm", package= "pacter"), tempdir())
    file.copy(system.file("onssm.tpl", package= "pacter"), tempdir())
    file.create(file.path(tempdir(), "onssm.dat"))
  } else {
    # copy the tpl
    file.copy(system.file("onssm.tpl", package= "pacter"), tempdir())
    file.create(file.path(tempdir(), "onssm.dat"))
    # modify the tpl
    the_tpl <- readLines(file.path(tempdir(), "onssm.tpl"))
    logbx <- paste("  logbx = ", params$logbx, ";", sep ="")
    logby <- paste("  logby = ", params$logby, ";", sep ="")
    logitcor <- paste("  logitcor = ", params$logitcor, ";", sep ="")
    logsigma <- paste("  logsigma = ", params$logsigma, ";", sep ="")
    the_tpl[352] <- logbx
    the_tpl[353] <- logby
    the_tpl[354] <- logitcor
    the_tpl[355] <- logsigma
    writeLines(the_tpl, file.path(tempdir(), "onssm.tpl"))
    # compile it
    # make the temporary directory the working directory while executing
    setwd(tempdir())
    compile_admb("onssm", 
                 safe = safe, 
                 re = FALSE,
                 verbose = verbose,
                 admb_errors = admb_errors)
  }

  
  
  
  # make the temporary directory the working directory while executing
  setwd(tempdir())
  # execute the state state space model
  admbOut <- 
    run_admb("onssm", 
             verbose = verbose, 
             mcmc = FALSE, 
             mcmc.opts = mcmc.control(),
             profile = FALSE,
             admb_errors = admb_errors)
  
  # read results back into R
  read_admb("onssm", verbose = verbose, admbOut = admbOut)
}