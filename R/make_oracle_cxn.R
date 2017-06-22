#' @title make_oracle_cxn
#' @description This function contains the logic for creating a connection to
#' Oracle.
#' @param usepkg default is \code{'rodbc'}. This indicates whether the connection to Oracle should
#' use \code{'rodbc'} or \code{'roracle'} to connect.  rodbc is slightly easier to setup, but
#' roracle will extract data ~ 5x faster.
#' @family data
#' @importFrom RODBC odbcConnect
#' @author  Mike McMahon, \email{Mike.McMahon@@dfo-mpo.gc.ca}
#' @export
#' @note This is stolen directly from the get_data() function of bio.datawrangling
#' (also created by me)
#' # make_oracle_cxn -------------------------------------------------------------
make_oracle_cxn <- function(usepkg='rodbc') {

  # use.roracle <-function(oracle.dsn, oracle.username, oracle.password){
  #   assign('oracle_cxn', ROracle::dbConnect( DBI::dbDriver("Oracle"), oracle.username,oracle.password,oracle.dsn))
  #   if (class(oracle_cxn)[1]=="OraConnection") {
  #     cat("\nSuccessfully connected to Oracle via ROracle\n")
  #     results = list('roracle', oracle_cxn)
  #     return(results)
  #   } else {
  #     cat("\nROracle attempt failed - trying RODBC\n")
  #     use.rodbc(oracle.dsn, oracle.username, oracle.password)
  #   }
  # }
  use.rodbc <-function(oracle.dsn, oracle.username, oracle.password){
    #warning supressed to hide password should the connection fail
    suppressWarnings(
      assign('oracle_cxn', odbcConnect(oracle.dsn, uid = oracle.username, pwd = oracle.password, believeNRows = F))
    )
    if (class(oracle_cxn) == "RODBC") {
      cat("\nSuccessfully connected to Oracle via RODBC")
      results = list('rodbc', oracle_cxn)
      return(results)
    } else {
      cat("\nRODBC connection attempt failed\n")
      return(-1)
    }
  }

  # if (class(oracle_cxn) == 'RODBC'){
  #   results = list('rodbc', oracle_cxn)
  #   return(results)
  # }else if (class(oracle_cxn)[1]=="OraConnection") {
  #   results = list('roracle', oracle_cxn)
  #   return(results)
  # } else {
    #get connection info - only prompt for values not in rprofile
    if (exists('oracle.username')){
      oracle.username <- oracle.username
      cat("\nUsing stored 'oracle.username'")
    }else{
      oracle.username <-
        readline(prompt = "Oracle Username: ")
      print(oracle.username)
    }

    if (exists('oracle.password')){
      oracle.password <- oracle.password
      cat("\nUsing stored 'oracle.password'")
    } else {
      oracle.password <-
        readline(prompt = "Oracle Password: ")
      print(oracle.password)
    }

    if (exists('oracle.dsn')){
      oracle.dsn <- oracle.dsn
      cat("\nUsing stored 'oracle.dsn'")
    }else{
      oracle.dsn <-
        readline(prompt = "Oracle DSN (e.g. PTRAN): ")
      print(oracle.dsn)
    }

    # if (usepkg=='roracle'){
    #   use.roracle(oracle.dsn, oracle.username, oracle.password)
    # }else{
      use.rodbc(oracle.dsn, oracle.username, oracle.password)
    # }
  }
# }
