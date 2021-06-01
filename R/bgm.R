#'  Bayesian global model, on encapsule un modèle global générique dans la fonction
#' @param tab_input mcmc list
#' @param n_proj nombre d'annee de projection
#' @examples
#'
#' #Les données d'entrée
#' data(tableau_pa)
#'
#' resultats_pa<-delta (tab=tableau_pa, esp="BOBO",list_param=c("mois","zone","type_pirogue","motorisation","zone2","saison","engin_peche2"), type_donnee="commercial", effort="auto", titre="PA", param_test=c("mois","zone","type_pirogue","motorisation","zone2","saison","engin_peche2"), espece_id_list='espece', var_eff_list= c("nb_jour_peche", "nb_sorties"), ope_id=c("presence","code_pays","annee", "mois", "zone", "zone2" , "saison", "type_pirogue" , "motorisation" , "engin_peche", "engin_peche2"), col_capture='captures', logtrans="auto", interactions="auto", facteur_flotille="engin_peche2", seuil=0.05)
#' pre_global_model_pa<-sqldf('select annee as Year,sum(i_ab) as C,avg("E.dens") as I from resultats_pa group by annee',drv='SQLite')

#' bgm(pre_global_model_pa)


#' @export

bgm<- function (tab_input,min_c,max_c,Er,CVr,CV_Cobs,CV_Iobs,n_proj=0){

  # Load the model
  # ----------------------------------------

  # Load the model, written in a .txt file

  model_file<-paste0(path.package("demerstem"), "/model_BiomProd_WithScenario_JAGS.txt")

  source(model_file)

  # Write the model specification into a virtual text file with name "model", to be found by JAGS
  # The file "model" will be used by functions "jags.model"
textConnection

model <- textConnection(model_JAGS)


  # Load and format data
  # -----------------------------------------


  n_obs <- length(tab_input$Year)

  Year <- tab_input$Year

  #Il manque 2014 à voir pourquoi ? (Dans les captures ?)


  # Needed to build the equilibrium curve
  B_e <- seq(from = 0, to = 1500, by = 50)
  n_equi <- length(B_e)

  # Format data as a list to be read in JAGS

  data <- list("I_obs" = tab_input$I, "C_obs" = tab_input$C, "n_obs" = n_obs, "n_proj" = n_proj,
               "B_e" = B_e, "n_equi" = n_equi, "Year"=Year,"min_c"=min_c,"max_c"=max_c,"Er"=Er,"CVr"=CVr,"CV_Cobs"=CV_Cobs,"CV_Iobs"=CV_Iobs)


  # MCMC options
  # ----------------------------------------------------------------------

  n.chains = 3

  # Adapt the MCMC samplers with n.adapt iterations

  n.adapt = 10000

  # Iteration after adapting

  n.burnin <- 1000
  n.stored <- 1000
  n.thin <- 10
  n.iter <- n.stored*n.thin

  # MANAGING NUMBER of ITERATIONS
  # (works for both para = T and para = F)
  # total number of REALIZED iterations = n.iter
  # total number of STORED iterations = n.iter/n.thin



  # Run the model to save MCMC results
  # ---------------------------------------------------------------------

  # Variable to store

  monitor <- c( "B", "h", "D", "C",
                "r", "r_p", "K", "K_p", "q", "sigma2p",
                "C_MSY", "C_MSY_p", "B_MSY", "h_MSY",
                "risk", "Over_C","Over_h",
                "C_e",
                "I_pred", "C_pred")

  # Compile the model, create a model object with the name "model.compiled.jags" and adapt the MCMC samplers
  # with n.adapt iterations

  print("adapting phase")
  model.compiled <- jags.model(file = model, data=data, n.chain=n.chains, n.adapt=n.adapt)

  # Iteration after adapting

  # Start to compute CPU time
  ptm <- proc.time()

  # Burnin period (not stored)

  print("burn-in")
  update(model.compiled, n.iter=n.burnin)

  # Store mcmc samples

  print("mcmc stored for results")
  mcmc <- coda.samples(model=model.compiled,variable.names=monitor,n.iter=n.iter,thin=n.thin)

  time.to.run.mcmc <- proc.time() - ptm
  print(time.to.run.mcmc)

  # ----------------------------------------------------------------------------------------------


  # Run the model to compute DIC
  # ---------------------------------------------------------------------

  # Start from a compiled model that has already been updated

  dic.pD <- dic.samples(model.compiled, n.iter, "pD")
  dic.pD		# Deviance Information Criterion

  # Alternative penalization of the Deviance
  # dic.popt <- dic.samples(model.compiled.jags, n.iter, "popt")
  # dic.popt
  # -----------------------------------------------------------------------



  # --------------------------------------------
  # Work with mcmc.list
  # --------------------------------------------

  # "mcmc" is an object of the class "mcmc.list" (see package library(coda)
  # to explore, plot ... mcmc objects

  is(mcmc)

  # Names of the variables stored in the mcmc list

  varnames(mcmc)

  return(mcmc)
  }

