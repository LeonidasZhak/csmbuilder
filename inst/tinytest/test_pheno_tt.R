library(tinytest)
library(csmbuilder)

###########################################
# Phenology Thermal Time Model Test Example
###########################################

generate_Tair <- function(doy = NULL){
  if(is.null(doy)) doy <- 1:365
  15*cos((doy-200)/365*2*pi)+15
}

# Initial conditions
state <- c(du = 0)

# parameters
parameters <-
  c(ko = 1,
    H = 136000,
    E = 95000,
    To = 20)

# Rate equation function
pheno_tt_ref_dt <- function(t, state, parms, wth){
  with(as.list(c(state, parms)), {
    Tt <- csm_get_at_t(wth[,2], wth[,1], t)
    list(
      csm_mod_arr(Tt, ko, H, E, To)
    )
  })
}

# Generate air temperature forcings
wth <-
  matrix(
    c(seq_along(c(300:365,1:180)) - 1,
      generate_Tair(c(300:365,1:180))),
         ncol = 2)

# Check mod_arr()
expect_equal(
  with(as.list(parameters), {
    csm_mod_arr(To, ko, H, E, To)
  }),
  parameters["ko"],
  info = "mod_arr",
  check.attributes = FALSE
)

# Check pheno_tt_ref_dt()
expect_equal(
  pheno_tt_ref_dt(0,
                      state = state,
                      parms = parameters,
                      wth = wth),
  with(as.list(parameters), {
    list(csm_mod_arr(wth[1,2], ko, H, E, To))
  }),
  info = "pheno_tt_ref_dt(); t=0"
)

expect_equal(
  pheno_tt_ref_dt(0.5,
                      state = state,
                      parms = parameters,
                      wth = wth),
  with(as.list(parameters), {
    list(csm_mod_arr(mean(wth[1:2,2]), ko, H, E, To))
  }),
  info = "pheno_tt_ref_dt(); t=0.5"
)

# Specify times at which to report output
times <- csm_time_vector(0, nrow(wth)-1, dt = 0.01)

# Create list of integration methods to test:
integ_list <- c("euler", "rk4")

# Run integration
pheno_tt_ref_out <-
  integ_list |>
  (\(.x) setNames(.x, .x))() |>
  lapply(\(.method){
    deSolve::ode(
      y = state,
      times = times,
      func = pheno_tt_ref_dt,
      parms = parameters,
      method = .method,
      wth = wth
    ) |>
      as.data.frame()
  })

#######################################################
# Create Phenology Thermal Time Model with csmbuilder
#######################################################

# Define state variables
sp_state <- csm_create_state(
  c("du"),
  definition = c("development units"),
  units = c("physiological days"),
  expression(~csm_mod_arr(Tair_t, ko, H, E, To)))

# Define parameters
sp_parameters <- csm_create_parameter(
  c("ko", "H", "E", "To"),
  definition = c("relative reaction rate",
                 "deactivation energy",
                 "activation energy",
                 "optimum temperature"),
  units = c("unitless", "Joules per mole",
            "Joules  per mole", "degrees Celsius"))

# Define weather inputs
sp_wth_inp <- csm_create_variable(
  c("wtime", "Tair"),
  definition = c("time of weather observation",
                 "air temperature"),
  units = c("days after planting", "degrees Celsius"))

sp_wth <- csm_create_data_structure(
  name = "wth",
  definition = "weather data",
  variables = c(sp_wth_inp),
  n_dim = 2
)

sp_wth_t <- csm_create_transform(
  "Tair_t",
  definition = "air temperature at t",
  units = "degrees Celsius",
  equation = ~csm_get_at_t(Tair, wtime, t))

# Define model
pheno_tt_model <-
  csm_create_model(
    name = "pheno_tt",
    state = sp_state,
    parms = sp_parameters,
    wth = sp_wth,
    wth_t = sp_wth_t)

# Create function for calculating rates
pheno_tt_dydt <-
  csm_render_model(
    model = pheno_tt_model,
    output_type = "function",
    language = "R")

pheno_tt_Rcpp <-
  csmbuilder::csm_render_model(
    model = pheno_tt_model,
    output_type = "function",
    language = "Rcpp",
    insert_functions = list(
      csm_get_at_t = c("double csm_get_at_t(NumericVector x, NumericVector t_ind,",
                       "NumericVector t) ",
                       "{",
                       "      int i = 0;",
                       "      int j = t_ind.size() - 1;",
                       "      while(i < j && (i + 1) != j){",
                       "        if(t_ind[i] > t[0]) i = 0;",
                       "        if(t_ind[j] < t[0]) j = t_ind.size() - 1;",
                       "        double m = (j - i)/(t_ind[j] - t_ind[i]);",
                       "        double b = j - m*t_ind[j];",
                       "        double ij = m*t[0]+b;",
                       "        i = std::floor(ij);",
                       "        if(t_ind[i] == t[0]) return x[i];",
                       "        j = std::ceil(ij);",
                       "      }",
                       "  if(t_ind[i] == t[0]){",
                       "    return x[i];",
                       "  }else{",
                       "    return x[i]*(t_ind[j]-t[0]) + x[j]*(t[0]-t_ind[i])/(t_ind[j] - t_ind[i]);",
                       "  }",
                       "}"),
      csm_mod_arr = c("double csm_mod_arr(double Tt, double ko, double H,",
                      "                   double E, double To)",
                      "{",
                      "    double R = 8.314;",
                      "    return ko * (H * exp(E/R * (1/(To + 273.15) - 1/(Tt + 273.15))))/(H - ",
                      "        E * (1 - exp(H/R * (1/(To + 273.15) - 1/(Tt + 273.15)))));",
                      "}")
    ))

# Run integration
pheno_tt_dydt_out <-
  integ_list |>
  (\(.x) setNames(.x, .x))() |>
  lapply(\(.method){
    csm_run_sim(
      model_function = pheno_tt_dydt,
      y_init = state,
      t = times,
      parms = parameters,
      wth = wth,
      method = .method
    )
  })

pheno_tt_Rcpp_out <-
  integ_list |>
  (\(.x) setNames(.x, .x))() |>
  lapply(\(.method){
    csm_run_sim(
      model_function = pheno_tt_Rcpp,
      y_init = state,
      t = times,
      parms = parameters,
      wth = wth,
      method = .method
    )
  })

for(integ_method in integ_list){
  expect_equal(
    current = pheno_tt_dydt_out[[integ_method]],
    target = pheno_tt_ref_out[[integ_method]],
    info = integ_method
  )
}

for(integ_method in integ_list){
  expect_equal(
    current = pheno_tt_Rcpp_out[[integ_method]],
    target = pheno_tt_ref_out[[integ_method]],
    info = integ_method
  )
}
