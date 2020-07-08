##################################################################################
###                Spatio-temporal distribution of Chikungunya in the          ###
###                         city of Rio de Janeiro, Brazil                     ###
##################################################################################

## Laís Picinini Freitas and Alexandra M. Schmidt

## Last update: 8 July 2020

## Preprint: https://www.medrxiv.org/content/10.1101/2020.06.08.20125757v1 

## Objective: 
# To study the spatio-temporal dynamics of the first chikungunya epidemic in Rio de 
# Janeiro city, exploring the effects of temperature, green area proportion and 
# sociodevelopment index.

rm(list=ls())

library(rstan)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

### Data ------------------------------------------------------------------------------------------

# Cases data can be downloaded at the source at the Rio de Janeiro city hall website: 
# http://www.rio.rj.gov.br/web/sms/exibeConteudo?id=4769664

# Population, sociodevelopment index and land use (green area) data were obtained from the Instituto Pereira Passos at www.data.rio. 
# Temperature data were obtained from the Brazilian National Institute of Meteorology (http://www.inmet.gov.br/portal/index.php?r=estacoes/estacoesAutomaticas), 
# the Brazilian Airspace Control Department (https://www.redemet.aer.mil.br/?i=produtos&p=consulta-de-mensagens-opmet), 
# the Rio de Janeiro State Environmental Institute (http://200.20.53.25/qualiar/home/index), 
# the Rio de Janeiro Municipal Environmental Secretariat (http://www.data.rio/datasets/dados-hor%C3%A1rios-do-monitoramento-da-qualidade-do-ar-monitorar?orderBy=Data) 
# and the Alerta Rio System (http://alertario.rio.rj.gov.br/download/dados-meteorologicos/).


# For the purpose of this script, we will use the number of cases by neighbourhood and week sampled from one chain of the posterior predicted values from the results of model4:
cases <- read.csv('simulated_cases.csv')
y <- matrix(cases$cases, nrow=160, byrow=TRUE)

# covariables: population, sociodevelopment index (sdi), green area proportion
covs_neigh <- read.csv('covariates_rio.csv') 

# covariables matrix
x <- matrix(ncol=2, nrow=length(covs_neigh$CodBairro))
x[,1] <- covs_neigh$sdi
x[,2] <- (covs_neigh$green_area^(1/3))

# temperature data
temperature <- read.csv('temperature_rio_2016.csv')
tmin <- matrix(temperature$min_temperature, nrow=160, byrow=TRUE)
tmincenter <- (tmin-mean(tmin))/sd(tmin)

N <- nrow(y)
TAM <- ncol(y)

# calculating the offset (expected cases)
E <- ((sum(y)/sum(covs_bairro$population))*covs_bairro$population) / TAM


### Neighbourhood matrix --------------------------------------------------------------------------

# functions to build the neighbourhood matrix in stan
library(devtools)
source_url("https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/mungeCARdata4stan.R?raw=TRUE")

# loading the neighbourhood matrix (number of neighbours and which are the neighbours)
numneigh <- scan("NumNeigh.csv")
neigh <- scan("Neigh.csv")

nbs <- mungeCARdata4stan(neigh, numneigh)
N <- nbs$N
node1 <- nbs$node1
node2 <- nbs$node2
N_edges <- nbs$N_edges


### Models ----------------------------------------------------------------------------------------


### Model 0 ####

# Model in stan
model0 <- stan("model0.stan", 
               data=list(N=N, T=TAM, m0=0, C0=5,
                         N_edges=N_edges, node1=node1, node2=node2,
                         y=y, E=E), 
               chains=4, iter=10000,
               control=list(adapt_delta=0.95, max_treedepth=15))

# Summary of the model, WAIC
print(model0, probs=c(0.05, 0.95), digits_summary=3)

require(loo)
loo(model0)
loglikGa <- extract_log_lik(model0)
waic0 <- waic(loglikGa)
waic0

check_hmc_diagnostics(model0)

# Traceplots of the chains

library(bayesplot)

traceplot(model0, pars=c("beta0","sigma"), inc_warmup=FALSE)


### Model 1 ####

# Model in stan
model1 <- stan("model1.stan", 
               data=list(N=N, T=TAM, m0=0, C0=5,
                         y=y, x=x, E=E, K=2), 
               chains=4, iter=10000,
               control=list(adapt_delta=0.95, max_treedepth=15))

# Summary of the model, WAIC
print(model1, probs=c(0.05, 0.95), digits_summary=3)

require(loo)
loo(model1)
loglikGa <- extract_log_lik(model1)
waic0 <- waic(loglikGa)
waic0

check_hmc_diagnostics(model1)

# Traceplots of the chains

library(bayesplot)

traceplot(model1, pars=c("beta0","W"), inc_warmup=FALSE)


### Model 2 ####

# Model in stan
model2 <- stan("model2.stan", 
               data=list(N=N, T=TAM, m0=0, C0=5,
                         N_edges=N_edges, node1=node1, node2=node2,
                         y=y, x=x, E=E, K=2), 
               chains=4, iter=10000,
               control=list(adapt_delta=0.95, max_treedepth=15))

# Summary of the model, WAIC
print(model2, probs=c(0.05, 0.95), digits_summary=3)

require(loo)
loo(model2)
loglikGa <- extract_log_lik(model2)
waic0 <- waic(loglikGa)
waic0

check_hmc_diagnostics(modelo)

# Traceplots of the chains

library(bayesplot)

traceplot(model2, pars=c("beta0","sigma"), inc_warmup=FALSE)
traceplot(model2, pars=c("W"), inc_warmup=FALSE)


### Model 3 ####

# Model in stan
model3 <- stan("model3.stan", 
               data=list(N=N, T=TAM, m0=0, C0=5,
                         N_edges=N_edges, node1=node1, node2=node2,
                         y=y, tmincenter=tmincenter, E=E), 
               chains=4, iter=10000,
               control=list(adapt_delta=0.95, max_treedepth=15))

# Summary of the model, WAIC
print(model3, probs=c(0.05, 0.95), digits_summary=3)

require(loo)
loo(model3)
loglikGa <- extract_log_lik(model3)
waic0 <- waic(loglikGa)
waic0

check_hmc_diagnostics(modelo)

# Traceplots of the chains

library(bayesplot)

traceplot(model3, pars=c("beta0","sigma"), inc_warmup=FALSE)


### Model 4 ---------------------------------------------------------------------------------------

# Model in stan
model4 <- stan("model4.stan", 
                         data=list(N=N, T=TAM, m0=0, C0=5,
                                   N_edges=N_edges, node1=node1, node2=node2,
                                   y=y, x=x, tmincenter=tmincenter, E=E, K=2), 
                         chains=4, iter=10000,
                         control=list(adapt_delta=0.95, max_treedepth=15))

# Summary of the model, WAIC
print(model4, probs=c(0.05, 0.95), digits_summary=3)

require(loo)
loo(model4)
loglikGa <- extract_log_lik(model4)
waic0 <- waic(loglikGa)
waic0

check_hmc_diagnostics(modelo)

# Traceplots of the chains

library(bayesplot)

traceplot(model4,pars=c("beta0","sigma"), inc_warmup=FALSE)
traceplot(model4, pars=c("W"), inc_warmup=FALSE)


