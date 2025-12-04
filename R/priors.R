# some common definitions around priors :-)
# now more like penalization functions

# moje super fast theta link :-) kterou jsem pouzival jeste u puvodnich modelů :-)

#'@import bayesmeta
#library(bayesmeta) # phalft,dhalft

inverse_halft_CDF <- function (x, ...) 1 - phalft(1/x, ...)
inverse_halft_PDF <- function (x, ...) -dhalft(1/x, ...)*(-1)/(x^2)

# definice priorů: log density (_lp) a derivace log density (_lpg) (notace podobna jako v GPstuff)
# nemusi byt nutne normalized (mit integral = 1) a toho hojne vyuzivam.
log_uniform_lp <- function (x, p = NULL) -log(x) # density: 1/x, log density, -log(x)
log_uniform_lpg <- function (x, p = NULL) -1/x
# !!! a ted jsem uplne nahodou zjistil, ze plot exp(log_uniform_lpg(x)), tj exp(-1/x), vypada uplne strasne podobne jako inverse_halft_CDF()!!!!

# tyjo a ted lepsi varianta toho prioru log_uniform... ! Dlouho jsem to hledal a je to jednoduche! Hyperbola, 
# ale fixuje to ten nedostatek, ze 1/0.0001 (to hodne maly minimum) da velky cislo i v logaritmu, a penalizace je pak prilis velka
# staci tam pridat do jmenovatele 1 a ono to pak ma v nule jasny peak => 1!
# zkousel jsem to resit ruznymi exp. variantami ale ty potom u tech vetsich cisel jdou dolu prilis rychle
# ta "2" je scale
# magn - vaha jak moc to ovlivni likelihood
inverse_w_peak_lp <- function (x, magn = 1, scale = 0.5) magn*log(1/(1+x/scale))
inverse_w_peak_lpg <- function (x, magn = 1, scale = 0.5) -magn/scale*1/(1+x/scale)


# definuji si prior, jehoz density je inverse_halft_CDF() :-) Prijde mi, ze je super pro length-scales of squared exponential
inverse_halft_CDF_lp <- function (x, p = NULL) log(inverse_halft_CDF(x, df = 4, scale = 1))
inverse_halft_CDF_lpg <- function (x, p = NULL) 1/inverse_halft_CDF(x, df = 4, scale = 1)*inverse_halft_PDF(x, df = 4, scale = 1)

# penalizacni fce podle Yi et al 2011!
# oni to maji definovane na w coz je 1/length-scale^2, so here x corresponds to length-scale^2
# zde pouzivam "LASSO" v tom spravnem vyznamu
# magn = 2.077589 je spoctena tak aby byla ekvivalentni magn = 1 u inverse_halft_CDF_lp na mych datech (vypocet viz priors_Vanhatalo.R)
#
# orig version for squared length-scales:
#lasso_ls <- function (x, magn = 2.077589) -magn*1/x
#lasso_lsg <- function (x, magn = 2.077589) magn*1/x^2
# version for when x corresponds to length-scale (not squared):
lasso_ls <- function (x, magn = 2.077589) -magn*1/x^2
lasso_lsg <- function (x, magn = 2.077589) magn*2/x^3

# uniform prior
uniform_lp <- function (x) rep(0, length(x))
uniform_lpg <- function (x) rep(0, length(x))

# my: matern-like prior for spatial length-scale(s): exp(-(x/scale)^something)*x
# dobra by byla i halft(x, nu = 4, scale = ...)*x jako ma Vanhatalo (aniz o tom vi :-)), ale nechce se mi pocitat ta derivace
# tak to delam podobnou, jednodussi funkci. Jde o to aby mela maximum a smerem do velkych lengths klesala pomalu.
# magn - vaha jak moc to ovlivni likelihood
exp_x_lp <- function (x, magn = 0.4, scale = 7) magn*(-(x/scale)^0.75 + log(x))
exp_x_lpg <- function (x, magn = 0.4, scale = 7) magn*(-0.75*(x/scale)^-0.25/scale + 1/x)

# Tento prior se pouziva v INLA, jako inla.spde2.pcmatern. Prior je popsan v:
# Fuglstad, G.-A., Simpson, D., Lindgren, F., & Rue, H. (2019). Constructing Priors that Penalize the Complexity of Gaussian Random Fields. Journal of the American Statistical Association, 114(525), 445–452.
# tam ho prezentuji jako joint s priorem pro sigmu, ale neni na tom nic, co by muselo byt joint; udelam z toho 2 priory
# pro kazdy parametr zvlast
# tam ho pouzivaji pro nu <= 1, mailoval jsem si o tom s Finnem Lindgrenem, ale jelikoz ten vzorec nezavisi na nu, mel by byt pouzitelny i pro moje 
# nu = 3/2
# x corresponds to length-scale (not squared)
# d je dimenze
# doporucena lambda dle formule ve Fuglstad, pro ro0 = 3, alpha = 0.05 vychazi jako lambda = 8.99 (!!! ale bacha, zalozeno na cislech vycucanych z prstu!!)
# a pri lambda = 8.99 mi vychazi prepoctova konstanta vuci exp_x_lp  0.1495346 (tj. kolik by melo byt lambda pro lambda = 1 toho exp_x_lp)
# zde je verze bez logaritmu (derivace neotestovana):
#ls_pcmatern_lp <- function (x, lambda, d = 2) d/2*lambda*x^(-d/2-1)*exp(-lambda*x^(-d/2))
#ls_pcmatern_lpg <- function (x, lambda, d = 2) d/2*lambda*  ((-d/2-1)*x^(-d/2-2)*exp(-lambda*x^(-d/2)) + 
#																	   x^(-d/2-1)*exp(-lambda*x^(-d/2))*lambda*d/2*x^(-d/2-1)
# zde s logaritmem:
ls_pcmatern_lp <- function (x, lambda, d = 2) log(d) - log(2) + log(lambda) + log(x)*(-d/2-1) - lambda*x^(-d/2)
ls_pcmatern_lpg <- function (x, lambda, d = 2) (-d/2-1)/x - lambda*(-d/2)*x^(-d/2-1)

# a zde je prior pro sigmu, ktery je soucasti "joint PC prioru" inla.spde2.pcmatern (viz ls_pcmatern_lp)
# sel by pouzit i zvlast pro jine sigma promenne
# u tohoto prioru vysla prepoctova konstanta vuci inverse_w_peak_lp 1.15 (tj. kolik by melo byt lambda pro lambda = 1 toho inverse.., aby 
# to na danych datech bylo stejne) - takze velmi podobne! (pak se ukazalo ze to je spis 1.29; a to se bude lisit per parametr a data... je to jen hruby odhad)
# zde je verze bez logaritmu:
#sigma2_pcmatern_lp <- function (sigma2, lambda) lambda*exp(-lambda*sqrt(sigma2))
#sigma2_pcmatern_lpg <- function (sigma2, lambda) lambda^2*exp(-lambda*sqrt(sigma2))/(2*sqrt(sigma2))
# zde s logaritmem:
sigma2_pcmatern_lp <- function (sigma2, lambda) sigma2_exp_lp(sigma2, lambda)
sigma2_pcmatern_lpg <- function (sigma2, lambda) sigma2_exp_lpg(sigma2, lambda)

sigma2_exp_lp <- function (sigma2, lambda) exp_lp(sqrt(sigma2), lambda)
sigma2_exp_lpg <- function (sigma2, lambda) exp_lpg(sqrt(sigma2), lambda)/(2*sqrt(sigma2))

exp_lp <- function (x, lambda) log(lambda) - lambda*x
exp_lpg <- function (x, lambda) -lambda