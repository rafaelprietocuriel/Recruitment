
require(bvpSolve); require(deSolve); require(ineq); require(viridis)
require(akima); require(fields); require(cmocean); require(phaseR)
require(knitr); require(phaseR); require(graphics); require(captioner)
require(latex2exp); require(igraph); require(nloptr); library(spatstat)

#### cartels have a power law distribution w beta = -1
#### conflict = (C + 1)/(A + 1)

#### steps are weeks
#### the system is considered for 520 steps (10 years)

#### load graph, initial members and conflict Matrix
#### the way the data is registered, it has directionality
#### so the first steps are to avoid direction
{
  Nod <- read.csv("BACRIM2020_Nodes.csv")
  Rivals <- read.csv("BACRIM2020_Rivals.csv")
  Allies <- read.csv("BACRIM2020_Alliances.csv")
  NCartels <- dim(Nod)[1]
  
  #### initial criminals
  Crime0 <- 115000*1 ### total criminals t0
  de <- (1:NCartels)^(-1) ### here, kappa for power law with coeff = -1
  epsF <- 1
  
  #### Rivals network
  names(Rivals)[2] <- "from"
  names(Rivals)[4] <- "to"
  names(Nod)[1] <- "name"
  G <- graph_from_data_frame(Rivals[, c(2,4)],
                             directed = T,
                             vertices = Nod)
  E(G)$weight <- Rivals$weight
  EM <- as_adjacency_matrix(G,
                            sparse = F,
                            attr = "weight")
  EM <- pmax(EM, t(EM)) ### so symmetric relation
  G <- graph_from_adjacency_matrix(EM,
                                   mode = "undirected",
                                   weighted = T)
  RivalsG <- as_data_frame(G, what = "edges")
  
  #### Allies network
  names(Allies)[2] <- "from"
  names(Allies)[4] <- "to"
  names(Nod)[1] <- "name"
  H <- graph_from_data_frame(Allies[, c(2,4)],
                             directed = T, 
                             vertices = Nod)
  E(H)$weight <- Allies$weight
  EM <- as_adjacency_matrix(H,
                            sparse = F,
                            attr = "weight")
  EM <- pmax(EM, t(EM)) ### so symmetric relation
  H <- graph_from_adjacency_matrix(EM,
                                   mode = "undirected",
                                   weighted = T)
  AlliesH <- as_data_frame(H, what = "edges")
  
  #### compute strength (like degree)
  Nod$strengthRivals <- strength(G)
  Nod$strengthAllies <- strength(H)
  Nod$strengthTotal <- Nod$strengthRivals + Nod$strengthAllies
  
  ####both netw and degree
  Nod <- Nod[order(Nod$strengthTotal, decreasing = T), ]
  G <- graph_from_data_frame(RivalsG,
                             directed = F, 
                             vertices = Nod)
  E(G)$weight <- RivalsG$weight
  H <- graph_from_data_frame(AlliesH,
                             directed = F, 
                             vertices = Nod)
  E(H)$weight <- AlliesH$weight
  
  #### weighted conflict matrix
  CMatrix <- as_adjacency_matrix(G, 
                                 attr = "weight",
                                 sparse = F,
                                 type = "both")
  
  ConfM <- (CMatrix + 1)/(Nod$strengthAllies + 1)
  diag(ConfM) <- 0
  
  #### assign number of criminals
  y0 <- rpois(NCartels, lambda = de/sum(de)*Crime0)
  Nod$EstCrime <-  y0
}

#### splines data
{
  ToIn <- read.csv("Trends2012_2021.csv")
  FH <- splinefun(ToIn$YEAR, ToIn$homicide/52, method = "natural")
  FM <- splinefun(ToIn$YEAR, ToIn$missings/52, method = "natural")
  FA <- splinefun(ToIn$YEAR, ToIn$arrests/52, method = "natural")
  tA <- seq(from = 2012, to = 2022, length.out = 520)
  WeeklyCasualties <- FH(tA) + FM(tA)
  WeeklyArrest <- FA(tA)
}

#### plot fun for OpPar
{
  PlOpt <- function(OpPar, DifEq, WeeklyCasualties){
    Crime0 <- OpPar$par[5]*1000
    NCartels = 150
    de <- (1:NCartels)^(-1)
    y0 <- rpois(NCartels, lambda = de/sum(de)*Crime0)
    yini <- c(y0, 
              D = 0, 
              R = 0)
    params <- c(rho =OpPar$par[1]/1000, #### recruitment
                eta =OpPar$par[2]*1,    #### incapacitation
                theta = OpPar$par[3]/100000000, ## conflict
                omega = OpPar$par[4]/1000000) ## fragmentation
    times <- seq(from = 1, to = Years*52, by = 1)
    out <- ode(times = times,
               y = yini, 
               func = DifEq, 
               ConfM = ConfM,
               St = 0*y0+1,
               parms = params)
    InstArrest <- params[2]*apply(out[, 2:(NCartels +1)], 1, sum)/apply(out[, 2:(NCartels +1)], 1, sum)
    InstDeath <- 0 * InstArrest
    InstRecruit <- 0 * InstArrest
    TotCrime <- apply(out[, 2:152], 1, sum)
    for(k in 1:length(InstDeath)){
      InstRecruit[k] <- params[1]*(0*(1:NCartels) + 1) %*% out[k, 2:(NCartels +1)]
      InstDeath[k] <- params[3]*out[k, 2:(NCartels +1)] %*% ConfM %*% out[k, 2:(NCartels +1)]
    } 
    par(mfrow = c(2,2), 
        mar = c(1,1,1,1)*1.8,
        oma = c(1,1,1,1)*0.2)
    plot(TotCrime, type = "h", 
         ylim = c(0, max(TotCrime)),
         col = 3)
    text(1, max(TotCrime)*0.9, "Total criminals", adj = 0)
    for(k in 1:10){points(c(1,520), c(1,1)*k*50000, type = "l", lty = 2)}
    plot(InstDeath, type = "h", col = 4,
         ylim = c(0, max(InstDeath)))
    points(WeeklyCasualties, type = "l", lwd = 2, col = "gold")
    text(1, max(InstDeath)*0.9, "weekly deaths", adj = 0)
    
    plot(InstRecruit, 
         ylim = c(0, max(InstRecruit)),
         type = "h", col = 2)
    text(1, max(InstRecruit)*0.9, "weekly recruits", adj = 0)
    
    plot(out[, 2], type = "l", 
         ylim = c(1, max(out[, 2])),
         log = "y",
         col = NA)
    for(k in 1:NCartels){points(out[, k+1], type = "l", 
                                col = rainbow(NCartels)[k])}
    text(1, max(out[, k+1])*0.9, "cartel size", adj = 0)
    text(1, max(out[, k+1])*0.09, 
         paste("Max C = ", round(1000*OpPar$par[1]/OpPar$par[4])),
         adj = 0)
  }
}

#### compute parameters from fixed values
{
  #### years to compute
  Years = 10
  #### theta lethality so there are 200 kills per week
  theta = 1.549231/100000000
  #### eta incapacitation
  eta = 110
  #### rho recruitment
  rho = 2.527474/1000
  #### omega fragmentation
  omega = 0.04/1000000
}

#### ALL CARTELS DYNAMICS. Diff eq and run dynamics once and keep Instants  
{
  # recruitment, incarceration, conflict, fragmentation
  DifEq <- function(t, y, params, ConfM, St){
    N <- 150
    TC <- sum(y[1:N]) ### total criminals
    dy1<- y[1]*(params[1] * St[1] - params[2]/TC - params[3]*ConfM[1,]%*%y[1:N] - params[4]*y[1])
    dy2<- y[2]*(params[1] * St[2] - params[2]/TC - params[3]*ConfM[2,]%*%y[1:N] - params[4]*y[2])
    dy3<- y[3]*(params[1] * St[3] - params[2]/TC - params[3]*ConfM[3,]%*%y[1:N] - params[4]*y[3])
    dy4<- y[4]*(params[1] * St[4] - params[2]/TC - params[3]*ConfM[4,]%*%y[1:N] - params[4]*y[4])
    dy5<- y[5]*(params[1] * St[5] - params[2]/TC - params[3]*ConfM[5,]%*%y[1:N] - params[4]*y[5])
    dy6<- y[6]*(params[1] * St[6] - params[2]/TC - params[3]*ConfM[6,]%*%y[1:N] - params[4]*y[6])
    dy7<- y[7]*(params[1] * St[7] - params[2]/TC - params[3]*ConfM[7,]%*%y[1:N] - params[4]*y[7])
    dy8<- y[8]*(params[1] * St[8] - params[2]/TC - params[3]*ConfM[8,]%*%y[1:N] - params[4]*y[8])
    dy9<- y[9]*(params[1] * St[9] - params[2]/TC - params[3]*ConfM[9,]%*%y[1:N] - params[4]*y[9])
    dy10<- y[10]*(params[1] * St[10] - params[2]/TC - params[3]*ConfM[10,]%*%y[1:N] - params[4]*y[10])
    dy11<- y[11]*(params[1] * St[11] - params[2]/TC - params[3]*ConfM[11,]%*%y[1:N] - params[4]*y[11])
    dy12<- y[12]*(params[1] * St[12] - params[2]/TC - params[3]*ConfM[12,]%*%y[1:N] - params[4]*y[12])
    dy13<- y[13]*(params[1] * St[13] - params[2]/TC - params[3]*ConfM[13,]%*%y[1:N] - params[4]*y[13])
    dy14<- y[14]*(params[1] * St[14] - params[2]/TC - params[3]*ConfM[14,]%*%y[1:N] - params[4]*y[14])
    dy15<- y[15]*(params[1] * St[15] - params[2]/TC - params[3]*ConfM[15,]%*%y[1:N] - params[4]*y[15])
    dy16<- y[16]*(params[1] * St[16] - params[2]/TC - params[3]*ConfM[16,]%*%y[1:N] - params[4]*y[16])
    dy17<- y[17]*(params[1] * St[17] - params[2]/TC - params[3]*ConfM[17,]%*%y[1:N] - params[4]*y[17])
    dy18<- y[18]*(params[1] * St[18] - params[2]/TC - params[3]*ConfM[18,]%*%y[1:N] - params[4]*y[18])
    dy19<- y[19]*(params[1] * St[19] - params[2]/TC - params[3]*ConfM[19,]%*%y[1:N] - params[4]*y[19])
    dy20<- y[20]*(params[1] * St[20] - params[2]/TC - params[3]*ConfM[20,]%*%y[1:N] - params[4]*y[20])
    dy21<- y[21]*(params[1] * St[21] - params[2]/TC - params[3]*ConfM[21,]%*%y[1:N] - params[4]*y[21])
    dy22<- y[22]*(params[1] * St[22] - params[2]/TC - params[3]*ConfM[22,]%*%y[1:N] - params[4]*y[22])
    dy23<- y[23]*(params[1] * St[23] - params[2]/TC - params[3]*ConfM[23,]%*%y[1:N] - params[4]*y[23])
    dy24<- y[24]*(params[1] * St[24] - params[2]/TC - params[3]*ConfM[24,]%*%y[1:N] - params[4]*y[24])
    dy25<- y[25]*(params[1] * St[25] - params[2]/TC - params[3]*ConfM[25,]%*%y[1:N] - params[4]*y[25])
    dy26<- y[26]*(params[1] * St[26] - params[2]/TC - params[3]*ConfM[26,]%*%y[1:N] - params[4]*y[26])
    dy27<- y[27]*(params[1] * St[27] - params[2]/TC - params[3]*ConfM[27,]%*%y[1:N] - params[4]*y[27])
    dy28<- y[28]*(params[1] * St[28] - params[2]/TC - params[3]*ConfM[28,]%*%y[1:N] - params[4]*y[28])
    dy29<- y[29]*(params[1] * St[29] - params[2]/TC - params[3]*ConfM[29,]%*%y[1:N] - params[4]*y[29])
    dy30<- y[30]*(params[1] * St[30] - params[2]/TC - params[3]*ConfM[30,]%*%y[1:N] - params[4]*y[30])
    dy31<- y[31]*(params[1] * St[31] - params[2]/TC - params[3]*ConfM[31,]%*%y[1:N] - params[4]*y[31])
    dy32<- y[32]*(params[1] * St[32] - params[2]/TC - params[3]*ConfM[32,]%*%y[1:N] - params[4]*y[32])
    dy33<- y[33]*(params[1] * St[33] - params[2]/TC - params[3]*ConfM[33,]%*%y[1:N] - params[4]*y[33])
    dy34<- y[34]*(params[1] * St[34] - params[2]/TC - params[3]*ConfM[34,]%*%y[1:N] - params[4]*y[34])
    dy35<- y[35]*(params[1] * St[35] - params[2]/TC - params[3]*ConfM[35,]%*%y[1:N] - params[4]*y[35])
    dy36<- y[36]*(params[1] * St[36] - params[2]/TC - params[3]*ConfM[36,]%*%y[1:N] - params[4]*y[36])
    dy37<- y[37]*(params[1] * St[37] - params[2]/TC - params[3]*ConfM[37,]%*%y[1:N] - params[4]*y[37])
    dy38<- y[38]*(params[1] * St[38] - params[2]/TC - params[3]*ConfM[38,]%*%y[1:N] - params[4]*y[38])
    dy39<- y[39]*(params[1] * St[39] - params[2]/TC - params[3]*ConfM[39,]%*%y[1:N] - params[4]*y[39])
    dy40<- y[40]*(params[1] * St[40] - params[2]/TC - params[3]*ConfM[40,]%*%y[1:N] - params[4]*y[40])
    dy41<- y[41]*(params[1] * St[41] - params[2]/TC - params[3]*ConfM[41,]%*%y[1:N] - params[4]*y[41])
    dy42<- y[42]*(params[1] * St[42] - params[2]/TC - params[3]*ConfM[42,]%*%y[1:N] - params[4]*y[42])
    dy43<- y[43]*(params[1] * St[43] - params[2]/TC - params[3]*ConfM[43,]%*%y[1:N] - params[4]*y[43])
    dy44<- y[44]*(params[1] * St[44] - params[2]/TC - params[3]*ConfM[44,]%*%y[1:N] - params[4]*y[44])
    dy45<- y[45]*(params[1] * St[45] - params[2]/TC - params[3]*ConfM[45,]%*%y[1:N] - params[4]*y[45])
    dy46<- y[46]*(params[1] * St[46] - params[2]/TC - params[3]*ConfM[46,]%*%y[1:N] - params[4]*y[46])
    dy47<- y[47]*(params[1] * St[47] - params[2]/TC - params[3]*ConfM[47,]%*%y[1:N] - params[4]*y[47])
    dy48<- y[48]*(params[1] * St[48] - params[2]/TC - params[3]*ConfM[48,]%*%y[1:N] - params[4]*y[48])
    dy49<- y[49]*(params[1] * St[49] - params[2]/TC - params[3]*ConfM[49,]%*%y[1:N] - params[4]*y[49])
    dy50<- y[50]*(params[1] * St[50] - params[2]/TC - params[3]*ConfM[50,]%*%y[1:N] - params[4]*y[50])
    dy51<- y[51]*(params[1] * St[51] - params[2]/TC - params[3]*ConfM[51,]%*%y[1:N] - params[4]*y[51])
    dy52<- y[52]*(params[1] * St[52] - params[2]/TC - params[3]*ConfM[52,]%*%y[1:N] - params[4]*y[52])
    dy53<- y[53]*(params[1] * St[53] - params[2]/TC - params[3]*ConfM[53,]%*%y[1:N] - params[4]*y[53])
    dy54<- y[54]*(params[1] * St[54] - params[2]/TC - params[3]*ConfM[54,]%*%y[1:N] - params[4]*y[54])
    dy55<- y[55]*(params[1] * St[55] - params[2]/TC - params[3]*ConfM[55,]%*%y[1:N] - params[4]*y[55])
    dy56<- y[56]*(params[1] * St[56] - params[2]/TC - params[3]*ConfM[56,]%*%y[1:N] - params[4]*y[56])
    dy57<- y[57]*(params[1] * St[57] - params[2]/TC - params[3]*ConfM[57,]%*%y[1:N] - params[4]*y[57])
    dy58<- y[58]*(params[1] * St[58] - params[2]/TC - params[3]*ConfM[58,]%*%y[1:N] - params[4]*y[58])
    dy59<- y[59]*(params[1] * St[59] - params[2]/TC - params[3]*ConfM[59,]%*%y[1:N] - params[4]*y[59])
    dy60<- y[60]*(params[1] * St[60] - params[2]/TC - params[3]*ConfM[60,]%*%y[1:N] - params[4]*y[60])
    dy61<- y[61]*(params[1] * St[61] - params[2]/TC - params[3]*ConfM[61,]%*%y[1:N] - params[4]*y[61])
    dy62<- y[62]*(params[1] * St[62] - params[2]/TC - params[3]*ConfM[62,]%*%y[1:N] - params[4]*y[62])
    dy63<- y[63]*(params[1] * St[63] - params[2]/TC - params[3]*ConfM[63,]%*%y[1:N] - params[4]*y[63])
    dy64<- y[64]*(params[1] * St[64] - params[2]/TC - params[3]*ConfM[64,]%*%y[1:N] - params[4]*y[64])
    dy65<- y[65]*(params[1] * St[65] - params[2]/TC - params[3]*ConfM[65,]%*%y[1:N] - params[4]*y[65])
    dy66<- y[66]*(params[1] * St[66] - params[2]/TC - params[3]*ConfM[66,]%*%y[1:N] - params[4]*y[66])
    dy67<- y[67]*(params[1] * St[67] - params[2]/TC - params[3]*ConfM[67,]%*%y[1:N] - params[4]*y[67])
    dy68<- y[68]*(params[1] * St[68] - params[2]/TC - params[3]*ConfM[68,]%*%y[1:N] - params[4]*y[68])
    dy69<- y[69]*(params[1] * St[69] - params[2]/TC - params[3]*ConfM[69,]%*%y[1:N] - params[4]*y[69])
    dy70<- y[70]*(params[1] * St[70] - params[2]/TC - params[3]*ConfM[70,]%*%y[1:N] - params[4]*y[70])
    dy71<- y[71]*(params[1] * St[71] - params[2]/TC - params[3]*ConfM[71,]%*%y[1:N] - params[4]*y[71])
    dy72<- y[72]*(params[1] * St[72] - params[2]/TC - params[3]*ConfM[72,]%*%y[1:N] - params[4]*y[72])
    dy73<- y[73]*(params[1] * St[73] - params[2]/TC - params[3]*ConfM[73,]%*%y[1:N] - params[4]*y[73])
    dy74<- y[74]*(params[1] * St[74] - params[2]/TC - params[3]*ConfM[74,]%*%y[1:N] - params[4]*y[74])
    dy75<- y[75]*(params[1] * St[75] - params[2]/TC - params[3]*ConfM[75,]%*%y[1:N] - params[4]*y[75])
    dy76<- y[76]*(params[1] * St[76] - params[2]/TC - params[3]*ConfM[76,]%*%y[1:N] - params[4]*y[76])
    dy77<- y[77]*(params[1] * St[77] - params[2]/TC - params[3]*ConfM[77,]%*%y[1:N] - params[4]*y[77])
    dy78<- y[78]*(params[1] * St[78] - params[2]/TC - params[3]*ConfM[78,]%*%y[1:N] - params[4]*y[78])
    dy79<- y[79]*(params[1] * St[79] - params[2]/TC - params[3]*ConfM[79,]%*%y[1:N] - params[4]*y[79])
    dy80<- y[80]*(params[1] * St[80] - params[2]/TC - params[3]*ConfM[80,]%*%y[1:N] - params[4]*y[80])
    dy81<- y[81]*(params[1] * St[81] - params[2]/TC - params[3]*ConfM[81,]%*%y[1:N] - params[4]*y[81])
    dy82<- y[82]*(params[1] * St[82] - params[2]/TC - params[3]*ConfM[82,]%*%y[1:N] - params[4]*y[82])
    dy83<- y[83]*(params[1] * St[83] - params[2]/TC - params[3]*ConfM[83,]%*%y[1:N] - params[4]*y[83])
    dy84<- y[84]*(params[1] * St[84] - params[2]/TC - params[3]*ConfM[84,]%*%y[1:N] - params[4]*y[84])
    dy85<- y[85]*(params[1] * St[85] - params[2]/TC - params[3]*ConfM[85,]%*%y[1:N] - params[4]*y[85])
    dy86<- y[86]*(params[1] * St[86] - params[2]/TC - params[3]*ConfM[86,]%*%y[1:N] - params[4]*y[86])
    dy87<- y[87]*(params[1] * St[87] - params[2]/TC - params[3]*ConfM[87,]%*%y[1:N] - params[4]*y[87])
    dy88<- y[88]*(params[1] * St[88] - params[2]/TC - params[3]*ConfM[88,]%*%y[1:N] - params[4]*y[88])
    dy89<- y[89]*(params[1] * St[89] - params[2]/TC - params[3]*ConfM[89,]%*%y[1:N] - params[4]*y[89])
    dy90<- y[90]*(params[1] * St[90] - params[2]/TC - params[3]*ConfM[90,]%*%y[1:N] - params[4]*y[90])
    dy91<- y[91]*(params[1] * St[91] - params[2]/TC - params[3]*ConfM[91,]%*%y[1:N] - params[4]*y[91])
    dy92<- y[92]*(params[1] * St[92] - params[2]/TC - params[3]*ConfM[92,]%*%y[1:N] - params[4]*y[92])
    dy93<- y[93]*(params[1] * St[93] - params[2]/TC - params[3]*ConfM[93,]%*%y[1:N] - params[4]*y[93])
    dy94<- y[94]*(params[1] * St[94] - params[2]/TC - params[3]*ConfM[94,]%*%y[1:N] - params[4]*y[94])
    dy95<- y[95]*(params[1] * St[95] - params[2]/TC - params[3]*ConfM[95,]%*%y[1:N] - params[4]*y[95])
    dy96<- y[96]*(params[1] * St[96] - params[2]/TC - params[3]*ConfM[96,]%*%y[1:N] - params[4]*y[96])
    dy97<- y[97]*(params[1] * St[97] - params[2]/TC - params[3]*ConfM[97,]%*%y[1:N] - params[4]*y[97])
    dy98<- y[98]*(params[1] * St[98] - params[2]/TC - params[3]*ConfM[98,]%*%y[1:N] - params[4]*y[98])
    dy99<- y[99]*(params[1] * St[99] - params[2]/TC - params[3]*ConfM[99,]%*%y[1:N] - params[4]*y[99])
    dy100<- y[100]*(params[1] * St[100] - params[2]/TC - params[3]*ConfM[100,]%*%y[1:N] - params[4]*y[100])
    dy101<- y[101]*(params[1] * St[101] - params[2]/TC - params[3]*ConfM[101,]%*%y[1:N] - params[4]*y[101])
    dy102<- y[102]*(params[1] * St[102] - params[2]/TC - params[3]*ConfM[102,]%*%y[1:N] - params[4]*y[102])
    dy103<- y[103]*(params[1] * St[103] - params[2]/TC - params[3]*ConfM[103,]%*%y[1:N] - params[4]*y[103])
    dy104<- y[104]*(params[1] * St[104] - params[2]/TC - params[3]*ConfM[104,]%*%y[1:N] - params[4]*y[104])
    dy105<- y[105]*(params[1] * St[105] - params[2]/TC - params[3]*ConfM[105,]%*%y[1:N] - params[4]*y[105])
    dy106<- y[106]*(params[1] * St[106] - params[2]/TC - params[3]*ConfM[106,]%*%y[1:N] - params[4]*y[106])
    dy107<- y[107]*(params[1] * St[107] - params[2]/TC - params[3]*ConfM[107,]%*%y[1:N] - params[4]*y[107])
    dy108<- y[108]*(params[1] * St[108] - params[2]/TC - params[3]*ConfM[108,]%*%y[1:N] - params[4]*y[108])
    dy109<- y[109]*(params[1] * St[109] - params[2]/TC - params[3]*ConfM[109,]%*%y[1:N] - params[4]*y[109])
    dy110<- y[110]*(params[1] * St[110] - params[2]/TC - params[3]*ConfM[110,]%*%y[1:N] - params[4]*y[110])
    dy111<- y[111]*(params[1] * St[111] - params[2]/TC - params[3]*ConfM[111,]%*%y[1:N] - params[4]*y[111])
    dy112<- y[112]*(params[1] * St[112] - params[2]/TC - params[3]*ConfM[112,]%*%y[1:N] - params[4]*y[112])
    dy113<- y[113]*(params[1] * St[113] - params[2]/TC - params[3]*ConfM[113,]%*%y[1:N] - params[4]*y[113])
    dy114<- y[114]*(params[1] * St[114] - params[2]/TC - params[3]*ConfM[114,]%*%y[1:N] - params[4]*y[114])
    dy115<- y[115]*(params[1] * St[115] - params[2]/TC - params[3]*ConfM[115,]%*%y[1:N] - params[4]*y[115])
    dy116<- y[116]*(params[1] * St[116] - params[2]/TC - params[3]*ConfM[116,]%*%y[1:N] - params[4]*y[116])
    dy117<- y[117]*(params[1] * St[117] - params[2]/TC - params[3]*ConfM[117,]%*%y[1:N] - params[4]*y[117])
    dy118<- y[118]*(params[1] * St[118] - params[2]/TC - params[3]*ConfM[118,]%*%y[1:N] - params[4]*y[118])
    dy119<- y[119]*(params[1] * St[119] - params[2]/TC - params[3]*ConfM[119,]%*%y[1:N] - params[4]*y[119])
    dy120<- y[120]*(params[1] * St[120] - params[2]/TC - params[3]*ConfM[120,]%*%y[1:N] - params[4]*y[120])
    dy121<- y[121]*(params[1] * St[121] - params[2]/TC - params[3]*ConfM[121,]%*%y[1:N] - params[4]*y[121])
    dy122<- y[122]*(params[1] * St[122] - params[2]/TC - params[3]*ConfM[122,]%*%y[1:N] - params[4]*y[122])
    dy123<- y[123]*(params[1] * St[123] - params[2]/TC - params[3]*ConfM[123,]%*%y[1:N] - params[4]*y[123])
    dy124<- y[124]*(params[1] * St[124] - params[2]/TC - params[3]*ConfM[124,]%*%y[1:N] - params[4]*y[124])
    dy125<- y[125]*(params[1] * St[125] - params[2]/TC - params[3]*ConfM[125,]%*%y[1:N] - params[4]*y[125])
    dy126<- y[126]*(params[1] * St[126] - params[2]/TC - params[3]*ConfM[126,]%*%y[1:N] - params[4]*y[126])
    dy127<- y[127]*(params[1] * St[127] - params[2]/TC - params[3]*ConfM[127,]%*%y[1:N] - params[4]*y[127])
    dy128<- y[128]*(params[1] * St[128] - params[2]/TC - params[3]*ConfM[128,]%*%y[1:N] - params[4]*y[128])
    dy129<- y[129]*(params[1] * St[129] - params[2]/TC - params[3]*ConfM[129,]%*%y[1:N] - params[4]*y[129])
    dy130<- y[130]*(params[1] * St[130] - params[2]/TC - params[3]*ConfM[130,]%*%y[1:N] - params[4]*y[130])
    dy131<- y[131]*(params[1] * St[131] - params[2]/TC - params[3]*ConfM[131,]%*%y[1:N] - params[4]*y[131])
    dy132<- y[132]*(params[1] * St[132] - params[2]/TC - params[3]*ConfM[132,]%*%y[1:N] - params[4]*y[132])
    dy133<- y[133]*(params[1] * St[133] - params[2]/TC - params[3]*ConfM[133,]%*%y[1:N] - params[4]*y[133])
    dy134<- y[134]*(params[1] * St[134] - params[2]/TC - params[3]*ConfM[134,]%*%y[1:N] - params[4]*y[134])
    dy135<- y[135]*(params[1] * St[135] - params[2]/TC - params[3]*ConfM[135,]%*%y[1:N] - params[4]*y[135])
    dy136<- y[136]*(params[1] * St[136] - params[2]/TC - params[3]*ConfM[136,]%*%y[1:N] - params[4]*y[136])
    dy137<- y[137]*(params[1] * St[137] - params[2]/TC - params[3]*ConfM[137,]%*%y[1:N] - params[4]*y[137])
    dy138<- y[138]*(params[1] * St[138] - params[2]/TC - params[3]*ConfM[138,]%*%y[1:N] - params[4]*y[138])
    dy139<- y[139]*(params[1] * St[139] - params[2]/TC - params[3]*ConfM[139,]%*%y[1:N] - params[4]*y[139])
    dy140<- y[140]*(params[1] * St[140] - params[2]/TC - params[3]*ConfM[140,]%*%y[1:N] - params[4]*y[140])
    dy141<- y[141]*(params[1] * St[141] - params[2]/TC - params[3]*ConfM[141,]%*%y[1:N] - params[4]*y[141])
    dy142<- y[142]*(params[1] * St[142] - params[2]/TC - params[3]*ConfM[142,]%*%y[1:N] - params[4]*y[142])
    dy143<- y[143]*(params[1] * St[143] - params[2]/TC - params[3]*ConfM[143,]%*%y[1:N] - params[4]*y[143])
    dy144<- y[144]*(params[1] * St[144] - params[2]/TC - params[3]*ConfM[144,]%*%y[1:N] - params[4]*y[144])
    dy145<- y[145]*(params[1] * St[145] - params[2]/TC - params[3]*ConfM[145,]%*%y[1:N] - params[4]*y[145])
    dy146<- y[146]*(params[1] * St[146] - params[2]/TC - params[3]*ConfM[146,]%*%y[1:N] - params[4]*y[146])
    dy147<- y[147]*(params[1] * St[147] - params[2]/TC - params[3]*ConfM[147,]%*%y[1:N] - params[4]*y[147])
    dy148<- y[148]*(params[1] * St[148] - params[2]/TC - params[3]*ConfM[148,]%*%y[1:N] - params[4]*y[148])
    dy149<- y[149]*(params[1] * St[149] - params[2]/TC - params[3]*ConfM[149,]%*%y[1:N] - params[4]*y[149])
    dy150<- y[150]*(params[1] * St[150] - params[2]/TC - params[3]*ConfM[150,]%*%y[1:N] - params[4]*y[150])
    dD <- params[3]*y[1:N] %*% ConfM %*% y[1:N]### deaths. Not by 2 since it is already inside!
    dR <- params[1] * St %*%y[1:N]
    list(c(dy1,
           dy2,
           dy3,
           dy4,
           dy5,
           dy6,
           dy7,
           dy8,
           dy9,
           dy10,
           dy11,
           dy12,
           dy13,
           dy14,
           dy15,
           dy16,
           dy17,
           dy18,
           dy19,
           dy20,
           dy21,
           dy22,
           dy23,
           dy24,
           dy25,
           dy26,
           dy27,
           dy28,
           dy29,
           dy30,
           dy31,
           dy32,
           dy33,
           dy34,
           dy35,
           dy36,
           dy37,
           dy38,
           dy39,
           dy40,
           dy41,
           dy42,
           dy43,
           dy44,
           dy45,
           dy46,
           dy47,
           dy48,
           dy49,
           dy50,
           dy51,
           dy52,
           dy53,
           dy54,
           dy55,
           dy56,
           dy57,
           dy58,
           dy59,
           dy60,
           dy61,
           dy62,
           dy63,
           dy64,
           dy65,
           dy66,
           dy67,
           dy68,
           dy69,
           dy70,
           dy71,
           dy72,
           dy73,
           dy74,
           dy75,
           dy76,
           dy77,
           dy78,
           dy79,
           dy80,
           dy81,
           dy82,
           dy83,
           dy84,
           dy85,
           dy86,
           dy87,
           dy88,
           dy89,
           dy90,
           dy91,
           dy92,
           dy93,
           dy94,
           dy95,
           dy96,
           dy97,
           dy98,
           dy99,
           dy100,
           dy101,
           dy102,
           dy103,
           dy104,
           dy105,
           dy106,
           dy107,
           dy108,
           dy109,
           dy110,
           dy111,
           dy112,
           dy113,
           dy114,
           dy115,
           dy116,
           dy117,
           dy118,
           dy119,
           dy120,
           dy121,
           dy122,
           dy123,
           dy124,
           dy125,
           dy126,
           dy127,
           dy128,
           dy129,
           dy130,
           dy131,
           dy132,
           dy133,
           dy134,
           dy135,
           dy136,
           dy137,
           dy138,
           dy139,
           dy140,
           dy141,
           dy142,
           dy143,
           dy144,
           dy145,
           dy146,
           dy147,
           dy148,
           dy149,
           dy150,
           dD,
           dR))
  }
  Years <- 10
  yini <- c(Nod$EstCrime, 
            D = 0, 
            R = 0)
  params <- c(rho = rho, #### recruitment
              eta = eta,    #### incapacitation
              theta = theta, ## conflict
              omega = omega) ## fragmentation
  times <- seq(from = 0, to = Years*52, by = 1)
  out <- ode(times = times,
             y = yini, 
             func = DifEq, 
             ConfM = ConfM,
             St = 0*y0+1,
             parms = params)
  InstArrest <- params[2]*apply(out[, 2:(NCartels +1)], 1, sum)/apply(out[, 2:(NCartels +1)], 1, sum)
  InstDeath <- 0 * InstArrest
  InstRecruit <- 0 * InstArrest
  InstFragm <- 0 * InstArrest
  CSize <- 0 * InstArrest
  for(k in 1:length(InstDeath)){
    InstRecruit[k] <- params[1]*(0*Nod$strengthAllies + 1) %*% out[k, 2:(NCartels +1)]
    InstDeath[k] <- params[3]*out[k, 2:(NCartels +1)] %*% ConfM %*% out[k, 2:(NCartels +1)]
    InstFragm[k] <- params[4]*sum(out[k, 2:(NCartels +1)]^2)
    CSize[k] <- sum(out[k,  2:(NCartels +1)])
    
  } 
  OutDF <- as.data.frame(out)
  OutDF$InstArrest <- InstArrest
  OutDF$InsDeath <- InstDeath
  OutDF$InstRecruit <- InstRecruit
  OutDF$InstFragm <- InstFragm
}

#### reduce recruitment to zero
{
  params2 <- params
  params2[1] <- 0
  out2 <-   ode(times = times,
                y = out[521, 2:153], 
                func = DifEq, 
                ConfM = ConfM,
                St = 0*y0+1,
                parms = params2)
  InstDeath2 <- 0 * InstArrest
  CSize2 <- 0*InstArrest
  for(k in 1:length(InstDeath2)){
    InstDeath2[k] <- params[3]*out2[k, 2:(NCartels +1)] %*% ConfM %*% out2[k, 2:(NCartels +1)]
    CSize2[k] <- sum(out2[k,  2:(NCartels +1)])
  } 
  
  params2 <- params
  params2[2] <- params2[2]*2
  out3 <-   ode(times = times,
                y = out[521, 2:153], 
                func = DifEq, 
                ConfM = ConfM,
                St = 0*y0+1,
                parms = params2)
  InstDeath3 <- 0 * InstArrest
  CSize3 <- 0*InstArrest
  for(k in 1:length(InstDeath2)){
    InstDeath3[k] <- params[3]*out3[k, 2:(NCartels +1)] %*% ConfM %*% out3[k, 2:(NCartels +1)]
    CSize3[k] <- sum(out3[k,  2:(NCartels +1)])
  } 
  
  params2 <- params
  out4 <-   ode(times = times,
                y = out[521, 2:153], 
                func = DifEq, 
                ConfM = ConfM,
                St = 0*y0+1,
                parms = params2)
  InstDeath4 <- 0 * InstArrest
  CSize4 <- 0*InstArrest
  for(k in 1:length(InstDeath2)){
    InstDeath4[k] <- params[3]*out4[k, 2:(NCartels +1)] %*% ConfM %*% out4[k, 2:(NCartels +1)]
    CSize4[k] <- sum(out4[k,  2:(NCartels +1)])
  } 
  
  params2 <- params
  params2[1] <- params2[1]*0.5
  out5 <-   ode(times = times,
                y = out[521, 2:153], 
                func = DifEq, 
                ConfM = ConfM,
                St = 0*y0+1,
                parms = params2)
  InstDeath5 <- 0 * InstArrest
  CSize5 <- 0*InstArrest
  for(k in 1:length(InstDeath2)){
    InstDeath5[k] <- params[3]*out5[k, 2:(NCartels +1)] %*% ConfM %*% out5[k, 2:(NCartels +1)]
    CSize5[k] <- sum(out5[k,  2:(NCartels +1)])
  } 
  
  
  plot(c(InstDeath, InstDeath4), type = "l",
       ylim = c(35, 210),
       xlim = c(0, 1000)) 
  points(c(InstDeath, InstDeath3), type = "l")
  points(c(InstDeath, InstDeath5), type = "l")
  points(c(InstDeath, InstDeath2), type = "l")

  for(k in 1:24){points(c(1,520*2), 10*c(k,k), type = "l", lty = 2)
    points(c(k,k)*52, c(1,200), type = "l", lty = 2, col = 2)
  }
  
  
  plot(c(CSize, CSize4), type = "l",
       #ylim = c(35, 210),
       xlim = c(0, 1000)) 
  points(c(CSize, CSize3), type = "l")
  points(c(CSize, CSize5), type = "l")
  points(c(CSize, CSize2), type = "l")
  
  for(k in 1:24){points(c(1,520*2), 100000+10000*c(k,k), type = "l", lty = 2)
    points(c(k,k)*52, c(100000,350000), type = "l", lty = 2, col = 2)
  }
  
}

#### error function
{
  ErrF <- function(x, 
                   ConfM, 
                   WeeklyCasualties, 
                   WeeklyArrest,
                   DifEq){
    rho = x[1]/1000; 
    eta = x[2]*1; 
    theta = x[3]/100000000; 
    omega = x[4]/1000000; 
    Crime0 = x[5]*1000
    NCartels <- 150
    de <- (1:150)^(-1)
    y0 <- de/sum(de)*Crime0
    yini <- c(y0, 
              D = 0, 
              R = 0)
    params <- c(rho = rho, #### recruitment
                eta = eta,    #### incapacitation
                theta = theta, ## conflict
                omega = omega) ## fragmentation
    times <- seq(from = 1, to = 520, by = 1)
    out <- ode(times = times,
               y = yini, 
               func = DifEq, 
               ConfM = ConfM,
               method = "euler",
               St = rep(1, 150),
               parms = params)
    InstDeath <- rep(0, dim(out)[1])
    InstArrest <- params[2]*apply(out[, 2:(NCartels +1)], 1, sum)/apply(out[, 2:(NCartels +1)], 1, sum)
    for(k in 1:length(InstDeath)){
      InstDeath[k] <- params[3]*out[k, 2:(NCartels +1)] %*% ConfM %*% out[k, 2:(NCartels +1)]
    } 
    return(sum((InstDeath-WeeklyCasualties)^2) + 0*sum((InstArrest-WeeklyArrest)^2))
  }
}

##call error fn
{
  ErrF(x = c(rho = 11.23, #e-3
             eta = 211,
             theta = 5.9,
             omega = 0.05,
             Crime0 = 160),
       ConfM = ConfM, 
       WeeklyCasualties = WeeklyCasualties*0.1,
       WeeklyArrest = WeeklyArrest*0.05,
       DifEq = DifEq)
}

#### optimise error function
{
  OpPar <- optim(par = c(rho = 2.68, #e-3
                         eta = 117.6,
                         theta = 1.59, #e-8
                         omega = .053, #e-6
                         Crime0 = 111.0), #e3
                 lower = c(2,110, 1, 0.04, 110),
                 upper = c(3, 120, 2, .1, 115),
                 method = "L-BFGS-B",
                 fn = ErrF,
                 ConfM = ConfM, 
                 WeeklyCasualties = WeeklyCasualties*0.1,
                 WeeklyArrest = WeeklyArrest*0.05,
                 DifEq = DifEq)
  OpParStar <- OpPar
  save(OpParStar, file = "OpParStar.RData")
}
