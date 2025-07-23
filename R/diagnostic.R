############################################################################
# DIAGNOSTIC FUNCTIONS
############################################################################
#' Takes output from deSolve and pulls out hospitalisations aggregated in the
#' age groups used for fitting for term first exposure, term subsequent exposures,
#' term combined, preterm first exposure, preterm subsequent exposures,
#' preterm combined, combined first exposure, combined second exposure, combined
#'
#' @param deSolveOut output matrix
#' @param times time vector
#' @param nAges number of age groups
#' @param numImm number of immunisation categories
#' @param old 1 = include older age group 24-<60 months,
#'            2 = include 26-<36 months, 36->48 months, 48-<60 months
#' @return list
#'
aggregate_output_immiPT <- function(deSolveOut, times, nAges, nState = 48, nRet = 55, old = 1)
{
  deSolveOut <- array(deSolveOut[,-1], dim = c(length(times), nAges, nRet))
  deSolveOut <- deSolveOut[which(times == round(times)),,]

  listNames <- c("deSolve", "hos_term", "hos_preterm", "hosp_comb",
                 "doses_matvacc", #total mat vacc doses
                 "doses_mAb_0", #total mAb doses at birth
                 "doses_mAb_c", #total mAb doses catchup
                 "doses_mAb_2") #total mAb doses 2nd season
                            #"doses_term", "doses_preterm", "doses_comb")
  modelOut <- sapply(listNames,function(x) NULL)

  colNames <- c("<3months", "3-<6months", "6-<12months", "12-<24months")
  if(old == 1) colNames <- c(colNames, "24-<60months")
  if(old == 2) colNames <- c(colNames, "24-<36months", "36-<48months", "48-<60months")

  for(i in (nState+1):nRet)
    deSolveOut[, , i] <- apply(as.data.frame(deSolveOut[, , i]), 2, function(x) diff(c(0, x)))
  modelOut[[listNames[1]]] <- deSolveOut

  for(i in 1:length(listNames[-1]))
  {
    detInc <- as.data.frame(deSolveOut[, , nState + i])
    DetInc_1_3 <- apply(detInc[, 1:3], 1, sum)
    DetInc_4_6 <- apply(detInc[, 4:6], 1, sum)
    DetInc_7_12 <- apply(detInc[, 7:12], 1, sum)
    DetInc_13_24 <- apply(detInc[, 13:24], 1, sum)
    if(old == 1) DetInc_25_60 <- apply(detInc[, 25:60], 1, sum)
    if(old == 2){
      DetInc_25_36 <- apply(detInc[, 25:36], 1, sum)
      DetInc_37_48 <- apply(detInc[, 37:48], 1, sum)
      DetInc_49_60 <- apply(detInc[, 49:60], 1, sum)
    }

    detIC <- cbind(DetInc_1_3, DetInc_4_6, DetInc_7_12, DetInc_13_24)
    if(old == 1) detIC <- cbind(detIC, DetInc_25_60)
    if(old == 2) detIC <- cbind(detIC, DetInc_25_36, DetInc_37_48, DetInc_49_60)
    colnames(detIC) <- colNames

    modelOut[[listNames[i+1]]] <- detIC
  }

  return(modelOut)
}
###############################################################################################
#' Plot age-to-risk function
#'
#' @param A the average maximum risk increase for a term infant (at age 0)
#' @param B the exponential function decay constant
#' @param C the average minimum risk over all ages for all infants
#' @param D the scaling factor (<1) modifying risk for the prior infection group
#' @param E the scaling factor (>1) modifying risk for those born preterm
#' @param maxAge maximum age in months to plot
#' @return ggplot
#'
plot_ageToRisk <- function(A, B, C, D, E, maxAge = 12)
{
  age <- seq(0, maxAge, by = 0.25) #months
  Y_T0 <- data.frame(Age = age, Risk = A * exp(-B * age) + C)
  Y_T1 <- data.frame(Age = age, Risk = (A * exp(-B * age) + C) * D)
  Y_PT0 <- data.frame(Age = age, Risk = (A * exp(-B * age) + C) * E)
  Y_PT1 <- data.frame(Age = age, Risk = (A * exp(-B * age) + C) * D * E)

  if(any(Y_T0$Risk > 1)) Y_T0$Risk[which(Y_T0$Risk > 1)] <- 1
  if(any(Y_T1$Risk > 1)) Y_T1$Risk[which(Y_T1$Risk > 1)] <- 1
  if(any(Y_PT0$Risk > 1)) Y_PT0$Risk[which(Y_PT0$Risk > 1)] <- 1
  if(any(Y_PT1$Risk > 1)) Y_PT1$Risk[which(Y_PT1$Risk > 1)] <- 1

  risk_df <- cbind(Y_T0, Y_T1$Risk, Y_PT0$Risk, Y_PT1$Risk)
  colnames(risk_df) <- c("Age", "T0", "T1", "PT0", "PT1")

  maxLimY <- plyr::round_any(((A + C)*E*1.05), 0.1, f = ceiling)
  if(maxLimY > 1) maxLimY <- 1.01

  g1 <- ggplot2::ggplot() +
    ggplot2::geom_line(data = Y_PT0, ggplot2::aes(x = Age, y = Risk, col = "red"), size = 1.1) +
    ggplot2::geom_line(data = Y_PT1, ggplot2::aes(x = Age, y = Risk, col = "red"), size = 1.1, linetype = "dashed") +
    ggplot2::geom_line(data = Y_T0, ggplot2::aes(x = Age, y = Risk, col = "blue"), size = 1.1) +
    ggplot2::geom_line(data = Y_T1, ggplot2::aes(x = Age, y = Risk, col = "blue"), size = 1.1, linetype = "dashed") +
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="none")+
    ggplot2::xlab("Age (months)")+
    ggplot2::ylab("Risk of hospitalisation")+
    ggplot2::ggtitle(label =   "Term (red), Preterm (blue), First exposure (line), Subsequent (dashed)",
            subtitle =  paste("A = ", A, " B = ", B, " C = ", C, " D = ", D, " E = ", E))+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12))+
    ggplot2::scale_x_continuous(limits = c(0, maxAge + 0.1), expand = c(0, 0), breaks = 0:maxAge)+
    ggplot2::scale_y_continuous(limits = c(0, maxLimY), expand = c(0,0),
                       breaks = seq(0, maxLimY, 0.01), labels = c(0, seq(0.01, maxLimY, 0.01)))


  return(list(g1, risk_df))
}

###############################################################################################
#' Plot fit to data for risk model (with preterms)
#'
#' @param term table of MODEL predicted monthly hospitalisations
#'             for term infants
#' @param pTerm table of MODEL predicted monthly hospitalisations
#'              for preterm infants
#' @param obs table of OBSERVED term and preterm population monthly
#'            hospitalisation data
#' @param years year labels for x-axis
#' @param lineCol default line colour is black
#' @param lineType default line type is dashed
#' @param lineSize default line size is 1
#' @param plotText if plotText is NA (as default), plotText is "RSV
#'                 hospitalisations observed (dots) and model (dashed)"
#' @param fileName if fileName is NA (as default), resulting plots are not
#'                  saved. Otherwise fileName is base for 3 resulting jpeg files
#' @return ggplot
#'
plot_riskFit <- function(term, pTerm, obs, years = 2010:2020, 
                         lineCol = "black", lineType = "dashed", lineSize = 1,
                         plotText = NA, fileName = NA)
{
  if(is.na(plotText)) plotText <- "RSV hospitalisations observed (dots) and model (dashed)"

  colNames <- c("<3months","3-<6months", "6-<12months", "12-<24months", "24-<60months")
  colnames(term) <- colNames
  colnames(pTerm) <- colNames

  melt_term <- reshape2::melt(as.matrix(term))
  colnames(melt_term) <- c("Month", "Age", "Hospitalisations")

  obs_term <- obs[, seq(1, ncol(obs), 2)]
  colnames(obs_term) <- colNames
  melt_obs_term <- reshape2::melt(as.matrix(obs_term))
  colnames(melt_obs_term) <- c("Month", "Age", "Hospitalisations")

  melt_pTerm <- reshape2::melt(as.matrix(pTerm))
  colnames(melt_pTerm) <- c("Month", "Age", "Hospitalisations")

  obs_pTerm <- obs[, seq(2, ncol(obs), 2)]
  colnames(obs_pTerm) <- colNames
  melt_obs_pTerm <- reshape2::melt(as.matrix(obs_pTerm))
  colnames(melt_obs_pTerm) <- c("Month", "Age", "Hospitalisations")

  g1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = melt_obs_term, ggplot2::aes(x = Month, y = Hospitalisations)) +
    ggplot2::geom_point(data = melt_obs_pTerm, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
    ggplot2::geom_line(data = melt_pTerm, ggplot2::aes(x = Month, y = Hospitalisations, col = Age), linetype = lineType) +
    ggplot2::geom_line(data = melt_term, ggplot2::aes(x = Month, y = Hospitalisations), col = lineCol, linetype = lineType) +
    ggplot2::scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="none")+
    ggplot2::ggtitle(label = plotText, subtitle = "Term (black) and preterm (colour)")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12))+
    ggplot2::facet_wrap(~Age, ncol = 1)

  g2 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = melt_obs_term, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
    ggplot2::geom_line(data = melt_term, ggplot2::aes(x = Month, y = Hospitalisations, col = Age), linetype = lineType) +
    ggplot2::scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="none")+
    ggplot2::ggtitle(label = plotText, subtitle = "Term")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12))+
    ggplot2::facet_wrap(~Age, ncol = 1)

  g3 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = melt_obs_pTerm, ggplot2::aes(x = Month, y = Hospitalisations, col = Age)) +
    ggplot2::geom_line(data = melt_pTerm, ggplot2::aes(x = Month, y = Hospitalisations, col = Age),  linetype = lineType) +
    ggplot2::scale_x_continuous(name = "Year", breaks = seq(0, nrow(obs), 12), labels = years)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="none")+
    ggplot2::ggtitle(label = plotText, subtitle =  "Preterm")+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12))+
    ggplot2::facet_wrap(~Age, ncol = 1)

  if(!is.na(fileName))
  {
    print(g1)
    ggplot2::ggsave(filename = paste(fileName, ".jpeg", sep = ''), dpi = 900)
    ggplot2::ggsave(filename = paste(fileName, ".pdf", sep = ''), dpi = 900)
    print(g2)
    ggplot2::ggsave(filename = paste(fileName, "_term.pdf", sep = ''), dpi = 900)
    print(g3)
    ggplot2::ggsave(filename = paste(fileName, "_preterm.pdf", sep = ''), dpi = 900)
  }

  return(list(g1, g2, g3))
}

