library(deSolve)
library(conmat)
library(dplyr)
library(pracma)
library(ggplot2)
library(reshape2)

load("data/data.southernWA.2016.rda")
load("data/data.contactDaily.ABSMicro.rda")

source("R/data.r")
source("R/functions.r")
source("R/parameters.r")
source("R/models.r")
source("R/immunisation.r")
source("R/diagnostic.r")

institute_cols <- c("#F1B434", "#F56B00", "#00A39C", "#00807A",
                    "#4A99DE", "#426EA8", "#1F3B73", "#565F5F", "#111921")
#Saffron, Pumpkin, Teal, Dark Teal, Celestial Blue, Azure Blue, 
#Midnight Blue, Cool Grey, Black
##############################################################################
#SET UP MODEL
age <- age_structure(year = "2016", area = 1)
mixing <- data.contactDaily.ABSMicro * 365/12

#create vectors
#reduced infectiousness - currently assuming no reduction
omega_vect <- as.vector(rep(1, age$nAges))
omega_vect[!(age$age_years < 10)] <- omega
#reduced susceptibility due to natural maternal immunity
sigma_vect <- as.vector(c(sigma, rep(1, (age$nAges-12))))

#Scenario parameters
source("R/scenarios.r")

#Initialise y - state variables 
#Model with no immunisation
y0_seir <- matrix(0, age$nAges, 16)
y0_seir[,1] <- (1 - alpha) * 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- (1 - alpha) * 0.01 * age$pop_groups #I0 term
y0_seir[,9] <- alpha * 0.99 * age$pop_groups #S0_bar preterm
y0_seir[,11] <- alpha * 0.01 * age$pop_groups #I0_bar preterm
y0_noImm <- cbind(y0_seir, matrix(0, age$nAges, 14))
#Model with maternal and mAb
y0_seir <- matrix(0, age$nAges, 48)
y0_seir[,1] <- (1 - alpha) * 0.99 * age$pop_groups #S0 term
y0_seir[,3] <- (1 - alpha) * 0.01 * age$pop_groups #I0 term
y0_seir[,9] <- alpha * 0.99 * age$pop_groups #S0_bar preterm
y0_seir[,11] <- alpha * 0.01 * age$pop_groups #I0_bar preterm
y0 <- cbind(y0_seir, matrix(0, age$nAges, 7))

##############################################################################
#MODEL FIT

#Plot age-to-risk function using fixed and fitted values
 plot_ageToRisk(A = A_fit, B = B_fit, C = C_fit, D = D_fit,
                E = E_fit, maxAge = 24)

##############################################################################
#SCENARIO 0 - Base, no immunisation
scenario0_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario0)

scenario0_out <- aggregate_output_immiPT (deSolveOut = scenario0_raw,
                                                 times = times,
                                                 nAges = age$nAges)

##################################################################################
#SCENARIO 1 - Monoclonal, year-round at birth, 50% - 70% - 90% coverage 
scenario1_50_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario1_50)

scenario1_50_out <- aggregate_output_immiPT (deSolveOut = scenario1_50_raw,
                                          times = times,
                                          nAges = age$nAges)


scenario1_70_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario1_70)

scenario1_70_out <- aggregate_output_immiPT (deSolveOut = scenario1_70_raw,
                                             times = times,
                                             nAges = age$nAges)

scenario1_90_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario1_90)

scenario1_90_out <- aggregate_output_immiPT (deSolveOut = scenario1_90_raw,
                                             times = times,
                                             nAges = age$nAges)
##################################################################################
#SCENARIO 2 - Monoclonal, Seasonal at birth (April-Sept), 50% - 70% - 90% coverage 
scenario2_50_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario2_50)

scenario2_50_out <- aggregate_output_immiPT (deSolveOut = scenario2_50_raw,
                                             times = times,
                                             nAges = age$nAges)


scenario2_70_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario2_70)

scenario2_70_out <- aggregate_output_immiPT (deSolveOut = scenario2_70_raw,
                                             times = times,
                                             nAges = age$nAges)

scenario2_90_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario2_90)

scenario2_90_out <- aggregate_output_immiPT (deSolveOut = scenario2_90_raw,
                                             times = times,
                                             nAges = age$nAges)
##################################################################################
#SCENARIO 3 - Monoclonal, Seasonal with catchup + extra dose at-risk (WA 2024 equivalent), 
#             Birth 90% Catchup 70% 2nd season 30% coverage 
scenario3_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario3)

scenario3_out <- aggregate_output_immiPT (deSolveOut = scenario3_raw,
                                             times = times,
                                             nAges = age$nAges)

##################################################################################
#SCENARIO 4 - Maternal Vacc, year-round, 50% - 70% - 90% coverage 
scenario4_50_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario4_50)

scenario4_50_out <- aggregate_output_immiPT (deSolveOut = scenario4_50_raw,
                                             times = times,
                                             nAges = age$nAges)


scenario4_70_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario4_70)

scenario4_70_out <- aggregate_output_immiPT (deSolveOut = scenario4_70_raw,
                                             times = times,
                                             nAges = age$nAges)

scenario4_90_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario4_90)

scenario4_90_out <- aggregate_output_immiPT (deSolveOut = scenario4_90_raw,
                                             times = times,
                                             nAges = age$nAges)

##################################################################################
#SCENARIO 5 - Maternal Vacc, seasonal, 50% - 70% - 90% coverage 
scenario5_50_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario5_50)

scenario5_50_out <- aggregate_output_immiPT (deSolveOut = scenario5_50_raw,
                                             times = times,
                                             nAges = age$nAges)


scenario5_70_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario5_70)

scenario5_70_out <- aggregate_output_immiPT (deSolveOut = scenario5_70_raw,
                                             times = times,
                                             nAges = age$nAges)

scenario5_90_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario5_90)

scenario5_90_out <- aggregate_output_immiPT (deSolveOut = scenario5_90_raw,
                                             times = times,
                                             nAges = age$nAges)

###############################################################################
#MONOCLONAL + MATERNAL VACCINE
#Scenario 6 - Maternal vacc year-round 50% - 70% - 90% coverage + 
#     monoclonal at birth seasonal (April - September)  for at-risk 100% coverage 
#     50% maternal vacc + 100% mAb preterm

scenario6_50_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario6_50)

scenario6_50_out <- aggregate_output_immiPT (deSolveOut = scenario6_50_raw,
                                             times = times,
                                             nAges = age$nAges)


scenario6_70_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario6_70)

scenario6_70_out <- aggregate_output_immiPT (deSolveOut = scenario6_70_raw,
                                             times = times,
                                             nAges = age$nAges)

scenario6_90_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario6_90)

scenario6_90_out <- aggregate_output_immiPT (deSolveOut = scenario6_90_raw,
                                             times = times,
                                             nAges = age$nAges)

###############################################################################
#MONOCLONAL + MATERNAL VACCINE
#Scenario 7 - Maternal vacc year-round OR mAb seasonal with combined coverage 80%
#    and 30% coverage of mAb dose in 2nd season for at-risk
scenario7_mat45mAb35_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario7_mat45mAb35)

scenario7_mat45mAb35_out <- aggregate_output_immiPT (deSolveOut = scenario7_mat45mAb35_raw,
                                             times = times,
                                             nAges = age$nAges)


scenario7_mat60mAb20_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario7_mat60mAb20)

scenario7_mat60mAb20_out <- aggregate_output_immiPT (deSolveOut = scenario7_mat60mAb20_raw,
                                             times = times,
                                             nAges = age$nAges)

################################################################################
#MONOCLONAL + MATERNAL VACCINE
#Scenario 8 - Maternal vacc year-round OR mAb seasonal with combined coverage 90%,
#    70% catchup for out of season and 30% coverage of mAb dose in 2nd season for at-risk
scenario8_mat45mAb35_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario8_mat45mAb35)

scenario8_mat45mAb35_out <- aggregate_output_immiPT (deSolveOut = scenario8_mat45mAb35_raw,
                                                     times = times,
                                                     nAges = age$nAges)


scenario8_mat60mAb20_raw <- deSolve::ode(y0, times, deSolve_risk_imm, scenario8_mat60mAb20)

scenario8_mat60mAb20_out <- aggregate_output_immiPT (deSolveOut = scenario8_mat60mAb20_raw,
                                                     times = times,
                                                     nAges = age$nAges)
############################################################################################
#Find number of hospitalisations each month for monthly 
#age groups up to 60 months for one year
#From deSolve matrix
# col 49: hos_term, #output - hospitalisations term
# col 50: hos_preterm, #hospitalisations preterm
# col 51: hos_comb, #hospitalisations combined
# col 52: doses_matvacc, #total mat vacc doses
# col 53: doses_mAb_0, #total mAb doses at birth
# col 54: doses_mAb_c, #total mAb doses catchup
# col 55: doses_mAb_2 #total mAb doses 2nd season


scenarioList <- c("scenario0", 
                  "scenario1_50", 
                  "scenario1_70",
                  "scenario1_90",
                  "scenario2_50",
                  "scenario2_70",
                  "scenario2_90",
                  "scenario3",
                  "scenario4_50",
                  "scenario4_70",
                  "scenario4_90",
                  "scenario5_50",
                  "scenario5_70",
                  "scenario5_90",
                  "scenario6_50",
                  "scenario6_70",
                  "scenario6_90",
                  "scenario7_mat45mAb35",
                  "scenario7_mat60mAb20",
                  "scenario8_mat45mAb35",
                  "scenario8_mat60mAb20")


yearStart <- run_in + model_time -11
yearEnd <- run_in + model_time
colHosComb <- 51


for(s in scenarioList)
  assign(paste(s, "_hosp", sep = "") , get(paste(s, "_out", sep = ""))$deSolve[yearStart:yearEnd, 1:60, colHosComb])

for(s in scenarioList[-1])
  assign(paste(s, "_avert", sep = ""), scenario0_hosp - get(paste(s, "_hosp", sep = "")) )

scenario1_50_avert[which(scenario1_50_avert < 0)] <- 0 
scenario1_70_avert[which(scenario1_70_avert < 0)] <- 0 
scenario1_90_avert[which(scenario1_90_avert < 0)] <- 0 
scenario2_50_avert[which(scenario2_50_avert < 0)] <- 0 
scenario2_70_avert[which(scenario2_70_avert < 0)] <- 0 
scenario2_90_avert[which(scenario2_90_avert < 0)] <- 0 
scenario3_avert[which(scenario3_avert < 0)] <- 0 
scenario4_50_avert[which(scenario4_50_avert < 0)] <- 0 
scenario4_70_avert[which(scenario4_70_avert < 0)] <- 0 
scenario4_90_avert[which(scenario4_90_avert < 0)] <- 0 
scenario5_50_avert[which(scenario5_50_avert < 0)] <- 0 
scenario5_70_avert[which(scenario5_70_avert < 0)] <- 0 
scenario5_90_avert[which(scenario5_90_avert < 0)] <- 0 
scenario6_50_avert[which(scenario6_50_avert < 0)] <- 0 
scenario6_70_avert[which(scenario6_70_avert < 0)] <- 0 
scenario6_90_avert[which(scenario6_90_avert < 0)] <- 0 
scenario7_mat45mAb35_avert[which(scenario7_mat45mAb35_avert < 0)] <- 0 
scenario7_mat60mAb20_avert[which(scenario7_mat60mAb20_avert < 0)] <- 0 
scenario8_mat45mAb35_avert[which(scenario8_mat45mAb35_avert < 0)] <- 0 
scenario8_mat60mAb20_avert[which(scenario8_mat60mAb20_avert < 0)] <- 0

numDoses <- matrix(NA, nrow = length(scenarioList), ncol = 3)
rownames(numDoses) <- scenarioList
colnames(numDoses) <- c("mAb", "vacc", "total")

i=0
for(s in scenarioList)
{  
  i = i + 1
  numDoses[i, "mAb"] <- sum(get(paste(s, "_out", sep = ""))$doses_mAb_0[yearStart:yearEnd,1:4],
                            get(paste(s, "_out", sep = ""))$doses_mAb_c[yearStart:yearEnd,1:4],
                            get(paste(s, "_out", sep = ""))$doses_mAb_2[yearStart:yearEnd,1:4])
  numDoses[i, "vacc"] <- sum(get(paste(s, "_out", sep = ""))$doses_matvacc[yearStart:yearEnd,1:4])
  numDoses[i, "total"] <- numDoses[i, "mAb"] + numDoses[i, "vacc"]
}

numNeedImm <- as.data.frame(rbind(c("Year-round birth mAb 50%", numDoses[2, "total"]/sum(scenario1_50_avert[,1:24])),
                                c("Year-round birth mAb 70%", numDoses[3, "total"]/sum(scenario1_70_avert[,1:24])),
                                c("Year-round birth mAb 90%", numDoses[4, "total"]/sum(scenario1_90_avert[,1:24])),
                                c("Seasonal birth mAb 50%", numDoses[5, "total"]/sum(scenario2_50_avert[,1:24])),
                                c("Seasonal birth mAb 70%", numDoses[6, "total"]/sum(scenario2_70_avert[,1:24])),
                                c("Seasonal birth mAb 90%", numDoses[7, "total"]/sum(scenario2_90_avert[,1:24])),
                                c("Seasonal birth mAb 90%, catchup 70% 2nd 30%", numDoses[8, "total"]/sum(scenario3_avert[,1:24])),
                                c("Year-round vac 50%", numDoses[9, "total"]/sum(scenario4_50_avert[,1:24])),
                                c("Year-round vac 70%", numDoses[10, "total"]/sum(scenario4_70_avert[,1:24])),
                                c("Year-round vac 90%", numDoses[11, "total"]/sum(scenario4_90_avert[,1:24])),
                                c("Seasonal vac 50%", numDoses[12, "total"]/sum(scenario5_50_avert[,1:24])),
                                c("Seasonal vac 70%", numDoses[13, "total"]/sum(scenario5_70_avert[,1:24])),
                                c("Seasonal vac 90%", numDoses[14, "total"]/sum(scenario5_90_avert[,1:24])),
                                c("Year-round vac 50%, seasonal birth mAb all pT", numDoses[15, "total"]/sum(scenario6_50_avert[,1:24])),
                                c("Year-round vac 70%, seasonal birth mAb all pT", numDoses[16, "total"]/sum(scenario6_70_avert[,1:24])),
                                c("Year-round vac 90%, seasonal birth mAb all pT", numDoses[17, "total"]/sum(scenario6_90_avert[,1:24])),
                                c("Year-round vac 45%, seasonal birth mAb 35%, 2nd season at-risk 30%", numDoses[18, "total"]/sum(scenario7_mat45mAb35_avert[,1:24])),
                                c("Year-round vac 60%, seasonal birth mAb 20%, 2nd season at-risk 30%", numDoses[19, "total"]/sum(scenario7_mat60mAb20_avert[,1:24])),
                                c("Year-round vac 45%, seasonal birth mAb 35%, catchup 70%, 2nd season at-risk 30%", numDoses[20, "total"]/sum(scenario8_mat45mAb35_avert[,1:24])),
                                c("Year-round vac 60%, seasonal birth mAb 20%, catchup 70%, 2nd season at-risk 30%", numDoses[21, "total"]/sum(scenario8_mat60mAb20_avert[,1:24]))))
colnames(numNeedImm) <- c("Scenario", "Efficiency")
numNeedImm$Efficiency <- as.numeric(numNeedImm$Efficiency)
################################################################################
#Create table
totalHospUnder3m <- sum(scenario0_hosp[,1:3])
totalHospUnder6m <- sum(scenario0_hosp[,4:6])
totalHospUnder12m <- sum(scenario0_hosp[,7:12])
totalHospUnder2 <- sum(scenario0_hosp[,13:24])
totalHospAllUnder2 <- sum(scenario0_hosp[,1:24])

tableRes <- matrix(NA, nrow = length(scenarioList)-1, ncol = 10)

i = 0
for(s in scenarioList[-1])
{
  i = i + 1
  tableRes[i,] <- c(sum(get(paste(s, "_avert", sep = ''))[,1:3]), 
                    sum(get(paste(s, "_avert", sep = ''))[,4:6]), 
                    sum(get(paste(s, "_avert", sep = ''))[,7:12]), 
                    sum(get(paste(s, "_avert", sep = ''))[,13:24]),
                    sum(get(paste(s, "_avert", sep = ''))[,1:24]),
                    sum(get(paste(s, "_avert", sep = ''))[,1:3])/totalHospUnder3m, 
                    sum(get(paste(s, "_avert", sep = ''))[,4:6])/totalHospUnder6m, 
                    sum(get(paste(s, "_avert", sep = ''))[,7:12])/totalHospUnder12m, 
                    sum(get(paste(s, "_avert", sep = ''))[,13:24])/totalHospUnder2,
                    sum(get(paste(s, "_avert", sep = ''))[,1:24])/totalHospAllUnder2)
}
colnames(tableRes) <- c("Hosp avert <3m", "Hosp avert 3-<6m", "Hosp avert 6-<12m", "Hosp avert 13-<24m", "Hosp avert <24m", 
                        "Prop avert <3m", "Prop avert 3-<6m", "Prop avert 6-<12m", "Prop avert 13-<24m", "Prop avert <24m")
tableRes <- cbind(tableRes, numDoses[-1,], numNeedImm$Efficiency)
colnames(tableRes)[length(colnames(tableRes))] <- "NNT" 

##FIGURE 2###############################################
#Stacked plot with age
plotData_Age <- tableRes[c("scenario1_70", "scenario2_70", "scenario3",
                           "scenario4_70", "scenario5_70", "scenario6_70",
                           "scenario7_mat45mAb35", "scenario8_mat45mAb35"),1:4]

rownames(plotData_Age) <- c("Year-round mAb", "Seasonal mAb","2024 program", "Year-round vac",
                            "Seasonal vac", "Birth at-risk", "Extra at-risk", "2025 program")
colnames(plotData_Age) <- c("<3m", "3-<6m", "6-<12m", "12-<24m")
plotData_Age <- plotData_Age/totalHospAllUnder2

ylabel <- "Proportion of annual hospitalisations averted < 2 y.o."

meltData_age <- melt(plotData_Age)
colnames(meltData_age) <- c("Scenario", "Age", "Proportions")
scenLabel <- c("1(b)", "2(b)", "3","4(b)", "5(b)", "6(b)", "7(a)", "8(a)")

sumProp <- apply(plotData_Age, 1, sum)

meltData_age$Scenario <- factor(meltData_age$Scenario, levels = c("Year-round mAb", "Seasonal mAb",
                                                                  "2024 program",
                                                                  "Year-round vac", "Seasonal vac",
                                                                  "Birth at-risk",
                                                                  "Extra at-risk",
                                                                  "2025 program"))

meltData_age$Age <- factor(meltData_age$Age, levels = c("12-<24m", "6-<12m", "3-<6m", "<3m"))

ggplot2::ggplot(data=meltData_age, ggplot2::aes(x = Scenario, y = Proportions, fill = Age)) +
  ggplot2::geom_bar(position="stack", stat="identity")+
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, 0.45),
                              breaks = seq(0, 0.45, 0.05), expand = c(0, 0),
                              labels = c(0, seq(0.05, 0.45, 0.05)))+
  scale_x_discrete(label = scenLabel, expand = c(0, 0)) +
  ggplot2::scale_fill_manual(values= c(institute_cols[8], institute_cols[6], institute_cols[3], institute_cols[1]))+
  ggplot2::theme_bw()+
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14), axis.title = element_text(size = 18),
        legend.position = c(0.8, 0.9),
        legend.background = element_rect(fill = "white"),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14))+
  geom_text(x = 1, y = sumProp[1]+0.01, label=round(sumProp[1],3), size=5, col = "black") +
  geom_text(x = 2, y = sumProp[2]+0.01, label=round(sumProp[2],3), size=5, col = "black") +
  geom_text(x = 3, y = sumProp[3]+0.01, label=round(sumProp[3],3), size=5, col = "black") +
  geom_text(x = 4, y = sumProp[4]+0.01, label=round(sumProp[4],3), size=5, col = "black") +
  geom_text(x = 5, y = sumProp[5]+0.01, label=round(sumProp[5],3), size=5, col = "black") +
  geom_text(x = 6, y = sumProp[6]+0.01, label=round(sumProp[6],3), size=5, col = "black") +
  geom_text(x = 7, y = sumProp[7]+0.01, label=round(sumProp[7],3), size=5, col = "black") +
  geom_text(x = 8, y = sumProp[8]+0.01, label=round(sumProp[8],3), size=5, col = "black") +
  geom_text(x = 1, y = 0.135, label = unique(meltData_age$Scenario)[1], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 2, y = 0.12, label = unique(meltData_age$Scenario)[2], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 3, y = 0.11, label = unique(meltData_age$Scenario)[3], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 4, y = 0.125, label = unique(meltData_age$Scenario)[4], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 5, y = 0.11, label = unique(meltData_age$Scenario)[5], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 6, y = 0.1, label = unique(meltData_age$Scenario)[6], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 7, y = 0.105, label = unique(meltData_age$Scenario)[7], vjust = 0.5, hjust = 1, angle = 90, size = 6)+
  geom_text(x = 8, y = 0.115, label = unique(meltData_age$Scenario)[8], vjust = 0.5, hjust = 1, angle = 90, size = 6)

##FIGURE 3#########################################################################
plotData <- as.data.frame(rbind(c("Year-round mAb", "50%", sum(scenario1_50_avert[,1:24])/totalHospAllUnder2),
                                c("Year-round mAb", "70%", sum(scenario1_70_avert[,1:24])/totalHospAllUnder2),
                                c("Year-round mAb", "90%", sum(scenario1_90_avert[,1:24])/totalHospAllUnder2),
                                c("Seasonal mAb", "50%", sum(scenario2_50_avert[,1:24])/totalHospAllUnder2),
                                c("Seasonal mAb", "70%", sum(scenario2_70_avert[,1:24])/totalHospAllUnder2),
                                c("Seasonal mAb", "90%", sum(scenario2_90_avert[,1:24])/totalHospAllUnder2),
                                c("Year-round vac", "50%", sum(scenario4_50_avert[,1:24])/totalHospAllUnder2),
                                c("Year-round vac", "70%", sum(scenario4_70_avert[,1:24])/totalHospAllUnder2),
                                c("Year-round vac", "90%", sum(scenario4_90_avert[,1:24])/totalHospAllUnder2),
                                c("Seasonal vac", "50%", sum(scenario5_50_avert[,1:24])/totalHospAllUnder2),
                                c("Seasonal vac", "70%", sum(scenario5_70_avert[,1:24])/totalHospAllUnder2),
                                c("Seasonal vac", "90%", sum(scenario5_90_avert[,1:24])/totalHospAllUnder2), 
                                c("Hybrid at-risk", "50%", sum(scenario6_50_avert[,1:24])/totalHospAllUnder2),
                                c("Hybrid at-risk", "70%", sum(scenario6_70_avert[,1:24])/totalHospAllUnder2),
                                c("Hybrid at-risk", "90%", sum(scenario6_90_avert[,1:24])/totalHospAllUnder2)))

plotData_tradeoff <- as.data.frame(rbind(c("Extra at-risk", "mat vac 45% mAb 35%", sum(scenario7_mat45mAb35_avert[,1:24])/totalHospAllUnder2),
                                         c("Extra at-risk", "mat vac 60% mAb 20%", sum(scenario7_mat60mAb20_avert[,1:24])/totalHospAllUnder2),
                                         c("2025 program", "mat vac 45% mAb 35%", sum(scenario8_mat45mAb35_avert[,1:24])/totalHospAllUnder2),
                                         c("2025 program", "mat vac 60% mAb 20%", sum(scenario8_mat60mAb20_avert[,1:24])/totalHospAllUnder2)))

colnames(plotData) <- c("Scenario", "Coverage", "Proportions")
colnames(plotData_tradeoff) <- colnames(plotData)
plotData$Proportions <- as.numeric(plotData$Proportions)
plotData_tradeoff$Proportions <- as.numeric(plotData_tradeoff$Proportions)

plotData$Scenario <- factor(plotData$Scenario, levels = c("Year-round mAb", "Seasonal mAb",
                                                          "Year-round vac", "Seasonal vac",
                                                          "Hybrid at-risk"))

plotData$Coverage <- factor(plotData$Coverage, levels = c("50%", "70%", "90%"))
scenLabel <- c(paste("Year-round", "\n", "mAb"),
               paste("Seasonal", "\n", "mAb"),
               paste("Year-round", "\n", "vac"), 
               paste("Seasonal", "\n", "vac"),
               paste("Birth mAb", "\n", "at-risk*", "\n", "hybrid"))

plotData_tradeoff$Scenario <- factor(plotData_tradeoff$Scenario, 
                                     levels = c("Extra at-risk", "2025 program")) 
scenLabel_tradeoff <- c(paste("Extra", "\n", "mAb at-risk", "\n", "hybrid"), 
                        paste("2025", "\n", "program", "\n", "hybrid"))


ggplot2::ggplot(data=plotData, ggplot2::aes(x = Scenario, y = Proportions, fill = Coverage)) +
  ggplot2::geom_bar(position="dodge", stat="identity", width = 0.9)+
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, 0.45),
                              breaks = seq(0, 0.45, 0.05), expand = c(0, 0),
                              labels = c(0, seq(0.05, 0.45, 0.05)))+
  ggplot2::xlab("")+
  scale_x_discrete(labels = scenLabel, expand = c(0, 0)) +
  ggplot2::scale_fill_manual(values= c(institute_cols[5], institute_cols[6], institute_cols[7]))+
  ggplot2::theme_bw()+
  theme(legend.position = c(0.8, 0.9),legend.background = element_rect(fill = "white"), 
        legend.title = element_text(size = 14), axis.title.y = element_text(size = 18), 
        axis.text.x = element_text(size = 14), legend.text = element_text(size = 14))

ggplot2::ggplot(data=plotData_tradeoff, ggplot2::aes(x = Scenario, y = Proportions, fill = Coverage)) +
  ggplot2::geom_bar(position="dodge", stat="identity", width = 0.9)+
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, 0.45),
                              breaks = seq(0, 0.45, 0.05), expand = c(0, 0),
                              labels = c(0, seq(0.05, 0.45, 0.05)))+
  ggplot2::xlab("")+
  scale_x_discrete(labels = scenLabel_tradeoff, expand = c(0, 0)) +
  ggplot2::scale_fill_manual(values= c(institute_cols[5], institute_cols[6], institute_cols[7]))+
  ggplot2::theme_bw()+
  theme(legend.title = element_text(size = 14), axis.title.y = element_text(size = 18), 
    axis.text.x = element_text(size = 14), legend.text = element_text(size = 14))

##FIGURE 4#######################################################################
#Num doses
numDoses <- matrix(NA, nrow = length(scenarioList), ncol = 3)
rownames(numDoses) <- scenarioList
colnames(numDoses) <- c("mAb", "vacc", "total")

i=0
for(s in scenarioList)
{  
  i = i + 1
  numDoses[i, "mAb"] <- sum(get(paste(s, "_out", sep = ""))$doses_mAb_0[yearStart:yearEnd,1:4],
                            get(paste(s, "_out", sep = ""))$doses_mAb_c[yearStart:yearEnd,1:4],
                            get(paste(s, "_out", sep = ""))$doses_mAb_2[yearStart:yearEnd,1:4])
  numDoses[i, "vacc"] <- sum(get(paste(s, "_out", sep = ""))$doses_matvacc[yearStart:yearEnd,1:4])
  numDoses[i, "total"] <- numDoses[i, "mAb"] + numDoses[i, "vacc"]
}

#Total dosage cost using CDC costs
CDC_Cost <- as.data.frame(rbind(c("mAb", "Year-round mAb"),
                                c("mAb", "Seasonal mAb"),
                                c("mAb", "2024 program"),
                                c("mat vac", "Year-round vac"),
                                c("mat vac", "Seasonal vac"),
                                c("hybrid", "Birth at-risk"),
                                c("hybrid", "Extra at-risk"),
                                c("hybrid", "2025 program"))) 

for(i in c(200, 300, 400))
{  
  CDC_Cost <- as.data.frame(cbind(CDC_Cost, 
                                  c(i*numDoses[3, "total"],
                                    i*numDoses[6, "total"],
                                    i*numDoses[8, "total"],
                                    200*numDoses[10, "total"],
                                    200*numDoses[13, "total"],
                                    (200*numDoses[16, "vacc"] + i*numDoses[16, "mAb"]),
                                    (200*numDoses[18, "vacc"] + i*numDoses[19, "mAb"]),
                                    (200*numDoses[20, "vacc"] + i*numDoses[21, "mAb"]))))
}
colnames(CDC_Cost) <- c("Program", "Scenario", "mAb = $200, vac = $200", "mAb = $300, vac = $200", "mAb = $400, vac = $200")

load("data/data.southernWA.2016.rda")
totalPop <- sum(data.southernWA.2016$population)
popUnder2 <-(2/5)*data.southernWA.2016$population[1]

#cost per 1000 people
CDC_cost1000 <- cbind(CDC_Cost[,1:2], (CDC_Cost[,3:5]/totalPop)*1000)
CDC_costUnder2 <- cbind(CDC_Cost[,1:2], (CDC_Cost[,3:5]/popUnder2))

dataCost <- melt(CDC_costUnder2)
colnames(dataCost)[3:4] <- c("mAbCost", "Cost")
dataCost$Cost <- as.numeric(dataCost$Cost)

dataCost$Scenario <- factor(dataCost$Scenario, levels = c("Year-round mAb", "Seasonal mAb",
                                                          "2024 program",
                                                          "Year-round vac", "Seasonal vac",
                                                          "Birth at-risk",
                                                          "Extra at-risk",
                                                          "2025 program"))

ylabel <- "Annual dosage cost per person < 2 y.o ($US)"

scenLabel <- c("1(b)", "2(b)", "3","4(b)", "5(b)", "6(b)", "7(a)", "8(a)")

ggplot2::ggplot(data=dataCost, ggplot2::aes(x = Scenario, y = Cost, fill = Program)) +
  ggplot2::geom_bar(position="dodge", stat="identity", width = 0.9)+
  ggplot2::scale_y_continuous(name = ylabel, limits = c(0, 150),
                              breaks = seq(0, 150, 20), expand = c(0, 0),
                              labels = c(0, seq(20, 150, 20)))+
  geom_text(
    aes(label = Scenario), 
    vjust = 0.5, hjust = 1, angle = 90, size = 5, nudge_y = -5
  )+
  scale_x_discrete(label = scenLabel, expand = c(0, 0)) +
  ggplot2::scale_fill_manual(values= c(institute_cols[1], institute_cols[3], institute_cols[6]))+
  ggplot2::theme_bw()+
  facet_wrap(~ mAbCost)+
  theme(strip.background = element_blank(), strip.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 14), axis.title = element_text(size = 18),
        legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white"),
        legend.title = element_text(size = 14), legend.text = element_text(size = 14))


