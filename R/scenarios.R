################################################################################
#SET UP SCENARIOS
################################################################################
#Fixed pars for all scenarios
pars_imm <- list(b0 = b0_fit,
                    b1 = b1_fit,
                    phi = phi_fit,
                    delta = delta,
                    gamma0 = gamma0,
                    gamma1 = gamma1,
                    nu = nu,
                    nuM = nuM, #durability of mAb
                    nuV = nuV, #durability of mat vacc
                    omega_vect = omega_vect,
                    omegaE = omegaE,
                    alpha = alpha,
                    A = A_fit,
                    B = B_fit,
                    C = C_fit,
                    D = D_fit,
                    E = E_fit,
                    rho = rho, #efficacy for terms
                    rho_PT = rho_PT, #efficacy for preterms
                    age_in = age$age_in,
                    age_out = age$age_out,
                    age_months = age$age_years*12,
                    sigma_vect = sigma_vect,
                    sigmaE = sigmaE,
                    mixing = mixing,
                    nAges = age$nAges)

#Zero coverage
kappa_zero <- rep(0, n_times)

kappaM_B <- c(0.2, 0.4, 0.5, 0.7, 0.9) #coverage of mAb at birth
kappaM_C <- 0.7 #coverage of mAb catchup 
kappaM_2 <- 0.3 #coverage of mAb 2nd dose   

kappaV <- c(0.5, 0.7, 0.9) #coverage of mat vacc

#Seasonal coverage
seasonal <- rep(0, 12) # Jan - Dec
seasonal[4:8] <- 1 #mAb at birth for those born in May - Sept, adjusting for t = 0
atBirth_April <- seasonal
atBirth_April[3] <- 1 #mAb at birth for those born in April - Sept
#kappaM_seasonal <-  c(rep(0, run_in), kappaM * rep(seasonal, model_time/12))
#kappaM_seasonalApril <- c(rep(0, run_in), kappaM * rep(atBirth_April, model_time/12))


#Catchup
catchUp <- rep(0, 12) # Jan - Dec
catchUp[3] <- 1 #assume optimal catchup occurring in April
#kappaM_catchUp <- c(rep(0, run_in), kappaM * rep(catchUp, model_time/12))

#First RSV season
age_catchUp <- 6 #Up to 7 month olds in first catchup month April [6,7)
#Second RSV season dose
age_catchUp_2seasons <- 18 #Up to and including 18 month olds [0,19) 

################################################################################
#For each scenario, define 
# Maternal coverage (monthly vector) 
#                     - kappaV (terms) 
#                     - kappaV_PT (preterms)
# Monoclonal coverage at birth (monthly vector)
#                     - kappaM_Birth (terms)
#                     - kappaM_BirthPT (preterms)
# Monoclonal coverage catchup (monthly vector)
#                     - kappaM_catchUp (terms)
#                     - kappaM_catchUpPT (preterms)
# Monoclonal catchup age, covers birth to this age group for all
#                     - catchAge 
# Monoclonal coverage dose 2 (monthly vector)
#                     - kappaM_dose2 (terms)
#                     - kappaM_dose2PT (preterms)
# Monoclonal dose 2 age, covers catchAge +1 to this age group for all
#                     - dose2Age
#################################################################################
#BASE
#Scenario 0 - No immunisation
scenario0 <- append(pars_imm, list(kappaV = kappa_zero,
                                  kappaV_PT = kappa_zero,
                                  kappaM_Birth = kappa_zero, 
                                  kappaM_BirthPT = kappa_zero, 
                                  kappaM_catchUp = kappa_zero,
                                  kappaM_catchUpPT = kappa_zero,
                                  catchAge = 0, 
                                  kappaM_dose2 = kappa_zero, 
                                  kappaM_dose2PT = kappa_zero, 
                                  dose2Age = 0)) 
#MONOCLONAL
#Scenario 1 - Monoclonal year-round
#50% 
scenario1_50 <- append(pars_imm, list(kappaV = kappa_zero,
                                   kappaV_PT = kappa_zero,
                                   kappaM_Birth = c(rep(0, run_in), rep(0.5, model_time)), 
                                   kappaM_BirthPT = c(rep(0, run_in), rep(0.5, model_time)), 
                                   kappaM_catchUp = kappa_zero,
                                   kappaM_catchUpPT = kappa_zero,
                                   catchAge = 0, 
                                   kappaM_dose2 = kappa_zero, 
                                   kappaM_dose2PT = kappa_zero, 
                                   dose2Age = 0)) 
#70% 
scenario1_70 <- append(pars_imm, list(kappaV = kappa_zero,
                                      kappaV_PT = kappa_zero,
                                      kappaM_Birth = c(rep(0, run_in), rep(0.7, model_time)), 
                                      kappaM_BirthPT = c(rep(0, run_in), rep(0.7, model_time)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0))
#90% 
scenario1_90 <- append(pars_imm, list(kappaV = kappa_zero,
                                      kappaV_PT = kappa_zero,
                                      kappaM_Birth = c(rep(0, run_in), rep(0.9, model_time)), 
                                      kappaM_BirthPT = c(rep(0, run_in), rep(0.9, model_time)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#Scenario 2 - Monoclonal seasonal at birth
#50% 
scenario2_50 <- append(pars_imm, list(kappaV = kappa_zero,
                                      kappaV_PT = kappa_zero,
                                      kappaM_Birth = c(rep(0, run_in), 0.5 * rep(atBirth_April, model_time/12)), 
                                      kappaM_BirthPT = c(rep(0, run_in), 0.5 * rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#70% 
scenario2_70 <- append(pars_imm, list(kappaV = kappa_zero,
                                      kappaV_PT = kappa_zero,
                                      kappaM_Birth = c(rep(0, run_in), 0.7 * rep(atBirth_April, model_time/12)), 
                                      kappaM_BirthPT = c(rep(0, run_in), 0.7 * rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0))
#90% 
scenario2_90 <- append(pars_imm, list(kappaV = kappa_zero,
                                      kappaV_PT = kappa_zero,
                                      kappaM_Birth = c(rep(0, run_in), 0.9 * rep(atBirth_April, model_time/12)), 
                                      kappaM_BirthPT = c(rep(0, run_in), 0.9 * rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#Scenario 3 - Seasonal with catchup + extra dose at-risk 79% newborn, 65% catchup, 30% at-risk
scenario3 <- append(pars_imm, list(kappaV = kappa_zero,
                                      kappaV_PT = kappa_zero,
                                      kappaM_Birth = c(rep(0, run_in), 0.79 * rep(atBirth_April, model_time/12)), 
                                      kappaM_BirthPT = c(rep(0, run_in), 0.79 * rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = c(rep(0, run_in), 0.65 * rep(catchUp, model_time/12)),
                                      kappaM_catchUpPT = c(rep(0, run_in), 0.65 * rep(catchUp, model_time/12)),
                                      catchAge = 6, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = c(rep(0, run_in), 0.3 * rep(catchUp, model_time/12)), 
                                      dose2Age = 18)) 


#MATERNAL VACCINE
#Scenario 4 - Maternal vacc year-round
#50% 
scenario4_50 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.5, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.5, model_time)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = kappa_zero, 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#70% 
scenario4_70 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.7, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.7, model_time)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = kappa_zero, 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0))
#90% 
scenario4_90 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.9, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.9, model_time)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = kappa_zero, 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#Scenario 5 - Maternal vacc seasonal
#50% 
scenario5_50 <- append(pars_imm, list(kappaV = c(rep(0, run_in), 0.5 * rep(atBirth_April, model_time/12)),
                                      kappaV_PT = c(rep(0, run_in), 0.5 * rep(atBirth_April, model_time/12)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = kappa_zero, 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#70% 
scenario5_70 <- append(pars_imm, list(kappaV = c(rep(0, run_in), 0.7 * rep(atBirth_April, model_time/12)),
                                      kappaV_PT = c(rep(0, run_in), 0.7 * rep(atBirth_April, model_time/12)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = kappa_zero, 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0))
#90% 
scenario5_90 <- append(pars_imm, list(kappaV = c(rep(0, run_in), 0.9 * rep(atBirth_April, model_time/12)),
                                      kappaV_PT = c(rep(0, run_in), 0.9 * rep(atBirth_April, model_time/12)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = kappa_zero, 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 

#MONOCLONAL + MATERNAL VACCINE
#Scenario 6 - Maternal year-round + seasonal monoclonal at birth for at-risk 
#             Overriding protection is mAb for at-risk from April-September
#50% maternal vacc + 100% mAb preterm
scenario6_50 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.5, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.5 * as.numeric(!atBirth_April), model_time/12)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = c(rep(0, run_in), rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#70% maternal vacc + 100% mAb preterm
scenario6_70 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.7, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.7 * as.numeric(!atBirth_April), model_time/12)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = c(rep(0, run_in), rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0))
#90% maternal vacc + 100% mAb preterm
scenario6_90 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.9, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.9 * as.numeric(!atBirth_April), model_time/12)),
                                      kappaM_Birth = kappa_zero, 
                                      kappaM_BirthPT = c(rep(0, run_in), rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 0, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = kappa_zero, 
                                      dose2Age = 0)) 
#Scenario 7 - Maternal year-round OR monoclonal seasonal at birth + extra dose at-risk
#Maternal vacc 45%, mAb 35%, 2nd at-risk 30%
scenario7_mat45mAb35 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.45, model_time)),
                                      kappaV_PT = c(rep(0, run_in), rep(0.45, model_time)),
                                      kappaM_Birth = c(rep(0, run_in), 0.35 * rep(atBirth_April, model_time/12)), 
                                      kappaM_BirthPT = c(rep(0, run_in), 0.35 * rep(atBirth_April, model_time/12)), 
                                      kappaM_catchUp = kappa_zero,
                                      kappaM_catchUpPT = kappa_zero,
                                      catchAge = 6, 
                                      kappaM_dose2 = kappa_zero, 
                                      kappaM_dose2PT = c(rep(0, run_in), 0.3 * rep(catchUp, model_time/12)), 
                                      dose2Age = 18)) 
#Maternal vacc 60%, mAb 20%, 2nd at-risk 30%
scenario7_mat60mAb20 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.6, model_time)),
                                              kappaV_PT = c(rep(0, run_in), rep(0.6, model_time)),
                                              kappaM_Birth = c(rep(0, run_in), 0.2 * rep(atBirth_April, model_time/12)), 
                                              kappaM_BirthPT = c(rep(0, run_in), 0.2 * rep(atBirth_April, model_time/12)), 
                                              kappaM_catchUp = kappa_zero,
                                              kappaM_catchUpPT = kappa_zero,
                                              catchAge = 6, 
                                              kappaM_dose2 = kappa_zero, 
                                              kappaM_dose2PT = c(rep(0, run_in), 0.3 * rep(catchUp, model_time/12)), 
                                              dose2Age = 18)) 
#Scenario 8 - Maternal year-round OR monoclonal seasonal at birth + catchup + extra dose at-risk
#Maternal vacc 45%, mAb 35%, catchup 70%, 2nd at-risk 30%
scenario8_mat45mAb35 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.45, model_time)),
                                              kappaV_PT = c(rep(0, run_in), rep(0.45, model_time)),
                                              kappaM_Birth = c(rep(0, run_in), 0.35 * rep(atBirth_April, model_time/12)), 
                                              kappaM_BirthPT = c(rep(0, run_in), 0.35 * rep(atBirth_April, model_time/12)), 
                                              kappaM_catchUp = c(rep(0, run_in), 0.7 * rep(catchUp, model_time/12)),
                                              kappaM_catchUpPT = c(rep(0, run_in), 0.7 * rep(catchUp, model_time/12)),
                                              catchAge = 6,  
                                              kappaM_dose2 = kappa_zero, 
                                              kappaM_dose2PT = c(rep(0, run_in), 0.3 * rep(catchUp, model_time/12)), 
                                              dose2Age = 18)) 
#Maternal vacc 60%, mAb 20%, catchup 70%, 2nd at-risk 30%
scenario8_mat60mAb20 <- append(pars_imm, list(kappaV = c(rep(0, run_in), rep(0.6, model_time)),
                                              kappaV_PT = c(rep(0, run_in), rep(0.6, model_time)),
                                              kappaM_Birth = c(rep(0, run_in), 0.2 * rep(atBirth_April, model_time/12)), 
                                              kappaM_BirthPT = c(rep(0, run_in), 0.2 * rep(atBirth_April, model_time/12)), 
                                              kappaM_catchUp = c(rep(0, run_in), 0.7 * rep(catchUp, model_time/12)),
                                              kappaM_catchUpPT = c(rep(0, run_in), 0.7 * rep(catchUp, model_time/12)),
                                              catchAge = 6, 
                                              kappaM_dose2 = kappa_zero, 
                                              kappaM_dose2PT = c(rep(0, run_in), 0.3 * rep(catchUp, model_time/12)), 
                                              dose2Age = 18)) 

