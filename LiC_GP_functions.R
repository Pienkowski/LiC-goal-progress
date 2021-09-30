############################################################################
######### Life in Conservation: Goal progress functions ####################
############################################################################

# Label re-coder 
label_recoder <- function(x){
  x <-  as.factor(ifelse(x == "years_cons", "Years in conservation",
                         ifelse(x == "DO", "Dispositional optimism", 
                                ifelse(x %in% c("gender_Male" , "genderMale"), "Male", 
                                       ifelse(x == "gender_Unknown", "Gender: unknown", 
                                              ifelse(x %in% c("position_simple_Practice" , "position_simplePractice"), "Practice/policy", 
                                                     ifelse(x == "position_simple_Unknown.other", "Position: unknown", 
                                                            ifelse(x == "years_cons_scaled", "Years in conservation",
                                                                   ifelse(x %in% c("Environment_CC","EnvironmentCC" ), "Cross-cutting",
                                                                          ifelse(x %in% c("Environment_Marine", "EnvironmentMarine") , "Marine",
                                                                                 ifelse(x == "Environment_Unknown", "Environment: unknown",
                                                                                        ifelse(x == "education_simple_University", "University",
                                                                                               ifelse(x %in%  c("SO_Region_CentralandSouthernAsia", "SO_RegionCentral and Southern Asia" ) ,"C. & S. Asia",
                                                                                                      ifelse(x  %in%  c("SO_Region_EasternandSouth_EasternAsia", "SO_RegionEastern and South-Eastern Asia" ),"E. & S.E. Asia", 
                                                                                                             ifelse(x  %in%  c("SO_Region_EuropeandNorthernAmerica"),"Europe & N. America", 
                                                                                                                    ifelse(x  %in%  c("SO_Region_LatinAmericaandtheCaribbean", "SO_RegionLatin America and the Caribbean"),"Lat. America & Carib.", 
                                                                                                                           ifelse(x  %in%  c("SO_Region_NorthernAfricaandWesternAsia", "SO_RegionNorthern Africa and Western Asia"),"N. Africa & W. Asia", 
                                                                                                                                  ifelse(x  %in%  c("SO_Region_Oceania", "SO_RegionOceania"),"Oceania", 
                                                                                                                                         ifelse(x %in%  c("SO_Region_Sub_SaharanAfrica", "SO_RegionSub-Saharan Africa"),"Sub-Saharan Africa", 
                                                                                                                                                ifelse(x == "rli.project","Projected Red List Index",
                                                                                                                                                       ifelse(x == "rli","Red List Index",
                                                                                                                                                              ifelse(x == "Prim.cover.20","Primary cover",
                                                                                                                                                                     ifelse(x == "I_VI_mean","Protected area cover", 
                                                                                                                                                                            ifelse(x == "spend.area","Spending", 
                                                                                                                                                                                   ifelse(x == "Gov.Total","Governance", x
                                                                                                                                                                                   )))))))))))))))))))))))))}

