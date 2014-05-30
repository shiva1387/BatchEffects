# plotting p values of features differing based on biochem traits

#Chlorella
Chlorella_p<-p_values_sig
png("day12_differential_highChlorellaVsOthers.png")
plot(seq(1:length(Chlorella_p)),Chlorella_p,main="Day 12 Chlorella")
dev.off()
#growth
growth_p<-p_values_sig
png("day4_differential_highGrowthVslow.png")
plot(seq(1:length(growth_p)),growth_p,main="Day 4 growth")
dev.off()
#biomass
biomass_p<-p_values_sig
png("day12_differential_highBiomassVslow.png")
plot(seq(1:length(biomass_p)),biomass_p,main="Day 12 biomass")
dev.off()
#biomassProd
biomassProd_p<-p_values_sig
png("day12_differential_highbiomassProdVslow.png")
plot(seq(1:length(biomassProd_p)),biomassProd_p,main="Day 12 biomassProd")
dev.off()
#palm
palm_p<-p_values_sig
png("day12_differential_highPalmVsOthers.png")
plot(seq(1:length(palm_p)),palm_p,main="Day 12 palm")
dev.off()
#totalLipid
totalLipid_p<-p_values_sig
png("day12_differential_hightotalLipidVsOthers.png")
plot(seq(1:length(totalLipid_p)),totalLipid_p,main="Day 12 totalLipid")
dev.off()
#lipidProd
lipidProd_p<-p_values_sig
png("day12_differential_highlipidProdVsOthers.png")
plot(seq(1:length(lipidProd_p)),lipidProd_p,main="Day 12 lipidProd")
dev.off()
#totalProtein
totalProtein_p<-p_values_sig
png("day12_differential_hightotalProteinVsOthers.png")
plot(seq(1:length(totalProtein_p)),totalProtein_p,main="Day 12 totalProtein")
dev.off()
#Carbon
Carbon_p<-p_values_sig
png("day12_differential_highCarbonVsOthers.png")
plot(seq(1:length(Carbon_p)),Carbon_p,main="Day 12 Carbon")
dev.off()




