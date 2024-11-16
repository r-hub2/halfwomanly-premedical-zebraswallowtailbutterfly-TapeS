n <- 1e3
TF <- logical(n)
for(i in 1:n){
  (spp <- sample(1:36, 2))
  (dbh <- sample(1:100, 2))
  (h <- sapply(seq(along=spp), function(a) fHoehenTarif(BaMap(spp[a], 6), dbh[a])))
  (tps <- biomass(BaMap(spp, 6), dbh, dbh*0.8, h))
  (bdt <- rBDATPRO::getBiomass(list(spp=spp, D1=dbh, D2=0.8*dbh, H2=0.3*h, H=h)))
  (TF[i] <- isTRUE(all.equal(tps, bdt, tolerance = 1e-3)))
  if(isFALSE(TF[i])) {print((tps-bdt)/bdt); print(spp); print(dbh); print(h)}
}
mean(TF)
