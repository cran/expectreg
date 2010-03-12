desc <- packageDescription("expectreg")
year <- sub(".*(2[[:digit:]]{3})-.*", "\\1", desc$Date)
vers <- paste("R package version", desc$Version)
citEntry(entry="Manual",
title = "expectreg: Expectile and Quantile Regression",
author = personList(as.person("Fabian Sobotka"),
as.person("Thomas Kneib"),
as.person("Sabine Schnabel"),
year = year,
note = vers,
textVersion =
paste("Fabian Sobotka, Thomas Kneib, Sabine Schnabel (",
year,
"). expectreg: Expectile and Quantile Regression. ",
vers, ".", sep=""))
