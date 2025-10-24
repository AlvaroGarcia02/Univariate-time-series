#==============================================================================
#===                         CÓDIGO TRABAJO I                               ===
#==============================================================================
#= Series de tiempo univ. - Álvaro García, María Paula Arévalo y Tomás Ávila  =
#===                              2025-2                                    ===
#==============================================================================

rm(list = ls())

data <- read.csv("Net_generation_for_all_sectors.csv", skip = 4)
data <- data[12, -c(1:3)]
data <- as.data.frame(t(data)); colnames(data) <- c("Generation")
series <- ts(as.numeric(data$Generation), start = c(2001, 1), frequency = 12)

plot(series, ylab = "Energía producida (en GWh)", xlab = "Tiempo")
abline(h = mean(series), lwd = 2, lty = 2)
lines(lowess(series), lwd = 2)
legend("topleft", legend = c("Media","Tendencia"),
       lwd = 2, lty = 2:1, bty = "n")

library(moments)
length(series)
mean(series)
sd(series)
sd(series)/mean(series)
skewness(series)
kurtosis(series)

boxplot(series,
        horizontal = TRUE,
        xlab = "Energía producida (en GWh)")
points(mean(series), 1, pch = 19, cex = 1.5)
summary(series)

boxplot(series ~ cycle(series),
        xlab = "Mes",
        ylab = "Energía producida (en GWh)")

par(mfrow = c(1, 2))
acf(series, main = "", xlab = "Rezago", lag.max = 36)
pacf(series, main = "", xlab = "Rezago", lag.max = 36)
par(mfrow = c(1, 1))

library(tseries)
adf.test(series)

d.series <- diff(series)
plot(d.series)
acf(d.series, lag.max = 36)
adf.test(d.series)

d.log.series <- diff(log(series))
plot(d.log.series, ylab = "Retorno", xlab = "Tiempo")
abline(h = mean(d.log.series), lwd = 2, lty = 2)
lines(lowess(d.log.series), lwd = 2)
legend("topright", legend = c("Media", "Tendencia"),
       lwd = 2, lty = 2:1, bty = "n")
par(mfrow = c(1, 2))
acf(d.log.series, main = "", xlab = "Rezago", lag.max = 36)
pacf(d.log.series, main = "", xlab = "Rezago", lag.max = 36)
par(mfrow = c(1, 1))

dd.log.series <- diff(d.log.series)
plot(dd.log.series, ylab = "Diferenciación del retorno", xlab = "Tiempo")
abline(h = mean(d.log.series), lwd = 2, lty = 2)
lines(lowess(d.log.series), lwd = 2)
legend("topright", legend = c("Media", "Tendencia"),
       lwd = 2, lty = 2:1, bty = "n")
par(mfrow = c(1, 2))
acf(dd.log.series, main = "", xlab = "Rezago", lag.max = 36)
pacf(dd.log.series, main = "", xlab = "Rezago", lag.max = 36)
par(mfrow = c(1, 1))

adf.test(d.log.series)
adf.test(dd.log.series)

# Ajustar un modelo ARMA para dd.log.series:

# AR(5)
m1 <- arima(dd.log.series, order = c(5, 0, 0))
m1
-0.8084/0.0567
-0.5354/0.0681
-0.4791/0.0696
-0.4897/0.0683
-0.2395/0.0570
library(FinTS)
library(polynom)
par(mfrow = c(1, 2))
pacf(dd.log.series, main = "", xlab = "Rezago", lag.max = 36)
plotArmaTrueacf(c(-.8084, -.5354, -.4791, -.4897, -.2395), pacf = TRUE, main = "", xlab = "Rezago", lag.max = 36)
par(mfrow = c(1, 1))

# MA(1)
m2 <- arima(dd.log.series, order = c(0, 0, 1))
m2
-0.9999/0.0085
library(FinTS)
library(polynom)
par(mfrow = c(1,2))
acf(dd.log.series, main = "", xlab = "Rezago", lag.max = 36)
plotArmaTrueacf(list(ma = -0.9999), main = "", xlab = "Rezago", lag.max = 36)
par(mfrow = c(1, 1))

# ARMA(5,1)
m3 <- arima(dd.log.series, order = c(5, 0, 1))
m3
-0.1145/0.0583
-0.0348/0.0568
-0.1898/0.0559
-0.2654/0.0569
-0.0671/0.0588
-1/0.0087

# ARMA(3,1)
m4 <- arima(dd.log.series, order = c(3, 0, 1))
m4
-0.0511/0.0576
-0.0155/0.0577
-0.1707/0.0579
-1/0.0086

# Selección automática de ARIMA
library(forecast)
m.auto <- auto.arima(dd.log.series, seasonal = FALSE)
m.auto
-0.5212/0.0497

par(mfrow = c(1,2))
acf(dd.log.series, main = "", xlab = "Rezago", lag.max = 36)
plotArmaTrueacf(list(ar=c(-0.5212)),
                main = "", xlab = "Rezago", lag.max = 36)

pacf(dd.log.series, main = "", xlab = "Rezago", lag.max = 36)
plotArmaTrueacf(list(ar=c(-0.5212)),
                pacf = TRUE, main = "", xlab = "Rezago", lag.max = 36)
par(mfrow = c(1,1))

# Residuales:
par(mfrow = c(1, 2))
acf(m1$residuals, main = "", xlab = "Rezago", lag.max = 36)
pacf(m1$residuals, main = "", xlab = "Rezago", lag.max = 36)
acf(m2$residuals, main = "", xlab = "Rezago", lag.max = 36)
pacf(m2$residuals, main = "", xlab = "Rezago", lag.max = 36)
acf(m3$residuals, main = "", xlab = "Rezago", lag.max = 36)
pacf(m3$residuals, main = "", xlab = "Rezago", lag.max = 36)
acf(m4$residuals, main = "", xlab = "Rezago", lag.max = 36)
pacf(m4$residuals, main = "", xlab = "Rezago", lag.max = 36)
acf(m.auto$residuals, main = "", xlab = "Rezago", lag.max = 36)
pacf(m.auto$residuals, main = "", xlab = "Rezago", lag.max = 36)

Box.test(m1$residuals, type = "Ljung-Box")
Box.test(m2$residuals, type = "Ljung-Box")
Box.test(m3$residuals, type = "Ljung-Box")
Box.test(m4$residuals, type = "Ljung-Box")
Box.test(m.auto$residuals, type = "Ljung-Box")

AIC(m1, m2, m3, m4, m.auto)
BIC(m1, m2, m3, m4, m.auto)



