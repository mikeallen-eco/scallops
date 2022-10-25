# https://bcss.org.my/tut/bayes-with-jags-a-tutorial-for-wildlife-researchers/occupancy-modelling/dynamic-occupancy-modelling/

# Get the crossbill data
cb <- read.csv("http://bcss.org.my/data/crossbills.csv", comment="#")
str(cb)
# Extract the detection histories
DH <- as.matrix(cb[, 7:42])
head(DH)

# convert to sites x occasions x year array
Y <- array(DH, dim=c(nrow(DH), 3, 12))
head(Y)

# Aggregate detection histories across occasions
y <- apply(Y, c(1, 3), sum, na.rm=TRUE)  # sites by years, number of detections
n <- apply(!is.na(Y), c(1, 3), sum)      # ... and no. of surveys

# Get known values of z
z <- (y > 0)*1
z[z == 0] <- NA
jagsData <- list(nSites = nrow(y), nYears = ncol(y), y = y,
                 n = n, z = z)
str(jagsData)

wanted <- c("psi1", "phi", "gamma", "p","N")
