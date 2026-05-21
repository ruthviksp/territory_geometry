rm(list = ls())

# Blackbuck Lekking Site Choice Model (Continuous-time Gillespie)
# Author: AR, June 2025
# Modified for reproducibility + CSV export

# ----------------------
# Scenario:
# Each male blackbuck selects a lek site (0 or 1).
# Choices are influenced by:
# (a) stochastic individual exploration (random switching), and
# (b) copying the choice of another individual (social influence).
# ----------------------

# Set seed
seed <- 12345
set.seed(seed)

# Parameters
N <- 100            # Number of male blackbucks
r <- 100             # Copying rate
s <- 0.5             # Spontaneous switching rate
Tint <- 0.5         # Recording interval
Tend <- 100         # Total simulation time
iter <- ceiling(Tend / Tint)   # Number of recording intervals

# Initial state
# All males start at lek site 1
x <- rep(0, N)

# State trackers
xa <- rep(0, iter)    # Number at site 1
xb <- rep(0, iter)    # Number at site 0
m <- rep(0, iter)     # Magnetization
time <- rep(0, iter)  # Recording times

# Event timers
# First N entries = spontaneous switching clocks
# Next N entries = copying clocks

t <- rep(0, 2 * N)

# Initialize spontaneous switching timers
for (j in 1:N) {
  t[j] <- rexp(1, rate = s)
}

# Initialize copying timers
for (j in (N + 1):(2 * N)) {
  t[j] <- rexp(1, rate = r)
}

# Simulation counters
T <- 0
Tprint <- 0
n <- 1

# Gillespie simulation loop
while (T < Tend) {
  ind <- which.min(t)      # Next event index
  T <- T + t[ind]          # Advance global time
  t <- t - t[ind]          # Decrement all timers
  
  if (ind <= N) {
    # Spontaneous flip event: individual switches site randomly
    x[ind] <- abs(x[ind] - 1)
    t[ind] <- rexp(1, rate = s)   # Reset timer for this individual
  } else {
    # Copying event: individual adopts another's site choice
    ind <- ind - N
    ind1 <- sample(1:N, 1)            # Sample random male
    x[ind] <- x[ind1]                 # Copy its site choice
    t[ind + N] <- rexp(1, rate = r)   # Reset copying timer
  }
  
  # Record state at each Tint interval
  while (T > Tprint && n <= iter) {
    xa[n] <- sum(x == 1)
    xb[n] <- sum(x == 0)
    m[n] <- (xa[n] - xb[n]) / N
    time[n] <- Tprint
    Tprint <- Tprint + Tint
    n <- n + 1
  }
}

# Create dataframe
df <- data.frame(time = time, site1 = xa, site2 = xb, 
                 proportion_site1 = xa / N, proportion_site2 = xb / N, magnetization = m)

# Save to CSV
write.csv(df, file = paste0("processed_data/abm_seed_", seed, ".csv"), row.names = FALSE)

# Plot
plot(time, m, type = "l", col = "cyan4", lwd = 2, xlab = "Time", ylab = "Magnetization", main = "Blackbuck Lekking Site Consensus Dynamics")
