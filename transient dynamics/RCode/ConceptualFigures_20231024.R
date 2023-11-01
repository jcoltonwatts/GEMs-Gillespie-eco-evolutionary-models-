#### packages needed ####
library(lattice)
library(ggplot2)

#### default parameter values ####
s = 2 # d_min = s*b_max^2
b11 = 0.0001
b12 = 0.00012
d11 = 0.0001
d12 = 0.00012

b22 = 0.0001
b21 = 0.00012
d22 = 0.0001
d21 = 0.00012

#### Figure 1 ####
# Contourplot of NEEAS for b_max1 for combinations of R1 and R2.
# NEEAs correspond to b_maxs where derivative of LRS (i.e., R0_i) with respect to 
# b_max = 0.
# For comparison to ESS, for s = 2, the ESS(b_max*1) = ESS(b_max2*) = 0.25.
# For default parameters, this gives single-species ecological states
# K1 = K2 = (0.25 - 2*0.25^2) / (b11 + d11) = 625.

# define LRS of competitor 1 (birth rate divided by death rate)
# first save as an expression
LRS_R1 = expression((b_max1 - b11*R1 - b12*R2) / (s*b_max1^2 + d11*R1 + d12*R2))
# take derivative with respect to b_max1
dLRS_R1_b_max1 = D(LRS_R1, "b_max1")
dLRS_R1_b_max1
# Write that into a function: 
# dLRS as response,  b_max1 as independent variable,
# specify all other params in "parms" list.
f.dLRS_R1_b_max1 <- function(b_max1.s, LRS, pars) {
   with (as.list(c(b_max1.s, pars)),{
    dLRS = 1/(s * b_max1.s^2 + d11 * R1 + d12 * R2) - (b_max1.s - b11 * R1 - b12 * R2) * 
      (s * (2 * b_max1.s))/(s * b_max1.s^2 + d11 * R1 + d12 * R2)^2
    dLRS
    })
}

# do the same for LRS of second population (swapping subscripts)
f.dLRS_R2_b_max2 <- function(b_max2.s, LRS, pars) {
  with (as.list(c(b_max2.s, pars)),{
    dLRS = 1/(s * b_max2.s^2 + d22 * R2 + d21 * R1) - (b_max2.s - b22 * R2 -  b21 * R1) * 
      (s * (2 * b_max2.s))/(s * b_max2.s^2 + d22 * R2 + d21 * R1)^2
    dLRS
  })
}

abund.dat = as.data.frame(expand.grid("R1" = seq(1, 1000, by = 10),
                                      "R2" = seq(1, 1000, by = 10),
                                      "b_max1.star" = NA,
                                      "b_max2.star" = NA))
for (i in 1:nrow(abund.dat)){
  # For each row in abund.dat,
  # update pars using globally defined parameter values above
  # and ith combinations of R1 and R2 values from abund.dat.
  pars <- c(s = s, b11 = b11, b12 = b12, d11 = d11, d12 = d12,
            b22 = b22, b21 = b21, d22 = d22, d21 = d21,
            R1 = abund.dat$R1[i], R2 = abund.dat$R2[i])
  # Find NEEA (i.e., bmax for which dLRS w/ respect to bmax is 0)
  # for each competitor.
  abund.dat$b_max1.star[i] = uniroot(f.dLRS_R1_b_max1, pars = pars, c(0, 100), extendInt = "yes")$root
  abund.dat$b_max2.star[i] = uniroot(f.dLRS_R2_b_max2, pars = pars, c(0, 100), extendInt = "yes")$root
}

## create the isoclines for the competitors to show relationship to NEEAs
K1 = (0.25 - s*0.25^2) / (b11 + d11)
K1
K2 = (0.25 - s*0.25^2) / (b22 + d22)
K2
# Calculate "alpha" values from original LV competition model
# to simplify plotting of isoclines
alpha12 = (b12 + d12) / (b11 + d11)
alpha21 = (b21 + d21) / (b22 + d22)

# Plot, use grey scale, put colored line where NEEA = ESS.
# Note, this will generate a warning, because the isoclines
# for one competitor become negative when the other competitor's
# abundance is above its single-species equilibrium value, but
# we set the plot bounds to exclude these negative values.
quartz(height = 4, width = 5) # does not work on windows!
ggplot(abund.dat, aes(x = R1, y = R2, z = b_max1.star)) + 
  geom_contour_filled(binwidth = 0.05)  +
  geom_contour(breaks = c(0.25), color = "orange", linewidth = 1.5) +
  scale_color_manual(aesthetics = "fill", values = colorRampPalette(c("black",'grey99'))(13),
                     guide = guide_legend(title = "NEEA", reverse = T),
                     labels = seq(0.00, 0.65, 0.05)) +
  geom_line(inherit.aes = FALSE,
            aes(x = R1, y = K2 - alpha21*R1), lty = 1, linewidth = 0.75) +
  geom_line(inherit.aes = FALSE,
            aes(x = K1 - alpha12*R2, y = R2), lty = 2, linewidth = 0.5) +
  ylim(0, 1000)+
  xlim(0, 1000)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = "right",
        legend.text = element_text(size = rel(0.65)),
        legend.title = element_text(size = rel(0.7))) +
  labs(x = expression(paste("Abundance of focal population ",
                            italic("i"),
                            ", ",
                            italic(N[i]))),
       y = expression(paste("Abundance of competitor population ",
                            italic("j"),
                            ", ",
                            italic(N[j]))))

#### Figure 2 ####
# Plot of the real part of the dominant eigenvalue of the ecological
# system at the unstable equilibrium as a function of competitors' b_max values.

# Set up a range of b_maxs to consider.
# Only go up to 0.5, 
# because for b_max >= 1/s, b_max - d_min is <= 0 (& K is <= 0).
b1s = seq(0.01, 0.5, 0.005) 
b2s = seq(0.01, 0.5, 0.005)

# create dataframe using expand.grid()
qual.dat.to = as.data.frame(expand.grid(b1 = b1s,
                                        b2 = b2s,
                                        Q = NA,
                                        eig = NA))

# loop over b1s and b2s to characterize ecological system for each combo
for (i in 1:nrow(qual.dat.to)){
  # Calculate K1, K2, K1/a12, K2/a12.
  # For birth-death growth eqn., K1 = (b1 - d1) / (b11 + d11),
  # and a12 = (b12 + d12) / (b11 + d11)
  # (replace "1"s with "2"s for species 2).
  # Thus, K1/a12 = (b1 - d1) / (b12 + d12)
  this.b1 = qual.dat.to$b1[i]
  this.b2 = qual.dat.to$b2[i]
  this.d1 = s*this.b1^2
  this.d2 = s*this.b2^2
  K1 = (this.b1 - this.d1) / (b11 + d11)
  K2 = (this.b2 - this.d2) / (b22 + d22)
  K1_div_a12 = (this.b1 - this.d1) / (b12 + d12)
  K2_div_a21 = (this.b2 - this.d2) / (b21 + d21)
  a12 = (b12 + d12) / (b11 + d11)
  a21 = (b21 + d21) / (b22 + d22)
  
  # now check invasibility conditions and fill in "Q" variable with corresponding result
  # "0" = exclusion of one competitor, "1" = unstable equil., "2" = stable equil.
  if ((K1_div_a12 > K2) & (K2_div_a21 < K1)){
    qual.dat.to$Q[i] = 0
  }
  if ((K1_div_a12 < K2) & (K2_div_a21 > K1)){
    qual.dat.to$Q[i] = 0
  }
  if ((K1_div_a12 < K2) & (K2_div_a21 < K1)){
    qual.dat.to$Q[i] = 1
    # in this case, check the eigenvalue of the comm matrix at the equil
    equil = c(K1 - a12*(((a21*K1) - K2) / ((a12*a21) - 1)), (((a21*K1) - K2) / ((a12*a21) - 1)))
    # set R1 and R2 to the equilibrium values
    R1 = equil[1]
    R2 = equil[2]
    # specify ODEs as expressions
    dR1dt = expression((this.b1 - b11*R1 - b12*R2)*R1 - (s*this.d1^2 + d11*R1 + d12*R2)*R1)
    dR2dt = expression((this.b2 - b22*R2 - b21*R1)*R2 - (s*this.d2^2 + d22*R2 + d21*R1)*R2)
    # set up community matrix, where each row is one of the ODEs, and each column is the partial derivative
    # of that ODE w/ respect to R1, R2, etc., evaluated at the equilibrium R1 and R2 above
    Mat = matrix(data = c(eval(D(dR1dt, "R1")), eval(D(dR1dt, "R2")),
                          eval(D(dR2dt, "R1")), eval(D(dR2dt, "R2"))),
                 nrow = 2, ncol = 2, byrow = T)
    # find the eigenvalue with the largest real part, and save it in qual.dat$eig
    qual.dat.to$eig[i] = max(Re(eigen(Mat)$values))
  }
  if ((K1_div_a12 > K2) & (K2_div_a21 > K1)){
    qual.dat.to$Q[i] = 2
  }
}

# Plot real part of dominant eigenvalue as a function of b1 * b2.
# Note that this generates a warning b/c many combinations produce
# an ecological system with no unstable equilibrium, and thus
# there is no corresponding dominant eigenvalue to plot.
quartz(height = 5, width = 5) # does not work on windows!
ggplot(qual.dat.to, aes(x = b1, y = b2, z = eig)) +
  geom_contour_filled(binwidth = 0.02)  +
  scale_color_manual(aesthetics = "fill",
                     values = colorRampPalette(c("black",'grey85'))(7),
                     guide = guide_legend(title = expression(paste(Re(lambda[1]))),
                                          reverse = T),
                     labels = seq(0.00, 0.14, 0.02)) +
  ylim(0, 0.5)+
  xlim(0, 0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.position = c(0.875, 0.5),
        legend.text = element_text(size = rel(0.65)),
        legend.title = element_text(size = rel(0.7))) +
  labs(x = expression(paste("Maximum birth rate of focal population ",
                            italic("i"),
                            ", ",
                            italic(b[max[i]]))),
       y = expression(paste("Maximum birth rate of competitor population ",
                            italic("j"),
                            ", ",
                            italic(b[max[j]]))))
