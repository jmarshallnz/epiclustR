nzmg_to_nzgd1949 <- function(n,e)
{
  # constants for NZGD1949
  a  <-  6378388.0;
  n0 <- 6023150.0;
  e0 <- 2510000.0;
  lt0 <- -41.0;
  ln0 <- 173.0;

  B <- c(0.7557853228 + 0i,
         0.249204646 + 0.003371507i,
        -0.001541739 + 0.041058560i,
        -0.10162907  + 0.01727609i,
        -0.26623489  - 0.36249218i,
        -0.6870983   - 1.1651967i);

  C <- c(1.3231270439,
        -0.577245789 - 0.007809598i,
         0.508307513 - 0.112208952i,
        -0.15094762 + 0.18200602i,
         1.01418179 + 1.64497696i,
         1.9660549  + 2.5127645i);

  D <- c(1.5627014243,
         0.5185406398,
        -0.03333098,
        -0.1052906,
        -0.0368594,
         0.007317,
         0.01220,
         0.00394,
        -0.0013)

  z <- complex(real=(n-n0)/a, imaginary=(e-e0)/a);

  theta <- sum(C*z^(1:6))

  for (i in 1:2)
    theta <- (z + sum((1:5)*B[2:6]*theta^(2:6))) / sum((1:6)*B*theta^(0:5))

  dpsi <- Re(theta)
  dphi <- sum(D*dpsi^(1:9))
  lat <- lt0 + 10^5/3600*dphi
  long <- ln0 + 180/pi*Im(theta)
  return(c(lat,long))
}

tocartesian <- function(lat,long, a, f)
{
  e2 <- 2*f - f^2
  v <- a / sqrt(1 - e2*(sin(lat))^2)
  x <- v*cos(lat)*cos(long);
  y <- v*cos(lat)*sin(long);
  z <- v*(1-e2)*sin(lat);
  return(c(x,y,z))
}

topolar <- function(x, y, z, a, f)
{
  e2 <- 2*f - f^2
  p <- sqrt(x^2 + y^2)
  r <- sqrt(p^2 + z^2)

  tanmu <- z/p*((1-f)+e2*a/r)

  mu <- atan(tanmu)

  long <- atan(y/x);

  lat <- atan((z + e2*a/(1-f)*(sin(mu))^3) / (p-e2*a*(cos(mu))^3))
  return(c(lat, long %% pi))
}

nzgd1949_to_nzgd2000 <- function(lat, long)
{
  # convert to cartesian
  xyz <- tocartesian(lat*pi/180, long*pi/180, 6378388, 1/297);
  # transform coordinates
  T <- c(59.47, -5.04, 187.44);
  arcsec_to_rad <- pi/(3600*180)
  Rx <- -0.470 * arcsec_to_rad
  Ry <-  0.100 * arcsec_to_rad
  Rz <- -1.024 * arcsec_to_rad
  R <- t(matrix(c(1, Rz, -Ry, -Rz, 1, Rx, Ry, -Rx, 1), nrow=3))
  dS <- -4.5993

  xyz <- T + (1 + dS*1e-6)*R %*% xyz

  # convert back to polar
  return(topolar(xyz[1], xyz[2], xyz[3], 6378137, 1/298.257222101)*180/pi);
}

nztm2000_to_nzgd2000 <- function(N, E)
{
  # based on formulae at:
  # 
  # http://www.linz.govt.nz/geodetic/conversion-coordinates/projection-conversions/transverse-mercator-preliminary-computations
  #
  #

  a <- 6378137;         # A
  f <- 1/298.257222101; # 1/RF        
  long0 <- 173/180*pi;  # CM/rad2deg    Meridian
  k0 <- 0.9996;         # SF            Scale F
  phi0 <- 0;            # OLAT/rad2deg  Origin latitude
  E0 <- 1600000;        # FE            False easting
  N0 <- 10000000;       # FN            False northing
# 1.0     # utom ??

  b <- a*(1-f);

  e2 <- 2*f - f^2;

  # m0 is zero for NZTM2000
  m0 <- a*((1 - e2/4 - 3/64*e2^2 - 5/256*e2^3)*phi0 - 3/8*(e2 + 1/4*e2^2 + 15/128*e2^3)*sin(2*phi0) + 15/256*(e2^2 + 3/4*e2^3)*sin(4*phi0) + 35/3072*e2^3*sin(6*phi0));

  mp <- m0 + (N - N0)/k0;  # cn1

  n <- (a - b)/(a + b);    # = f/(2-f);

  G <- a*(1 - n)*(1 - n^2)*(1 + 9/4*n^2 + 225/64*n^4);

  sigma <- mp/G;

  phi_p <- sigma + (3/2*n - 27/32*n^3)*sin(2*sigma) 
                 + (21/16*n^2 - 55/32*n^4)*sin(4*sigma)
                 + 151/96*n^3*sin(6*sigma)
                 + 1097/512*n^4*sin(8*sigma);  #fphi

  eslt <- 1 - e2*sin(phi_p)^2;
  eta <- a / sqrt(eslt);
  rho <- eta * (1 - e2) / eslt;
  psi <- eta / rho;

  Ep <- E - E0;
  x <- Ep / (k0*eta);

  t <- tan(phi_p);

  t1 <- Ep * x/2;
  t2 <- Ep * x^3/24 * (-4*psi^2 + 9*psi*(1 - t^2) + 12*t^2);
  t3 <- Ep * x^5/720 * (8*psi^4*(11 - 24*t^2) - 12*psi^3*(21 - 71*t^2) + 15*psi^2*(15 - 98*t^2 + 15*t^4) + 180*psi*(5*t^2 - 3*t^4) + 360*t^4);
  t4 <- Ep * x^7/40320 * (1385 - 3633*t^2 + 4095*t^4 + 1575*t^6);

  lat <- phi_p + t/(k0*rho)*(-t1 + t2 - t3 + t4);

  t1 <- x;
  t2 <- x^3/6 * (psi + 2*t^2);
  t3 <- x^5/120 * (-4*psi^3*(1 - t^2) + psi^2*(9 - 68*t^2) + 72*psi*t^2 + 24*t^4);
  t4 <- x^7/5040 * (61 + 662*t^2 + 1320*t^4 + 720*t^6);

  long <- long0 + 1/cos(phi_p)*(t1 - t2 + t3 - t4);
  return(c(long %% pi, lat));
}

#print_minsec(nztm_to_wgs84(5380181.71, 1558376.32))
#print_minsec(nztm_to_wgs84(5040771.40, 1197666.98))

nztm2000_to_wgs84 <- function(n, e)
{
  # wgs84 is essentially the same as nzgd2000 at sealevel...
  out <- data.frame(long=0, lat=0)
  for (i in 1:length(n))
  {
    out[i,] <- nztm2000_to_nzgd2000(n[i], e[i]) * 180/pi
  }
  return(out)
}

nzmg_to_wgs84 <- function(n, e)
{
  # wgs84 is essentially the same as nzgd2000 at sealevel...
  out <- data.frame(long=0, lat=0)
  for (i in 1:length(n))
  {
    latlong <- nzmg_to_nzgd1949(n[i], e[i])
    out[i,] <- nzgd1949_to_nzgd2000(latlong[1], latlong[2])
  }
  return(out)
}

print_minsec <- function(v)
{
  for (i in 1:length(v))
  {
    t <- abs(v[i]*3600) # convert to seconds
    deg <- floor(t/3600) # degrees
    min <- floor((t - deg*3600)/60)
    sec <- t - deg*3600 - min*60
    cat(deg, min, sec, ' ')
  }
  cat('\n')
}
