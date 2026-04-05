
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
from sqlalchemy import label



# Galah data
fits_file = fits.open("galah_dr4_allstar_240705.fits")
galah_data = fits_file[1].data # extract the data
galah_df = pd.DataFrame(galah_data) # convert data to a dataframe

# Gaia data
fits_file1 = fits.open("galah_dr4_vac_wise_tmass_gaiadr3_240705.fits")
gaia_data = fits_file1[1].data # extract the data
gaia_df = pd.DataFrame(gaia_data) # convert data to a dataframe

# Merge galah and gaia data
merged_df = pd.merge(
    galah_df,
    gaia_df,
    left_on="gaiadr3_source_id",
    right_on="source_id",
    how="inner"
)
# Filters out stars with flags indicating potential issues with the data
clean_df = merged_df[
    (merged_df["flag_sp"] == 0) & 
    (merged_df["flag_sp_fit"] == 0)
].copy()




# metal poor stars with [Fe/H] < -0.5 are often considered to be part of the halo population of the Milky Way, 
# which is thought to contain older stars that formed early in the galaxy's history. 
halo_stars = clean_df[clean_df["fe_h"] < -0.5]
# metal rich stars with [Fe/H] > -0.5 are often considered to be part of the disk population of the Milky Way, 
# which is thought to contain younger stars that formed more recently in the galaxy's history. 
disk_stars = clean_df[clean_df["fe_h"] > -0.5]



# massive stars with masses greater than 2 solar masses are often associated with younger stellar populations
high_mass = clean_df[clean_df["mass"] > 2.0]
# low mass stars with masses less than 1 solar mass are often associated with older stellar populations
low_mass = clean_df[clean_df["mass"] < 1.0]



# calculates the [alpha/Fe] ratio for each star by averaging the abundances of magnesium, silicon, 
# calcium, and titanium relative to iron
clean_df["alpha_fe"] = (clean_df["mg_fe"] + clean_df["si_fe"] + clean_df["ca_fe"] + clean_df["ti_fe"]) / 4

# alpha elements are produced in large quantities by massive stars that end their lives as supernovae. 
# Stars with high [alpha/Fe] ratios are often associated with older stellar populations
alpha_stars = clean_df[clean_df["alpha_fe"] > 0.2]




# convert the equatorial coordinates (RA, Dec) to galactic coordinates (l, b)
coords = SkyCoord(ra=clean_df["ra_y"].values * u.degree, dec=clean_df["dec_y"].values * u.degree)
# put into dataframe
clean_df["l"] = coords.galactic.l.degree
clean_df["b"] = coords.galactic.b.degree

plt.scatter(clean_df["l"], clean_df["b"], s=1)
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("Spatial Map of Stars in the Milky Way")
plt.show()




# create a scatter plot of the galactic coordinates (l, b) of the stars
# colored by their radial velocity
# limit is set due to outliers that would skew colour scale
plt.scatter(
    clean_df["l"],
    clean_df["b"],
    c=clean_df["radial_velocity"],
    s=1,
    vmin=-100,
    vmax=100
)
plt.colorbar(label="Radial Velocity (km/s)")
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("Radial Velocity of Stars in the Milky Way")
plt.show()

# --------

# spatial distribution of disk and halo stars
# disk stars (metal-rich)
plt.scatter(
    disk_stars["l"],
    disk_stars["b"],
    s=1,
    color="blue",
    label="Disk",
    alpha=0.5
)
# halo stars (metal-poor)
plt.scatter(
    halo_stars["l"],
    halo_stars["b"],
    s=1,
    color="red",
    label="Halo",
    alpha=0.5
)
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("Spatial Distribution of Disk and Halo Stars")
plt.legend()
plt.show()


# create a scatter plot of the galactic coordinates (l, b) of the halo stars
# colored by their radial velocity
plt.scatter(
    halo_stars["l"],
    halo_stars["b"],
    c=halo_stars["radial_velocity"],
    s=1,
    vmin=-100,
    vmax=100
)
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("Halo Stars (Metal-poor)")
plt.colorbar(label="Radial Velocity")
plt.show()
# create a scatter plot of the galactic coordinates (l, b) of the disk stars
# colored by their radial velocity
plt.scatter(
    disk_stars["l"],
    disk_stars["b"],
    c=disk_stars["radial_velocity"],
    s=1,
    vmin=-100,
    vmax=100
)
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("Disk Stars (Metal-rich)")
plt.colorbar(label="Radial Velocity")
plt.show()



# create a histogram of the radial velocities of the halo and disk stars
plt.hist(disk_stars["radial_velocity"], bins=50, alpha=0.5, label="Disk")
plt.hist(halo_stars["radial_velocity"], bins=50, alpha=0.5, label="Halo")
plt.legend()
plt.xlabel("Radial Velocity")
plt.ylabel("Number of Stars")
plt.title("Radial Velocity Distribution")
plt.show()


# ----------


# low-mass stars
plt.scatter(
    low_mass["l"],
    low_mass["b"],
    c=low_mass["radial_velocity"],
    s=1,
    vmin=-100,
    vmax=100
)
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("Low-Mass Stars")
plt.colorbar(label="Radial Velocity")
plt.show()
# high-mass stars
plt.scatter(
    high_mass["l"],
    high_mass["b"],
    c=high_mass["radial_velocity"],
    s=1,
    vmin=-100,
    vmax=100
)
plt.xlabel("Galactic Longitude (l)")
plt.ylabel("Galactic Latitude (b)")
plt.title("High-Mass Stars")
plt.colorbar(label="Radial Velocity")
plt.show()



# create a histogram of the radial velocities of the low-mass and high-mass stars
plt.hist(low_mass["radial_velocity"], bins=50, alpha=0.5, label="Low mass")
plt.hist(high_mass["radial_velocity"], bins=50, alpha=0.5, label="High mass")
plt.legend()
plt.xlabel("Radial Velocity")
plt.ylabel("Number of Stars")
plt.title("Radial Velocity Distribution (Mass Groups)")
plt.show()




