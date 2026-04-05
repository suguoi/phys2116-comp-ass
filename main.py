
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
print("HELLO")
import os
print(os.getcwd())

# open the FITS file and reads the data
# galah stuff
fits_file = fits.open("phys2116 computational assignment/galah_dr4_allstar_240705.fits")
star_data = fits_file[1].data # extract the data
stars_df = pd.DataFrame(star_data) # convert data to a dataframe

# gaia stuff, i think it has the galah stuff too? already cross matched :O wow
fits_file1 = fits.open("phys2116 computational assignment/galah_dr4_vac_wise_tmass_gaiadr3_240705.fits")
gaia_data = fits_file1[1].data # extract the data
gaia_df = pd.DataFrame(gaia_data) # convert data to a dataframe

# merge data
merged_df = pd.merge(
    stars_df,
    gaia_df,
    left_on="gaiadr3_source_id",
    right_on="source_id",
    how="inner"
)
# filters out stars with flags indicating potential issues with the data
clean_df = merged_df[
    (merged_df["flag_sp"] == 0) & 
    (merged_df["flag_sp_fit"] == 0)
].copy()
# calculates the [alpha/Fe] ratio for each star by averaging the abundances of magnesium, silicon, calcium, and titanium relative to iron
clean_df["alpha_fe"] = (clean_df["mg_fe"] + clean_df["si_fe"] + clean_df["ca_fe"] + clean_df["ti_fe"]) / 4

# why does ai give detailed comments
# metal poor stars with [Fe/H] < -1.0 are often considered to be part of the halo population of the Milky Way, 
# which is thought to contain older stars that formed early in the galaxy's history. 
halo_stars = clean_df[clean_df["fe_h"] < -1.0]
# metal rich stars with [Fe/H] > -0.5 are often considered to be part of the disk population of the Milky Way, 
# which is thought to contain younger stars that formed more recently in the galaxy's history. 
disk_stars = clean_df[clean_df["fe_h"] > -0.5]
# alpha elements are produced in large quantities by massive stars that end their lives as supernovae. 
# Stars with high [alpha/Fe] ratios are often associated with older stellar populations
alpha_stars = clean_df[clean_df["alpha_fe"] > 0.2]

# massive stars with masses greater than 2 solar masses are often associated with younger stellar populations,
high_mass = clean_df[clean_df["mass"] > 2.0]
# low mass stars with masses less than 1 solar mass are often associated with older stellar populations
low_mass = clean_df[clean_df["mass"] < 1.0]

# visualise the relationship between [Fe/H] and [Alpha/Fe] for the clean dataset using a scatter plot
plt.scatter(clean_df["fe_h"], clean_df["alpha_fe"], s=1)
plt.xlabel("[Fe/H]")
plt.ylabel("[Alpha/Fe]")
plt.show()


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


# create a scatter plot of the galactic coordinates (l, b) of the stars
# colored by their proper motion, which is calculated as the total proper motion from the pmra
clean_df["pm_total"] = np.sqrt(clean_df["pmra"]**2 + clean_df["pmdec"]**2)
plt.scatter(
    clean_df["l"],
    clean_df["b"],
    c=clean_df["pm_total"],
    s=1,
    vmin=0,
    vmax=50
)
plt.colorbar(label="Proper Motion")
plt.show()