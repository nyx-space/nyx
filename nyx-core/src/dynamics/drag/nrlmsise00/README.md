# NRLMSISE-00 empirical atmosphere model

This implementation is cloned from the [\`tobari\` crate](https://github.com/sksat/orts/tree/main/tobari),
developed by sksat. All credits for the Rust port of NRLMSISE-00 go to sksat and contributors of the \`tobari\` crate.

The implementation is a clean-room rewrite based on the following references:
- Picone, J.M. et al. (2002), "NRLMSISE-00 empirical model of the atmosphere:
  Statistical comparisons and scientific issues", J. Geophys. Res., 107(A12), 1468,
  doi:10.1029/2002JA009430
- Hedin, A.E. (1991), "Extension of the MSIS thermosphere model into the middle
  and lower atmosphere", J. Geophys. Res., 96(A2), 1159-1172.
- Hedin, A.E. (1987), "MSIS-86 thermospheric model",
  J. Geophys. Res., 92(A5), 4649-4662.

Coefficient values are from the official NRL distribution, treated as published data.
NRLMSISE-00 is believed to be in the public domain as a U.S. Government work
(17 U.S.C. § 105), though no explicit license was provided by NRL.

Note: MSIS is a registered trademark. This module uses the name "NRLMSISE-00"
for nominative fair use (identifying compatibility with the NRL model).

Validated against \`pymsis\` (official NRL Fortran wrapper, \`version=0\`).
