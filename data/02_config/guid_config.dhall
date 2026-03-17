< Activity :
    { disabled : Bool
    , guidance :
        Optional
          { disable_prop_mass : Bool
          , law :
              < FiniteBurn :
                  { end : Text
                  , frame : < Inertial | RCN | RIC | VNC >
                  , representation :
                      < Angles :
                          { azimuth :
                              < Constant : { a : Double }
                              | Linear : { a : Double, b : Double }
                              | Quadratic :
                                  { a : Double, b : Double, c : Double }
                              >
                          , elevation :
                              < Constant : { a : Double }
                              | Linear : { a : Double, b : Double }
                              | Quadratic :
                                  { a : Double, b : Double, c : Double }
                              >
                          }
                      | Vector : { _1 : Double, _2 : Double, _3 : Double }
                      >
                  , start : Text
                  , thrust_prct : Double
                  }
              | Kluever :
                  { max_eclipse_prct : Optional Double
                  , objectives :
                      List
                        { objective :
                            { additive_factor : Double
                            , desired_value : Double
                            , multiplicative_factor : Double
                            , parameter :
                                < BLTOF
                                | BdotR
                                | BdotT
                                | Cd
                                | Cr
                                | DryMass
                                | Element :
                                    < AoL
                                    | AoP
                                    | ApoapsisAltitude
                                    | ApoapsisRadius
                                    | BrouwerMeanShortAoP
                                    | BrouwerMeanShortEccentricity
                                    | BrouwerMeanShortInclination
                                    | BrouwerMeanShortMeanAnomaly
                                    | BrouwerMeanShortRAAN
                                    | BrouwerMeanShortSemiMajorAxis
                                    | C3
                                    | Custom
                                    | Declination
                                    | EccentricAnomaly
                                    | Eccentricity
                                    | Energy
                                    | EquinoctialH
                                    | EquinoctialK
                                    | EquinoctialLambda
                                    | EquinoctialP
                                    | EquinoctialQ
                                    | FlightPathAngle
                                    | HX
                                    | HY
                                    | HZ
                                    | Height
                                    | Hmag
                                    | HyperbolicAnomaly
                                    | Inclination
                                    | Latitude
                                    | Longitude
                                    | MeanAnomaly
                                    | PeriapsisAltitude
                                    | PeriapsisRadius
                                    | Period
                                    | RAAN
                                    | RightAscension
                                    | Rmag
                                    | SemiMajorAxis
                                    | SemiMinorAxis
                                    | SemiParameter
                                    | TrueAnomaly
                                    | TrueLongitude
                                    | VX
                                    | VY
                                    | VZ
                                    | VelocityDeclination
                                    | Vmag
                                    | X
                                    | Y
                                    | Z
                                    >
                                | Epoch
                                | GuidanceMode
                                | Isp
                                | PropMass
                                | Thrust
                                | TotalMass
                                >
                            , tolerance : Double
                            }
                        , weight : Double
                        }
                  }
              | Ruggiero :
                  { max_eclipse_prct : Optional Double
                  , objectives :
                      List
                        { efficiency : Double
                        , objective :
                            { additive_factor : Double
                            , desired_value : Double
                            , multiplicative_factor : Double
                            , parameter :
                                < BLTOF
                                | BdotR
                                | BdotT
                                | Cd
                                | Cr
                                | DryMass
                                | Element :
                                    < AoL
                                    | AoP
                                    | ApoapsisAltitude
                                    | ApoapsisRadius
                                    | BrouwerMeanShortAoP
                                    | BrouwerMeanShortEccentricity
                                    | BrouwerMeanShortInclination
                                    | BrouwerMeanShortMeanAnomaly
                                    | BrouwerMeanShortRAAN
                                    | BrouwerMeanShortSemiMajorAxis
                                    | C3
                                    | Custom
                                    | Declination
                                    | EccentricAnomaly
                                    | Eccentricity
                                    | Energy
                                    | EquinoctialH
                                    | EquinoctialK
                                    | EquinoctialLambda
                                    | EquinoctialP
                                    | EquinoctialQ
                                    | FlightPathAngle
                                    | HX
                                    | HY
                                    | HZ
                                    | Height
                                    | Hmag
                                    | HyperbolicAnomaly
                                    | Inclination
                                    | Latitude
                                    | Longitude
                                    | MeanAnomaly
                                    | PeriapsisAltitude
                                    | PeriapsisRadius
                                    | Period
                                    | RAAN
                                    | RightAscension
                                    | Rmag
                                    | SemiMajorAxis
                                    | SemiMinorAxis
                                    | SemiParameter
                                    | TrueAnomaly
                                    | TrueLongitude
                                    | VX
                                    | VY
                                    | VZ
                                    | VelocityDeclination
                                    | Vmag
                                    | X
                                    | Y
                                    | Z
                                    >
                                | Epoch
                                | GuidanceMode
                                | Isp
                                | PropMass
                                | Thrust
                                | TotalMass
                                >
                            , tolerance : Double
                            }
                        }
                  }
              >
          , thruster_model : Text
          }
    , name : Text
    , on_entry :
        Optional
          < Docking :
              { impulsive_maneuver :
                  Optional
                    { dv_x_km_s : Double
                    , dv_y_km_s : Double
                    , dv_z_km_s : Double
                    , local_frame : < Inertial | RCN | RIC | VNC >
                    }
              , increment_properties :
                  Optional
                    { drag : Optional { area_m2 : Double, coeff_drag : Double }
                    , mass :
                        Optional
                          { dry_mass_kg : Double
                          , extra_mass_kg : Double
                          , prop_mass_kg : Double
                          }
                    , srp :
                        Optional
                          { area_m2 : Double, coeff_reflectivity : Double }
                    }
              }
          | FrameSwap :
              { new_frame :
                  { ephemeris_id : Integer
                  , mu_km3_s2 : Optional Double
                  , orientation_id : Integer
                  , shape :
                      Optional
                        { polar_radius_km : Double
                        , semi_major_equatorial_radius_km : Double
                        , semi_minor_equatorial_radius_km : Double
                        }
                  }
              }
          | Staging :
              { decrement_properties :
                  Optional
                    { drag : Optional { area_m2 : Double, coeff_drag : Double }
                    , mass :
                        Optional
                          { dry_mass_kg : Double
                          , extra_mass_kg : Double
                          , prop_mass_kg : Double
                          }
                    , srp :
                        Optional
                          { area_m2 : Double, coeff_reflectivity : Double }
                    }
              , impulsive_maneuver :
                  Optional
                    { dv_x_km_s : Double
                    , dv_y_km_s : Double
                    , dv_z_km_s : Double
                    , local_frame : < Inertial | RCN | RIC | VNC >
                    }
              }
          >
    , propagator : Text
    }
| Terminate
>.Activity
  { disabled = False
  , guidance = Some
    { disable_prop_mass = False
    , law =
        < FiniteBurn :
            { end : Text
            , frame : < Inertial | RCN | RIC | VNC >
            , representation :
                < Angles :
                    { azimuth :
                        < Constant : { a : Double }
                        | Linear : { a : Double, b : Double }
                        | Quadratic : { a : Double, b : Double, c : Double }
                        >
                    , elevation :
                        < Constant : { a : Double }
                        | Linear : { a : Double, b : Double }
                        | Quadratic : { a : Double, b : Double, c : Double }
                        >
                    }
                | Vector : { _1 : Double, _2 : Double, _3 : Double }
                >
            , start : Text
            , thrust_prct : Double
            }
        | Kluever :
            { max_eclipse_prct : Optional Double
            , objectives :
                List
                  { objective :
                      { additive_factor : Double
                      , desired_value : Double
                      , multiplicative_factor : Double
                      , parameter :
                          < BLTOF
                          | BdotR
                          | BdotT
                          | Cd
                          | Cr
                          | DryMass
                          | Element :
                              < AoL
                              | AoP
                              | ApoapsisAltitude
                              | ApoapsisRadius
                              | BrouwerMeanShortAoP
                              | BrouwerMeanShortEccentricity
                              | BrouwerMeanShortInclination
                              | BrouwerMeanShortMeanAnomaly
                              | BrouwerMeanShortRAAN
                              | BrouwerMeanShortSemiMajorAxis
                              | C3
                              | Custom
                              | Declination
                              | EccentricAnomaly
                              | Eccentricity
                              | Energy
                              | EquinoctialH
                              | EquinoctialK
                              | EquinoctialLambda
                              | EquinoctialP
                              | EquinoctialQ
                              | FlightPathAngle
                              | HX
                              | HY
                              | HZ
                              | Height
                              | Hmag
                              | HyperbolicAnomaly
                              | Inclination
                              | Latitude
                              | Longitude
                              | MeanAnomaly
                              | PeriapsisAltitude
                              | PeriapsisRadius
                              | Period
                              | RAAN
                              | RightAscension
                              | Rmag
                              | SemiMajorAxis
                              | SemiMinorAxis
                              | SemiParameter
                              | TrueAnomaly
                              | TrueLongitude
                              | VX
                              | VY
                              | VZ
                              | VelocityDeclination
                              | Vmag
                              | X
                              | Y
                              | Z
                              >
                          | Epoch
                          | GuidanceMode
                          | Isp
                          | PropMass
                          | Thrust
                          | TotalMass
                          >
                      , tolerance : Double
                      }
                  , weight : Double
                  }
            }
        | Ruggiero :
            { max_eclipse_prct : Optional Double
            , objectives :
                List
                  { efficiency : Double
                  , objective :
                      { additive_factor : Double
                      , desired_value : Double
                      , multiplicative_factor : Double
                      , parameter :
                          < BLTOF
                          | BdotR
                          | BdotT
                          | Cd
                          | Cr
                          | DryMass
                          | Element :
                              < AoL
                              | AoP
                              | ApoapsisAltitude
                              | ApoapsisRadius
                              | BrouwerMeanShortAoP
                              | BrouwerMeanShortEccentricity
                              | BrouwerMeanShortInclination
                              | BrouwerMeanShortMeanAnomaly
                              | BrouwerMeanShortRAAN
                              | BrouwerMeanShortSemiMajorAxis
                              | C3
                              | Custom
                              | Declination
                              | EccentricAnomaly
                              | Eccentricity
                              | Energy
                              | EquinoctialH
                              | EquinoctialK
                              | EquinoctialLambda
                              | EquinoctialP
                              | EquinoctialQ
                              | FlightPathAngle
                              | HX
                              | HY
                              | HZ
                              | Height
                              | Hmag
                              | HyperbolicAnomaly
                              | Inclination
                              | Latitude
                              | Longitude
                              | MeanAnomaly
                              | PeriapsisAltitude
                              | PeriapsisRadius
                              | Period
                              | RAAN
                              | RightAscension
                              | Rmag
                              | SemiMajorAxis
                              | SemiMinorAxis
                              | SemiParameter
                              | TrueAnomaly
                              | TrueLongitude
                              | VX
                              | VY
                              | VZ
                              | VelocityDeclination
                              | Vmag
                              | X
                              | Y
                              | Z
                              >
                          | Epoch
                          | GuidanceMode
                          | Isp
                          | PropMass
                          | Thrust
                          | TotalMass
                          >
                      , tolerance : Double
                      }
                  }
            }
        >.FiniteBurn
          { end = "2010-12-22T01:30:45 UTC"
          , frame = < Inertial | RCN | RIC | VNC >.VNC
          , representation =
              < Angles :
                  { azimuth :
                      < Constant : { a : Double }
                      | Linear : { a : Double, b : Double }
                      | Quadratic : { a : Double, b : Double, c : Double }
                      >
                  , elevation :
                      < Constant : { a : Double }
                      | Linear : { a : Double, b : Double }
                      | Quadratic : { a : Double, b : Double, c : Double }
                      >
                  }
              | Vector : { _1 : Double, _2 : Double, _3 : Double }
              >.Vector
                { _1 = 1.0, _2 = 0.0, _3 = 0.0 }
          , start = "2010-12-22T01:30:00 UTC"
          , thrust_prct = 1.0
          }
    , thruster_model = "BiProp"
    }
  , name = "Finite Maneuver"
  , on_entry =
      None
        < Docking :
            { impulsive_maneuver :
                Optional
                  { dv_x_km_s : Double
                  , dv_y_km_s : Double
                  , dv_z_km_s : Double
                  , local_frame : < Inertial | RCN | RIC | VNC >
                  }
            , increment_properties :
                Optional
                  { drag : Optional { area_m2 : Double, coeff_drag : Double }
                  , mass :
                      Optional
                        { dry_mass_kg : Double
                        , extra_mass_kg : Double
                        , prop_mass_kg : Double
                        }
                  , srp :
                      Optional { area_m2 : Double, coeff_reflectivity : Double }
                  }
            }
        | FrameSwap :
            { new_frame :
                { ephemeris_id : Integer
                , mu_km3_s2 : Optional Double
                , orientation_id : Integer
                , shape :
                    Optional
                      { polar_radius_km : Double
                      , semi_major_equatorial_radius_km : Double
                      , semi_minor_equatorial_radius_km : Double
                      }
                }
            }
        | Staging :
            { decrement_properties :
                Optional
                  { drag : Optional { area_m2 : Double, coeff_drag : Double }
                  , mass :
                      Optional
                        { dry_mass_kg : Double
                        , extra_mass_kg : Double
                        , prop_mass_kg : Double
                        }
                  , srp :
                      Optional { area_m2 : Double, coeff_reflectivity : Double }
                  }
            , impulsive_maneuver :
                Optional
                  { dv_x_km_s : Double
                  , dv_y_km_s : Double
                  , dv_z_km_s : Double
                  , local_frame : < Inertial | RCN | RIC | VNC >
                  }
            }
        >
  , propagator = "Cislunar"
  }
