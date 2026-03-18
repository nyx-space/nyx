{ propagators =
  [ { _1 = "Near Earth"
    , _2 =
      { accel_models =
        { gravity_field = Some
          { _1 =
            { degree = 21
            , filepath = "data/01_planetary/EGM2008_to2190_TideFree.gz"
            , gunzipped = True
            , order = 21
            }
          , _2 = { ephemeris_id = +399, orientation_id = +399 }
          }
        , point_masses = Some
          { celestial_objects = [ +399, +301 ]
          , correction =
              None { converged : Bool, stellar : Bool, transmit_mode : Bool }
          }
        }
      , force_models =
        { drag = Some
          { density =
              < Constant : Double
              | Exponential : { r0 : Double, ref_alt_m : Double, rho0 : Double }
              | StdAtm : { max_alt_m : Double }
              >.StdAtm
                { max_alt_m = 1000000.0 }
          , drag_frame =
            { ephemeris_id = +399
            , mu_km3_s2 = Some 398600.435436096
            , orientation_id = +399
            , shape = Some
              { polar_radius_km = 6356.75
              , semi_major_equatorial_radius_km = 6378.14
              , semi_minor_equatorial_radius_km = 6378.14
              }
            }
          , estimate = False
          }
        , solar_pressure =
            None
              { estimate : Bool
              , phi : Double
              , shadow_model :
                  { light_source :
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
                  , shadow_bodies :
                      List
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
              }
        }
      , method =
          < CashKarp45
          | DormandPrince45
          | DormandPrince78
          | RungeKutta4
          | RungeKutta89
          | Verner56
          >.RungeKutta89
      , options =
        { attempts = 50
        , error_ctrl =
            < LargestError
            | LargestState
            | LargestStep
            | RSSCartesianState
            | RSSCartesianStep
            | RSSState
            | RSSStep
            >.RSSCartesianStep
        , fixed_step = False
        , init_step = "1 min"
        , integration_frame =
            None
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
        , max_step = "45 min"
        , min_step = "1 ms"
        , tolerance = 0.000000000001
        }
      }
    }
  , { _1 = "Cislunar"
    , _2 =
      { accel_models =
        { gravity_field = Some
          { _1 =
            { degree = 8
            , filepath = "data/01_planetary/EGM2008_to2190_TideFree.gz"
            , gunzipped = True
            , order = 8
            }
          , _2 = { ephemeris_id = +399, orientation_id = +399 }
          }
        , point_masses = Some
          { celestial_objects = [ +399, +301 ]
          , correction =
              None { converged : Bool, stellar : Bool, transmit_mode : Bool }
          }
        }
      , force_models =
        { drag =
            None
              { density :
                  < Constant : Double
                  | Exponential :
                      { r0 : Double, ref_alt_m : Double, rho0 : Double }
                  | StdAtm : { max_alt_m : Double }
                  >
              , drag_frame :
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
              , estimate : Bool
              }
        , solar_pressure = Some
          { estimate = False
          , phi = 1367.0
          , shadow_model =
            { light_source =
              { ephemeris_id = +10
              , mu_km3_s2 = Some 132712440041.93938
              , orientation_id = +1
              , shape = Some
                { polar_radius_km = 696000.0
                , semi_major_equatorial_radius_km = 696000.0
                , semi_minor_equatorial_radius_km = 696000.0
                }
              }
            , shadow_bodies =
              [ { ephemeris_id = +399
                , mu_km3_s2 = Some 398600.435436096
                , orientation_id = +1
                , shape = Some
                  { polar_radius_km = 6356.75
                  , semi_major_equatorial_radius_km = 6378.14
                  , semi_minor_equatorial_radius_km = 6378.14
                  }
                }
              , { ephemeris_id = +301
                , mu_km3_s2 = Some 4902.800066163796
                , orientation_id = +1
                , shape = Some
                  { polar_radius_km = 1737.4
                  , semi_major_equatorial_radius_km = 1737.4
                  , semi_minor_equatorial_radius_km = 1737.4
                  }
                }
              ]
            }
          }
        }
      , method =
          < CashKarp45
          | DormandPrince45
          | DormandPrince78
          | RungeKutta4
          | RungeKutta89
          | Verner56
          >.RungeKutta89
      , options =
        { attempts = 50
        , error_ctrl =
            < LargestError
            | LargestState
            | LargestStep
            | RSSCartesianState
            | RSSCartesianStep
            | RSSState
            | RSSStep
            >.RSSCartesianStep
        , fixed_step = False
        , init_step = "1 min"
        , integration_frame =
            None
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
        , max_step = "45 min"
        , min_step = "1 ms"
        , tolerance = 0.000000000001
        }
      }
    }
  ]
, seq =
  [ { _1 = "2010-12-21T00:00:00 UTC"
    , _2 =
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
                              | Vector :
                                  { _1 : Double, _2 : Double, _3 : Double }
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
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      , increment_properties :
                          Optional
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
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
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
                            }
                      , impulsive_maneuver :
                          Optional
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      }
                  >
            , propagator : Text
            }
        | Terminate
        >.Activity
          { disabled = False
          , guidance =
              None
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
          , name = "Parking orbit checkout"
          , on_entry =
              None
                < Docking :
                    { impulsive_maneuver :
                        Optional
                          { dv_km_s : { _1 : Double, _2 : Double, _3 : Double }
                          , local_frame : < Inertial | RCN | RIC | VNC >
                          }
                    , increment_properties :
                        Optional
                          { drag :
                              Optional { area_m2 : Double, coeff_drag : Double }
                          , mass :
                              Optional
                                { dry_mass_kg : Double
                                , extra_mass_kg : Double
                                , prop_mass_kg : Double
                                }
                          , srp :
                              Optional
                                { area_m2 : Double
                                , coeff_reflectivity : Double
                                }
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
                          { drag :
                              Optional { area_m2 : Double, coeff_drag : Double }
                          , mass :
                              Optional
                                { dry_mass_kg : Double
                                , extra_mass_kg : Double
                                , prop_mass_kg : Double
                                }
                          , srp :
                              Optional
                                { area_m2 : Double
                                , coeff_reflectivity : Double
                                }
                          }
                    , impulsive_maneuver :
                        Optional
                          { dv_km_s : { _1 : Double, _2 : Double, _3 : Double }
                          , local_frame : < Inertial | RCN | RIC | VNC >
                          }
                    }
                >
          , propagator = "Near Earth"
          }
    }
  , { _1 = "2010-12-21T01:30:00 UTC"
    , _2 =
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
                              | Vector :
                                  { _1 : Double, _2 : Double, _3 : Double }
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
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      , increment_properties :
                          Optional
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
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
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
                            }
                      , impulsive_maneuver :
                          Optional
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      }
                  >
            , propagator : Text
            }
        | Terminate
        >.Activity
          { disabled = False
          , guidance =
              None
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
          , name = "Separation and vehicle checkout"
          , on_entry = Some
              ( < Docking :
                    { impulsive_maneuver :
                        Optional
                          { dv_km_s : { _1 : Double, _2 : Double, _3 : Double }
                          , local_frame : < Inertial | RCN | RIC | VNC >
                          }
                    , increment_properties :
                        Optional
                          { drag :
                              Optional { area_m2 : Double, coeff_drag : Double }
                          , mass :
                              Optional
                                { dry_mass_kg : Double
                                , extra_mass_kg : Double
                                , prop_mass_kg : Double
                                }
                          , srp :
                              Optional
                                { area_m2 : Double
                                , coeff_reflectivity : Double
                                }
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
                          { drag :
                              Optional { area_m2 : Double, coeff_drag : Double }
                          , mass :
                              Optional
                                { dry_mass_kg : Double
                                , extra_mass_kg : Double
                                , prop_mass_kg : Double
                                }
                          , srp :
                              Optional
                                { area_m2 : Double
                                , coeff_reflectivity : Double
                                }
                          }
                    , impulsive_maneuver :
                        Optional
                          { dv_km_s : { _1 : Double, _2 : Double, _3 : Double }
                          , local_frame : < Inertial | RCN | RIC | VNC >
                          }
                    }
                >.Staging
                  { decrement_properties =
                      None
                        { drag :
                            Optional { area_m2 : Double, coeff_drag : Double }
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
                  , impulsive_maneuver = Some
                    { dv_km_s = { _1 = 0.000025, _2 = 0.0, _3 = 0.0 }
                    , local_frame = < Inertial | RCN | RIC | VNC >.VNC
                    }
                  }
              )
          , propagator = "Near Earth"
          }
    }
  , { _1 = "2010-12-22T01:30:00 UTC"
    , _2 =
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
                              | Vector :
                                  { _1 : Double, _2 : Double, _3 : Double }
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
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      , increment_properties :
                          Optional
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
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
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
                            }
                      , impulsive_maneuver :
                          Optional
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
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
                >.FiniteBurn
                  { end = "2010-12-22T01:30:45 UTC"
                  , frame = < Inertial | RCN | RIC | VNC >.VNC
                  , representation =
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
                          { dv_km_s : { _1 : Double, _2 : Double, _3 : Double }
                          , local_frame : < Inertial | RCN | RIC | VNC >
                          }
                    , increment_properties :
                        Optional
                          { drag :
                              Optional { area_m2 : Double, coeff_drag : Double }
                          , mass :
                              Optional
                                { dry_mass_kg : Double
                                , extra_mass_kg : Double
                                , prop_mass_kg : Double
                                }
                          , srp :
                              Optional
                                { area_m2 : Double
                                , coeff_reflectivity : Double
                                }
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
                          { drag :
                              Optional { area_m2 : Double, coeff_drag : Double }
                          , mass :
                              Optional
                                { dry_mass_kg : Double
                                , extra_mass_kg : Double
                                , prop_mass_kg : Double
                                }
                          , srp :
                              Optional
                                { area_m2 : Double
                                , coeff_reflectivity : Double
                                }
                          }
                    , impulsive_maneuver :
                        Optional
                          { dv_km_s : { _1 : Double, _2 : Double, _3 : Double }
                          , local_frame : < Inertial | RCN | RIC | VNC >
                          }
                    }
                >
          , propagator = "Cislunar"
          }
    }
  , { _1 = "2011-01-20T00:00:00 UTC"
    , _2 =
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
                              | Vector :
                                  { _1 : Double, _2 : Double, _3 : Double }
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
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      , increment_properties :
                          Optional
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
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
                            { drag :
                                Optional
                                  { area_m2 : Double, coeff_drag : Double }
                            , mass :
                                Optional
                                  { dry_mass_kg : Double
                                  , extra_mass_kg : Double
                                  , prop_mass_kg : Double
                                  }
                            , srp :
                                Optional
                                  { area_m2 : Double
                                  , coeff_reflectivity : Double
                                  }
                            }
                      , impulsive_maneuver :
                          Optional
                            { dv_km_s :
                                { _1 : Double, _2 : Double, _3 : Double }
                            , local_frame : < Inertial | RCN | RIC | VNC >
                            }
                      }
                  >
            , propagator : Text
            }
        | Terminate
        >.Terminate
    }
  ]
, thruster_sets = [ { _1 = "BiProp", _2 = { isp_s = 300.0, thrust_N = 25.0 } } ]
}
