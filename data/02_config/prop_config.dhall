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
