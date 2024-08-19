-- Default Almanac
{ files =
  [ { crc32 = Some 1551416339 -- 0x6f2456aa
    , uri = "http://public-data.nyxspace.com/anise/de421.bsp"
    }
  , { crc32 = Some 2182330089
    , uri = "http://public-data.nyxspace.com/anise/v0.4/pck11.pca"
    }
  , { crc32 = Some 3454388861
    , uri = "http://public-data.nyxspace.com/anise/moon_pa_de440_200625.bpc"
    } 
    -- Download the latest (of time of writing) LRO definitive ephemeris from the public Nyx Space cloud.
    -- Note that the original file is in _big endian_ format, and my machine is little endian, so I've used the
    -- `bingo` tool from https://naif.jpl.nasa.gov/naif/utilities_PC_Linux_64bit.html to convert the original file
    -- to little endian and upload it to the cloud.
    -- Refer to https://naif.jpl.nasa.gov/pub/naif/pds/data/lro-l-spice-6-v1.0/lrosp_1000/data/spk/?C=M;O=D for original file.
  , { crc32 = Some 3882673077
    , uri =
        "http://public-data.nyxspace.com/nyx/examples/lrorg_2023349_2024075_v01_LE.bsp"
    }
  ]
}
