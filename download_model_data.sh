#!/bin/bash
mkdir -p data

# Define a function to download only if missing or if it's a Git LFS pointer
download_if_missing() {
    url="$1"
    # Ensure we strip any existing paths from the second argument to enforce the data/ directory
    filename=$(basename "$2")
    file="data/01_planetary/$filename"

    should_download=false

    if [ ! -f "$file" ]; then
        should_download=true
    else
        # File exists, check if it is an LFS pointer.
        # 1. Check size: LFS pointers are text files ~130 bytes.
        #    Real SPICE kernels are binary and much larger.
        #    We use 300 bytes as a safe upper bound for a pointer.
        fsize=$(wc -c < "$file" | tr -d ' ')

        if [ "$fsize" -lt 300 ]; then
            # 2. Check content: Look for the LFS signature url
            if grep -q "version" "$file"; then
                echo "Detected Git LFS pointer for $filename ($fsize bytes). Deleting..."
                rm "$file"
                should_download=true
            fi
        fi
    fi

    if [ "$should_download" = true ]; then
        echo "Downloading $filename..."
        wget -q -O "$file" "$url"
    else
        echo "Found $file (valid), skipping download."
    fi
}

# Download ANISE and SPICE data
download_if_missing "http://public-data.nyxspace.com/anise/de440s.bsp" "de440s.bsp"
download_if_missing "http://public-data.nyxspace.com/anise/v0.10/pck08.pca" "pck08.pca"
download_if_missing "http://public-data.nyxspace.com/anise/v0.10/pck11.pca" "pck11.pca"
download_if_missing "http://public-data.nyxspace.com/anise/v0.10/moon_fk.epa" "moon_fk.epa"
download_if_missing "http://public-data.nyxspace.com/anise/v0.10/moon_fk_de440.epa" "moon_fk_de440.epa"
download_if_missing "http://public-data.nyxspace.com/anise/moon_pa_de440_200625.bpc" "moon_pa_de440_200625.bpc"
download_if_missing "http://public-data.nyxspace.com/anise/ci/earth_latest_high_prec-2023-09-08.bpc" "earth_latest_high_prec.bpc"
download_if_missing "http://public-data.nyxspace.com/anise/ci/earth_longterm_000101_251211_250915.bpc" "earth_longterm_000101_251211_250915.bpc"
# Download Nyx models
download_if_missing "http://public-data.nyxspace.com/nyx/models/EGM2008_to2190_TideFree_sha.gz" "EGM2008_to2190_TideFree_sha.gz"
download_if_missing "http://public-data.nyxspace.com/nyx/models/EGM2008_to2190_TideFree.gz" "EGM2008_to2190_TideFree.gz"
download_if_missing "http://public-data.nyxspace.com/nyx/models/SpaceWeather-2021-01-01_2026-09-06.csv.gz" "SpaceWeather-2021-01-01_2026-09-06.csv.gz"
download_if_missing "http://public-data.nyxspace.com/nyx/models/Luna_jggrx_1500e_sha.tab.gz" "Luna_jggrx_1500e_sha.tab.gz"
download_if_missing "http://public-data.nyxspace.com/nyx/models/JGM3.cof.gz" "JGM3.cof.gz"
