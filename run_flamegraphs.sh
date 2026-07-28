#!/bin/bash

if ! command -v cargo-flamegraph &> /dev/null; then
    echo "cargo-flamegraph could not be found. Please install it using 'cargo install flamegraph'."
    exit 1
fi

echo "Running flamegraph for 'single translation'..."
CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=cc CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS="" cargo flamegraph -p nyx-space --bench anise_bench -o flamegraph_translation.svg -- "single translation" --bench

echo "Running flamegraph for 'single rotation'..."
CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=cc CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS="" cargo flamegraph -p nyx-space --bench anise_bench -o flamegraph_rotation.svg -- "single rotation" --bench

echo "Running flamegraph for 'full transformation'..."
CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=cc CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS="" cargo flamegraph -p nyx-space --bench anise_bench -o flamegraph_transformation.svg -- "full transformation" --bench

echo "Flamegraphs generated: flamegraph_translation.svg, flamegraph_rotation.svg, flamegraph_transformation.svg"
