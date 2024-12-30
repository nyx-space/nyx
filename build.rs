use shadow_rs::ShadowBuilder;

fn main() {
    ShadowBuilder::builder()
        .build()
        .expect("shadow init for nyx_space failed");
}
