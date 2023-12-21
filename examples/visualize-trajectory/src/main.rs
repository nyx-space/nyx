use eframe;
use pretty_env_logger;

mod gui;
mod scenario;

use gui::NyxGui;

fn main() -> eframe::Result<()> {
    pretty_env_logger::init(); // log events to stdout (export RUST_LOG=info)

    // set default options for the UI window when it spawns
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default(),
        ..Default::default()
    };

    // run the app
    eframe::run_native(
        "Nyx Space",
        native_options,
        Box::new(|cc| Box::new(NyxGui::new(cc))),
    )
}
